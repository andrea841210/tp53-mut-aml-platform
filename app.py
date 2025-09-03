import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os

st.set_page_config(page_title="TP53 Node-1 Demo • LAML IC50 Explorer", layout="wide")
st.title("TP53 Node-1 Demo • LAML IC50 Explorer")
st.set_page_config(page_title="Node-1 Demo • LAML IC50 Explorer", layout="wide")
st.title("Node-1 Demo • LAML IC50 Explorer")
st.caption("Andrea × GPT Mirror Project — SO-06a debug build")

# -------------------------------
# Helpers
# -------------------------------
def normalize_columns(df: pd.DataFrame):
    mapping = {}
    new_cols = []
    for c in df.columns:
        norm = (
            str(c)
            .strip()
            .replace("\u00A0", " ")
            .replace(" ", "_")
            .replace("-", "_")
            .replace("/", "_")
            .replace("(", "")
            .replace(")", "")
        )
        norm = "".join(ch for ch in norm if ch.isalnum() or ch == "_").lower()
        mapping[c] = norm
        new_cols.append(norm)
    out = df.copy()
    out.columns = new_cols
    return out, mapping

def guess_col(df: pd.DataFrame, candidates):
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    no_us = {c: c.replace("_", "") for c in df.columns}
    for cand in candidates:
        k = cand.replace("_", "")
        for orig, transformed in no_us.items():
            if transformed == k:
                return orig
    return None

def summarize_columns(label, df, mapping):
    with st.expander(f"🔎 {label} — column inspection", expanded=False):
        st.write("**Original → Normalized** (first 40 shown):")
        view = pd.DataFrame({"original": list(mapping.keys()), "normalized": list(mapping.values())})
        st.dataframe(view.head(40), use_container_width=True)
        st.write("**All normalized columns:**", list(df.columns))

def build_mutation_status(mut_df: pd.DataFrame, gene_name: str):
    gene_col = guess_col(mut_df, ["gene", "gene_name", "hugo_symbol", "symbol"]) or "gene"
    depmap_col = guess_col(mut_df, ["depmap_id", "depmapid", "modelid", "model_id"]) or "depmap_id"
    status_col = guess_col(mut_df, ["mut_wt", "mutwt", "mutation_status", "status"])
    missing = []
    for needed in [gene_col, depmap_col]:
        if needed not in mut_df.columns:
            missing.append(needed)
    if missing:
        raise ValueError(f"Mutation file is missing required columns: {missing}.")
    if gene_col in mut_df.columns:
        tmp = mut_df[mut_df[gene_col].astype(str).str.upper() == gene_name.upper()].copy()
    else:
        tmp = mut_df.copy()
    if tmp.empty:
        tmp = mut_df[[depmap_col]].drop_duplicates().copy()
        tmp["mutwt"] = "WT"
    else:
        if status_col and status_col in tmp.columns:
            s = tmp[[depmap_col, status_col]].copy()
            s.rename(columns={status_col: "mutwt", depmap_col: "depmap_id"}, inplace=True)
            s["mutwt"] = s["mutwt"].astype(str).str.strip().str.upper().replace({"MUTANT": "MUT", "WILDTYPE": "WT"})
            s.loc[~s["mutwt"].isin(["MUT", "WT"]), "mutwt"] = "MUT"
            s = s.drop_duplicates("depmap_id")
        else:
            s = tmp[[depmap_col]].copy()
            s = s.drop_duplicates()
            s.rename(columns={depmap_col: "depmap_id"}, inplace=True)
            s["mutwt"] = "MUT"
    return s

# -------------------------------
# Mode selection
# -------------------------------
mode = st.sidebar.radio("Data source", ["Lite (built-in)", "Advanced (upload)"], index=0)

# Sidebar Inputs
gdsc_file = st.sidebar.file_uploader("GDSC IC50 table (CSV or Excel)", type=["csv", "xlsx", "xls"])
mut_file = st.sidebar.file_uploader("Mutation table (CSV)", type=["csv"])
map_file = st.sidebar.file_uploader("DepMap mapper (optional, CSV)", type=["csv"])

col_a, col_b = st.sidebar.columns(2)
with col_a:
    drug_name = st.text_input("Drug name", value="Rapamycin")
with col_b:
    st.markdown("**Gene selection (single-gene test)**")

# Independent dropdowns for TP53 / FLT3 / NPM1
opt_tp53 = st.sidebar.selectbox("TP53", ["Analyze", "Skip"], index=0)
opt_flt3 = st.sidebar.selectbox("FLT3", ["Skip", "Analyze"], index=0)
opt_npm1 = st.sidebar.selectbox("NPM1", ["Skip", "Analyze"], index=0)

# Determine selected gene (single active)
selected = [g for g, o in [("TP53", opt_tp53), ("FLT3", opt_flt3), ("NPM1", opt_npm1)] if o == "Analyze"]
if len(selected) == 0:
    gene_name = "TP53"  # default fallback
elif len(selected) > 1:
    st.warning("Multiple genes selected; using the first one only for this test run.")
    gene_name = selected[0]
else:
    gene_name = selected[0]

# --- helpers for reading file or path ---
def read_any_table(obj):
    if isinstance(obj, str):
        name = obj.lower()
        if name.endswith((".xlsx", ".xls")):
            return pd.read_excel(obj, engine="openpyxl")
        else:
            return pd.read_csv(obj)
    else:
        name = obj.name.lower()
        if name.endswith((".xlsx", ".xls")):
            return pd.read_excel(obj, engine="openpyxl")
        else:
            return pd.read_csv(obj)

DATA_DIR = "data"
DEFAULT_GDSC_PATH = os.path.join(DATA_DIR, "GDSC1_LAML.xlsx")
DEFAULT_MUT_PATH  = os.path.join(DATA_DIR, "CCLE_TP53xNPMxFLT3.csv")
DEFAULT_MAP_PATH  = os.path.join(DATA_DIR, "Model_DepMap.csv")

if mode == "Lite (built-in)":
    gdsc_source = DEFAULT_GDSC_PATH
    mut_source  = DEFAULT_MUT_PATH
    map_source  = DEFAULT_MAP_PATH
    cancer_type = "LAML"
    st.caption("Lite mode: built-in LAML data. You can change Drug/Gene; Cancer type is fixed to LAML.")
else:
    gdsc_source = gdsc_file
    mut_source  = mut_file
    map_source  = map_file
    cancer_type = st.text_input("Cancer type (TCGA class)", value="LAML")

# -------------------------------
# Main workflow
# -------------------------------
run = st.sidebar.button("Run")
if run:
    try:
        if gdsc_source is None or mut_source is None:
            st.warning("Please provide GDSC and Mutation data.")
            st.stop()

        df_gdsc_raw = read_any_table(gdsc_source)
        df_mut_raw  = read_any_table(mut_source)
        df_gdsc, map_gdsc = normalize_columns(df_gdsc_raw)
        df_mut,  map_mut  = normalize_columns(df_mut_raw)

        if map_source is not None:
            df_map_raw = read_any_table(map_source)
            df_map, map_map = normalize_columns(df_map_raw)
        else:
            df_map, map_map = None, None

        summarize_columns("GDSC table", df_gdsc, map_gdsc)
        summarize_columns("Mutation table", df_mut, map_mut)
        if df_map is not None:
            summarize_columns("Mapper table", df_map, map_map)

        col_drug = guess_col(df_gdsc, ["drug_name", "drug", "compound"])
        col_cell = guess_col(df_gdsc, ["cell_line_name", "cell_line", "model_name"])
        col_ic50 = guess_col(df_gdsc, ["ic50_um", "ic50", "ln_ic50"])
        col_tcga = guess_col(df_gdsc, ["tcga_classification", "tcga", "cancer_type"])
        col_depmap_in_gdsc = guess_col(df_gdsc, ["depmap_id", "depmapid", "modelid"])

        df_g = df_gdsc[df_gdsc[col_drug].astype(str).str.upper() == drug_name.upper()].copy()
        if col_tcga:
            df_g = df_g[df_g[col_tcga].astype(str).str.upper() == cancer_type.upper()].copy()
        if df_g.empty:
            st.error("No GDSC rows found.")
            st.stop()

        mut_status = build_mutation_status(df_mut, gene_name)

        if not col_depmap_in_gdsc:
            if df_map is None:
                st.warning("DepMap_ID not found; please upload mapper.")
            else:
                map_cell = guess_col(df_map, ["cell_line_name", "cell", "model_name"]) or "cell_line_name"
                map_depmap = guess_col(df_map, ["depmap_id", "depmapid", "modelid"]) or "depmap_id"
                df_g = df_g.merge(df_map[[map_cell, map_depmap]].drop_duplicates(), left_on=col_cell, right_on=map_cell, how="left")
                if "depmap_id" not in df_g.columns and map_depmap in df_g.columns:
                    df_g.rename(columns={map_depmap: "depmap_id"}, inplace=True)
        else:
            if col_depmap_in_gdsc != "depmap_id":
                df_g.rename(columns={col_depmap_in_gdsc: "depmap_id"}, inplace=True)

        df_join = df_g.merge(mut_status, on="depmap_id", how="left")
        df_join["mutwt"].fillna("WT", inplace=True)

        out_cols = [c for c in [col_cell, "depmap_id", col_ic50, "mutwt"] if c in df_join.columns]
        st.subheader("Filtered IC50 table")
        st.dataframe(df_join[out_cols].sort_values(col_ic50, ascending=True), use_container_width=True)

        st.subheader("IC50 distribution by mutation status")
        fig, ax = plt.subplots(figsize=(10, 5))
        plot_df = df_join[[col_cell, col_ic50, "mutwt"]].dropna().sort_values(col_ic50)
        colors = plot_df["mutwt"].map({"MUT": "red", "WT": "blue"}).fillna("blue")
        ax.bar(plot_df[col_cell].astype(str), plot_df[col_ic50], color=colors)
        ax.set_xlabel("Cell line")
        ax.set_ylabel("IC50 (uM)")
        ax.set_title(f"{drug_name} in {cancer_type} cell lines — {gene_name} Mut (red) vs WT (blue)")
        ax.tick_params(axis='x', rotation=75)
        st.pyplot(fig)

        with st.expander("📄 Debug notes"):
            st.write("Rows plotted:", len(plot_df))
            st.write("Unique DepMap IDs:", plot_df.shape[0])

        st.success("Done. Inspect columns above if something looks off.")

    except Exception as e:
        st.error(f"Runtime error: {e}")
else:
    st.info("Upload files on the left (Advanced mode) or use built-in data (Lite mode), set parameters, and click **Run**.")

st.markdown("---")
st.caption("Debug strategy enabled: auto-inspection + fallback logic for DepMap_ID alignment and Mut/WT derivation.")
