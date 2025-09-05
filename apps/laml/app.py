import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os

st.set_page_config(page_title="Node-1 Demo • LAML IC50 Explorer (SO-06c)", layout="wide")
st.title("Node-1 Demo • LAML IC50 Explorer (SO-06c)")
st.caption("Andrea × GPT Mirror Project — Multi-variant AND subset build • baseline = SO-06b single-gene")

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
    # soft match removing underscores
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

def build_mutation_status_single(mut_df: pd.DataFrame, gene_name: str):
    """Return depmap_id + mutwt for a single gene (MUT/WT)."""
    gene_col = guess_col(mut_df, ["gene", "gene_name", "hugo_symbol", "symbol"]) or "gene"
    depmap_col = guess_col(mut_df, ["depmap_id", "depmapid", "modelid", "model_id"]) or "depmap_id"
    status_col = guess_col(mut_df, ["mut_wt", "mutwt", "mutation_status", "status"])

    missing = []
    for needed in [gene_col, depmap_col]:
        if needed not in mut_df.columns:
            missing.append(needed)
    if missing:
        raise ValueError(f"Mutation file is missing required columns: {missing}.")

    # Filter to the requested gene if present, else pass-through
    if gene_col in mut_df.columns:
        tmp = mut_df[mut_df[gene_col].astype(str).str.upper() == gene_name.upper()].copy()
    else:
        tmp = mut_df.copy()

    if tmp.empty:
        # no rows for the gene -> assume everyone WT
        out = mut_df[[depmap_col]].drop_duplicates().copy()
        out.rename(columns={depmap_col: "depmap_id"}, inplace=True)
        out["mutwt"] = "WT"
        return out

    if status_col and status_col in tmp.columns:
        s = tmp[[depmap_col, status_col]].copy()
        s.rename(columns={status_col: "mutwt", depmap_col: "depmap_id"}, inplace=True)
        s["mutwt"] = (
            s["mutwt"].astype(str).str.strip().str.upper().replace({"MUTANT": "MUT", "WILDTYPE": "WT"})
        )
        s.loc[~s["mutwt"].isin(["MUT", "WT"]), "mutwt"] = "MUT"
        s = s.drop_duplicates("depmap_id")
        return s
    else:
        # presence implies MUT
        s = tmp[[depmap_col]].drop_duplicates().copy()
        s.rename(columns={depmap_col: "depmap_id"}, inplace=True)
        s["mutwt"] = "MUT"
        return s

def build_mutation_status_multi(mut_df: pd.DataFrame, genes: list):
    """Return depmap_id + per-gene MUT/WT + AND flag across selected genes."""
    parts = []
    for g in genes:
        sg = build_mutation_status_single(mut_df, g)
        sg = sg.rename(columns={"mutwt": f"status_{g.upper()}"})
        parts.append(sg)
    # outer merge to keep full universe
    out = None
    for p in parts:
        out = p if out is None else out.merge(p, on="depmap_id", how="outer")
    # any missing status -> treat as WT
    for g in genes:
        col = f"status_{g.upper()}"
        if col not in out.columns:
            out[col] = "WT"
        else:
            out[col] = out[col].fillna("WT").astype(str).str.upper()
    # AND across requested genes
    out["mut_and"] = (out[[f"status_{g.upper()}" for g in genes]] == "MUT").all(axis=1)
    return out

# -------------------------------
# Mode selection
# -------------------------------
mode = st.sidebar.radio("Data source", ["Lite (built-in)", "Advanced (upload)"], index=0)

# Sidebar Inputs
gdsc_file = st.sidebar.file_uploader("GDSC IC50 table (CSV or Excel)", type=["csv", "xlsx", "xls"])
mut_file = st.sidebar.file_uploader("Mutation table (CSV)", type=["csv"])
map_file = st.sidebar.file_uploader("DepMap mapper (optional, CSV)", type=["csv"])

# Subset controls
subset_mode = st.sidebar.selectbox("Subset mode", ["Single gene", "Multi-gene (AND)"])

# Presets for AND combos
preset = None
selected_genes = []
if subset_mode == "Single gene":
    gene_name = st.sidebar.selectbox("Gene name", ["TP53", "FLT3", "NPM1"], index=0)
else:
    preset = st.sidebar.selectbox(
        "Gene combo (AND)",
        ["TP53 + FLT3", "TP53 + NPM1", "TP53 + FLT3 + NPM1", "Custom…"],
        index=0,
    )
    if preset == "TP53 + FLT3":
        selected_genes = ["TP53", "FLT3"]
    elif preset == "TP53 + NPM1":
        selected_genes = ["TP53", "NPM1"]
    elif preset == "TP53 + FLT3 + NPM1":
        selected_genes = ["TP53", "FLT3", "NPM1"]
    else:
        selected_genes = st.sidebar.multiselect("Pick genes (AND)", ["TP53", "FLT3", "NPM1"], default=["TP53", "FLT3"]) or ["TP53", "FLT3"]

# Shared inputs
col_a, col_b = st.sidebar.columns(2)
with col_a:
    drug_name = st.text_input("Drug name", value="Rapamycin")
with col_b:
    only_and = st.checkbox("Only show AND-mutant lines", value=False)

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
    st.caption(
        "Lite mode: built-in LAML data. Sources — GDSC1 (downloaded 2025-09), "
        "CCLE mutation (DepMap 22Q2), Model mapper (DepMap 25Q2)."
    )
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

        # Load & normalize
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

        # Identify key columns
        col_drug = guess_col(df_gdsc, ["drug_name", "drug", "compound"]) or "drug_name"
        col_cell = guess_col(df_gdsc, ["cell_line_name", "cell_line", "model_name"]) or "cell_line_name"
        col_ic50 = guess_col(df_gdsc, ["ic50_um", "ic50", "ln_ic50"]) or "ic50_um"
        col_tcga = guess_col(df_gdsc, ["tcga_classification", "tcga", "cancer_type"]) or "tcga_classification"
        col_depmap_in_gdsc = guess_col(df_gdsc, ["depmap_id", "depmapid", "modelid"])  # may be None

        # Filter GDSC by drug & cancer type
        df_g = df_gdsc[df_gdsc[col_drug].astype(str).str.upper() == drug_name.upper()].copy()
        if col_tcga in df_g.columns and pd.notnull(cancer_type):
            df_g = df_g[df_g[col_tcga].astype(str).str.upper() == str(cancer_type).upper()].copy()
        if df_g.empty:
            st.error("No GDSC rows found after filtering by Drug/Cancer type.")
            st.stop()

        # Ensure depmap_id exists in df_g
        if not col_depmap_in_gdsc:
            if df_map is None:
                st.warning("DepMap_ID not found in GDSC; please upload mapper. Proceeding without join may reduce accuracy.")
            else:
                map_cell = guess_col(df_map, ["cell_line_name", "cell", "model_name"]) or "cell_line_name"
                map_depmap = guess_col(df_map, ["depmap_id", "depmapid", "modelid"]) or "depmap_id"
                df_g = df_g.merge(
                    df_map[[map_cell, map_depmap]].drop_duplicates(),
                    left_on=col_cell,
                    right_on=map_cell,
                    how="left",
                )
                if "depmap_id" not in df_g.columns and map_depmap in df_g.columns:
                    df_g.rename(columns={map_depmap: "depmap_id"}, inplace=True)
        else:
            if col_depmap_in_gdsc != "depmap_id":
                df_g.rename(columns={col_depmap_in_gdsc: "depmap_id"}, inplace=True)

        if "depmap_id" not in df_g.columns:
            st.error("DepMap ID could not be resolved in filtered GDSC table. Upload a mapper with cell line → DepMap_ID.")
            st.stop()

        # Build mutation status
        if subset_mode == "Single gene":
            mut_status = build_mutation_status_single(df_mut, gene_name)
            df_join = df_g.merge(mut_status, on="depmap_id", how="left")
            df_join["mutwt"] = df_join["mutwt"].fillna("WT")
            group_col = "mutwt"
            title_suffix = f"{gene_name} Mut (red) vs WT (blue)"
            color_map = {"MUT": "red", "WT": "blue"}
            if only_and:
                # In single-gene mode, only_and behaves as show only MUT
                df_join = df_join[df_join[group_col] == "MUT"].copy()
        else:
            # Multi-gene AND
            mut_multi = build_mutation_status_multi(df_mut, selected_genes)
            df_join = df_g.merge(mut_multi, on="depmap_id", how="left")
            df_join["mut_and"] = df_join["mut_and"].fillna(False)
            group_col = "mut_and"
            title_suffix = f"{' + '.join(selected_genes)} AND-MUT (red) vs others (blue)"
            color_map = {True: "red", False: "blue"}
            if only_and:
                df_join = df_join[df_join[group_col] == True].copy()

        # Output table
        out_cols = [c for c in [col_cell, "depmap_id", col_ic50, group_col] if c in df_join.columns]
        st.subheader("Filtered IC50 table")
        st.dataframe(df_join[out_cols].sort_values(col_ic50, ascending=True), use_container_width=True)

        # Plot
        st.subheader("IC50 distribution by subset")
        fig, ax = plt.subplots(figsize=(10, 5))
        plot_df = df_join[[col_cell, col_ic50, group_col]].dropna().sort_values(col_ic50)
        colors = plot_df[group_col].map(color_map).fillna("blue")
        ax.bar(plot_df[col_cell].astype(str), plot_df[col_ic50], color=colors)
        ax.set_xlabel("Cell line")
        ax.set_ylabel("IC50 (uM)")
        ax.set_title(f"{drug_name} in {cancer_type} cell lines — {title_suffix}")
        ax.tick_params(axis='x', rotation=75)
        st.pyplot(fig)

        # Debug section
        with st.expander("📄 Debug notes"):
            st.write("Rows plotted:", len(plot_df))
            st.write("Unique DepMap IDs:", plot_df.shape[0])
            if subset_mode == "Multi-gene (AND)":
                st.write("Selected genes:", selected_genes)
                cols_present = [c for c in df_join.columns if c.startswith("status_")]
                st.write("Per-gene status columns present:", cols_present)

        st.success("Done. If something looks off, expand the inspection panels above to verify column mappings.")

    except Exception as e:
        st.error(f"Runtime error: {e}")
else:
    st.info("Upload files on the left (Advanced mode) or use built-in data (Lite mode), set parameters, choose subset mode, and click **Run**.")

st.markdown("---")
st.caption("SO-06c debug strategy: column auto-inspection + robust DepMap_ID alignment + single/multi-gene subset logic (AND). Versions — GDSC1 (downloaded 2025-09), CCLE mutation (DepMap 22Q2), Model mapper (DepMap 25Q2).")

