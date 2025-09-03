import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

st.set_page_config(page_title="TP53 Node-1 Demo • LAML IC50 Explorer", layout="wide")
st.title("TP53 Node-1 Demo • LAML IC50 Explorer")
st.caption("Andrea × GPT Mirror Project — SO-06a debug build")

# -------------------------------
# Helpers
# -------------------------------

def normalize_columns(df: pd.DataFrame):
    """Return (df_copy, mapping) with normalized lowercase underscore cols and a map original->normalized."""
    mapping = {}
    new_cols = []
    for c in df.columns:
        norm = (
            str(c)
            .strip()
            .replace("\u00A0", " ")  # non-breaking space
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
    """Guess a column by a list of normalized candidate names."""
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    # fuzzy: remove underscores
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
    """
    Return a DataFrame with columns: [depmap_id, mutwt]
    Strategy:
      1) Prefer explicit Mut/WT column if present (accept many spellings)
      2) Else: derive Mut if any record exists for the gene (non-empty rows), otherwise WT
    """
    gene_col = guess_col(mut_df, ["gene", "gene_name", "hugo_symbol", "symbol"]) or "gene"
    depmap_col = guess_col(mut_df, ["depmap_id", "depmapid", "modelid", "model_id"]) or "depmap_id"
    status_col = guess_col(mut_df, ["mut_wt", "mutwt", "mutation_status", "status"])  # may be None

    missing = []
    for needed in [gene_col, depmap_col]:
        if needed not in mut_df.columns:
            missing.append(needed)
    if missing:
        raise ValueError(
            f"Mutation file is missing required columns after normalization: {missing}. "
            "Please ensure a DepMap/Model ID column and a Gene column exist."
        )

    # filter for target gene if present
    if gene_col in mut_df.columns:
        tmp = mut_df[mut_df[gene_col].astype(str).str.upper() == gene_name.upper()].copy()
    else:
        tmp = mut_df.copy()

    if tmp.empty:
        # fall back to whole table; mark all as WT later
        tmp = mut_df[[depmap_col]].drop_duplicates().copy()
        tmp["mutwt"] = "WT"
    else:
        if status_col and status_col in tmp.columns:
            s = tmp[[depmap_col, status_col]].copy()
            s.rename(columns={status_col: "mutwt", depmap_col: "depmap_id"}, inplace=True)
            # normalize values
            s["mutwt"] = s["mutwt"].astype(str).str.strip().str.upper().replace({"MUTANT": "MUT", "WILDTYPE": "WT"})
            s.loc[~s["mutwt"].isin(["MUT", "WT"]), "mutwt"] = "MUT"  # anything else treated as Mut
            s = s.drop_duplicates("depmap_id")
        else:
            # derive: presence of a row for the gene => Mut
            s = tmp[[depmap_col]].copy()
            s = s.drop_duplicates()
            s.rename(columns={depmap_col: "depmap_id"}, inplace=True)
            s["mutwt"] = "MUT"
    return s


# -------------------------------
# Sidebar: uploads & inputs
# -------------------------------
st.sidebar.header("Inputs")

gdsc_file = st.sidebar.file_uploader("GDSC IC50 table (CSV or Excel)", type=["csv", "xlsx", "xls"])  
mut_file = st.sidebar.file_uploader("Mutation table (CSV)", type=["csv"])  
map_file = st.sidebar.file_uploader("DepMap mapper (optional, CSV)", type=["csv"])  

col_a, col_b = st.sidebar.columns(2)
with col_a:
    drug_name = st.text_input("Drug name", value="Rapamycin")
with col_b:
    gene_name = st.text_input("Gene name", value="TP53")

cancer_type = st.text_input("Cancer type (TCGA class)", value="LAML")
run = st.sidebar.button("Run")

# -------------------------------
# Main workflow
# -------------------------------
if run:
    try:
        if gdsc_file is None or mut_file is None:
            st.warning("Please upload both a GDSC IC50 file and a Mutation file, then click Run.")
            st.stop()

        # Read GDSC
        if gdsc_file.name.lower().endswith((".xlsx", ".xls")):
            df_gdsc_raw = pd.read_excel(gdsc_file, engine="openpyxl")
        else:
            df_gdsc_raw = pd.read_csv(gdsc_file)
        df_gdsc, map_gdsc = normalize_columns(df_gdsc_raw)

        # Read mutation
        df_mut_raw = pd.read_csv(mut_file)
        df_mut, map_mut = normalize_columns(df_mut_raw)

        # Optional mapper
        if map_file is not None:
            df_map_raw = pd.read_csv(map_file)
            df_map, map_map = normalize_columns(df_map_raw)
        else:
            df_map, map_map = None, None

        # Inspect columns
        summarize_columns("GDSC table", df_gdsc, map_gdsc)
        summarize_columns("Mutation table", df_mut, map_mut)
        if df_map is not None:
            summarize_columns("Mapper table", df_map, map_map)

        # Locate key columns in GDSC
        col_drug = guess_col(df_gdsc, ["drug_name", "drug", "compound"])
        col_cell = guess_col(df_gdsc, ["cell_line_name", "cell_line", "model_name"])  # for label
        col_ic50 = guess_col(df_gdsc, ["ic50_um", "ic50", "ln_ic50", "ic50_uM"])  # accept variants
        col_tcga = guess_col(df_gdsc, ["tcga_classification", "tcga", "cancer_type"])  # for filter
        col_depmap_in_gdsc = guess_col(df_gdsc, ["depmap_id", "depmapid", "modelid"])  # may be None

        missing_gdsc = [x for x in [col_drug, col_cell, col_ic50] if x is None]
        if missing_gdsc:
            st.error(f"GDSC file is missing key columns: {missing_gdsc}. Please check your input table.")
            st.stop()

        # Filter GDSC by drug & cancer type if available
        df_g = df_gdsc[df_gdsc[col_drug].astype(str).str.upper() == drug_name.upper()].copy()
        if col_tcga:
            df_g = df_g[df_g[col_tcga].astype(str).str.upper() == cancer_type.upper()].copy()
        if df_g.empty:
            st.error("No GDSC rows found for the specified Drug and/or Cancer type.")
            st.stop()

        # Mutation status per DepMap
        mut_status = build_mutation_status(df_mut, gene_name)

        # If GDSC has no DepMap ID, try mapper join via cell line name
        if not col_depmap_in_gdsc:
            if df_map is None:
                st.warning(
                    "DepMap_ID not found in GDSC table. Upload a mapper CSV with at least columns "
                    "[cell_line_name, depmap_id] to enable joining.")
                # create dummy mapping from cell name if mutation table already has depmap ids mapped by cell line elsewhere
                # Proceed with left join on cell line if mutation file also has cell line names
            else:
                # guess mapper columns
                map_cell = guess_col(df_map, ["cell_line_name", "cell", "model_name"]) or "cell_line_name"
                map_depmap = guess_col(df_map, ["depmap_id", "depmapid", "modelid"]) or "depmap_id"
                if map_cell not in df_map.columns or map_depmap not in df_map.columns:
                    st.error("Mapper must contain cell_line_name and depmap_id columns (or equivalents).")
                    st.stop()
                df_g = df_g.merge(
                    df_map[[map_cell, map_depmap]].drop_duplicates(),
                    left_on=col_cell,
                    right_on=map_cell,
                    how="left",
                )
                if "depmap_id" not in df_g.columns and map_depmap in df_g.columns:
                    df_g.rename(columns={map_depmap: "depmap_id"}, inplace=True)
                if "depmap_id" not in df_g.columns:
                    st.error("Failed to create depmap_id column from mapper. Check column names.")
                    st.stop()
        else:
            # ensure the column is named depmap_id
            if col_depmap_in_gdsc != "depmap_id":
                df_g.rename(columns={col_depmap_in_gdsc: "depmap_id"}, inplace=True)

        # Join mutation status
        df_join = df_g.merge(mut_status, on="depmap_id", how="left")
        df_join["mutwt"].fillna("WT", inplace=True)

        # Final display table
        out_cols = [c for c in [col_cell, "depmap_id", col_ic50, "mutwt"] if c in df_join.columns]
        st.subheader("Filtered IC50 table")
        st.dataframe(df_join[out_cols].sort_values(col_ic50, ascending=True), use_container_width=True)

        # Plot
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

        st.success("Done. If something still looks off, expand the inspection panels above to verify column names.")

    except Exception as e:
        st.error(f"Runtime error: {e}")
        st.info(
            "Tip: This build auto-normalizes column names and will display the detected schema for each file.\n"
            "If your mutation file lacks an explicit Mut/WT column, the app marks samples as MUT when any row exists for the target gene.")
else:
    st.info("Upload files on the left, set parameters, and click **Run**.")

st.markdown("---")
st.caption("Debug strategy enabled: auto-inspection + fallback logic for DepMap_ID alignment and Mut/WT derivation.")

