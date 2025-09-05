# Streamlit app skeleton v0.11 — Switch to TCGA_Classification direct filter + right-side TCGA code reference

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional

# --- Defaults ---
DEFAULT_GENES = ["KRAS", "BRAF", "PIK3CA", "TP53", "STK11", "KEAP1", "SMAD4", "CDKN2A"]

# --- Normalization helper ---
def normalize_cell_name(name: str) -> str:
    if not isinstance(name, str):
        return ""
    return name.strip().replace(" ", "").replace("-", "").upper()

# --- Load data ---
@st.cache_data(show_spinner=True)
def load_gdsc(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)
    df = df.rename(columns=lambda x: x.strip())
    if "IC50 (uM)" in df.columns:
        df = df.rename(columns={"IC50 (uM)": "IC50_uM"})
    return df

@st.cache_data(show_spinner=True)
def load_models(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns=lambda x: x.strip())
    depmap_col = None
    for c in ["DepMap_ID", "depmap_id", "ModelID", "modelid"]:
        if c in df.columns:
            depmap_col = c
            break
    if depmap_col is None:
        st.error("Model_DepMap.csv missing DepMap_ID/ModelID column.")
        st.stop()
    cell_candidates = [
        "CellLineName","Cell Line Name","StrippedCellLineName","stripped_cell_line_name","Model Name","CCLEName"
    ]
    cell_col = None
    for c in cell_candidates:
        if c in df.columns:
            cell_col = c
            break
    tcga_candidates = ["TCGA_Classification","tcga_classification","primary_disease","lineage","OncotreePrimaryDisease"]
    tcga_col = None
    for c in tcga_candidates:
        if c in df.columns:
            tcga_col = c
            break
    rename_map = {depmap_col: "DepMap_ID"}
    if cell_col:
        rename_map[cell_col] = "CellLine"
    if tcga_col:
        rename_map[tcga_col] = "TCGA_Classification"
    df = df.rename(columns=rename_map)
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    return df

@st.cache_data(show_spinner=True)
def load_ccle(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns=lambda x: x.strip())
    return df

@st.cache_data(show_spinner=True)
def build_gdsc_joined(gdsc: pd.DataFrame, models: pd.DataFrame) -> pd.DataFrame:
    df = gdsc.copy()
    if "Cell Line Name" in df.columns:
        df["__norm_cell"] = df["Cell Line Name"].map(normalize_cell_name)
    if "__norm_cell" in df.columns and "__norm_cell" in models.columns:
        df = df.merge(models[["DepMap_ID","__norm_cell","TCGA_Classification"]], on="__norm_cell", how="left")
    df["TumorGroup"] = df.get("TCGA_Classification")
    return df

# --- Sidebar ---
with st.sidebar:
    st.header("Filters")
    gene = st.selectbox("Gene", DEFAULT_GENES, index=0)
    # Build TCGA class list from models
    model_path = Path("data/Model_DepMap.csv")
    models_tmp = pd.read_csv(model_path, low_memory=False)
    tcga_classes = sorted([c for c in models_tmp.get("TCGA_Classification", models_tmp.get("tcga_classification", pd.Series())).dropna().unique()])
    tumors = st.multiselect("TCGA Classification (direct)", tcga_classes, default=[])
    drug_query = st.text_input("Drug name contains (optional)", value="")
    top_n = st.number_input("Top N bars to show", min_value=10, max_value=200, value=200, step=10)
    run_button = st.button("Run")
    st.caption("Duplicates per (Drug, Cell line) are collapsed by **median IC50** for stability.")

# --- Load datasets ---
gdsc = load_gdsc(Path("data/GDSC1_IC50.csv"))
models = load_models(Path("data/Model_DepMap.csv"))
ccle = load_ccle(Path("data/CCLE_KRAS geneset.csv"))

# --- Main logic ---
if run_button:
    joined = build_gdsc_joined(gdsc, models)
    if tumors:
        joined = joined[joined["TCGA_Classification"].isin(tumors)]
    if gene in ccle.columns:
        mut_df = ccle[["DepMap_ID", gene]].rename(columns={gene: "Mut"})
        joined = joined.merge(mut_df, on="DepMap_ID", how="left")
    ic_col = "IC50_uM"
    plot_df = joined.dropna(subset=[ic_col]).copy()

    # Collapse duplicates by median
    if "DepMap_ID" in plot_df.columns:
        grp_keys = ["Drug_Name","DepMap_ID"] if "Drug_Name" in plot_df.columns else ["DepMap_ID"]
    else:
        grp_keys = ["Drug_Name","CellLine"] if "Drug_Name" in plot_df.columns else ["CellLine"]

    def _agg_block(g: pd.DataFrame) -> pd.Series:
        ic = g[ic_col].astype(float).median()
        out = {
            ic_col: ic,
            "Mut": bool(g.get("Mut", pd.Series([False])).any()),
            "CellLine": g.get("CellLine", pd.Series([None])).iloc[0],
            "TCGA_Classification": g.get("TCGA_Classification", pd.Series([None])).iloc[0],
            "TumorGroup": g.get("TumorGroup", pd.Series([None])).iloc[0],
            "DepMap_ID": g.get("DepMap_ID", pd.Series([None])).iloc[0],
            "Drug_Name": g.get("Drug_Name", pd.Series([None])).iloc[0],
            "_dup_count": len(g),
        }
        return pd.Series(out)

    plot_df = plot_df.groupby(grp_keys, as_index=False).apply(_agg_block).reset_index(drop=True)

    dup_before = len(joined)
    dup_after = len(plot_df)
    collapsed = dup_before - dup_after

    # Diagnostics
    st.subheader("Diagnostics")
    total = len(plot_df)
    depmap_matched = plot_df["DepMap_ID"].notna().sum() if "DepMap_ID" in plot_df.columns else 0
    mut_count = int(plot_df.get("Mut", pd.Series(dtype=bool)).sum()) if "Mut" in plot_df.columns else 0
    st.markdown(f"**Rows:** {total} | **DepMap_ID matched:** {depmap_matched} | **Mutants:** {mut_count} | **Collapsed duplicates:** {collapsed}")

    # Filtered table
    st.subheader("Filtered table (Median-collapsed)")
    show_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","TumorGroup",ic_col,"Mut"] if c in plot_df.columns]
    st.dataframe(plot_df[show_cols], use_container_width=True, hide_index=True)

    # Raw replicates
    with st.expander("Show raw replicates (before collapsing)"):
        raw_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","TumorGroup",ic_col,"Mut"] if c in joined.columns]
        st.dataframe(joined[raw_cols], use_container_width=True, hide_index=True)

    # Plot
    st.subheader("IC50 distribution (lower is more sensitive)")
    if plot_df.empty:
        st.warning("No rows to plot. Try broadening filters or check dataset mappings.")
    else:
        fig, ax = plt.subplots(figsize=(12,5))
        x_labels = plot_df.get("CellLine", plot_df.get("DepMap_ID", pd.Series(range(len(plot_df)))))
        y = plot_df[ic_col].astype(float)
        mut_series = plot_df.get("Mut", pd.Series(False, index=plot_df.index))
        colors = ["red" if bool(m) else "blue" for m in mut_series]
        ax.bar(x_labels.astype(str), y, color=colors)
        from matplotlib.patches import Patch
        legend_handles = [Patch(color="red", label="Mutant"), Patch(color="blue", label="WT")]
        ax.legend(handles=legend_handles, loc="best")
        ax.set_ylabel("IC50 (µM)")
        ax.set_xlabel("Cell line")
        title_drug = drug_query if drug_query else "(all drugs)"
        title_tumor = ", ".join(tumors) if tumors else "All"
        ax.set_title(f"{title_drug} — {title_tumor} — Gene: {gene}")
        ax.tick_params(axis='x', labelrotation=90)
        st.pyplot(fig, clear_figure=True)

    # TCGA reference
    with st.expander("TCGA code reference (common groups)", expanded=False):
        st.markdown("""
        - **LUAD** = Lung adenocarcinoma (subset of NSCLC)  
        - **LUSC** = Lung squamous cell carcinoma (subset of NSCLC)  
        - **COREAD** = Colorectal adenocarcinoma (CRC)  
        - **BRCA** = Breast invasive carcinoma  
        - **SKCM** = Skin cutaneous melanoma  
        - **GBM** = Glioblastoma multiforme  
        - **PAAD** = Pancreatic adenocarcinoma (PDAC)
        """)
