# Node-1 IC50 Explorer • Solid Tumors (KRAS track)
# Streamlit app v1.0 — TCGA-class filter in sidebar, robust joins, median collapse, diagnostics

from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, Set
import re

import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

# ---------------------------
# Page meta
# ---------------------------
st.set_page_config(page_title="Node‑1 IC50 Explorer — KRAS track", layout="wide")
st.title("Node‑1 IC50 Explorer • Solid Tumors (KRAS track)")

DATA_DIR = Path("data")
FILES = {
    "CCLE": DATA_DIR / "CCLE_mutation_22Q2.csv",   # placeholder may be KRAS geneset
    "MODEL": DATA_DIR / "Model_DepMap.csv",
    "GDSC1": DATA_DIR / "GDSC1_IC50.csv",
}

def _fmt_exists(p: Path) -> str:
    return f"✅ {p.name}" if p.exists() else f"❌ {p.name}"

st.caption(
    "Datasets → "
    f"CCLE {_fmt_exists(FILES['CCLE'])} • "
    f"Model {_fmt_exists(FILES['MODEL'])} • "
    f"GDSC1 {_fmt_exists(FILES['GDSC1'])} | IC50 aggregation: **median**"
)

DEFAULT_GENES = ["KRAS", "BRAF", "PIK3CA", "TP53", "STK11", "KEAP1", "SMAD4", "CDKN2A"]

# ---------------------------
# Helpers
# ---------------------------

def find_col(df: pd.DataFrame, candidates: Iterable[str], *, required: bool = True) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        st.error(f"Required column not found. Tried: {list(candidates)}")
        st.stop()
    return None

_norm_cache: dict[str, str] = {}

def normalize_cell_name(s: str) -> str:
    if s in _norm_cache:
        return _norm_cache[s]
    n = re.sub(r"[^A-Za-z0-9]", "", str(s)).upper()
    _norm_cache[s] = n
    return n

# ---------------------------
# Loaders (cached)
# ---------------------------
@st.cache_data(show_spinner=True)
def load_models(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns=lambda x: x.strip())
    depmap_col = find_col(df, ["DepMap_ID", "depmap_id", "ModelID", "modelid"])  # ACH-xxxxx
    # Cell-line name candidates
    cell_candidates = [
        "CellLineName", "Cell Line Name", "stripped_cell_line_name", "StrippedCellLineName",
        "CCLEName", "Model Name"
    ]
    cell_col = find_col(df, cell_candidates, required=False)
    tcga_col = find_col(df, ["TCGA_Classification", "TCGA Classification", "tcga_classification"], required=False)
    rename_map = {depmap_col: "DepMap_ID"}
    if cell_col:
        rename_map[cell_col] = "CellLine"
    if tcga_col:
        rename_map[tcga_col] = "TCGA_Classification"
    df = df.rename(columns=rename_map)
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    # Also build extra normalized keys if available
    for alt in [c for c in ["CCLEName", "StrippedCellLineName", "stripped_cell_line_name"] if c in df.columns]:
        df[f"__norm_{alt}"] = df[alt].map(normalize_cell_name)
    return df

@st.cache_data(show_spinner=True)
def load_gdsc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns=lambda x: x.strip())
    rename_map = {}
    if "IC50 (uM)" in df.columns:
        rename_map["IC50 (uM)"] = "IC50_uM"
    if "Drug Name" in df.columns:
        rename_map["Drug Name"] = "Drug_Name"
    if "Cell Line Name" in df.columns:
        rename_map["Cell Line Name"] = "CellLine"
    if "TCGA Classification" in df.columns:
        rename_map["TCGA Classification"] = "TCGA_Classification"
    df = df.rename(columns=rename_map)
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    return df

@st.cache_data(show_spinner=True)
def load_ccle(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns=lambda x: x.strip())
    # Expect long format: DepMap_ID + Hugo_Symbol
    depmap_col = find_col(df, ["DepMap_ID", "depmap_id"], required=False)
    gene_col = find_col(df, ["Hugo_Symbol", "Gene", "gene_symbol"], required=False)
    if depmap_col and gene_col:
        return df.rename(columns={depmap_col: "DepMap_ID", gene_col: "Hugo_Symbol"})[["DepMap_ID", "Hugo_Symbol"]]
    else:
        # If a wide geneset was provided, try to melt it (boolean columns per gene)
        return df

@st.cache_data(show_spinner=True)
def build_join(gdsc: pd.DataFrame, models: pd.DataFrame) -> pd.DataFrame:
    df = gdsc.copy()
    # Cascade merges to enrich DepMap_ID and TCGA_Classification
    keys_to_try = [k for k in ["__norm_cell", "__norm_CCLEName", "__norm_StrippedCellLineName", "__norm_stripped_cell_line_name"] if k in models.columns]
    if "__norm_cell" not in df.columns and "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    for k in keys_to_try:
        if k in df.columns and k in models.columns:
            df = df.merge(models[["DepMap_ID", "TCGA_Classification", k]], on=k, how="left", suffixes=("", "_m"))
            if "DepMap_ID_m" in df.columns:
                df["DepMap_ID"] = df["DepMap_ID"].fillna(df["DepMap_ID_m"])
                df = df.drop(columns=["DepMap_ID_m"])
            if "TCGA_Classification_m" in df.columns and "TCGA_Classification" not in df.columns:
                df = df.rename(columns={"TCGA_Classification_m": "TCGA_Classification"})
    return df

@st.cache_data(show_spinner=True)
def mutated_depmap_ids(ccle: pd.DataFrame, gene: str) -> Set[str]:
    if {"DepMap_ID", "Hugo_Symbol"}.issubset(ccle.columns):
        sub = ccle[ccle["Hugo_Symbol"].str.upper() == gene.upper()]
        return set(sub["DepMap_ID"].astype(str).unique())
    # Wide format fallback: treat non-empty truthy cells as mutation flag
    if gene in ccle.columns and "DepMap_ID" in ccle.columns:
        sub = ccle[ccle[gene].astype(str).str.upper().isin(["1", "TRUE", "YES", "Y"])]["DepMap_ID"]
        return set(sub.astype(str).unique())
    return set()

# ---------------------------
# Sidebar (build TCGA list from models)
# ---------------------------
with st.sidebar:
    st.header("Filters")
    gene = st.selectbox("Gene", DEFAULT_GENES, index=0)
    models_for_list = load_models(FILES["MODEL"])  # cached
    if "TCGA_Classification" in models_for_list.columns:
        tcga_classes = sorted(models_for_list["TCGA_Classification"].dropna().astype(str).unique().tolist())
    else:
        tcga_classes = []
    tumors = st.multiselect("TCGA Classification (direct)", tcga_classes, default=[])
    drug_query = st.text_input("Drug name contains (optional)", value="")
    top_n = st.number_input("Top N bars to show", min_value=10, max_value=200, value=200, step=10)
    run_button = st.button("Run")
    st.caption("Duplicates per (Drug, Cell line) are collapsed by **median IC50** for stability.")

# Always show a compact TCGA reference (right side)
with st.expander("TCGA code reference (common groups)", expanded=False):
    st.markdown(
        "- **LUAD** = Lung adenocarcinoma (NSCLC)  \n"
        "- **LUSC** = Lung squamous cell carcinoma (NSCLC)  \n"
        "- **COREAD** = Colorectal adenocarcinoma (CRC)  \n"
        "- **BRCA** = Breast invasive carcinoma  \n"
        "- **SKCM** = Skin cutaneous melanoma  \n"
        "- **GBM** = Glioblastoma multiforme  \n"
        "- **PAAD** = Pancreatic adenocarcinoma (PDAC)"
    )

# ---------------------------
# Main run
# ---------------------------
if not run_button:
    st.stop()

# Load datasets (cached)
gdsc = load_gdsc(FILES["GDSC1"])  # expects Drug_Name, CellLine, IC50_uM
models = models_for_list
ccle = load_ccle(FILES["CCLE"])   # expects DepMap_ID, Hugo_Symbol (long) or wide

joined = build_join(gdsc, models)

# Apply filters
if tumors:
    if "TCGA_Classification" in joined.columns:
        joined = joined[joined["TCGA_Classification"].astype(str).isin(tumors)]
if drug_query:
    if "Drug_Name" in joined.columns:
        joined = joined[joined["Drug_Name"].str.contains(drug_query, case=False, na=False)]

# Mut/WT assignment via CCLE mutations
mut_ids = mutated_depmap_ids(ccle, gene)
if "DepMap_ID" in joined.columns:
    joined["Mut"] = joined["DepMap_ID"].astype(str).isin(mut_ids)
else:
    joined["Mut"] = False

# Require IC50
if "IC50_uM" not in joined.columns:
    st.error("IC50 column not found. Ensure GDSC1 has 'IC50 (uM)' header.")
    st.stop()

# Drop NA and collapse duplicates by median per (Drug, DepMap_ID) or (Drug, CellLine)
work = joined.dropna(subset=["IC50_uM"]).copy()
if "DepMap_ID" in work.columns and work["DepMap_ID"].notna().any():
    grp_keys = ["Drug_Name", "DepMap_ID"] if "Drug_Name" in work.columns else ["DepMap_ID"]
else:
    grp_keys = ["Drug_Name", "CellLine"] if "Drug_Name" in work.columns else ["CellLine"]

before = len(work)

def _agg(g: pd.DataFrame) -> pd.Series:
    s = pd.Series({
        "IC50_uM": g["IC50_uM"].astype(float).median(),
        "Mut": bool(g.get("Mut", pd.Series([False])).any()),
        "CellLine": g.get("CellLine", pd.Series([None])).iloc[0],
        "TCGA_Classification": g.get("TCGA_Classification", pd.Series([None])).iloc[0],
        "DepMap_ID": g.get("DepMap_ID", pd.Series([None])).iloc[0],
        "Drug_Name": g.get("Drug_Name", pd.Series([None])).iloc[0],
        "_dup_count": len(g),
    })
    return s

plot_df = work.groupby(grp_keys, as_index=False).apply(_agg).reset_index(drop=True)
after = len(plot_df)
collapsed = max(0, before - after)

# Sort & trim
plot_df = plot_df.sort_values(by=["IC50_uM"], ascending=True).head(int(top_n))

# Diagnostics
st.subheader("Diagnostics")
try:
    total = int(len(plot_df))
    depmap_matched = int(plot_df["DepMap_ID"].notna().sum()) if "DepMap_ID" in plot_df.columns else 0
    mut_count = int(plot_df.get("Mut", pd.Series(dtype=bool)).sum()) if "Mut" in plot_df.columns else 0
    st.markdown(f"**Rows:** {total} | **DepMap_ID matched:** {depmap_matched} | **Mutants:** {mut_count} | **Collapsed duplicates:** {collapsed}")
except Exception:
    pass

# Filtered table (median-collapsed)
st.subheader("Filtered table (Median‑collapsed)")
show_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","IC50_uM","Mut"] if c in plot_df.columns]
st.dataframe(plot_df[show_cols], use_container_width=True, hide_index=True)

# Raw replicates (optional)
with st.expander("Show raw replicates (before collapsing)"):
    raw_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","IC50_uM","Mut"] if c in joined.columns]
    st.dataframe(joined[raw_cols], use_container_width=True, hide_index=True)

# Plot
st.subheader("IC50 distribution (lower is more sensitive)")
if plot_df.empty:
    st.warning("No rows to plot. Try broadening filters or check dataset mappings.")
else:
    fig, ax = plt.subplots(figsize=(12, 5))
    x_labels = plot_df.get("CellLine", plot_df.get("DepMap_ID", pd.Series(range(len(plot_df)))))
    y = plot_df["IC50_uM"].astype(float)
    colors = ["red" if bool(m) else "blue" for m in plot_df.get("Mut", pd.Series(False, index=plot_df.index))]
    ax.bar(x_labels.astype(str), y, color=colors)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="red", label="Mutant"), Patch(color="blue", label="WT")])
    ax.set_ylabel("IC50 (µM)")
    ax.set_xlabel("Cell line")
    title_drug = drug_query if drug_query else "(all drugs)"
    title_tcga = ", ".join(tumors) if tumors else "All TCGA"
    ax.set_title(f"{title_drug} — {title_tcga} — Gene: {gene}")
    ax.tick_params(axis='x', labelrotation=90)
    st.pyplot(fig, clear_figure=True)
