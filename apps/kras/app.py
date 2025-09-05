# Node-1 IC50 Explorer • Solid Tumors (KRAS track)
# Streamlit app skeleton v0.9 — Median collapse (fixed) + duplicate diagnostics + optional raw replicates view

from __future__ import annotations
import os
import re
from pathlib import Path
from typing import Iterable, Optional, Set

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

st.set_page_config(page_title="Node-1 IC50 Explorer — KRAS track", layout="wide")
st.title("Node-1 IC50 Explorer • Solid Tumors (KRAS track)")

DATA_DIR = Path(os.environ.get("DATA_DIR", "data"))
FILES = {
    "CCLE": (DATA_DIR / "CCLE_mutation_22Q2.csv", "22Q2"),
    "MODEL": (DATA_DIR / "Model_DepMap.csv", "25Q2"),
    "GDSC1": (DATA_DIR / "GDSC1_IC50.csv", "as-downloaded"),
}

def _fmt_exists(p: Path) -> str:
    return f"✅ {p}" if p.exists() else f"❌ {p}"

st.caption(
    "Datasets → "
    f"CCLE {_fmt_exists(FILES['CCLE'][0])} • "
    f"Model {_fmt_exists(FILES['MODEL'][0])} • "
    f"GDSC1 {_fmt_exists(FILES['GDSC1'][0])} | "
    f"versions: CCLE {FILES['CCLE'][1]}, Model {FILES['MODEL'][1]}, GDSC1 {FILES['GDSC1'][1]}"
)

DEFAULT_GENES = ["KRAS", "BRAF", "PIK3CA", "TP53", "STK11", "KEAP1", "SMAD4", "CDKN2A"]
DEFAULT_TUMORS = ["CRC", "NSCLC", "PDAC", "BRCA", "SKCM", "GBM"]

def find_col(df: pd.DataFrame, candidates: Iterable[str], *, required: bool = True) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        st.error(f"Required column not found. Tried: {list(candidates)}")
        st.stop()
    return None

_norm_cache = {}
def normalize_cell_name(s: str) -> str:
    if s in _norm_cache:
        return _norm_cache[s]
    n = re.sub(r"[^a-z0-9]", "", str(s).lower())
    _norm_cache[s] = n
    return n

@st.cache_data(show_spinner=True)
def load_ccle_mutations(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    depmap_col = find_col(df, ["DepMap_ID", "depmap_id", "depMapID"])
    gene_col = find_col(df, ["Hugo_Symbol", "Gene", "gene_symbol"])
    df = df.rename(columns={depmap_col: "DepMap_ID", gene_col: "Hugo_Symbol"})
    keep = ["DepMap_ID", "Hugo_Symbol"]
    keep = [c for c in keep if c in df.columns]
    return df[keep].dropna()

@st.cache_data(show_spinner=True)
def load_models(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns=lambda x: x.strip())
    # Prefer explicit DepMap_ID if present; else fall back to ModelID (ACH-xxxxx)
    depmap_col = None
    for c in ["DepMap_ID", "depmap_id", "ModelID", "modelid"]:
        if c in df.columns:
            depmap_col = c
            break
    if depmap_col is None:
        st.error("Model_DepMap.csv missing DepMap_ID/ModelID column.")
        st.stop()
    # Cell line naming candidates
    cell_candidates = [
        "CellLineName","Cell Line Name","stripped_cell_line_name","StrippedCellLineName",
        "Model Name","cell_line_name","model_name","model_name_s","CCLEName"
    ]
    cell_col = None
    for c in cell_candidates:
        if c in df.columns:
            cell_col = c
            break
    # TCGA/lineage candidates
    tcga_candidates = [
        "TCGA_Classification","TCGA Classification","tcga_classification","primary_disease","lineage","OncotreePrimaryDisease"
    ]
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
    # Build multiple normalization keys to improve match rate
    src_name_cols = [c for c in ["CellLine","CCLEName","StrippedCellLineName","stripped_cell_line_name"] if c in df.columns]
    if not src_name_cols:
        st.warning("Model_DepMap.csv has no cell-line name columns; only DepMap_ID joins will work.")
    # Primary normalized key from the first available name column
    if src_name_cols:
        df["__norm_cell"] = df[src_name_cols[0]].map(normalize_cell_name)
    # Secondary keys for broader matching
    for c in src_name_cols[1:]:
        key = f"__norm_{normalize_cell_name(c)}"
        df[key] = df[c].map(normalize_cell_name)
    return df

@st.cache_data(show_spinner=True)
def load_gdsc1(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns=lambda x: x.strip())
    rename_map = {}
    if "DepMap ID" in df.columns and "DepMap_ID" not in df.columns:
        rename_map["DepMap ID"] = "DepMap_ID"
    if "IC50 (uM)" in df.columns and "IC50_uM" not in df.columns:
        rename_map["IC50 (uM)"] = "IC50_uM"
    if "Drug Name" in df.columns and "Drug_Name" not in df.columns:
        rename_map["Drug Name"] = "Drug_Name"
    if "Cell Line Name" in df.columns and "CellLine" not in df.columns:
        rename_map["Cell Line Name"] = "CellLine"
    if "TCGA Classification" in df.columns and "TCGA_Classification" not in df.columns:
        rename_map["TCGA Classification"] = "TCGA_Classification"
    df = df.rename(columns=rename_map)
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    return df

def infer_onco_group(tcga_str: str) -> Optional[str]:
    if not isinstance(tcga_str, str):
        return None
    s = tcga_str.upper()
    if any(k in s for k in ["COAD","READ","COADREAD","COLORECT","COLON","RECTUM"]):
        return "CRC"
    if "LUNG" in s and "SMALL" not in s:
        return "NSCLC"
    if any(k in s for k in ["PAAD","PANCREAS","PANCREATIC"]):
        return "PDAC"
    if any(k in s for k in ["BRCA","BREAST"]):
        return "BRCA"
    if any(k in s for k in ["SKCM","MELANOMA"]):
        return "SKCM"
    if any(k in s for k in ["GBM","GLIOBLASTOMA","BRAIN"]):
        return "GBM"
    return None

@st.cache_data(show_spinner=True)
def build_gdsc_joined(gdsc: pd.DataFrame, models: pd.DataFrame) -> pd.DataFrame:
    df = gdsc.copy()
    # Try a cascade of merges to obtain DepMap_ID and TCGA from models
    name_keys = [c for c in models.columns if c.startswith("__norm_")]
    keys_to_try = ["__norm_cell"] + name_keys if "__norm_cell" in df.columns else name_keys
    gained_depmap = 0
    for k in keys_to_try:
        if k in df.columns and k in models.columns:
            before = df["DepMap_ID"].notna().sum() if "DepMap_ID" in df.columns else 0
            df = df.merge(models[["DepMap_ID", k]], on=k, how="left", suffixes=("","_m2"))
            # Fill only where missing
            if "DepMap_ID_m2" in df.columns:
                df["DepMap_ID"] = df["DepMap_ID"].astype(object)
                df.loc[df["DepMap_ID"].isna(), "DepMap_ID"] = df.loc[df["DepMap_ID"].isna(), "DepMap_ID_m2"]
                df = df.drop(columns=["DepMap_ID_m2"])
            after = df["DepMap_ID"].notna().sum()
            gained_depmap += max(0, after - before)
    if "TCGA_Classification" not in df.columns and "__norm_cell" in df.columns and "__norm_cell" in models.columns:
        df = df.merge(models[["__norm_cell","TCGA_Classification"]], on="__norm_cell", how="left", suffixes=("","_m"))
        if "TCGA_Classification" not in df.columns and "TCGA_Classification_m" in df.columns:
            df = df.rename(columns={"TCGA_Classification_m":"TCGA_Classification"})
    df["TumorGroup"] = df.get("TCGA_Classification").map(infer_onco_group) if "TCGA_Classification" in df.columns else None
    return df

@st.cache_data(show_spinner=True)
def mutated_depmap_ids(ccle: pd.DataFrame, gene: str) -> Set[str]:
    sub = ccle[ccle["Hugo_Symbol"].str.upper() == gene.upper()]
    return set(sub["DepMap_ID"].dropna().astype(str).unique())

with st.sidebar:
    st.header("Filters")
    gene = st.selectbox("Gene", DEFAULT_GENES, index=0)
    tumors = st.multiselect("Cancer types", DEFAULT_TUMORS, default=[DEFAULT_TUMORS[0]])
    drug_query = st.text_input("Drug name contains (optional)", value="")
    top_n = st.number_input("Top N bars to show", min_value=10, max_value=200, value=200, step=10)
    run_button = st.button("Run")
    st.caption("Duplicates per (Drug, Cell line) are collapsed by **median IC50** for stability.")

if not run_button:
    st.stop()

try:
    ccle = load_ccle_mutations(FILES["CCLE"][0])
    models = load_models(FILES["MODEL"][0])
    gdsc = load_gdsc1(FILES["GDSC1"][0])
except Exception as e:
    st.error(f"Failed to load datasets: {e}")
    st.stop()

joined = build_gdsc_joined(gdsc, models)
if tumors:
    joined = joined[joined["TumorGroup"].isin(tumors)]
if drug_query:
    if "Drug_Name" in joined.columns:
        joined = joined[joined["Drug_Name"].str.contains(drug_query, case=False, na=False)]
    else:
        st.info("'Drug_Name' column not found in GDSC1; skipping drug filter.")

mut_ids = mutated_depmap_ids(ccle, gene)
if "DepMap_ID" in joined.columns:
    joined["Mut"] = joined["DepMap_ID"].astype(str).isin(mut_ids)
else:
    joined["Mut"] = False

ic_col = "IC50_uM" if "IC50_uM" in joined.columns else None
if ic_col is None:
    st.error("IC50 column not found (expected 'IC50 (uM)' → 'IC50_uM'). Please check GDSC1 file.")
    st.stop()

plot_df = joined.dropna(subset=[ic_col]).copy()

# --- Collapse duplicates per (Drug, DepMap_ID or CellLine) using MEDIAN ---
# Prefer grouping by DepMap_ID when available; else fall back to CellLine
if "DepMap_ID" in plot_df.columns:
    grp_keys = ["Drug_Name", "DepMap_ID"] if "Drug_Name" in plot_df.columns else ["DepMap_ID"]
else:
    grp_keys = ["Drug_Name", "CellLine"] if "Drug_Name" in plot_df.columns else ["CellLine"]

dup_before = len(plot_df)

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

dup_after = len(plot_df)
collapsed = dup_before - dup_after
# Global sort by IC50 ascending (no Mut-first grouping) to mirror Google Sheet behavior
plot_df = plot_df.sort_values(by=[ic_col], ascending=True)
plot_df = plot_df.head(int(top_n))

# Diagnostics summary
st.subheader("Diagnostics")
try:
    total = int(len(plot_df))
    depmap_matched = int(plot_df["DepMap_ID"].notna().sum()) if "DepMap_ID" in plot_df.columns else 0
    mut_count = int(plot_df.get("Mut", pd.Series(dtype=bool)).sum()) if "Mut" in plot_df.columns else 0
    st.markdown(f"**Rows:** {total} | **DepMap_ID matched:** {depmap_matched} | **Mutants:** {mut_count} | **Collapsed duplicates:** {collapsed}")
except Exception:
    pass

st.subheader("Filtered table (Median-collapsed)")
show_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","TumorGroup",ic_col,"Mut"] if c in plot_df.columns]
st.dataframe(plot_df[show_cols], use_container_width=True, hide_index=True)

# --- Raw replicates table (optional) ---
with st.expander("Show raw replicates (before collapsing)"):
    raw_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","TumorGroup",ic_col,"Mut"] if c in joined.columns]
    st.dataframe(joined[raw_cols], use_container_width=True, hide_index=True)

st.subheader("IC50 distribution (lower is more sensitive)")
if plot_df.empty:
    st.warning("No rows to plot. Try broadening filters or check dataset mappings.")
else:
    fig, ax = plt.subplots(figsize=(12, 5))
    x_labels = plot_df.get("CellLine", plot_df.get("DepMap_ID", pd.Series(range(len(plot_df)))))
    y = plot_df[ic_col].astype(float)
    mut_series = plot_df.get("Mut", pd.Series(False, index=plot_df.index))
    colors = ["red" if bool(m) else "blue" for m in mut_series]
    ax.bar(x_labels.astype(str), y, color=colors)
    # Legend
    from matplotlib.patches import Patch
    legend_handles = [Patch(color="red", label="Mutant"), Patch(color="blue", label="WT")]
    ax.legend(handles=legend_handles, loc="best")
    ax.set_ylabel("IC50 (µM)")
    ax.set_xlabel("Cell line")
    title_drug = drug_query if drug_query else "(all drugs)"
    ax.set_title(f"{title_drug} — {', '.join(tumors) if tumors else 'All tumor groups'} — Gene: {gene}")
    ax.tick_params(axis='x', labelrotation=90)
    st.pyplot(fig, clear_figure=True)

st.caption("Notes: Mut=any variant recorded in CCLE 22Q2 for the selected gene; TumorGroup is inferred. If DepMap_ID missing in GDSC1, normalized cell-line name join is used.")
