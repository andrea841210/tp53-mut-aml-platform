# Node-1 IC50 Explorer • Solid Tumors (KRAS track)
# Streamlit app skeleton v0.3 — added Run button in sidebar

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
    depmap_col = find_col(df, ["DepMap_ID", "depmap_id", "ModelID"])
    cell_candidates = [
        "stripped_cell_line_name","Model Name","Cell Line Name","cell_line_name","model_name","model_name_s"
    ]
    cell_col = find_col(df, cell_candidates, required=False)
    tcga_col = find_col(df,["TCGA_Classification","TCGA Classification","tcga_classification","primary_disease","lineage","OncotreePrimaryDisease"],required=False)
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
    if "DepMap_ID" not in df.columns and "__norm_cell" in df.columns and "__norm_cell" in models.columns:
        df = df.merge(models[["DepMap_ID","__norm_cell"]], on="__norm_cell", how="left")
    if "TCGA_Classification" not in df.columns and "__norm_cell" in df.columns and "__norm_cell" in models.columns:
        df = df.merge(models[["__norm_cell","TCGA_Classification"]], on="__norm_cell", how="left", suffixes=("","_m"))
        if "TCGA_Classification" not in df.columns and "TCGA_Classification_m" in df.columns:
            df = df.rename(columns={"TCGA_Classification_m":"TCGA_Classification"})
    if "TCGA_Classification" in df.columns:
        df["TumorGroup"] = df["TCGA_Classification"].map(infer_onco_group)
    else:
        df["TumorGroup"] = None
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
    top_n = st.number_input("Top N bars to show", min_value=10, max_value=200, value=50, step=10)
    run_button = st.button("Run")

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
plot_df = plot_df.sort_values(by=["Mut", ic_col], ascending=[False, True])
plot_df = plot_df.head(int(top_n))

st.subheader("Filtered table")
show_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","TumorGroup",ic_col,"Mut"] if c in plot_df.columns]
st.dataframe(plot_df[show_cols], use_container_width=True, hide_index=True)

st.subheader("IC50 distribution (lower is more sensitive)")
if plot_df.empty:
    st.warning("No rows to plot. Try broadening filters or check dataset mappings.")
else:
    fig, ax = plt.subplots(figsize=(12, 5))
    x = plot_df.get("CellLine", plot_df.get("DepMap_ID", pd.Series(range(len(plot_df)))))
    y = plot_df[ic_col]
    ax.bar(x.astype(str), y)
    if plot_df["Mut"].any():
        ax.bar(x[plot_df["Mut"]].astype(str), y[plot_df["Mut"]], label="Mut")
    ax.set_ylabel("IC50 (µM)")
    ax.set_xlabel("Cell line")
    title_drug = drug_query if drug_query else "(all drugs)"
    ax.set_title(f"{title_drug} — {', '.join(tumors) if tumors else 'All tumor groups'} — Gene: {gene}")
    ax.tick_params(axis='x', labelrotation=90)
    st.pyplot(fig, clear_figure=True)

st.caption("Notes: Mut=any variant recorded in CCLE 22Q2 for the selected gene; TumorGroup is inferred. If DepMap_ID missing in GDSC1, normalized cell-line name join is used.")
