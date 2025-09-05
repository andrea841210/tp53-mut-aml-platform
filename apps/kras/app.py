# Node-1 IC50 Explorer • Solid Tumors (KRAS track)
# Streamlit app v1.4 — Reorder sections, hide Debug by default

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
VERSION = "v1.4"
st.set_page_config(page_title=f"Node-1 IC50 Explorer — KRAS track {VERSION}", layout="wide")
st.title(f"Node-1 IC50 Explorer • Solid Tumors (KRAS track) — {VERSION}")

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
    f"GDSC1 {_fmt_exists(FILES['GDSC1'])} | IC50 aggregation: **median** | "
    f"Build: {VERSION}"
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
    # Keep potential name/ontology columns
    cell_candidates = [
        "CellLineName", "Cell Line Name", "stripped_cell_line_name", "StrippedCellLineName",
        "CCLEName", "Model Name"
    ]
    cell_col = find_col(df, cell_candidates, required=False)
    tcga_col = find_col(df, ["TCGA_Classification", "TCGA Classification", "tcga_classification"], required=False)
    # Prepare base rename
    rename_map = {depmap_col: "DepMap_ID"}
    if cell_col:
        rename_map[cell_col] = "CellLine"
    if tcga_col:
        rename_map[tcga_col] = "TCGA_Classification"
    df = df.rename(columns=rename_map)
    # Build normalization keys
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    for alt in [c for c in ["CCLEName", "StrippedCellLineName", "stripped_cell_line_name"] if c in df.columns]:
        df[f"__norm_{alt}"] = df[alt].map(normalize_cell_name)
    # --- NEW: Oncotree → TCGA mapping fallback ---
    if "TCGA_Classification" not in df.columns and "OncotreePrimaryDisease" in df.columns:
        def oncotree_to_tcga(x: str) -> Optional[str]:
            if not isinstance(x, str):
                return None
            s = x.strip().upper()
            mapping = {
                "LUNG ADENOCARCINOMA": "LUAD",
                "LUNG SQUAMOUS CELL CARCINOMA": "LUSC",
                "NON-SMALL CELL LUNG CANCER": "NSCLC",
                "SMALL CELL LUNG CANCER": "SCLC",
                "COLORECTAL ADENOCARCINOMA": "COREAD",
                "BREAST CARCINOMA": "BRCA",
                "SKIN CUTANEOUS MELANOMA": "SKCM",
                "GLIOBLASTOMA": "GBM",
                "PANCREATIC ADENOCARCINOMA": "PAAD",
                "ESOPHAGEAL CARCINOMA": "ESCA",
                "HEAD AND NECK SQUAMOUS CELL CARCINOMA": "HNSC",
                "OVARIAN SEROUS CYSTADENOCARCINOMA": "OV",
            }
            if s in mapping:
                return mapping[s]
            if ("ADENOCARCINOMA" in s) and ("COLON" in s or "RECT" in s):
                return "COREAD"
            if "MELANOMA" in s:
                return "SKCM"
            if "BREAST" in s:
                return "BRCA"
            if "GLIOBLAST" in s:
                return "GBM"
            if "PANCREA" in s:
                return "PAAD"
            if "LUNG" in s and ("SMALL" not in s and "SCLC" not in s):
                return "NSCLC"
            if "LUNG" in s and ("SMALL" in s or "SCLC" in s):
                return "SCLC"
            return None
        df["TCGA_Classification"] = df["OncotreePrimaryDisease"].map(oncotree_to_tcga)
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
    return df

@st.cache_data(show_spinner=True)
def build_join(gdsc: pd.DataFrame, models: pd.DataFrame) -> pd.DataFrame:
    df = gdsc.copy()
    if "CellLine" in df.columns:
        df["__norm_cell"] = df["CellLine"].map(normalize_cell_name)
    if "__norm_cell" in df.columns and "__norm_cell" in models.columns:
        df = df.merge(models[["DepMap_ID","__norm_cell","TCGA_Classification"]], on="__norm_cell", how="left")
    return df

@st.cache_data(show_spinner=True)
def mutated_depmap_ids(ccle: pd.DataFrame, gene: str) -> Set[str]:
    if {"DepMap_ID", "Hugo_Symbol"}.issubset(ccle.columns):
        sub = ccle[ccle["Hugo_Symbol"].str.upper() == gene.upper()]
        return set(sub["DepMap_ID"].astype(str).unique())
    if gene in ccle.columns and "DepMap_ID" in ccle.columns:
        sub = ccle[ccle[gene].astype(str).str.upper().isin(["1","TRUE","YES","Y"])]
        return set(sub["DepMap_ID"].astype(str).unique())
    return set()

# ---------------------------
# Sidebar
# ---------------------------
with st.sidebar:
    st.header("Filters")
    gene = st.selectbox("Gene", DEFAULT_GENES, index=0)
    models_for_list = load_models(FILES["MODEL"])
    gdsc_for_list = load_gdsc(FILES["GDSC1"])
    tcga_opts = set()
    if "TCGA_Classification" in models_for_list.columns:
        tcga_opts.update(models_for_list["TCGA_Classification"].dropna().astype(str).unique().tolist())
    if "TCGA_Classification" in gdsc_for_list.columns:
        tcga_opts.update(gdsc_for_list["TCGA_Classification"].dropna().astype(str).unique().tolist())
    tcga_classes = sorted(tcga_opts)
    tumors = st.multiselect("TCGA Classification (direct)", tcga_classes, default=[])
    drug_query = st.text_input("Drug name contains (optional)", value="")
    top_n = st.number_input("Top N bars to show", min_value=10, max_value=1000, value=200, step=10)
    run_button = st.button("Run")
    st.caption("Duplicates per (Drug, Cell line) are collapsed by **median IC50** for stability. **Note:** SCLC (small cell lung cancer) is excluded from LUAD/LUSC/NSCLC selections.")

# ---------------------------
# Main run
# ---------------------------
if not run_button:
    st.stop()

gdsc = load_gdsc(FILES["GDSC1"])
models = load_models(FILES["MODEL"])
ccle = load_ccle(FILES["CCLE"])

joined = build_join(gdsc, models)

if tumors:
    if "TCGA_Classification" in joined.columns:
        joined = joined[joined["TCGA_Classification"].astype(str).isin(tumors)]
    else:
        st.warning("TCGA_Classification not found in data; skipping TCGA filter for this dataset.")joined["TCGA_Classification"].astype(str).isin(tumors)]
if drug_query and "Drug_Name" in joined.columns:
    joined = joined[joined["Drug_Name"].str.contains(drug_query, case=False, na=False)]

mut_ids = mutated_depmap_ids(ccle, gene)
if "DepMap_ID" in joined.columns:
    joined["Mut"] = joined["DepMap_ID"].astype(str).isin(mut_ids)

if "IC50_uM" not in joined.columns:
    st.error("IC50 column not found. Ensure GDSC1 has 'IC50 (uM)' header.")
    st.stop()

work = joined.dropna(subset=["IC50_uM"]).copy()
if "DepMap_ID" in work.columns and work["DepMap_ID"].notna().any():
    grp_keys = ["Drug_Name", "DepMap_ID"] if "Drug_Name" in work.columns else ["DepMap_ID"]
else:
    grp_keys = ["Drug_Name", "CellLine"] if "Drug_Name" in work.columns else ["CellLine"]

before = len(work)

def _agg(g: pd.DataFrame) -> pd.Series:
    return pd.Series({
        "IC50_uM": g["IC50_uM"].astype(float).median(),
        "Mut": bool(g.get("Mut", pd.Series([False])).any()),
        "CellLine": g.get("CellLine", pd.Series([None])).iloc[0],
        "TCGA_Classification": g.get("TCGA_Classification", pd.Series([None])).iloc[0],
        "DepMap_ID": g.get("DepMap_ID", pd.Series([None])).iloc[0],
        "Drug_Name": g.get("Drug_Name", pd.Series([None])).iloc[0],
        "_dup_count": len(g),
    })

plot_df = work.groupby(grp_keys, as_index=False).apply(_agg).reset_index(drop=True)
after = len(plot_df)
collapsed = max(0, before - after)

sorted_df = plot_df.sort_values(by=["IC50_uM"], ascending=True)
available_rows = len(sorted_df)
plot_df = sorted_df.head(int(top_n))

# --- Ordered sections ---
st.subheader("Diagnostics")

try:
    total = int(len(plot_df))
    depmap_matched = int(plot_df["DepMap_ID"].notna().sum()) if "DepMap_ID" in plot_df.columns else 0
    mut_count = int(plot_df.get("Mut", pd.Series(dtype=bool)).sum()) if "Mut" in plot_df.columns else 0
    st.markdown(f"**Rows (after Top N):** {total} | **Available (before Top N):** {available_rows} | **DepMap_ID matched:** {depmap_matched} | **Mutants:** {mut_count} | **Collapsed duplicates:** {collapsed}")
except Exception:
    pass

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

st.subheader("Filtered table (Median-collapsed)")
show_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","IC50_uM","Mut"] if c in plot_df.columns]
st.dataframe(plot_df[show_cols], use_container_width=True, hide_index=True)

with st.expander("Show raw replicates (before collapsing)", expanded=False):
    raw_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification","IC50_uM","Mut"] if c in joined.columns]
    st.dataframe(joined[raw_cols], use_container_width=True, hide_index=True)

with st.expander("TCGA code reference (common groups)", expanded=False):
    st.markdown("""
    - **LUAD** = Lung adenocarcinoma (NSCLC)  
    - **LUSC** = Lung squamous cell carcinoma (NSCLC)  
    - **COREAD** = Colorectal adenocarcinoma (CRC)  
    - **BRCA** = Breast invasive carcinoma  
    - **SKCM** = Skin cutaneous melanoma  
    - **GBM** = Glioblastoma multiforme  
    - **PAAD** = Pancreatic adenocarcinoma (PDAC)
    """)

# --- Debug (kept folded at bottom) ---
with st.expander("Debug: TCGA dropdown sources", expanded=False):
    try:
        src_cols_models = [c for c in ["TCGA_Classification", "TCGA Classification", "primary_disease", "lineage", "OncotreePrimaryDisease"] if c in models_for_list.columns]
        src_cols_gdsc = [c for c in ["TCGA_Classification", "TCGA Classification"] if c in gdsc_for_list.columns]
        st.write({"models_has_cols": src_cols_models, "gdsc_has_cols": src_cols_gdsc, "num_options": len(tcga_classes)})
        if tcga_classes:
            st.write("First 20 TCGA options:")
            st.write(pd.Series(tcga_classes).head(20))
    except Exception:
        pass
