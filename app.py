"""
Streamlit Lite module for gene-subset filtering (TP53 / FLT3 / NPM1)
- Reads long table: data/CCLE_TP53_FLT3_NPM1_long.csv  (DepMap_ID, Gene, Status)
- Joins to your IC50/AUC master table (GDSC1_LAML.xlsx) with DepMap mapping (Model_DepMap.csv)

Usage:
1) Put this file's content into your app.py (or merge the parts you need).
2) Ensure the three data files exist in repo/data: CCLE_TP53_FLT3_NPM1_long.csv, GDSC1_LAML.xlsx, Model_DepMap.csv
3) Run: streamlit run app.py
"""

from pathlib import Path
import pandas as pd
import streamlit as st

# -------------------
# Config & file paths
# -------------------
DATA_DIR = Path("data")
LONG_PATH = DATA_DIR / "CCLE_TP53_FLT3_NPM1_long.csv"
IC50_PATH = DATA_DIR / "GDSC1_LAML.xlsx"     # IC50/AUC master table
MODEL_PATH = DATA_DIR / "Model_DepMap.csv"   # Mapping table (ModelID/CCLEName/COSMICID → DepMap_ID)

# -------------------
# Loaders (cached)
# -------------------
@st.cache_data
def load_long() -> pd.DataFrame:
    df = pd.read_csv(LONG_PATH)
    df.columns = [c.strip() for c in df.columns]
    df["DepMap_ID"] = df["DepMap_ID"].astype(str)
    df["Gene"] = df["Gene"].astype(str).str.upper()
    df["Status"] = df["Status"].astype(str).str.upper()
    return df

@st.cache_data
def load_ic50_with_depmap() -> pd.DataFrame:
    """Load GDSC1_LAML.xlsx and attach DepMap_ID via Model_DepMap.csv (by COSMIC ID)."""
    # Load IC50/AUC table
    if IC50_PATH.suffix.lower() in [".xls", ".xlsx"]:
        ic50 = pd.read_excel(IC50_PATH)
    else:
        ic50 = pd.read_csv(IC50_PATH)

    # Load mapping table
    model = pd.read_csv(MODEL_PATH)

    # Normalize COSMIC IDs
    ic50["Cosmic ID"] = pd.to_numeric(ic50["Cosmic ID"], errors="coerce")
    model["COSMICID"] = pd.to_numeric(model["COSMICID"], errors="coerce")

    # Merge on COSMIC → bring in ModelID (DepMap)
    merged = ic50.merge(
        model[["ModelID", "COSMICID", "CCLEName", "CellLineName", "StrippedCellLineName"]],
        left_on="Cosmic ID",
        right_on="COSMICID",
        how="left",
    )

    # Standardize: ModelID → DepMap_ID
    merged = merged.rename(columns={"ModelID": "DepMap_ID"})
    merged["DepMap_ID"] = merged["DepMap_ID"].astype(str)

    # Optional: clean IC50 column name
    if "IC50_uM" not in merged.columns and "IC50" in merged.columns:
        merged = merged.rename(columns={"IC50": "IC50_uM"})

    return merged

# -------------------
# Subset logic
# -------------------
GENES = ["TP53", "FLT3", "NPM1"]

def get_ids(df_long: pd.DataFrame, gene: str) -> set:
    return set(df_long[(df_long["Gene"] == gene) & (df_long["Status"] == "MUT")]["DepMap_ID"])

@st.cache_data
def compute_gene_sets(df_long: pd.DataFrame):
    return {g: get_ids(df_long, g) for g in GENES}

def subset_ids_from_selection(gsets: dict, selection: str) -> set:
    tp53, flt3, npm1 = gsets.get("TP53", set()), gsets.get("FLT3", set()), gsets.get("NPM1", set())
    mapping = {
        "TP53 only": tp53,
        "TP53 + FLT3": tp53 & flt3,
        "TP53 + FLT3 + NPM1": tp53 & flt3 & npm1,
        # Advanced (hidden by default)
        "FLT3 only": flt3,
        "NPM1 only": npm1,
        "TP53 + NPM1": tp53 & npm1,
        "FLT3 + NPM1": flt3 & npm1,
    }
    return mapping.get(selection, set())

# -------------------
# UI
# -------------------
st.subheader("Gene Subset (Lite Mode)")

try:
    df_long = load_long()
except FileNotFoundError:
    st.error(f"Long table not found: {LONG_PATH}.")
    st.stop()

try:
    df_ic50 = load_ic50_with_depmap()
except Exception as e:
    st.error(str(e))
    st.stop()

# Options
lite_options = ["TP53 only", "TP53 + FLT3", "TP53 + FLT3 + NPM1"]
advanced = st.toggle("Advanced subsets", value=False, help="Enable all 7 combinations for demo")
all_options = [
    "TP53 only", "FLT3 only", "NPM1 only",
    "TP53 + FLT3", "TP53 + NPM1", "FLT3 + NPM1",
    "TP53 + FLT3 + NPM1",
]
options = lite_options if not advanced else all_options

selection = st.selectbox("Choose subset", options, index=0)

# Compute sets & pick IDs
gene_sets = compute_gene_sets(df_long)
ids = subset_ids_from_selection(gene_sets, selection)

st.caption(f"Matched cell lines: {len(ids)}")

if selection.endswith("NPM1") and len(ids) == 0:
    st.info("NPM1 mutations are rare in DepMap/CCLE cell lines (common in primary AML cohorts). 0 is expected.")

if not ids:
    st.warning("No samples in this subset. Try another selection or disable Advanced.")
    st.stop()

# Join and show
view = df_ic50[df_ic50["DepMap_ID"].isin(ids)].copy()

rename_map = {
    "Drug Name": "Drug",
    "drug": "Drug",
    "IC50 (uM)": "IC50_uM",
    "IC50": "IC50_uM",
    "Cell Line Name": "CellLine",
    "Cell line name": "CellLine",
}
for k, v in list(rename_map.items()):
    if k in view.columns and v not in view.columns:
        view = view.rename(columns={k: v})

st.dataframe(view, use_container_width=True)

# (Optional) plot hook
# if "IC50_uM" in view.columns:
#     st.bar_chart(view.sort_values("IC50_uM")["IC50_uM"])

