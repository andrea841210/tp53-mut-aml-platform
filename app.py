"""
App.py — Node-1 (Lite) with Gene Subset + Drug & Cancer-Type Filters

Data files expected in repo/data:
- CCLE_TP53_FLT3_NPM1_long.csv  (DepMap_ID, Gene, Status)
- GDSC1_LAML.xlsx               (IC50/AUC master; has 'Cosmic ID')
- Model_DepMap.csv              (mapping: COSMICID → ModelID[=DepMap_ID])

What this app shows (left sidebar → main panel):
1) Gene Subset filters (TP53 only / TP53+FLT3 / TP53+FLT3+NPM1; Advanced toggle)
2) Drug selector + Cancer-type filters (TCGA, Tissue, Tissue Sub-type)
3) IC50 table (filtered) and bar chart per cell line for the chosen drug

Run:
  streamlit run app.py
"""

from pathlib import Path
import pandas as pd
import streamlit as st

# -------------------
# Config & paths
# -------------------
DATA_DIR = Path("data")
LONG_PATH = DATA_DIR / "CCLE_TP53_FLT3_NPM1_long.csv"
IC50_PATH = DATA_DIR / "GDSC1_LAML.xlsx"
MODEL_PATH = DATA_DIR / "Model_DepMap.csv"

# -------------------
# Loaders
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
    # IC50/AUC master
    if IC50_PATH.suffix.lower() in [".xls", ".xlsx"]:
        ic50 = pd.read_excel(IC50_PATH)
    else:
        ic50 = pd.read_csv(IC50_PATH)

    # Mapping: COSMIC → ModelID (=DepMap)
    model = pd.read_csv(MODEL_PATH)
    ic50["Cosmic ID"] = pd.to_numeric(ic50["Cosmic ID"], errors="coerce")
    model["COSMICID"] = pd.to_numeric(model["COSMICID"], errors="coerce")

    merged = ic50.merge(
        model[["ModelID", "COSMICID", "CCLEName", "CellLineName", "StrippedCellLineName"]],
        left_on="Cosmic ID",
        right_on="COSMICID",
        how="left",
    )
    merged = merged.rename(columns={"ModelID": "DepMap_ID"})
    merged["DepMap_ID"] = merged["DepMap_ID"].astype(str)

    # Make a friendly CellLine label if available
    if "CellLine" not in merged.columns:
        # Prefer CCLEName → CellLineName → StrippedCellLineName → original 'Cell Line Name'
        for c in ["CCLEName", "CellLineName", "StrippedCellLineName", "Cell Line Name"]:
            if c in merged.columns:
                merged = merged.rename(columns={c: "CellLine"})
                break

    # Normalize IC50 name
    if "IC50_uM" not in merged.columns and "IC50" in merged.columns:
        merged = merged.rename(columns={"IC50": "IC50_uM"})

    # Normalize Drug name
    if "Drug" not in merged.columns and "Drug Name" in merged.columns:
        merged = merged.rename(columns={"Drug Name": "Drug"})

    return merged

# -------------------
# Gene subset helpers
# -------------------
GENES = ["TP53", "FLT3", "NPM1"]

def ids_for(df_long: pd.DataFrame, gene: str) -> set:
    return set(df_long[(df_long["Gene"] == gene) & (df_long["Status"] == "MUT")]["DepMap_ID"])  # type: ignore

@st.cache_data
def compute_gene_sets(df_long: pd.DataFrame):
    return {g: ids_for(df_long, g) for g in GENES}

def subset_ids_from(selection: str, gsets: dict) -> set:
    tp53, flt3, npm1 = gsets.get("TP53", set()), gsets.get("FLT3", set()), gsets.get("NPM1", set())
    mapping = {
        "TP53 only": tp53,
        "TP53 + FLT3": tp53 & flt3,
        "TP53 + FLT3 + NPM1": tp53 & flt3 & npm1,
        # Advanced
        "FLT3 only": flt3,
        "NPM1 only": npm1,
        "TP53 + NPM1": tp53 & npm1,
        "FLT3 + NPM1": flt3 & npm1,
    }
    return mapping.get(selection, set())

# -------------------
# App UI
# -------------------
st.set_page_config(page_title="Node-1 LAML • IC50 Explorer", layout="wide")

# Load data
try:
    df_long = load_long()
    df_ic50 = load_ic50_with_depmap()
except Exception as e:
    st.error(f"Data loading error: {e}")
    st.stop()

# Sidebar — Gene subset first
with st.sidebar:
    st.header("Filters")
    st.subheader("Gene Subset (Lite)")
    lite_options = ["TP53 only", "TP53 + FLT3", "TP53 + FLT3 + NPM1"]
    advanced = st.toggle("Advanced subsets", value=False)
    all_options = [
        "TP53 only", "FLT3 only", "NPM1 only",
        "TP53 + FLT3", "TP53 + NPM1", "FLT3 + NPM1",
        "TP53 + FLT3 + NPM1",
    ]
    options = lite_options if not advanced else all_options
    subset_choice = st.selectbox("Choose subset", options, index=0)

# Compute subset ids and filter IC50 table early
gene_sets = compute_gene_sets(df_long)
subset_ids = subset_ids_from(subset_choice, gene_sets)

if subset_choice.endswith("NPM1") and len(subset_ids) == 0:
    st.sidebar.info("NPM1 mutations are rare in DepMap/CCLE cell lines (0 is expected).")

if not subset_ids:
    st.warning("No samples in this subset. Try another selection or disable Advanced.")
    st.stop()

base = df_ic50[df_ic50["DepMap_ID"].isin(subset_ids)].copy()

# Sidebar — Drug & cancer-type filters
with st.sidebar:
    st.subheader("Drug & Cancer Type")
    drugs = sorted(base["Drug"].dropna().unique().tolist())
    drug = st.selectbox("Drug", drugs, index=0 if drugs else None)

    # Cancer-type filters (optional)
    tcga_vals = sorted(base.get("TCGA Classification", pd.Series(dtype=str)).dropna().unique().tolist())
    tissue_vals = sorted(base.get("Tissue", pd.Series(dtype=str)).dropna().unique().tolist())
    subtype_vals = sorted(base.get("Tissue Sub-type", pd.Series(dtype=str)).dropna().unique().tolist())

    sel_tcga = st.multiselect("TCGA Classification", tcga_vals, default=tcga_vals)
    sel_tissue = st.multiselect("Tissue", tissue_vals, default=tissue_vals)
    sel_subtype = st.multiselect("Tissue Sub-type", subtype_vals, default=subtype_vals)

# Apply filters
df = base.copy()
df = df[df["Drug"] == drug]
if "TCGA Classification" in df.columns:
    df = df[df["TCGA Classification"].isin(sel_tcga)]
if "Tissue" in df.columns:
    df = df[df["Tissue"].isin(sel_tissue)]
if "Tissue Sub-type" in df.columns:
    df = df[df["Tissue Sub-type"].isin(sel_subtype)]

# Main — Header & metrics
st.title("Gene Subset (Lite Mode)")
st.caption(f"Matched cell lines after subset: {len(subset_ids)} | Rows after all filters: {len(df)}")

# Show table
st.dataframe(
    df[[c for c in ["Drug", "Drug ID", "CellLine", "Cosmic ID", "TCGA Classification", "Tissue", "Tissue Sub-type", "IC50_uM", "AUC", "DepMap_ID"] if c in df.columns]],
    use_container_width=True,
)

# Plot — IC50 bar chart per cell line
if "IC50_uM" in df.columns and not df.empty:
    chart_df = df[["CellLine", "IC50_uM"]].dropna().copy()
    chart_df = chart_df.groupby("CellLine", as_index=False)["IC50_uM"].median()
    chart_df = chart_df.sort_values("IC50_uM", ascending=True)
    st.subheader("IC50 (µM) by Cell Line — median per cell line")
    st.bar_chart(chart_df.set_index("CellLine"))
else:
    st.info("No IC50 data available for the current selection.")

