import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

# Title
st.title("Drug Sensitivity Explorer • LAML Subset")

# Load data
@st.cache_data

def load_data():
    gdsc = pd.read_excel("data/GDSC1_LAML.xlsx")
    depmap = pd.read_csv("data/Model_DepMap.csv")
    ccle_long = pd.read_csv("data/CCLE_TP53_FLT3_NPM1_long.csv")
    return gdsc, depmap, ccle_long

gdsc, depmap, ccle_long = load_data()

# Sidebar - subset selection
subset_option = st.sidebar.selectbox("Select Genetic Subset", [
    "TP53 only", "TP53 + FLT3", "TP53 + FLT3 + NPM1"
])

# Apply subset logic
if subset_option == "TP53 only":
    subset_df = ccle_long[(ccle_long["TP53"] == 1) & (ccle_long["FLT3"] == 0) & (ccle_long["NPM1"] == 0)]
elif subset_option == "TP53 + FLT3":
    subset_df = ccle_long[(ccle_long["TP53"] == 1) & (ccle_long["FLT3"] == 1) & (ccle_long["NPM1"] == 0)]
elif subset_option == "TP53 + FLT3 + NPM1":
    subset_df = ccle_long[(ccle_long["TP53"] == 1) & (ccle_long["FLT3"] == 1) & (ccle_long["NPM1"] == 1)]

# Sidebar inputs for drug and tissue
st.sidebar.write("### Filter options")
drug = st.sidebar.selectbox("Select Drug", sorted(gdsc["Drug"].dropna().unique().tolist()))
cancer = st.sidebar.selectbox("Select Cancer Type", sorted(gdsc["TCGA Classification"].dropna().unique().tolist()))

# Filter base
base = gdsc[(gdsc["Drug"] == drug) & (gdsc["TCGA Classification"] == cancer)]

# Join with subset
base = base.merge(subset_df, how="inner", left_on="DepMap_ID", right_on="DepMap_ID")

# Display data if available
if not base.empty:
    st.write("### Filtered Drug Sensitivity Table")
    st.dataframe(base[["Drug", "Cell Line Name", "IC50 (uM)", "AUC"]])

    st.write("### IC50 Distribution")
    fig, ax = plt.subplots()
    base[["Cell Line Name", "IC50 (uM)"]].set_index("Cell Line Name").plot(kind='bar', ax=ax)
    st.pyplot(fig)
else:
    st.warning("No matching data found. Try another combination.")
