import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("TP53 Node-1 Demo • LAML IC50 Explorer")

# Upload and Load Data
@st.cache_data
def load_data():
    return pd.read_csv("/mnt/data/CCLE_TP53_FLT3_NPM1_wide.csv")

ccle_wide = load_data()

# Sidebar - Select Subset
genetic_subset = st.sidebar.selectbox(
    "Select Genetic Subset",
    ["TP53 only", "TP53 + FLT3", "TP53 + FLT3 + NPM1"]
)

# Filter Based on Subset
if genetic_subset == "TP53 only":
    subset_df = ccle_wide[(ccle_wide["TP53"] == 1) & (ccle_wide["FLT3"] == 0) & (ccle_wide["NPM1"] == 0)]
elif genetic_subset == "TP53 + FLT3":
    subset_df = ccle_wide[(ccle_wide["TP53"] == 1) & (ccle_wide["FLT3"] == 1) & (ccle_wide["NPM1"] == 0)]
elif genetic_subset == "TP53 + FLT3 + NPM1":
    subset_df = ccle_wide[(ccle_wide["TP53"] == 1) & (ccle_wide["FLT3"] == 1) & (ccle_wide["NPM1"] == 1)]

# Sidebar - Drug and Cancer Type Filter
base = subset_df.copy()

drug_name = st.sidebar.selectbox(
    "Select Drug",
    sorted(base["Drug"].dropna().unique())
)
cancer_type = st.sidebar.selectbox(
    "Select Cancer Type",
    sorted(base["TCGA Classification"].dropna().unique())
)

# Filter for Selection
data = base[(base["Drug"] == drug_name) & (base["TCGA Classification"] == cancer_type)]

# Title
st.markdown(f"### Drug Sensitivity for `{drug_name}` in `{cancer_type}`")

# Plot IC50
st.markdown("#### IC50 Distribution")
fig, ax = plt.subplots()
colors = data["Mut/WT"].map({"Mut": "red", "WT": "blue"})
ax.bar(data["Cell Line Name"], data["IC50 (uM)"], color=colors)
ax.set_ylabel("IC50 (uM)")
ax.set_xlabel("Cell Line")
ax.set_xticklabels(data["Cell Line Name"], rotation=90)
ax.set_title("IC50 by Cell Line")
st.pyplot(fig)

# Plot AUC
st.markdown("#### AUC Distribution")
fig2, ax2 = plt.subplots()
ax2.bar(data["Cell Line Name"], data["AUC"], color=colors)
ax2.set_ylabel("AUC")
ax2.set_xlabel("Cell Line")
ax2.set_xticklabels(data["Cell Line Name"], rotation=90)
ax2.set_title("AUC by Cell Line")
st.pyplot(fig2)

# Show Raw Table
display_table = st.checkbox("Show Raw Data Table")
if display_table:
    st.dataframe(data.reset_index(drop=True))

# Footer
st.markdown("---")
st.markdown("**Data Source:** GDSC / CCLE with mutation filtering")
