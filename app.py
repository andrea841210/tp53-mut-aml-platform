import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.title("TP53 Drug Sensitivity Explorer")

# Sidebar Inputs
drug_name = st.sidebar.text_input("Enter Drug Name", "Rapamycin")
gene_name = st.sidebar.text_input("Enter Gene Symbol", "TP53")

# Upload CSVs
drug_data_file = st.file_uploader("Upload GDSC Drug Sensitivity CSV", type=["csv", "xlsx"])
mutation_data_file = st.file_uploader("Upload CCLE Mutation CSV", type=["csv"])

if drug_data_file and mutation_data_file:
    # Load files
    if drug_data_file.name.endswith(".csv"):
        drug_data = pd.read_csv(drug_data_file)
    else:
        drug_data = pd.read_excel(drug_data_file)

    mutation_data = pd.read_csv(mutation_data_file)
    # 標準化欄位名稱，避免 DepMap ID / DepMap_ID 不一致
    mutation_data.columns = mutation_data.columns.str.replace(" ", "_")

    # Filter for drug
    drug_filtered = drug_data[drug_data["Drug Name"].str.lower() == drug_name.lower()].copy()

    # Merge on DepMap_ID
    merged = pd.merge(drug_filtered, mutation_data, on="DepMap_ID", how="left")

    # Label mutation status
    merged["Mutation Status"] = merged[gene_name].apply(lambda x: "Mut" if pd.notna(x) and x != "WT" else "WT")

    # Sort by IC50
    merged.sort_values("IC50 (uM)", inplace=True)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    colors = merged["Mutation Status"].map({"WT": "blue", "Mut": "red"})
    ax.bar(merged["Cell Line Name"], merged["IC50 (uM)"], color=colors)
    ax.set_ylabel("IC50 (uM)")
    ax.set_title(f"{drug_name} Sensitivity in Cell Lines ({gene_name} Mutation Status)")
    ax.tick_params(axis="x", rotation=90)
    st.pyplot(fig)

    # Optional CSV export
    st.download_button("Download Merged CSV", merged.to_csv(index=False), file_name="merged_output.csv", mime="text/csv")
else:
    st.info("Please upload both GDSC and Mutation CSV files.")
