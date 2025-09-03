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
    # Load and standardize column names
    if drug_data_file.name.endswith(".csv"):
        drug_data = pd.read_csv(drug_data_file)
    else:
        drug_data = pd.read_excel(drug_data_file)
    mutation_data = pd.read_csv(mutation_data_file)

    drug_data.columns = drug_data.columns.str.strip().str.replace(" ", "_")
    mutation_data.columns = mutation_data.columns.str.strip().str.replace(" ", "_")

    # Confirm merge key exists
    if "DepMap_ID" not in drug_data.columns or "DepMap_ID" not in mutation_data.columns:
        st.error("Error: 'DepMap_ID' column not found in both files. Please check your inputs.")
    elif gene_name not in mutation_data.columns:
        st.error(f"Error: Gene '{gene_name}' not found in mutation file columns.")
    else:
        # Filter for drug
        drug_filtered = drug_data[drug_data["Drug_Name"].str.lower() == drug_name.lower()].copy()

        # Merge
        merged = pd.merge(drug_filtered, mutation_data, on="DepMap_ID", how="left")

        # Label mutation status
        merged["Mutation_Status"] = merged[gene_name].apply(
            lambda x: "Mut" if pd.notna(x) and x != "WT" else "WT"
        )

        # Sort
        merged.sort_values("IC50_(uM)", inplace=True)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 4))
        colors = merged["Mutation_Status"].map({"WT": "blue", "Mut": "red"})
        ax.bar(merged["Cell_Line_Name"], merged["IC50_(uM)"], color=colors)
        ax.set_ylabel("IC50 (uM)")
        ax.set_title(f"{drug_name} Sensitivity by {gene_name} Mutation Status")
        ax.tick_params(axis="x", rotation=90)
        st.pyplot(fig)

        # Download merged file
        st.download_button(
            label="Download Merged CSV",
            data=merged.to_csv(index=False),
            file_name="merged_output.csv",
            mime="text/csv"
        )
else:
    st.info("Please upload both GDSC and Mutation CSV files.")
