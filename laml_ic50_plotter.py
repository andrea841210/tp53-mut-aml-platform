
import pandas as pd
import matplotlib.pyplot as plt
import io
import base64

def load_data(gdsc_path, ccle_path, depmap_path):
    gdsc_df = pd.read_excel(gdsc_path)
    ccle_df = pd.read_csv(ccle_path)
    depmap_df = pd.read_csv(depmap_path)
    return gdsc_df, ccle_df, depmap_df

def filter_laml_data(gdsc_df, drug_name):
    laml_df = gdsc_df[gdsc_df['TCGA Classification'] == 'LAML']
    drug_df = laml_df[laml_df['Drug Name'].str.lower() == drug_name.lower()]
    return drug_df

def merge_with_mutation(drug_df, ccle_df, depmap_df, gene_name):
    merged_df = pd.merge(drug_df, depmap_df, on='DepMap ID', how='left')
    merged_df = pd.merge(merged_df, ccle_df, on='ModelID', how='left')
    merged_df['Mut/WT'] = merged_df['Hugo_Symbol'].apply(
        lambda x: 'Mut' if pd.notnull(x) and x.lower() == gene_name.lower() else 'WT'
    )
    return merged_df

def plot_ic50(merged_df):
    merged_df_sorted = merged_df.sort_values(by='IC50 (uM)', ascending=True)
    colors = ['red' if mut == 'Mut' else 'blue' for mut in merged_df_sorted['Mut/WT']]
    plt.figure(figsize=(10, 6))
    bars = plt.bar(merged_df_sorted['Cell Line Name'], merged_df_sorted['IC50 (uM)'], color=colors)
    plt.xticks(rotation=90)
    plt.ylabel('IC50 (uM)')
    plt.title('IC50 values of LAML Cell Lines (Mut vs WT)')
    plt.tight_layout()

    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png')
    img_buf.seek(0)
    img_base64 = base64.b64encode(img_buf.read()).decode('utf-8')
    plt.close()
    return img_base64
