TP53 Node-1 Demo App

This Streamlit app visualizes IC50 values of selected drugs in AML (LAML) cell lines,
highlighting mutant (Mut) vs wild-type (WT) status for a given gene (e.g., TP53).

🚀 Usage
🔹 Lite Mode (default)

Built-in dataset (data/ folder in this repo)

GDSC1_LAML.xlsx

CCLE_TP53.csv

Model_DepMap.csv

Cancer type is fixed to LAML

Drug and Gene are free text inputs

👉 Easiest way to try: just type a drug (e.g., Rapamycin) and a gene (TP53), then click Run

🔹 Advanced Mode

Upload your own datasets:

GDSC IC50 table (CSV or Excel, full version)

Mutation table (CSV, e.g., CCLE mutation subset)

DepMap metadata / mapper (CSV)

Cancer type / Drug / Gene all customizable

Useful if you want to test genes other than TP53, or explore cancer types beyond LAML.

📦 Installation (Local)

Clone this repo and install requirements:

pip install -r requirements.txt
streamlit run app.py

🌐 Deployment

This app can be deployed to:

Streamlit Cloud

HuggingFace Spaces

Or run locally with Python

📝 Notes

Lite mode uses small built-in files for quick demo and teaching.

Advanced mode allows flexible exploration with full-size GDSC/CCLE datasets.

Visual output:

Filtered IC50 table

Bar chart: red = Mut, blue = WT

Created for demo purposes by Andrea × GPT Mirror Project.
