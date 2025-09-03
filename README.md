# TP53 Node-1 Demo App

This interactive Python app visualizes IC50 values of selected drugs on LAML cell lines,
highlighting mutant (Mut) vs wild-type (WT) status for a given gene (e.g., TP53).

## 💻 Usage

1. Upload the following files:
   - GDSC1_LAML.xlsx
   - CCLE_mutation_subset.csv (e.g., TP53)
   - DepMap metadata (ModelID / DepMapID mapper)

2. Input three parameters:
   - Drug name (e.g., Rapamycin)
   - Gene name (e.g., TP53)
   - Cancer type (e.g., LAML)

3. Output:
   - Filtered IC50 table for matching cell lines
   - Bar chart: red (Mut), blue (WT)

## 📦 Installation (Local)

```bash
pip install -r requirements.txt
python app.py
```

## 🌐 Deployment (Optional)
You may deploy this Streamlit app on platforms like Streamlit Cloud, HuggingFace Spaces, or locally via Python.

---

Created for demo purposes by Andrea × GPT Mirror Project.