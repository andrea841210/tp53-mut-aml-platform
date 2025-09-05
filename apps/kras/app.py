# --- Diagnostics ---
st.subheader("Diagnostics")

try:
    total = int(len(plot_df))
    depmap_matched = int(plot_df["DepMap_ID"].notna().sum()) if "DepMap_ID" in plot_df.columns else 0
    mut_count = int(plot_df.get("Mut", pd.Series(dtype=bool)).sum()) if "Mut" in plot_df.columns else 0
    st.markdown(f"**Rows (after Top N):** {total} | **Available (before Top N):** {available_rows} | **DepMap_ID matched:** {depmap_matched} | **Mutants:** {mut_count} | **Collapsed duplicates:** {collapsed}")
except Exception:
    pass

# --- Aggregation audit (median vs min/max, per cell line) ---
with st.expander("Aggregation audit (median vs min/max, per cell line)", expanded=False):
    try:
        rep = joined.copy()
        rep = rep[rep["IC50_uM"].notna()]
        key_cols = [c for c in ["Drug_Name","CellLine","DepMap_ID","TCGA_Classification"] if c in rep.columns]
        agg = rep.groupby(key_cols, as_index=False).agg(
            n_reps=("IC50_uM","size"),
            min_ic50=("IC50_uM","min"),
            median_ic50=("IC50_uM","median"),
            max_ic50=("IC50_uM","max")
        ).sort_values("median_ic50")
        st.dataframe(agg.head(200), use_container_width=True, hide_index=True)
        st.caption("Audit view shows replicate count and min/median/max per (Drug, Cell line). Plot uses median; Sheets may resemble min.")
    except Exception as e:
        st.write("audit error:", e)

# --- TCGA reference codes ---
with st.expander("TCGA code reference (common groups)", expanded=False):
    st.markdown("- CRC → COREAD  \n- NSCLC → LUAD+LUSC  \n- PDAC → PAAD  \n- BRCA → BRCA  \n- SKCM → SKCM  \n- GBM → GBM")

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
