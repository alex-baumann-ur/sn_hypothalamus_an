### Differential expression between treatment and control groups

celltypes = ["GLU_Neurons", "Oligodendrocytes", "GABA_Neurons", "Astrocytes", "OPC", "Microglia", "Ependymal", "Vascular", "ParsTuber", "Immune", "all_markers"]

# Load unfiltered differential-expression tables for each cell type.
comparison_dir = (Path(de_analysis_path) /
                   f"de_markers_{control_label.lower()}_vs_{treatment_label.lower()}")
de_by_cell_type = {}
for celltype in celltypes:
    path = (comparison_dir / "unfiltered" /
            f"{dataset_id}_{celltype}_de_markers_{treatment_label}_vs_{control_label}.csv")
    de_by_cell_type[celltype] = pd.read_csv(path)

# Filter DE tables and write the files consumed by the downstream R scripts.
celltypes = ["GLU_Neurons", "Oligodendrocytes", "GABA_Neurons", "Astrocytes", "OPC", "Microglia", "Ependymal", "Vascular", "ParsTuber", "Immune", "all_markers"]
celltypes_log2fc_cutoff_de = {celltype: (0.5 if celltype in {"GLU_Neurons", "Oligodendrocytes",
                           "GABA_Neurons", "Astrocytes", "OPC", "Microglia", "all_markers"} else None)
                           for celltype in celltypes}
output_dir = comparison_dir / "filtered"
output_dir.mkdir(parents=True, exist_ok=True)
columns = ["gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
de_dict_filtered_all = {}


for celltype in celltypes:
    cutoff = celltypes_log2fc_cutoff_de[celltype]
    if cutoff is None:
        continue
    de = de_by_cell_type[celltype].copy()
    detected = (de["pct.1"] > 0.1) | (de["pct.2"] > 0.1)
    significant = de["p_val_adj"] < 0.05
    positive = de.loc[significant & (de["avg_log2FC"] > cutoff) & detected, columns].sort_values("avg_log2FC", ascending=False)
    negative = de.loc[significant & (de["avg_log2FC"] < -cutoff) & detected, columns].sort_values("avg_log2FC")
    prefix = f"{dataset_id}_{celltype}_de_filtered_"
    suffix = "_log2fc0_5_pval005_pct01_"
    positive.to_csv(output_dir / f"{prefix}pos{suffix}{treatment_label}.csv", index=False)
    negative.to_csv(output_dir / f"{prefix}neg{suffix}{control_label}.csv", index=False)
    combined = de.loc[significant & (de["avg_log2FC"].abs() > cutoff) & detected, columns].copy()
    combined = combined.sort_values("avg_log2FC", ascending=False)
    combined.to_csv(output_dir / f"{prefix}pos{suffix}ALL.csv", index=False)
    de_dict_filtered_all[celltype] = combined

### Cell-type-specific analyses

# Load unfiltered differential-expression tables for the selected neuronal populations.
specific_celltypes = ["agrp_neurons", "pomc_neurons"]
specific_dir = comparison_dir / "other_explicit_de" / "unfiltered"
de_dict_spec_cells = {}
for celltype in specific_celltypes:
    path = specific_dir / f"de_markers_{celltype}_{treatment_label}_vs_{control_label}.csv"
    de_dict_spec_cells[celltype] = pd.read_csv(path)

# Load unfiltered differential-expression tables for the selected neuronal populations.
specific_celltypes = ["agrp_neurons", "pomc_neurons"]
specific_dir = comparison_dir / "other_explicit_de" / "unfiltered"
de_dict_spec_cells = {}
for celltype in specific_celltypes:
    path = specific_dir / f"de_markers_{celltype}_{treatment_label}_vs_{control_label}.csv"
    de_dict_spec_cells[celltype] = pd.read_csv(path)

# Filter the selected neuronal-population tables and save combined and direction-specific results.
specific_output_dir = comparison_dir / "other_explicit_de" / "filtered"
specific_output_dir.mkdir(parents=True, exist_ok=True)
de_dict_filtered_all_spec_cells = {}
for celltype, de in de_dict_spec_cells.items():
    detected = (de["pct.1"] > 0.1) | (de["pct.2"] > 0.1)
    significant = de["p_val_adj"] < 0.05
    positive = de.loc[significant & (de["avg_log2FC"] > 0.5) & detected, columns].sort_values("avg_log2FC", ascending=False)
    negative = de.loc[significant & (de["avg_log2FC"] < -0.5) & detected, columns].sort_values("avg_log2FC")
    combined = de.loc[significant & (de["avg_log2FC"].abs() > 0.5) & detected, columns].sort_values("avg_log2FC", ascending=False)
    prefix = f"de_markers_{celltype}_filtered_"
    positive.to_csv(specific_output_dir / f"{prefix}pos_log2fc0_5_pval005_pct01_{treatment_label}.csv", index=False)
    negative.to_csv(specific_output_dir / f"{prefix}neg_log2fc0_5_pval005_pct01_{control_label}.csv", index=False)
    combined.to_csv(specific_output_dir / f"{prefix}log2fc0_5_pval005_pct01_ALL.csv", index=False)
    de_dict_filtered_all_spec_cells[celltype] = combined

