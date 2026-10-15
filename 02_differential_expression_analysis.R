# Differential expression analysis
# Set PROJECT_DIR, DATASET_ID, CONTROL_LABEL, and TREATMENT_LABEL before running.
# Group labels must match the values in orig.ident and Condition; set the cluster IDs when
# running the cell-type-specific Agrp and Pomc analyses.

# Configuration
project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/project")
dataset_id <- Sys.getenv("DATASET_ID", unset = "dataset")
control_label <- Sys.getenv("CONTROL_LABEL", unset = "control")
treatment_label <- Sys.getenv("TREATMENT_LABEL", unset = "treatment")
agrp_cluster_id <- Sys.getenv("AGRP_CLUSTER_ID", unset = "your_agrp_cluster")
pomc_cluster_id <- Sys.getenv("POMC_CLUSTER_ID", unset = "your_pomc_cluster")
object_dir <- file.path(project_dir, "results", "seurat_objects", "")
de_dir <- file.path(project_dir, "results", "differential_expression", "")
figure_dir <- file.path(project_dir, "results", "figures", "")
overview_dir <- file.path(figure_dir, "overview_tables", "")
volcano_dir <- file.path(figure_dir, "volcano_plots", "")
plots_dir <- file.path(project_dir, "results", "marker_gene_plots", "")
module_scores_dir <- file.path(project_dir, "results", "module_scores", "")

# Load required packages
library(dplyr)
library(Seurat) # Single cell analysis package install.packages('Seurat')
library(ggplot2) # Needed for plotting
library(ggrepel)
library(devtools) # To be able to install tools from github
library(harmony) # https://github.com/immunogenomics/harmony or install.packages("harmony") # package install necessary
library(Rcpp)
library(data.table)
library(DropletUtils) # package install necessary
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(zeallot) # For %<-% that unpacks lists in the Python manner # package install necessary
library(magrittr)
library("biomaRt") # Used for transcript annotation # package install necessary
library("psych") # needed for the dfOrder function
library("gridExtra")
library("cowplot")
library("patchwork")
library(presto)


# Load the integrated Seurat object produced by script 01
seurat_object <- readRDS(file.path(object_dir, paste0(dataset_id, "_integrated_cca_singlets.rds")))

# really only select singlets here
seurat_object <- subset(seurat_object, subset = scDblFinder.class == "singlet")

head(seurat_object@meta.data)
seurat_object
write.csv(table(seurat_object$predicted), paste0(
  overview_dir,
  "cell_numbers_clusters_all.csv"
))
write.csv(table(seurat_object$predicted, seurat_object$orig.ident), paste0(
  overview_dir,
  "cell_numbers_clusters_by_sample.csv"
))
write.csv(table(seurat_object$predicted, seurat_object$predicted_renamed), paste0(
  overview_dir,
  "cell_numbers_clusters_by_big_cell_types.csv"
))


### Here begins the differential expression analysis ###

pdf(file.path(figure_dir, paste0(dataset_id, "_integrated_cca_singlets_UMAP.pdf")), width = 10, height = 7)
plot(DimPlot(seurat_object, reduction = "integrated.cca", group.by = "orig.ident", label = TRUE))
plot(DimPlot(seurat_object, reduction = "umap_scvi", group.by = "predicted_renamed", label = TRUE))
plot(DimPlot(seurat_object, reduction = "umap_scvi", group.by = "orig.ident", label = TRUE))
plot(DimPlot(seurat_object, reduction = "umap_cca", group.by = "predicted_renamed", label = TRUE))
plot(DimPlot(seurat_object, reduction = "umap_cca", group.by = "orig.ident", label = TRUE))
dev.off()

pdf(file.path(figure_dir, paste0(dataset_id, "_integrated_cca_singlets_UMAP_large.pdf")), width = 30, height = 20)
plot(DimPlot(seurat_object, reduction = "umap_scvi", group.by = "predicted", label = TRUE))
plot(DimPlot(seurat_object, reduction = "umap_cca", group.by = "predicted", label = TRUE))
dev.off()

### Identify conserved cell type markers ###
Idents(seurat_object) <- "predicted_renamed"

head(seurat_object@meta.data)

conserved_markers_list <- list()
for (cell_type in unique(seurat_object$predicted_renamed)) {
  de <- FindConservedMarkers(seurat_object, ident.1 = cell_type, grouping.var = "orig.ident")
  de$gene <- rownames(de)
  conserved_markers_list[[cell_type]] <- de
}

# unfiltered markers to CSV
for (cell_type in names(conserved_markers_list)) {
  write.csv(conserved_markers_list[[cell_type]], paste0(
    de_dir, "conserved_markers/unfiltered/", paste0(dataset_id, "_"), cell_type,
    "_conserved_markers.csv"
  ))
}

# filter the markers
# -> with notebook 02_filter_de_genes.ipynb

# Load the filtered conserved markers
celltypes_log2fc_cutoff_cons <- c(
  "GLU_Neurons" = "1", "Oligodendrocytes" = "2", "GABA_Neurons" = "1", "Astrocytes" = "2", "OPC" = "2",
  "Microglia" = "2", "Ependymal" = "2", "Vascular" = "2", "ParsTuber" = "2", "Immune" = "2"
)

conserved_markers_list_filtered <- list()
for (cell_type in names(celltypes_log2fc_cutoff_cons)) {
  conserved_markers_list_filtered[[cell_type]] <- read.csv(paste0(
    de_dir, "conserved_markers/filtered/", paste0(dataset_id, "_"), cell_type,
    "_conserved_markers_filtered_pos_log2fc", celltypes_log2fc_cutoff_cons[[cell_type]], "_pval005_pct01.csv"
  ))
}


# DotPlots with conserved markers
Idents(seurat_object) <- "predicted_renamed"

pdf(paste0(de_dir, "conserved_markers/", paste0(dataset_id, "_top_10_conserved_markers_dotplot.pdf")), width = 10, height = 5)
for (cell_type in names(conserved_markers_list_filtered)) {
  plot(DotPlot(seurat_object, features = conserved_markers_list_filtered[[cell_type]]$gene[1:10], group.by = "predicted_renamed") + ggtitle(paste0("Top 10 conserved markers (upregulated) for ", cell_type, " (by log2FC)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  plot(DotPlot(seurat_object,
    features = conserved_markers_list_filtered[[cell_type]][order(conserved_markers_list_filtered[[cell_type]]$max_pval), ]$gene[1:10],
    group.by = "predicted_renamed"
  ) + ggtitle(
    paste0("Top 10 conserved markers (upregulated) for ", cell_type, " (by pval)")
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
dev.off()


# Compare treatment and control groups within each cell type
# Aggregate expression for each cell type and condition
aggregate_seurat <- AggregateExpression(seurat_object, group.by = c("predicted_renamed", "orig.ident"), return.seurat = TRUE)
aggregate_seurat$celltype_condition <- rownames(aggregate_seurat@meta.data)
aggregate_seurat$condition <- sapply(strsplit(aggregate_seurat$celltype_condition, "_"), function(x) x[2])
aggregate_seurat$celltype <- sapply(strsplit(aggregate_seurat$celltype_condition, "_"), function(x) x[1])

head(aggregate_seurat@meta.data)


# Unfiltered DE markers and cell-scatter plots for each cell type
Idents(seurat_object) <- "celltype_condition"

de_markers_list <- list()
for (cell_type in unique(seurat_object$predicted_renamed)) {
  markers <- FindMarkers(seurat_object, ident.1 = paste0(paste0(treatment_label, "_"), cell_type), ident.2 = paste0(paste0(control_label, "_"), cell_type))
  markers$gene <- rownames(markers)
  de_markers_list[[cell_type]] <- markers
  write.csv(markers, paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/unfiltered/"), paste0(dataset_id, "_"), cell_type, paste0("_de_markers_", treatment_label, paste0("_vs_", control_label, ".csv"))))

  p1 <- CellScatter(aggregate_seurat, paste0(gsub("_", "-", cell_type), paste0("_", treatment_label)), paste0(gsub("_", "-", cell_type), paste0("_", control_label)),
    highlight = markers$gene[1:10]
  )
  p2 <- LabelPoints(plot = p1, points = markers$gene[1:10], repel = TRUE, xnudge = 0, ynudge = 0, color = "darkred", max.overlaps = 100)
  pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), paste0("/cellscatter_plots/", dataset_id, "_")), cell_type, paste0("_de_markers_", treatment_label, paste0("_vs_", control_label, "_cellscatter.pdf"))))
  plot(p2)
  dev.off()
}

# Markers across all cell types
Idents(seurat_object) <- "orig.ident"
de_markers_list[["all_markers"]] <- FindMarkers(seurat_object, ident.1 = treatment_label, ident.2 = control_label)
de_markers_list[["all_markers"]]$gene <- rownames(de_markers_list[["all_markers"]])
write.csv(de_markers_list[["all_markers"]], paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/unfiltered/"),
  paste0(dataset_id, paste0("_all_markers_de_markers_", treatment_label, paste0("_vs_", control_label, ".csv")))
))

aggregate_orig_seurat <- AggregateExpression(seurat_object, group.by = "orig.ident", return.seurat = TRUE)

p1 <- CellScatter(aggregate_orig_seurat, treatment_label, control_label, highlight = de_markers_list[["all_markers"]]$gene[1:10])
p2 <- LabelPoints(
  plot = p1, points = de_markers_list[["all_markers"]]$gene[1:10], repel = TRUE, xnudge = 0, ynudge = 0, color = "darkred",
  max.overlaps = 100
)
pdf(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/cellscatter_plots/"),
  paste0(dataset_id, paste0("_all_markers_de_markers_", treatment_label, paste0("_vs_", control_label, "_cellscatter.pdf")))
))
plot(p2)
dev.off()

unique(seurat_object$predicted_renamed)

# Read the DE marker lists
samples_de <- c("GLU_Neurons", "GABA_Neurons", "Oligodendrocytes", "Astrocytes", "OPC", "Microglia", "all_markers")

de_markers_list_unfiltered <- list()
for (cell_type in samples_de) {
  de_markers_list_unfiltered[[cell_type]] <- read.csv(paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/unfiltered/"),
    paste0(dataset_id, "_"), cell_type, paste0("_de_markers_", treatment_label, paste0("_vs_", control_label, ".csv"))
  ))
}
de_markers_list_unfiltered[["all_markers"]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/unfiltered/"),
  paste0(dataset_id, paste0("_all_markers_de_markers_", treatment_label, paste0("_vs_", control_label, ".csv")))
))

de_markers_list_unfiltered[["Agrp"]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/unfiltered/"),
  paste0("de_markers_agrp_neurons_", treatment_label, paste0("_vs_", control_label, ".csv"))
))
de_markers_list_unfiltered[["Pomc"]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/unfiltered/"),
  paste0("de_markers_pomc_neurons_", treatment_label, paste0("_vs_", control_label, ".csv"))
))

names(de_markers_list_unfiltered)

# Filter the DE markers
# -> with notebook 03_filter_de_genes.ipynb

# Filtered DE markers and cell-scatter plots for each cell type
samples_de <- c("GLU_Neurons", "GABA_Neurons", "Oligodendrocytes", "Astrocytes", "OPC", "Microglia", "all_markers")

de_markers_filtered_list <- list()

for (cell_type in samples_de) {
  de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]] <- read.csv(paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/filtered/"),
    paste0(dataset_id, "_"), cell_type, paste0("_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
  ))
  de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]] <- read.csv(paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/filtered/"),
    paste0(dataset_id, "_"), cell_type, paste0("_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
  ))
}

de_markers_filtered_list[[paste0("Agrp_", treatment_label, "")]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/"),
  paste0("de_markers_agrp_neurons_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
de_markers_filtered_list[[paste0("Agrp_", control_label, "")]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/"),
  paste0("de_markers_agrp_neurons_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
de_markers_filtered_list[[paste0("Pomc_", treatment_label, "")]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/"),
  paste0("de_markers_pomc_neurons_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
de_markers_filtered_list[[paste0("Pomc_", control_label, "")]] <- read.csv(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/"),
  paste0("de_markers_pomc_neurons_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))

names(de_markers_filtered_list)

head(de_markers_filtered_list[[paste0("GLU_Neurons_", treatment_label, "")]])
head(de_markers_filtered_list[[paste0("GLU_Neurons_", control_label, "")]])

head(de_markers_filtered_list[[paste0("Pomc_", treatment_label, "")]])
head(de_markers_filtered_list[[paste0("Pomc_", control_label, "")]])

head(Idents(seurat_object))

# DotPlots with Top 10 genes of positive and negative log2FC -> only log2FC of 0.5
Idents(seurat_object) <- "orig.ident"
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/"), paste0(dataset_id, paste0("_top_10_de_markers_", treatment_label, paste0("_vs_", control_label, "_dotplot.pdf")))),
  width = 10, height = 5
)
for (cell_type in samples_de) {
  plot(DotPlot(seurat_object,
    features = de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]]$gene[1:10],
    group.by = "orig.ident"
  ) + ggtitle(paste0(
    "Top 10 upregulated genes for ", cell_type,
    paste0(" in ", treatment_label, " (by log2FC)")
  )) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  plot(DotPlot(seurat_object,
    features = de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]]$gene[1:10],
    group.by = "orig.ident"
  ) + ggtitle(paste0(
    "Top 10 downregulated genes for ", cell_type,
    paste0(" in ", treatment_label, paste0(" (by log2FC, up in ", control_label, ")"))
  )) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  plot(DotPlot(seurat_object,
    features = de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]][
      order(de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]]$p_val_adj),
    ]$gene[1:10],
    group.by = "orig.ident"
  ) + ggtitle(paste0(
    "Top 10 upregulated genes for ", cell_type,
    paste0(" in ", treatment_label, " (by pval)")
  )) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  plot(DotPlot(seurat_object,
    features = de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]][
      order(de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]]$p_val_adj),
    ]$gene[1:10],
    group.by = "orig.ident"
  ) + ggtitle(paste0(
    "Top 10 downregulated genes for ", cell_type,
    paste0(" in ", treatment_label, paste0(" (by pval, up in ", control_label, ")"))
  )) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
dev.off()


### Volcano plots ###
samples_de <- c(
  "GLU_Neurons", "GABA_Neurons", "Oligodendrocytes", "Astrocytes", "OPC", "Microglia", "all_markers",
  "Agrp", "Pomc"
)

volcano_plot_files <- list()
prepare_data_for_volcano_plots <- function(de_genes, pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1) {
  de_genes_inside <- de_genes
  de_genes_inside$significance <- "Not significant"
  de_genes_inside$significance[de_genes_inside$p_val_adj < pval_threshold &
    de_genes_inside$avg_log2FC > log2FC_threshold &
    (de_genes_inside$pct.1 > pct_threshold | de_genes_inside$pct.2 > pct_threshold)] <- treatment_label
  de_genes_inside$significance[de_genes_inside$p_val_adj < pval_threshold &
    de_genes_inside$avg_log2FC < -log2FC_threshold &
    (de_genes_inside$pct.1 > pct_threshold | de_genes_inside$pct.2 > pct_threshold)] <- control_label
  return(de_genes_inside)
}

for (i in samples_de) {
  volcano_plot_files[[i]] <- prepare_data_for_volcano_plots(de_markers_list_unfiltered[[i]],
    pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
  )
}

names(volcano_plot_files)
table(volcano_plot_files[["Agrp"]]$significance)

head(volcano_plot_files[["GLU_Neurons"]])

plot_volcano_plot <- function(cell_type_id, cell_type_name, volcano_plot_tables, pval_threshold = 0.05, log2FC_threshold = 1.5, pct_threshold = 0.1) {
  de_genes_inside <- volcano_plot_tables[[cell_type_id]]
  genes_1 <- de_genes_inside[de_genes_inside$significance == treatment_label, ]
  genes_2 <- de_genes_inside[de_genes_inside$significance == control_label, ]
  top_genes_top1_log2fc <- genes_1[order(-genes_1$avg_log2FC), ][1:10, ]
  top_genes_top2_log2fc <- genes_2[order(genes_2$avg_log2FC), ][1:10, ]
  top_genes_top1_pval <- genes_1[order(genes_1$p_val_adj), ][1:10, ]
  top_genes_top2_pval <- genes_2[order(genes_2$p_val_adj), ][1:10, ]

  top_genes_all <- unique(rbind(top_genes_top1_log2fc, top_genes_top2_log2fc, top_genes_top1_pval, top_genes_top2_pval))

  tiff(
    file.path(volcano_dir, paste0("VolcanoPlot_", cell_type_id, "_", treatment_label, "_vs_", control_label, ".tiff")),
    units = "cm", width = 20, height = 16, res = 600
  )
  print(paste0("Plotting volcano plot for ", cell_type_id))
  plot(ggplot(de_genes_inside, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(data = top_genes_all, aes(label = gene), color = "black", max.overlaps = Inf) +
    scale_color_manual(values = setNames(c("red", "blue", "grey"), c(treatment_label, control_label, "Not significant"))) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    labs(
      x = "Average log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme_minimal() +
    ggtitle(paste0("Volcano plot for ", cell_type_name, " (", treatment_label, " vs. ", control_label, ")")))
  dev.off()
}


plot_volcano_plot(
  cell_type_id = "GLU_Neurons", cell_type_name = "GLU neurons", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "Oligodendrocytes", cell_type_name = "Oligodendrocytes", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "GABA_Neurons", cell_type_name = "GABA neurons", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "all_markers", cell_type_name = "All clusters", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "OPC", cell_type_name = "OPC", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "Astrocytes", cell_type_name = "Astrocytes", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "Agrp", cell_type_name = "Agrp neurons", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)
plot_volcano_plot(
  cell_type_id = "Pomc", cell_type_name = "Pomc neurons", volcano_plot_tables = volcano_plot_files,
  pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
)

for (i in names(volcano_plot_files)) {
  plot_volcano_plot(
    cell_type_id = i, volcano_plot_tables = volcano_plot_files,
    pval_threshold = 0.05, log2FC_threshold = 0.5, pct_threshold = 0.1
  )
}



Idents(seurat_object) <- "celltype_condition"
head(Idents(seurat_object))

# VlnPlot with first 10 genes of positive log2FC and first 10 genes of negative log2FC
for (cell_type in samples_de) {
  if (cell_type != "all_markers") {
    pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), "Violinplot_", cell_type, paste0(paste0(paste0("_de_markers_", treatment_label, "_vs_"), control_label, "_log2fc_"), treatment_label, ".pdf")),
      width = 20, height = 25
    )
    plot(VlnPlot(seurat_object,
      features = de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]]$gene[1:10],
      group.by = "orig.ident", idents = c(paste0(paste0(treatment_label, "_"), cell_type), paste0(paste0(control_label, "_"), cell_type)), ncol = 3
    ) + ggtitle(
      paste0("Top upregulated genes for ", cell_type, paste0(" in ", treatment_label, " (by log2fc)"))
    ))
    dev.off()
    pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), "Violinplot_", cell_type, paste0("_de_markers_", treatment_label, paste0(paste0("_vs_", control_label, "_log2fc_"), control_label, ".pdf"))),
      width = 20, height = 25
    )
    plot(VlnPlot(seurat_object,
      features = de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]]$gene[1:10],
      group.by = "orig.ident", idents = c(paste0(paste0(treatment_label, "_"), cell_type), paste0(paste0(control_label, "_"), cell_type)), ncol = 3
    ) + ggtitle(
      paste0("Top downregulated genes for ", cell_type, paste0(" in ", treatment_label, paste0(" (by log2fc, up in ", control_label, ")")))
    ))
    dev.off()
    pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), "Violinplot_", cell_type, paste0(paste0(paste0("_de_markers_", treatment_label, "_vs_"), control_label, "_adj_p_val_"), treatment_label, ".pdf")),
      width = 20, height = 25
    )
    plot(VlnPlot(seurat_object,
      features = de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]][order(
        de_markers_filtered_list[[paste0(cell_type, paste0("_", treatment_label))]]$p_val_adj
      ), ]$gene[1:10],
      group.by = "orig.ident", idents = c(paste0(paste0(treatment_label, "_"), cell_type), paste0(paste0(control_label, "_"), cell_type)), ncol = 3
    ) + ggtitle(
      paste0("Top upregulated genes for ", cell_type, paste0(" in ", treatment_label, " (by adj_p_val)"))
    ))
    dev.off()
    pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), "Violinplot_", cell_type, paste0("_de_markers_", treatment_label, paste0(paste0("_vs_", control_label, "_adj_p_val_"), control_label, ".pdf"))),
      width = 20, height = 25
    )
    plot(VlnPlot(seurat_object,
      features = de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]][order(
        de_markers_filtered_list[[paste0(cell_type, paste0("_", control_label))]]$p_val_adj
      ), ]$gene[1:10],
      group.by = "orig.ident", idents = c(paste0(paste0(treatment_label, "_"), cell_type), paste0(paste0(control_label, "_"), cell_type)), ncol = 3
    ) + ggtitle(
      paste0("Top downregulated genes for ", cell_type, paste0(" in ", treatment_label, paste0(" (by adj_p_val, up in ", control_label, ")")))
    ))
    dev.off()
  }
}


# VlnPlot for all markers
Idents(seurat_object) <- "orig.ident"

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), paste0(paste0(paste0("Violinplot_all_markers_de_markers_", treatment_label, "_vs_"), control_label, "_log2fc_"), treatment_label, ".pdf")),
  width = 20, height = 25
)
plot(VlnPlot(seurat_object,
  features = de_markers_filtered_list[[paste0("all_markers_", treatment_label, "")]]$gene[1:10],
  group.by = "orig.ident", ncol = 3
) + ggtitle(
  paste0(paste0("Top upregulated genes for all markers in ", treatment_label, " (by log2fc)"))
))
dev.off()
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), paste0("Violinplot_all_markers_de_markers_", treatment_label, paste0(paste0("_vs_", control_label, "_log2fc_"), control_label, ".pdf"))),
  width = 20, height = 25
)
plot(VlnPlot(seurat_object,
  features = de_markers_filtered_list[[paste0("all_markers_", control_label, "")]]$gene[1:10],
  group.by = "orig.ident", ncol = 3
) + ggtitle(
  paste0(paste0("Top downregulated genes for all markers in ", treatment_label, paste0(" (by log2fc, up in ", control_label, ")")))
))
dev.off()
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), paste0(paste0(paste0("Violinplot_all_markers_de_markers_", treatment_label, "_vs_"), control_label, "_adj_p_val_"), treatment_label, ".pdf")),
  width = 20, height = 25
)
plot(VlnPlot(seurat_object,
  features = de_markers_filtered_list[[paste0("all_markers_", treatment_label, "")]][order(
    de_markers_filtered_list[[paste0("all_markers_", treatment_label, "")]]$p_val_adj
  ), ]$gene[1:10],
  group.by = "orig.ident", ncol = 3
) + ggtitle(
  paste0(paste0("Top upregulated genes for all markers in ", treatment_label, " (by adj_p_val)"))
))
dev.off()
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), paste0("Violinplot_all_markers_de_markers_", treatment_label, paste0(paste0("_vs_", control_label, "_adj_p_val_"), control_label, ".pdf"))),
  width = 20, height = 25
)
plot(VlnPlot(seurat_object,
  features = de_markers_filtered_list[[paste0("all_markers_", control_label, "")]][order(
    de_markers_filtered_list[[paste0("all_markers_", control_label, "")]]$p_val_adj
  ), ]$gene[1:10],
  group.by = "orig.ident", ncol = 3
) + ggtitle(
  paste0(paste0("Top downregulated genes for all markers in ", treatment_label, paste0(" (by adj_p_val, up in ", control_label, ")")))
))
dev.off()


##### DE analysis of specific cell types #####
head(seurat_object@meta.data)
seurat_object$predicted_condition <- paste0(seurat_object$predicted, "_", seurat_object$orig.ident)
Idents(seurat_object) <- "predicted_condition"

unique(Idents(seurat_object))

unique(seurat_object$predicted_condition)

de_agrp_neurons <- FindMarkers(
  seurat_object,
  ident.1 = paste0(agrp_cluster_id, "_", treatment_label),
  ident.2 = paste0(agrp_cluster_id, "_", control_label)
)
de_agrp_neurons$gene <- rownames(de_agrp_neurons)
write.csv(de_agrp_neurons, paste0(de_dir, paste0(paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/de_markers_agrp_neurons_"), treatment_label, paste0("_vs_", control_label, ".csv"))), row.names = FALSE)

de_pomc_neurons <- FindMarkers(
  seurat_object,
  ident.1 = paste0(pomc_cluster_id, "_", treatment_label),
  ident.2 = paste0(pomc_cluster_id, "_", control_label)
)
de_pomc_neurons$gene <- rownames(de_pomc_neurons)
write.csv(de_pomc_neurons, paste0(de_dir, paste0(paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/de_markers_pomc_neurons_"), treatment_label, paste0("_vs_", control_label, ".csv"))), row.names = FALSE)


de_dict_specific <- list()
for (celltype in c("agrp_neurons", "pomc_neurons")) {
  de_file <- paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/unfiltered/de_markers_"), celltype, paste0("_", treatment_label, paste0("_vs_", control_label, ".csv")))
  de <- read.csv(de_file)
  de_dict_specific[[celltype]] <- de
}

head(de_dict_specific[["agrp_neurons"]])
head(de_dict_specific[["pomc_neurons"]])

de_dict_specific_filtered <- list()
for (celltype in names(de_dict_specific)) {
  de_file1 <- paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/de_markers_"), celltype, paste0("_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv"))
  de_file2 <- paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/filtered/de_markers_"), celltype, paste0("_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv"))
  de1 <- read.csv(de_file1)
  de2 <- read.csv(de_file2)
  de <- rbind(de1, de2)
  de_dict_specific_filtered[[paste0(celltype, paste0("_", treatment_label))]] <- de1
  de_dict_specific_filtered[[paste0(celltype, paste0("_", control_label))]] <- de2
}

head(de_dict_specific_filtered[[paste0("agrp_neurons_", treatment_label, "")]])
head(de_dict_specific_filtered[[paste0("agrp_neurons_", control_label, "")]])

Idents(seurat_object) <- "predicted"

DotPlot(seurat_object, features = de_dict_specific_filtered[[paste0("agrp_neurons_", treatment_label, "")]]$gene[1:10], idents = agrp_cluster_id, group.by = "orig.ident") + ggtitle(
  paste0(paste0("Top upregulated genes for Agrp neurons in ", treatment_label, " (by log2fc)"))
)

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), paste0("DotPlot_agrp_neurons_de_markers_", treatment_label, paste0("_vs_", control_label, ".pdf"))), width = 10, height = 5)
plot(DotPlot(seurat_object, features = de_dict_specific_filtered[[paste0("agrp_neurons_", treatment_label, "")]]$gene[1:10], idents = agrp_cluster_id, group.by = "orig.ident") + ggtitle(
  paste0(paste0("Top upregulated genes for Agrp neurons in ", treatment_label, " (by log2fc)"))
))
plot(DotPlot(seurat_object, features = de_dict_specific_filtered[[paste0("agrp_neurons_", control_label, "")]]$gene[1:10], idents = agrp_cluster_id, group.by = "orig.ident") + ggtitle(
  paste0(paste0("Top downregulated genes for Agrp neurons in ", treatment_label, paste0(" (by log2fc, up in ", control_label, ")")))
))
dev.off()
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), paste0("VlnPlot_agrp_neurons_de_markers_", treatment_label, paste0("_vs_", control_label, ".pdf"))), width = 10, height = 5)
plot(VlnPlot(seurat_object, features = de_dict_specific_filtered[[paste0("agrp_neurons_", treatment_label, "")]]$gene[1:10], idents = agrp_cluster_id, group.by = "orig.ident"))
plot(VlnPlot(seurat_object, features = de_dict_specific_filtered[[paste0("agrp_neurons_", control_label, "")]]$gene[1:10], idents = agrp_cluster_id, group.by = "orig.ident"))
dev.off()


pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), paste0("DotPlot_pomc_neurons_de_markers_", treatment_label, paste0("_vs_", control_label, ".pdf"))), width = 10, height = 5)
plot(DotPlot(seurat_object, features = de_dict_specific_filtered[[paste0("pomc_neurons_", treatment_label, "")]]$gene[1:10], idents = pomc_cluster_id, group.by = "orig.ident") + ggtitle(
  paste0(paste0("Top upregulated genes for Pomc neurons in ", treatment_label, " (by log2fc)"))
))
plot(DotPlot(seurat_object, features = de_dict_specific_filtered[[paste0("pomc_neurons_", control_label, "")]]$gene[1:10], idents = pomc_cluster_id, group.by = "orig.ident") + ggtitle(
  paste0(paste0("Top downregulated genes for Pomc neurons in ", treatment_label, paste0(" (by log2fc, up in ", control_label, ")")))
))
dev.off()
pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), paste0("VlnPlot_pomc_neurons_de_markers_", treatment_label, paste0("_vs_", control_label, ".pdf"))), width = 10, height = 5)
plot(VlnPlot(seurat_object, features = de_dict_specific_filtered[[paste0("pomc_neurons_", treatment_label, "")]]$gene[1:10], idents = pomc_cluster_id, group.by = "orig.ident"))
plot(VlnPlot(seurat_object, features = de_dict_specific_filtered[[paste0("pomc_neurons_", control_label, "")]]$gene[1:10], idents = pomc_cluster_id, group.by = "orig.ident"))
dev.off()


### DE analysis for specific marker genes ###

leptin_markers <- c("Pomc", "Npy", "Agrp", "Lepr")

Idents(seurat_object) <- "predicted_renamed"

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "FeaturePlots_leptin_markers.pdf"), width = 10, height = 5)
plot(FeaturePlot(seurat_object, features = "Lepr", label = TRUE, split.by = "orig.ident"))
plot(FeaturePlot(seurat_object, features = "Lepr", label = TRUE, split.by = "orig.ident", idents = c("GLU_Neurons", "GABA_Neurons")))
plot(FeaturePlot(seurat_object, features = "Pomc", label = TRUE, split.by = "orig.ident"))
plot(FeaturePlot(seurat_object, features = "Pomc", label = TRUE, split.by = "orig.ident", idents = c("GLU_Neurons", "GABA_Neurons")))
plot(FeaturePlot(seurat_object, features = "Npy", label = TRUE, split.by = "orig.ident"))
plot(FeaturePlot(seurat_object, features = "Npy", label = TRUE, split.by = "orig.ident", idents = c("GLU_Neurons", "GABA_Neurons")))
plot(FeaturePlot(seurat_object, features = "Agrp", label = TRUE, split.by = "orig.ident"))
plot(FeaturePlot(seurat_object, features = "Agrp", label = TRUE, split.by = "orig.ident", idents = c("GLU_Neurons", "GABA_Neurons")))
dev.off()


pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "DotPlots_leptin_markers_big_clusters.pdf"), width = 10, height = 5)
plot(DotPlot(seurat_object, features = leptin_markers, group.by = "predicted_renamed") + ggtitle(
  paste0("Leptin markers in all clusters")
))
plot(DotPlot(seurat_object, features = leptin_markers, group.by = "predicted_renamed", split.by = "orig.ident") + ggtitle(
  paste0("Leptin markers in all clusters, split by sample")
))
plot(DotPlot(seurat_object, features = leptin_markers, idents = c("GLU_Neurons", "GABA_Neurons"), group.by = "predicted_renamed") + ggtitle(
  paste0("Leptin markers in GLU and GABA neurons")
))
plot(DotPlot(seurat_object, features = leptin_markers, idents = c("GLU_Neurons", "GABA_Neurons"), group.by = "predicted_renamed", split.by = "orig.ident") + ggtitle(
  paste0("Leptin markers in GLU and GABA neurons, split by sample")
))
dev.off()

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "VlnPlots_leptin_markers_big_clusters.pdf"), width = 10, height = 7)
# plot(VlnPlot(seurat_object, features = leptin_markers, idents = c(agrp_cluster_id, pomc_cluster_id), group.by = "predicted"))
for (leptin_marker in leptin_markers) {
  plot(VlnPlot(seurat_object, features = leptin_marker, group.by = "predicted", split.by = "orig.ident"))
}
dev.off()

Idents(seurat_object) <- "predicted"

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "DotPlots_leptin_markers_small_clusters.pdf"), width = 10, height = 5)
plot(DotPlot(seurat_object, features = leptin_markers, idents = c(agrp_cluster_id, pomc_cluster_id), group.by = "predicted") + ggtitle(
  paste0("Leptin markers in Agrp and Pomc neurons")
))
plot(DotPlot(seurat_object, features = leptin_markers, idents = c(agrp_cluster_id, pomc_cluster_id), group.by = "predicted", split.by = "orig.ident") + ggtitle(
  paste0("Leptin markers in Agrp and Pomc neurons, split by sample")
))
dev.off()

pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "VlnPlots_leptin_markers_small_clusters.pdf"), width = 10, height = 7)
# plot(VlnPlot(seurat_object, features = leptin_markers, idents = c(agrp_cluster_id, pomc_cluster_id), group.by = "predicted"))
for (leptin_marker in leptin_markers) {
  plot(VlnPlot(seurat_object, features = leptin_marker, idents = c(agrp_cluster_id, pomc_cluster_id), group.by = "predicted", split.by = "orig.ident"))
}
dev.off()


pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/vln_plots/"), paste0(paste0(paste0("Violinplot_agrp_neurons_de_markers_", treatment_label, "_vs_"), control_label, "_log2fc_"), treatment_label, ".pdf")),
  width = 20, height = 25
)
plot(VlnPlot(subset_agrp_neurons, features = de_agrp_neurons$gene[1:10], group.by = "orig.ident", ncol = 3) + ggtitle(
  paste0(paste0("Top upregulated genes for Agrp neurons in ", treatment_label, " (by log2fc)"))
))
dev.off()
subset_agrp_neurons <- subset(seurat_object, subset = predicted %in% c(agrp_cluster_id))
DimPlot(subset_agrp_neurons, group.by = "orig.ident", label = TRUE)

subset_pomc_neurons <- subset(seurat_object, subset = predicted %in% c(pomc_cluster_id))
DimPlot(subset_pomc_neurons, group.by = "orig.ident", label = TRUE)


### Marker genes
genes_astrocytes <- c("Gfap", "Aldh1l1", "Aqp4")
genes_microglia <- c("Aif1", "Tmem119", "Cx3cr1")
genes_oligodendrocytes <- c("Olig2", "Olig1", "Sox10", "Mbp", "Mog", "Plp1")
genes_neurons <- c("Rbfox3", "Map2") # Rbfox3 also known as NeuN
genes_neuroendocrine <- c("Pomc", "Npy", "Agrp", "Cartpt", "Lepr", "Ghrl", "Ghsr", "Hcrt", "Bdnf", "Esr1", "Esr2", "Gabra1", "Gabra2")
genes_hypothalamus <- c("Rgs16", "Crh", "Trh", "Sim1", "Nr5a1")
genes_func_upr_isr <- c("Atf4", "Ddit3", "Atf3", "Hspa5") # Ddit3 is also known as Chop, Hspa5 is also known as BiP
genes_func_autophagy <- c("Map1lc3b", "Becn1", "Sqstm1", "Atg5", "Atg7") # Map1lc3b is also known as LC3B, Sqstm1 is also known as p62
genes_func_stress <- c("Hmox1", "Nfe2l2", "Tnf", "Il1b", "Il6", "Serpina3n", "C1qa", "Fos", "Ddit4") # Nfe2l2 is also known as Nrf2
genes_func_apoptosis <- c("Casp3", "Casp9", "Bax", "Bcl2", "Trp53", "H2ax") # Trp53 is also known as p53, H2afx replaced with H2ax
genes_func_mitochondrial <- c("Dnm1l", "Mfn1", "Mfn2", "Opa1", "Ucp2", "Slc25a27", "Ppargc1a", "Tfam") # Dnm1l is also known as Drp1, Ucp4 replaced with Slc25a27
genes_func_prolif <- c("Mki67", "Pcna", "Ccnd1", "Cdkn1a") # Cdkn1a is also known as p21
genes_func_synaptic <- c("Tnip1", "Nefm", "Nefh", "Syn1", "Dlg4", "Shank3", "Syp", "Snap25", "Nrxn1", "Nlgn1", "Nlgn3", "Gria1") # Dlg4 is also known as PSD-95, Nef replaced with Tnip1, Nrnx1 replaced with Nrxn1
genes_func_leptin <- c("Jak2", "Stat3", "Socs3", "Ptpn1", "Foxo1", "Cish", "Ptpn2", "Sh2b1", "Mc4r")

# Marker genes for module scores (high vs. low expression expected)
genes_astrocytes_high <- c("Gfap", "Aldh1l1", "Aqp4")
genes_microglia_high <- c("Aif1", "Cx3cr1")
genes_microglia_low <- "Tmem119"
genes_oligodendrocytes_low <- c("Olig2", "Olig1", "Sox10", "Mbp", "Mog", "Plp1")
genes_neurons_low <- c("Rbfox3", "Map2") # Rbfox3 also known as NeuN
genes_neuroendocrine_low <- c("Pomc", "Cartpt", "Bdnf", "Esr1", "Esr2")
genes_neuroendocrine_high <- c("Npy", "Agrp", "Lepr", "Ghsr", "Hcrt")
genes_neuroendocrine_unknown <- c("Ghrl", "Gabra1", "Gabra2")
genes_hypothalamus_low <- "Trh"
genes_hypothalamus_high <- "Crh"
genes_hypothalamus_unknown <- c("Rgs16", "Sim1", "Nr5a1")
genes_func_upr_isr_high <- c("Atf4", "Ddit3", "Atf3", "Hspa5") # Ddit3 is also known as Chop, Hspa5 is also known as BiP
genes_func_autophagy_high <- c("Map1lc3b", "Becn1", "Atg5", "Atg7") # Map1lc3b is also known as LC3B
genes_func_autophagy_low <- "Sqstm1" # Sqstm1 is also known as p62
genes_func_stress_high <- c("Hmox1", "Nfe2l2", "Tnf", "Il1b", "Il6", "Serpina3n", "C1qa", "Fos", "Ddit4") # Nfe2l2 is also known as Nrf2
genes_func_apoptosis_high <- c("Casp3", "Casp9", "Bax", "Trp53", "H2ax") # Trp53 is also known as p53, H2afx replaced with H2ax
genes_func_apoptosis_low <- "Bcl2"
genes_func_mitochondrial_low <- c("Mfn1", "Mfn2", "Opa1", "Ppargc1a", "Tfam") # Dnm1l is also known as Drp1, Ucp4 replaced with Slc25a27
genes_func_mitochondrial_high <- c("Dnm1l", "Ucp2", "Slc25a27") # Dnm1l is also known as Drp1, Ucp4 replaced with Slc25a27
genes_func_prolif_low <- c("Mki67", "Pcna", "Ccnd1")
genes_func_prolif_high <- "Cdkn1a" # Cdkn1a is also known as p21
genes_func_synaptic_low <- c("Tnip1", "Nefm", "Nefh", "Syn1", "Dlg4", "Shank3", "Syp", "Snap25", "Nrxn1", "Nlgn1", "Nlgn3", "Gria1") # Dlg4 is also known as PSD-95, Nef replaced with Tnip1, Nrnx1 replaced with Nrxn1
genes_func_leptin_low <- c("Socs3", "Cish")
genes_func_leptin_high <- c("Foxo1", "Ptpn1")
genes_func_leptin_unknown <- c("Jak2", "Stat3", "Ptpn2", "Sh2b1", "Mc4r")

marker_gene_dict <- list(
  "cell_type" = list(
    "Astrocytes" = genes_astrocytes,
    "Astrocytes_high" = genes_astrocytes_high,
    "Microglia" = genes_microglia,
    "Microglia_high" = genes_microglia_high,
    "Microglia_low" = genes_microglia_low,
    "Oligodendrocytes" = genes_oligodendrocytes,
    "Oligodendrocytes_low" = genes_oligodendrocytes_low,
    "Neurons" = genes_neurons,
    "Neurons_low" = genes_neurons_low,
    "Neuroendocrine" = genes_neuroendocrine,
    "Neuroendocrine_low" = genes_neuroendocrine_low,
    "Neuroendocrine_high" = genes_neuroendocrine_high,
    "Neuroendocrine_unknown" = genes_neuroendocrine_unknown,
    "Hypothalamus" = genes_hypothalamus,
    "Hypothalamus_low" = genes_hypothalamus_low,
    "Hypothalamus_high" = genes_hypothalamus_high,
    "Hypothalamus_unknown" = genes_hypothalamus_unknown
  ),
  "functional" = list(
    "UPR_ISR" = genes_func_upr_isr,
    "UPR_ISR_high" = genes_func_upr_isr_high,
    "Autophagy" = genes_func_autophagy,
    "Autophagy_high" = genes_func_autophagy_high,
    "Autophagy_low" = genes_func_autophagy_low,
    "Stress" = genes_func_stress,
    "Stress_high" = genes_func_stress_high,
    "Apoptosis" = genes_func_apoptosis,
    "Apoptosis_high" = genes_func_apoptosis_high,
    "Apoptosis_low" = genes_func_apoptosis_low,
    "Mitochondrial" = genes_func_mitochondrial,
    "Mitochondrial_high" = genes_func_mitochondrial_high,
    "Mitochondrial_low" = genes_func_mitochondrial_low,
    "Proliferation" = genes_func_prolif,
    "Proliferation_high" = genes_func_prolif_high,
    "Proliferation_low" = genes_func_prolif_low,
    "Synaptic" = genes_func_synaptic,
    "Synaptic_low" = genes_func_synaptic_low,
    "Leptin" = genes_func_leptin,
    "Leptin_high" = genes_func_leptin_high,
    "Leptin_low" = genes_func_leptin_low,
    "Leptin_unknown" = genes_func_leptin_unknown
  )
)


# ### Specific neuron clusters ###
# subset_agrp_neurons <- subset(seurat_object, subset = predicted %in% c(agrp_cluster_id))
# DimPlot(subset_agrp_neurons, group.by = "orig.ident", label = TRUE)

# subset_pomc_neurons <- subset(seurat_object, subset = predicted %in% c(pomc_cluster_id))
# DimPlot(subset_pomc_neurons, group.by = "orig.ident", label = TRUE)

Idents(seurat_object) <- "predicted_renamed"


pdf(paste0(de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"), "FeaturePlots_agrp_pomc_neurons_four_markers_whole_dataset.pdf"),
  width = 5, height = 5
)
plot(FeaturePlot(seurat_object, features = "Lepr", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Lepr", label = TRUE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Agrp", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Agrp", label = TRUE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Npy", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Npy", label = TRUE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Pomc", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object, features = "Pomc", label = TRUE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()


svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    "FeaturePlots_agrp_pomc_neurons_four_markers_whole_dataset_Lepr_no_label.svg"
  ),
  width = 5, height = 5
)
plot(FeaturePlot(seurat_object, features = "Lepr", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Lepr", label = TRUE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    "FeaturePlots_agrp_pomc_neurons_four_markers_whole_dataset_Agrp_no_label.svg"
  ),
  width = 5, height = 5
)
plot(FeaturePlot(seurat_object, features = "Agrp", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Agrp", label = TRUE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    "FeaturePlots_agrp_pomc_neurons_four_markers_whole_dataset_Npy_no_label.svg"
  ),
  width = 5, height = 5
)
plot(FeaturePlot(seurat_object, features = "Npy", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Npy", label = TRUE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    "FeaturePlots_agrp_pomc_neurons_four_markers_whole_dataset_Pomc_no_label.svg"
  ),
  width = 5, height = 5
)
plot(FeaturePlot(seurat_object, features = "Pomc", label = FALSE) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
# plot(FeaturePlot(seurat_object, features = "Pomc", label = TRUE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
dev.off()

pdf(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_whole_dataset.pdf"))
  ),
  width = 10, height = 5
)
plot(FeaturePlot(seurat_object,
  features = "Lepr", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Lepr", label = TRUE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Agrp", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Agrp", label = TRUE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Npy", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Npy", label = TRUE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Pomc", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(seurat_object,
  features = "Pomc", label = TRUE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_whole_dataset_Lepr_no_label.svg"))
  ),
  width = 10, height = 5
)
plot(FeaturePlot(seurat_object,
  features = "Lepr", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Lepr", label = TRUE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_whole_dataset_Agrp_no_label.svg"))
  ),
  width = 10, height = 5
)
plot(FeaturePlot(seurat_object,
  features = "Agrp", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Agrp", label = TRUE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_whole_dataset_Npy_no_label.svg"))
  ),
  width = 10, height = 5
)
plot(FeaturePlot(seurat_object,
  features = "Npy", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Npy", label = TRUE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
svg(
  paste0(
    de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
    paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_whole_dataset_Pomc_no_label.svg"))
  ),
  width = 10, height = 5
)
plot(FeaturePlot(seurat_object,
  features = "Pomc", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
# plot(FeaturePlot(seurat_object, features = "Pomc", label = TRUE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))


pdf(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
  paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_subclusters.pdf"))
), width = 10, height = 5)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Lepr", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Agrp", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Npy", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
plot(FeaturePlot(subset(seurat_object, subset = predicted == pomc_cluster_id),
  features = "Pomc", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()


svg(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
  paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_subclusters_Lepr.svg"))
), width = 10, height = 5)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Lepr", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
svg(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
  paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_subclusters_Agrp.svg"))
), width = 10, height = 5)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Agrp", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
svg(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
  paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_subclusters_Npy.svg"))
), width = 10, height = 5)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Npy", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()
svg(paste0(
  de_dir, paste0(paste0("de_markers_", tolower(control_label), "_vs_"), tolower(treatment_label), "/other_explicit_de/plots/"),
  paste0("FeaturePlots_agrp_pomc_neurons_four_markers_", treatment_label, paste0("_vs_", control_label, "_subclusters_Pomc.svg"))
), width = 10, height = 5)
plot(FeaturePlot(subset(seurat_object, subset = predicted == pomc_cluster_id),
  features = "Pomc", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

### -> Enlarge the Agrp cluster view and display only UMAP1 and UMAP2 axes.
### -> Histo-Schnitte und Färbungen?
### -> konserviert in Menschen? transfer learning, translationale Modelle...
### -> Consider classification models for the treatment groups.
### -> CyberSort für bulk RNA-Seq?


### Marker gene plots
# DotPlots for cell type-specific markers
for (cells in names(marker_gene_dict$cell_type)) {
  pdf(paste0(plots_dir, paste0("cell_type/DotPlot/", dataset_id, "_renamed_cell_type_"), cells, "_DotPlot.pdf"), width = 8, height = 15)
  plot(DotPlot(seurat_object, features = marker_gene_dict$cell_type[[cells]], group.by = "predicted_renamed") + ggtitle(paste0("Cell type: ", cells)))
  plot(DotPlot(seurat_object,
    features = marker_gene_dict$cell_type[[cells]], idents = c("GLU_Neurons", "GABA_Neurons"),
    group.by = "predicted"
  ) + ggtitle(paste0("Cell type: ", cells, " (in more defined neuron clusters)")))
  dev.off()
}

# DotPlots for functional markers
for (functional_marker in names(marker_gene_dict$functional)) {
  pdf(paste0(plots_dir, paste0("functional/DotPlot/", dataset_id, "_renamed_functional_"), functional_marker, "_DotPlot.pdf"), width = 8, height = 15)
  plot(DotPlot(seurat_object, features = marker_gene_dict$functional[[functional_marker]], group.by = "predicted_renamed") +
    ggtitle(paste0("Functional: ", functional_marker)))
  plot(DotPlot(seurat_object,
    features = marker_gene_dict$functional[[functional_marker]], idents = c("GLU_Neurons", "GABA_Neurons"),
    group.by = "predicted"
  ) + ggtitle(paste0("Functional: ", functional_marker, " (in more defined neuron clusters)")))
  plot(DotPlot(seurat_object,
    features = marker_gene_dict$functional[[functional_marker]], group.by = "predicted_renamed",
    split.by = "orig.ident", cols = c("darkgreen", "blue")
  ) + ggtitle(paste0("Functional: ", functional_marker, " (split by sample)")))
  plot(DotPlot(seurat_object,
    features = marker_gene_dict$functional[[functional_marker]], idents = c("GLU_Neurons", "GABA_Neurons"),
    group.by = "predicted", split.by = "orig.ident", cols = c("darkgreen", "blue")
  ) + ggtitle(
    paste0("Functional: ", functional_marker, " (in more defined neuron clusters, split by sample)")
  ))
  dev.off()
}


# VlnPlots for cell type-specific markers
for (cells in names(marker_gene_dict$cell_type)) {
  pdf(paste0(plots_dir, paste0("cell_type/VlnPlot/", dataset_id, "_renamed_cell_type_"), cells, "_VlnPlot.pdf"), width = 10, height = 5)
  for (gene in marker_gene_dict$cell_type[[cells]]) {
    plot(VlnPlot(seurat_object, features = gene, group.by = "predicted_renamed", split.by = "orig.ident"))
  }
  dev.off()
}

# VlnPlots for functional markers
for (functional_marker in names(marker_gene_dict$functional)) {
  pdf(paste0(plots_dir, paste0("functional/VlnPlot/", dataset_id, "_renamed_functional_"), functional_marker, "_VlnPlot.pdf"), width = 10, height = 5)
  for (gene in marker_gene_dict$functional[[functional_marker]]) {
    plot(VlnPlot(seurat_object, features = gene, group.by = "predicted_renamed", split.by = "orig.ident"))
  }
  dev.off()
}


# FeaturePlots for cell type-specific markers
for (cells in names(marker_gene_dict$cell_type)) {
  pdf(paste0(plots_dir, paste0("cell_type/FeaturePlot/", dataset_id, "_renamed_cell_type_"), cells, "_FeaturePlot.pdf"), width = 10, height = 5)
  for (gene in marker_gene_dict$cell_type[[cells]]) {
    plot(FeaturePlot(seurat_object, features = gene, label = TRUE))
    plot(FeaturePlot(seurat_object, features = gene, label = TRUE, split.by = "orig.ident"))
  }
  dev.off()
}

# FeaturePlots for functional markers
for (functional_marker in names(marker_gene_dict$functional)) {
  pdf(paste0(plots_dir, paste0("functional/FeaturePlot/", dataset_id, "_renamed_functional_"), functional_marker, "_FeaturePlot.pdf"), width = 10, height = 5)
  for (gene in marker_gene_dict$functional[[functional_marker]]) {
    plot(FeaturePlot(seurat_object, features = gene, label = TRUE))
    plot(FeaturePlot(seurat_object, features = gene, label = TRUE, split.by = "orig.ident"))
  }
  dev.off()
}


marker_gene_dict <- list(
  "cell_type" = list(
    "Astrocytes" = genes_astrocytes,
    "Microglia" = genes_microglia,
    "Oligodendrocytes" = genes_oligodendrocytes,
    "Neurons" = genes_neurons,
    "Neuroendocrine" = genes_neuroendocrine,
    "Hypothalamus" = genes_hypothalamus
  ),
  "functional" = list(
    "UPR_ISR" = genes_func_upr_isr,
    "Autophagy" = genes_func_autophagy,
    "Stress" = genes_func_stress,
    "Apoptosis" = genes_func_apoptosis,
    "Mitochondrial" = genes_func_mitochondrial,
    "Proliferation" = genes_func_prolif,
    "Synaptic" = genes_func_synaptic,
    "Leptin" = genes_func_leptin
  )
)

## AUCell analysis ##
library(AUCell)

counts_AUCell <- GetAssayData(seurat_object, slot = "counts")

cell_rankings_AUCell <- AUCell_buildRankings(counts_AUCell)

cells_AUCell_unfolded_prot_response <- AUCell_calcAUC(marker_gene_dict$functional[["UPR_ISR"]], cell_rankings_AUCell)

cells_assignment_unfolded_prot_response <- AUCell_exploreThresholds(cells_AUCell_unfolded_prot_response, plotHist = TRUE, assign = TRUE)

new_cells_unfolded_prot_response <- names(which(getAUC(cells_AUCell_unfolded_prot_response)["geneSet", ] > 0.15))


seurat_object$unfolded_prot_response <- ifelse(colnames(seurat_object) %in% new_cells_unfolded_prot_response, "unfolded_prot_response", "not_unfolded_prot_response")

DimPlot(seurat_object, group.by = "unfolded_prot_response", split.by = "orig.ident", label = TRUE)

head(seurat_object@meta.data)

head(cells_AUCell_unfolded_prot_response)


FeaturePlot(seurat_object, features = "Ddit3", label = TRUE)


## AddModuleScore ##
for (functional_marker in names(marker_gene_dict$functional)) {
  seurat_object <- AddModuleScore(seurat_object,
    features = marker_gene_dict$functional[[functional_marker]],
    name = paste0(functional_marker, "_addmodulescore")
  )
}

for (functional_marker in names(marker_gene_dict$functional)) {
  pdf(
    paste0(
      module_scores_dir,
      functional_marker, "_addmodulescore_FeaturePlot.pdf"
    ),
    width = 10, height = 5
  )
  plot(FeaturePlot(seurat_object,
    features = paste0(functional_marker, "_addmodulescore1"), label = TRUE,
    cols = c("lightgrey", "lightgrey", "red")
  ))
  plot(FeaturePlot(seurat_object,
    features = paste0(functional_marker, "_addmodulescore1"), label = TRUE,
    split.by = "orig.ident", cols = c("lightgrey", "lightgrey", "red")
  ))
  dev.off()
}

for (functional_marker in names(marker_gene_dict$functional)) {
  pdf(
    paste0(
      module_scores_dir,
      functional_marker, "_addmodulescore_VlnPlot.pdf"
    ),
    width = 10, height = 5
  )
  plot(VlnPlot(seurat_object,
    features = paste0(functional_marker, "_addmodulescore1"), group.by = "predicted_renamed",
    split.by = "orig.ident"
  ))
  plot(VlnPlot(seurat_object,
    features = paste0(functional_marker, "_addmodulescore1"),
    idents = c("GLU_Neurons", "GABA_Neurons"), group.by = "predicted", split.by = "orig.ident"
  ))
  dev.off()
}
