# Publication figure generation
# Set PROJECT_DIR, DATASET_ID, CONTROL_LABEL, and TREATMENT_LABEL before running.
# The neuron cluster IDs and marker panels below are study-specific inputs.

# Configuration
project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/project")
dataset_id <- Sys.getenv("DATASET_ID", unset = "dataset")
control_label <- Sys.getenv("CONTROL_LABEL", unset = "control")
treatment_label <- Sys.getenv("TREATMENT_LABEL", unset = "treatment")
agrp_cluster_id <- Sys.getenv("AGRP_CLUSTER_ID", unset = "your_agrp_cluster")
pomc_cluster_id <- Sys.getenv("POMC_CLUSTER_ID", unset = "your_pomc_cluster")
pmch_cluster_id <- Sys.getenv("PMCH_CLUSTER_ID", unset = "your_pmch_cluster")
rfx4_cluster_id <- Sys.getenv("RFX4_CLUSTER_ID", unset = "your_rfx4_cluster")
tac2_cluster_id <- Sys.getenv("TAC2_CLUSTER_ID", unset = "your_tac2_cluster")
object_dir <- file.path(project_dir, "results", "seurat_objects")
figure_dir <- file.path(project_dir, "results", "figures", "")

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
library(ggpubr)


# Load the singlet-only object and the full object for doublet summaries
seurat_object <- readRDS(file.path(object_dir, paste0(dataset_id, "_integrated_cca_singlets.rds")))
seurat_object_doublets <- readRDS(file.path(object_dir, paste0(dataset_id, "_integrated_cca.rds")))

table(seurat_object_doublets$scDblFinder.class)

# really only select singlets here
seurat_object <- subset(seurat_object, subset = scDblFinder.class == "singlet")

head(seurat_object@meta.data)
seurat_object

# create a new column in the meta.data with the specific neurons predicted names (predicted) and assign NA to all other cells
seurat_object$specific_neurons <- seurat_object$predicted
seurat_object$specific_neurons[!seurat_object$predicted_renamed %in% c("GLU_Neurons", "GABA_Neurons")] <- NA

seurat_object$rename_labels <- seurat_object$predicted_renamed
seurat_object$rename_labels[seurat_object$predicted_renamed == "GLU_Neurons"] <- "GLU-Neurons"
seurat_object$rename_labels[seurat_object$predicted_renamed == "GABA_Neurons"] <- "GABA-Neurons"

head(seurat_object@meta.data)
table(seurat_object$scDblFinder.class)


### Tables with cell numbers per cluster
write.csv(table(seurat_object$rename_labels), file = paste0(
  figure_dir,
  "Fig1_UMAP_clusters_all_rename_labels.csv"
), row.names = FALSE)

head(seurat_object@meta.data)
write.csv(table(seurat_object$predicted), file = paste0(
  figure_dir,
  "Fig1_UMAP_clusters_all_specific_labels.csv"
), row.names = FALSE)

### Figure 1 ###
# UMAP of all clusters and specific neurons
# pdf(paste0(figure_dir,
#     "Fig1_UMAP_clusters_all.pdf"), width = 15, height = 15)
# p1 <- DimPlot(seurat_object, reduction = "umap_scvi", group.by = "predicted_renamed", label = FALSE) + NoLegend() + ggtitle(NULL) + NoAxes()
# plot(LabelClusters(p1, id = "predicted_renamed", repel = TRUE, max.overlaps = 100, size = 7, face = "bold"))
# p2 <- DimPlot(seurat_object, reduction = "umap_scvi", group.by = "specific_neurons", label = FALSE) + NoLegend() + ggtitle(NULL) + NoAxes()
# plot(LabelClusters(p2, id = "specific_neurons", repel = TRUE, max.overlaps = 100))
# dev.off()

tiff(paste0(
  figure_dir,
  "Fig1_UMAP_clusters_all_big_clusters.tiff"
), units = "cm", width = 35, height = 35, res = 600)
p1 <- DimPlot(seurat_object, reduction = "umap_scvi", group.by = "rename_labels", label = FALSE) + NoLegend() + ggtitle(NULL) + NoAxes()
plot(LabelClusters(p1, id = "rename_labels", repel = TRUE, max.overlaps = 100, size = 7, face = "bold"))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig1_UMAP_clusters_all_specific_neurons.tiff"
), units = "cm", width = 35, height = 35, res = 600)
p2 <- DimPlot(seurat_object, reduction = "umap_scvi", group.by = "specific_neurons", label = FALSE) + NoLegend() + ggtitle(NULL) + NoAxes()
plot(LabelClusters(p2, id = "specific_neurons", repel = TRUE, max.overlaps = 100))
dev.off()

### Figure 2 ###
# Hypothalamus-specific neurons and cell type-specific genes (= positive controls)

# cell type-specific genes
seurat_object$predicted_renamed <- factor(seurat_object$predicted_renamed,
  levels = c("Astrocytes", "Ependymal", "GABA_Neurons", "GLU_Neurons", "Immune", "Microglia", "Oligodendrocytes", "OPC", "ParsTuber", "Vascular")
)

# pdf(paste0(figure_dir,
#     "Fig2_DotPlot_cell_types_positive_controls.pdf"), width = 10, height = 10)
# p1 <- DotPlot(subset(seurat_object, subset = predicted_renamed %in% c("Astrocytes", "GABA_Neurons", "GLU_Neurons", "Microglia", "Oligodendrocytes", "OPC")),
#     features = c("Aqp4", "Rbfox3", "Cx3cr1", "Mbp", "Olig1"), group.by = "predicted_renamed") + xlab(NULL) + ylab(NULL)
# plot(p1)
# dev.off()

# cell type-specific genes in reverse order
seurat_object$predicted_renamed <- factor(seurat_object$predicted_renamed,
  levels = rev(c("Astrocytes", "Ependymal", "GABA_Neurons", "GLU_Neurons", "Immune", "Microglia", "Oligodendrocytes", "OPC", "ParsTuber", "Vascular"))
)

# pdf(paste0(figure_dir,
#     "Fig2_DotPlot_cell_types_positive_controls_rev.pdf"), width = 10, height = 10)
# p1 <- DotPlot(subset(seurat_object, subset = predicted_renamed %in% c("Astrocytes", "GABA_Neurons", "GLU_Neurons", "Microglia", "Oligodendrocytes", "OPC")),
#     features = rev(c("Aqp4", "Rbfox3", "Cx3cr1", "Mbp", "Olig1")), group.by = "predicted_renamed") + xlab(NULL) + ylab(NULL)
# plot(p1)
# dev.off()

# cell type-specific genes in reverse order
seurat_object$predicted_renamed <- factor(seurat_object$predicted_renamed,
  levels = rev(c("GABA_Neurons", "GLU_Neurons", "Oligodendrocytes", "Astrocytes", "OPC", "Microglia", "Ependymal", "Immune", "ParsTuber", "Vascular"))
)

tiff(paste0(
  figure_dir,
  "Fig2_DotPlot_cell_types_positive_controls_rev.tiff"
), units = "cm", width = 18, height = 10, res = 600)
p1 <- DotPlot(subset(seurat_object, subset = predicted_renamed %in% c("GABA_Neurons", "GLU_Neurons", "Oligodendrocytes", "Astrocytes", "OPC", "Microglia")),
  features = rev(c("Rbfox3", "Mbp", "Aqp4", "Olig1", "Cx3cr1")), group.by = "predicted_renamed"
) + xlab(NULL) + ylab(NULL)
plot(p1)
dev.off()

# hypothalamus-specific neurons
hypo_seurat_object <- subset(seurat_object, subset = predicted %in% c(
  agrp_cluster_id, pomc_cluster_id, pmch_cluster_id, rfx4_cluster_id,
  tac2_cluster_id
))
hypo_seurat_object$predicted <- factor(hypo_seurat_object$predicted,
  levels = rev(c(rfx4_cluster_id, pomc_cluster_id, tac2_cluster_id, pmch_cluster_id, agrp_cluster_id))
)

# pdf(paste0(figure_dir,
#     "Fig2_DotPlot_hypothalamus_specific_neurons.pdf"), width = 10, height = 10)
# p1 <- DotPlot(hypo_seurat_object, features = rev(c("Agrp", "Npy", "Lepr", "Pomc", "Cartpt", "Hcrt", "Esr1")), group.by = "predicted") + xlab(NULL) + ylab(NULL)
# plot(p1)
# dev.off()

tiff(paste0(
  figure_dir,
  "Fig2_DotPlot_hypothalamus_specific_neurons.tiff"
), units = "cm", width = 20, height = 10, res = 600)
p1 <- DotPlot(hypo_seurat_object, features = rev(c("Npy", "Lepr", "Agrp", "Cartpt", "Esr1", "Pomc", "Hcrt")), group.by = "predicted") + xlab(NULL) + ylab(NULL)
plot(p1)
dev.off()

tiff(paste0(
  figure_dir,
  "Fig2_DotPlot_hypothalamus_specific_neurons_grid.tiff"
), units = "cm", width = 20, height = 10, res = 600)
p1 <- DotPlot(hypo_seurat_object, features = rev(c("Npy", "Lepr", "Agrp", "Cartpt", "Esr1", "Pomc", "Hcrt")), group.by = "predicted") +
  xlab(NULL) + ylab(NULL) + theme(panel.grid = element_line(
    color = "lightgrey",
    size = 0.75,
    linetype = 2
  ))
plot(p1)
dev.off()

# hypothalamus-specific neurons in reverse order
# hypo_seurat_object <- subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id, rfx4_cluster_id,
#         tac2_cluster_id))
# hypo_seurat_object$predicted <- factor(hypo_seurat_object$predicted,
#     levels = c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id, rfx4_cluster_id, tac2_cluster_id))

# pdf(paste0(figure_dir,
#     "Fig2_DotPlot_hypothalamus_specific_neurons_rev.pdf"), width = 10, height = 10)
# p1 <- DotPlot(hypo_seurat_object, features = c("Agrp", "Npy", "Lepr", "Pomc", "Cartpt", "Hcrt", "Esr1"), group.by = "predicted") + xlab(NULL) + ylab(NULL)
# plot(p1)
# dev.off()

tiff(paste0(
  figure_dir,
  "Fig2_DotPlot_hypothalamus_specific_neurons_rev.tiff"
), units = "cm", width = 20, height = 10, res = 600)
p1 <- DotPlot(hypo_seurat_object, features = c("Agrp", "Lepr", "Npy", "Cartpt", "Esr1", "Pomc", "Hcrt"), group.by = "predicted") + xlab(NULL) + ylab(NULL)
plot(p1)
dev.off()

### Figure 3 ###
# Controls for SIH model: food regulation (SIH vs. Ctrl)
Idents(seurat_object) <- "predicted_renamed"

# over all clusters and conditions together
# pdf(paste0(figure_dir,
#     "Fig3_FeaturePlots_controls_sih_model_all_clusters.pdf"), width = 10, height = 10)
# plot(FeaturePlot(seurat_object, features = "Lepr", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(seurat_object, features = "Pomc", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(seurat_object, features = "Npy", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(seurat_object, features = "Agrp", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(seurat_object, features = "Cartpt", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# dev.off()

# over all clusters and conditions together
tiff(paste0(
  figure_dir,
  "Fig3_FeaturePlots_controls_sih_model_all_clusters_Lepr.tiff"
), units = "cm", width = 20, height = 20, res = 600)
plot(FeaturePlot(seurat_object, features = "Lepr", label = FALSE, label.size = 10) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2), plot.title = element_text(size = 32),
  legend.key.size = unit(1.15, "cm"), legend.text = element_text(size = 20)
))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig3_FeaturePlots_controls_sih_model_all_clusters_Pomc.tiff"
), units = "cm", width = 20, height = 20, res = 600)
plot(FeaturePlot(seurat_object, features = "Pomc", label = FALSE, label.size = 10) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2), plot.title = element_text(size = 32),
  legend.key.size = unit(1.15, "cm"), legend.text = element_text(size = 20)
))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig3_FeaturePlots_controls_sih_model_all_clusters_Npy.tiff"
), units = "cm", width = 20, height = 20, res = 600)
plot(FeaturePlot(seurat_object, features = "Npy", label = FALSE, label.size = 10) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2), plot.title = element_text(size = 32),
  legend.key.size = unit(1.15, "cm"), legend.text = element_text(size = 20)
))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig3_FeaturePlots_controls_sih_model_all_clusters_Agrp.tiff"
), units = "cm", width = 20, height = 20, res = 600)
plot(FeaturePlot(seurat_object, features = "Agrp", label = FALSE, label.size = 10) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2), plot.title = element_text(size = 32),
  legend.key.size = unit(1.15, "cm"), legend.text = element_text(size = 20)
))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig3_FeaturePlots_controls_sih_model_all_clusters_Cartpt.tiff"
), units = "cm", width = 20, height = 20, res = 600)
plot(FeaturePlot(seurat_object, features = "Cartpt", label = FALSE, label.size = 10) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2), plot.title = element_text(size = 32),
  legend.key.size = unit(1.15, "cm"), legend.text = element_text(size = 20)
))
dev.off()

# specific neurons only and separated by condition
# pdf(paste0(figure_dir,
#     paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), ".pdf"))), width = 10, height = 10)
# plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id), features = "Lepr", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.line=element_blank(), axis.ticks=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id), features = "Npy", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id), features = "Agrp", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(subset(seurat_object, subset = predicted == pomc_cluster_id), features = "Pomc", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# plot(FeaturePlot(subset(seurat_object, subset = predicted == pmch_cluster_id), features = "Cartpt", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# dev.off()


# png(paste0(figure_dir,
#     paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Cartpt.png"))), width = 1000, height = 500)
# plot(FeaturePlot(subset(seurat_object, subset = predicted == pmch_cluster_id), features = "Cartpt", label = FALSE,
#     split.by = "orig.ident")&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# dev.off()

# png(paste0(figure_dir,
#     paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Cartpt_all.png"))), width = 1000, height = 1000)
# plot(FeaturePlot(seurat_object, features = "Cartpt", label = FALSE)&theme(axis.title=element_blank(), axis.text=element_blank(),
#     axis.ticks=element_blank(), axis.line=element_blank(), legend.position=c(0.05, 0.2)))
# dev.off()


tiff(paste0(
  figure_dir,
  paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Lepr.tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Lepr", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Npy.tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Npy", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp.tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(FeaturePlot(subset(seurat_object, subset = predicted == agrp_cluster_id),
  features = "Agrp", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Pomc.tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(FeaturePlot(subset(seurat_object, subset = predicted == pomc_cluster_id),
  features = "Pomc", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_FeaturePlots_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Cartpt.tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(FeaturePlot(subset(seurat_object, subset = predicted == pmch_cluster_id),
  features = "Cartpt", label = FALSE,
  split.by = "orig.ident"
) & theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), axis.line = element_blank(), legend.position = c(0.05, 0.2)
))
dev.off()


# DotPlot of specific neurons and separated by condition
# pdf(paste0(figure_dir,
#     paste0("Fig3_DotPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), ".pdf"))), width = 10, height = 10)
# plot(DotPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id)),
#     features = c("Agrp", "Npy", "Lepr", "Cartpt", "Pomc"), group.by = "predicted", split.by = "orig.ident") + xlab(NULL) + ylab(NULL))
# dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_DotPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), ".tiff"))
), units = "cm", width = 20, height = 10, res = 600)
plot(DotPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id)),
  features = c("Agrp", "Npy", "Lepr", "Cartpt", "Pomc"), group.by = "predicted", split.by = "orig.ident"
) + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_line(
    color = "lightgrey",
    size = 0.75,
    linetype = 2
  )))
dev.off()

tiff(paste0(
  figure_dir,
  "Fig3_DotPlot_test_Agrp.tiff"
), units = "cm", width = 20, height = 10, res = 600)
plot(DotPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id)),
  features = c("Agrp"), group.by = "predicted", split.by = "orig.ident"
) + xlab(NULL) + ylab(NULL))
dev.off()

DotPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id)),
  features = c("Agrp"), group.by = "predicted", split.by = "orig.ident", cols = c("blue", "red")
) + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_line(
    color = "lightgrey",
    size = 0.75,
    linetype = 2
  ))

# Factor for orig.ident
seurat_object$orig.ident <- factor(seurat_object$orig.ident,
  levels = c(control_label, treatment_label)
)

tail(seurat_object@meta.data)
seurat_object$Condition <- ifelse(seurat_object$orig.ident == control_label, "Control", "Starvation")
seurat_object$Condition <- factor(seurat_object$Condition,
  levels = c("Control", "Starvation")
)

# tiff(paste0(figure_dir,
#     paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_all.tiff"))), units = "cm", width = 10, height = 10, res = 600)
# plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id, pomc_cluster_id, pmch_cluster_id)),
#     features = c("Agrp", "Npy", "Lepr", "Cartpt", "Pomc"), group.by = "predicted", split.by = "Condition",
#     cols = c("lightgrey", "darkviolet")) & xlab(NULL) & theme(axis.text.x = element_text(angle = 45, hjust = 1)
#     ) & theme(legend.position = "right"))
# dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_all_Agrp.tiff"))
), units = "cm", width = 22, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Agrp", "Npy", "Lepr"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet"), y.max = 5.5
) & xlab(NULL) & theme(axis.text.x = element_text(angle = 0, hjust = 0.5))) # theme(axis.text.x = element_text(angle = 45, hjust = 1)))
# ) & theme(legend.position = "right"))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Agrp.tiff"))
), units = "cm", width = 7, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Agrp"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet"), y.max = 5.5
) + geom_segment(aes(x = 0.78, xend = 1.23, y = 5.5, yend = 5.5)) + geom_segment(aes(x = 0.78, xend = 0.78, y = 5.4, yend = 5.5)) + geom_segment(aes(x = 1.23, xend = 1.23, y = 5.4, yend = 5.5)) + geom_text(x = 1, y = 5.6, label = "***") & xlab(NULL) & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5),
  legend.position = "none"
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Npy.tiff"))
), units = "cm", width = 7, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Npy"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet"), y.max = 5.5
) + geom_segment(aes(x = 0.78, xend = 1.23, y = 5.5, yend = 5.5)) + geom_segment(aes(x = 0.78, xend = 0.78, y = 5.4, yend = 5.5)) + geom_segment(aes(x = 1.23, xend = 1.23, y = 5.4, yend = 5.5)) + geom_text(x = 1, y = 5.6, label = "***") & xlab(NULL) & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5),
  legend.position = "none"
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Lepr.tiff"))
), units = "cm", width = 7, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Lepr"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet"), y.max = 5.5
) + geom_segment(aes(x = 0.78, xend = 1.23, y = 5.5, yend = 5.5)) + geom_segment(aes(x = 0.78, xend = 0.78, y = 5.4, yend = 5.5)) + geom_segment(aes(x = 1.23, xend = 1.23, y = 5.4, yend = 5.5)) + geom_text(x = 1, y = 5.6, label = "***") & xlab(NULL) & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5),
  legend.position = "none"
))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_all_Pomc.tiff"))
), units = "cm", width = 10.5, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(pomc_cluster_id)),
  features = c("Pomc"), group.by = "predicted", split.by = "Condition", same.y.lims = TRUE,
  cols = c("lightgrey", "darkviolet"), y.max = 5.5
) + geom_segment(aes(x = 0.78, xend = 1.23, y = 5.5, yend = 5.5)) + geom_segment(aes(x = 0.78, xend = 0.78, y = 5.4, yend = 5.5)) + geom_segment(aes(x = 1.23, xend = 1.23, y = 5.4, yend = 5.5)) + geom_text(x = 1, y = 5.7, label = "n.s.") & xlab(NULL) & theme(axis.text.x = element_text(
  angle = 0, hjust = 0.5
)) & theme(legend.position = "right"))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Agrp.tiff"))
), units = "cm", width = 10, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Agrp"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet")
) & xlab(NULL) & theme(legend.position = "right") & theme(axis.text.x = element_blank()))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Npy.tiff"))
), units = "cm", width = 10, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Npy"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet")
) & xlab(NULL) & theme(legend.position = "right") & theme(axis.text.x = element_blank()))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Agrp_neurons_Lepr.tiff"))
), units = "cm", width = 10, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(agrp_cluster_id)),
  features = c("Lepr"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet")
) & xlab(NULL) & theme(legend.position = "right") & theme(axis.text.x = element_blank()))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Pomc_neurons_Pomc.tiff"))
), units = "cm", width = 10, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(pomc_cluster_id)),
  features = c("Pomc"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet")
) & xlab(NULL) & theme(legend.position = "right") & theme(axis.text.x = element_blank()))
dev.off()

tiff(paste0(
  figure_dir,
  paste0("Fig3_VlnPlot_controls_sih_model_specific_neurons_", tolower(treatment_label), paste0("_vs_", tolower(control_label), "_Pmch_neurons_Cartpt.tiff"))
), units = "cm", width = 10, height = 10, res = 600)
plot(VlnPlot(subset(seurat_object, subset = predicted %in% c(pmch_cluster_id)),
  features = c("Cartpt"), group.by = "predicted", split.by = "Condition",
  cols = c("lightgrey", "darkviolet")
) & xlab(NULL) & theme(legend.position = "right") & theme(axis.text.x = element_blank()))
dev.off()

### Cell counts by condition ###
