# Seurat preprocessing, integration, and doublet detection
# Set PROJECT_DIR, SAMPLE_ID, and DATASET_ID before running. The merged renamed
# Seurat object expected below must be created by the preceding integration step.

# Configuration
project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/project")
input_dir <- file.path(project_dir, "data")
object_dir <- file.path(project_dir, "results", "seurat_objects")
figure_dir <- file.path(project_dir, "results", "figures")
de_dir <- file.path(project_dir, "results", "differential_expression")
sample_id <- Sys.getenv("SAMPLE_ID", unset = "sample")
dataset_id <- Sys.getenv("DATASET_ID", unset = "dataset")
cluster_map_file <- file.path(project_dir, "metadata", "cluster_renaming.csv")

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
library(scDblFinder)

# Read the sample-level Seurat object
seurat_object <- readRDS(file.path(input_dir, paste0(sample_id, "_raw_data.rds")))

# Map the original cluster labels to publication-ready cell-type names
file_renaming <- read.csv(cluster_map_file)

head(file_renaming)

seurat_object@meta.data$predicted_renamed <- file_renaming$Clusters[match(seurat_object@meta.data$predicted, file_renaming$Var1)]

head(seurat_object@meta.data)

Idents(seurat_object) <- "predicted_renamed"
head(Idents(seurat_object))

pdf(file.path(figure_dir, paste0(sample_id, "_renamed_UMAP.pdf")), width = 10, height = 5)
plot(DimPlot(seurat_object, group.by = "predicted_renamed", label = TRUE))
plot(DimPlot(seurat_object, group.by = "orig.ident", label = TRUE, alpha = 0.3))
plot(DimPlot(seurat_object, group.by = "predicted_renamed", split.by = "orig.ident", label = TRUE))
dev.off()

head(seurat_object@meta.data)
seurat_object$celltype_condition <- paste(seurat_object$orig.ident, seurat_object$predicted_renamed, sep = "_")
Idents(seurat_object) <- "celltype_condition"
head(Idents(seurat_object))
table(Idents(seurat_object))

# Save the renamed sample-level object if needed
# saveRDS(seurat_object, file.path(object_dir, paste0(sample_id, "_renamed.rds")))


# Load the merged, renamed object for integration
seurat_object <- readRDS(file.path(
  object_dir,
  paste0(dataset_id, "_sc_data_seurat_object_hypomap_run_renamed.rds")
))


### See whether integration of the data is necessary
seurat_object_integration <- seurat_object

seurat_object_integration <- NormalizeData(seurat_object_integration)
seurat_object_integration <- FindVariableFeatures(seurat_object_integration)
seurat_object_integration <- ScaleData(seurat_object_integration)
seurat_object_integration <- RunPCA(seurat_object_integration)

seurat_object_integration <- FindNeighbors(seurat_object_integration, dims = 1:30, reduction = "pca")
seurat_object_integration <- FindClusters(seurat_object_integration, resolution = 0.8, cluster.name = "unintegrated_clusters")
seurat_object_integration <- RunUMAP(seurat_object_integration, dims = 1:30, reduction = "pca", reduction.name = "umap_unintegrated", reduction.key = "ccaUMAP_")

DimPlot(seurat_object_integration, reduction = "pca", label = TRUE, group.by = "orig.ident")

DimPlot(seurat_object_integration,
  reduction = "umap_unintegrated", label = TRUE,
  group.by = c("orig.ident", "predicted_renamed", "unintegrated_clusters"), ncol = 1
)

# -> integration with CCA would be necessary

# Integration with CCA
seurat_object_integration <- seurat_object

seurat_object_integration[["RNA"]] <- split(seurat_object_integration[["RNA"]], f = seurat_object_integration$orig.ident)
seurat_object_integration

seurat_object_integration <- NormalizeData(seurat_object_integration)
seurat_object_integration <- FindVariableFeatures(seurat_object_integration)
seurat_object_integration <- ScaleData(seurat_object_integration)
seurat_object_integration <- RunPCA(seurat_object_integration)

seurat_object_integration <- IntegrateLayers(
  object = seurat_object_integration, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

seurat_object_integration <- JoinLayers(seurat_object_integration)

seurat_object_integration <- FindNeighbors(seurat_object_integration, dims = 1:30, reduction = "integrated.cca")
seurat_object_integration <- FindClusters(seurat_object_integration, resolution = 0.8, cluster.name = "integrated_ccs_clusters")

seurat_object_integration <- RunUMAP(seurat_object_integration,
  dims = 1:30, reduction = "integrated.cca", reduction.name = "umap_cca",
  reduction.key = "ccaUMAP_"
)

seurat_object_integration

pdf(file.path(figure_dir, paste0(dataset_id, "_renamed_UMAP_integration_cca.pdf")), width = 10, height = 5)
plot(DimPlot(seurat_object_integration, reduction = "integrated.cca", label = TRUE, group.by = "orig.ident"))
plot(DimPlot(seurat_object_integration, reduction = "umap_cca", label = TRUE, group.by = "orig.ident", ncol = 1))
plot(DimPlot(seurat_object_integration, reduction = "umap_cca", label = TRUE, group.by = "predicted_renamed", ncol = 1))
plot(DimPlot(seurat_object_integration, reduction = "umap_cca", label = TRUE, group.by = "integrated_ccs_clusters", ncol = 1))
dev.off()

seurat_object_integration

# Save the integrated Seurat object
saveRDS(
  seurat_object_integration,
  file.path(object_dir, paste0(dataset_id, "_integrated_cca.rds"))
)


### Doublet detection ###
seurat_object_dbl <- as.SingleCellExperiment(seurat_object_integration)

seurat_object_dbl <- scDblFinder(seurat_object_dbl, samples = "orig.ident")
head(seurat_object_dbl)
table(seurat_object_dbl$scDblFinder.class)

seurat_object_dbl@colData@listData %>%
  as.data.frame() %>%
  head()

meta_scdblfinder <- seurat_object_dbl@colData@listData %>%
  as.data.frame() %>%
  dplyr::select(starts_with("scDblFinder")) # 'scDblFinder.class')
head(meta_scdblfinder)

rownames(meta_scdblfinder) <- seurat_object_dbl@colData@rownames
head(meta_scdblfinder)
seurat_object <- AddMetaData(
  object = seurat_object_integration,
  metadata = meta_scdblfinder %>% dplyr::select("scDblFinder.class")
)
saveRDS(seurat_object, file.path(object_dir, paste0(dataset_id, "_integrated_cca.rds")))

head(seurat_object@meta.data)

table(seurat_object$scDblFinder.class)
rm(list = c("meta_scdblfinder", "seurat_object_dbl"))

pdf(file.path(de_dir, paste0(dataset_id, "_doublet_detection_violin_plot.pdf")))
plot(VlnPlot(seurat_object,
  group.by = "orig.ident", split.by = "scDblFinder.class",
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2, pt.size = 0
) + theme(legend.position = "right"))
dev.off()


# Table of cell types and their counts in each sample orig.ident
table(seurat_object$predicted_renamed, seurat_object$orig.ident)

# write table to csv
write.csv(
  table(seurat_object$predicted_renamed, seurat_object$scDblFinder.class, seurat_object$orig.ident),
  file.path(de_dir, "table_cell_types_doublet_finder.csv")
)


seurat_object_no_dbl <- subset(seurat_object, subset = scDblFinder.class == "singlet")

SaveSeuratRds(
  seurat_object_no_dbl,
  file.path(object_dir, paste0(dataset_id, "_integrated_cca_singlets.rds"))
)
