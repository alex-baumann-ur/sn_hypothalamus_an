# Import Cell Ranger outputs, map cells to the HypoMap reference, and save
# per-sample and merged Seurat objects.
#
# Configure PROJECT_DIR, SAMPLE_IDS, and DATASET_ID in the environment.
# Provide a metadata/sample_conditions.csv file with columns sample_id and
# condition, plus metadata/cluster_renaming.csv with columns Var1 and Clusters.
# Sample IDs must match the Cell Ranger output folder names.

# Configuration
project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/project")
cellranger_dir <- Sys.getenv(
  "CELLRANGER_DIR",
  unset = file.path(project_dir, "data", "cellranger")
)
dataset_id <- Sys.getenv("DATASET_ID", unset = "dataset")
sample_ids_prev <- trimws(strsplit(
  Sys.getenv("SAMPLE_IDS", unset = "sample_1,sample_2"),
  ",",
  fixed = TRUE
)[[1]])
sample_id <- "SAMPLE_ID_MERGED"
sample_metadata_file <- file.path(project_dir, "metadata", "sample_conditions.csv")
cluster_map_file <- file.path(project_dir, "metadata", "cluster_renaming.csv")

object_dir <- file.path(project_dir, "results", "seurat_objects")
raw_object_dir <- file.path(project_dir, "data")
predicted_cluster_dir <- file.path(project_dir, "results", "hypomap_predictions")
figure_dir <- file.path(project_dir, "results", "figures", "hypomap")

# Load required packages. Install mapscvi separately if it is not available
devtools::install_github("lsteuernagel/mapscvi")
library(mapscvi)
library(Seurat)

# Import the filtered feature-barcode matrix and save one Seurat object per sample
Read10xDataIn <- function(sample) {
  comb_dirs <- paste0(cellranger_dir, sample, "/outs/filtered_feature_bc_matrix")
  sc_data.data <- Read10X(data.dir = comb_dirs)
  seur_obj <- CreateSeuratObject(counts = sc_data.data)
  saveRDS(seur_obj, file = paste0(raw_object_dir, sample, "_sc_data_seurat_object.rds"))
  }

for (sample in sample_ids) {
    Read10xDataIn(sample)
}

# Read in raw SeuratObjects
crtl_sc_seurat <- readRDS(paste0(raw_object_dir, sample_ids_prev[1], "_sc_data_seurat_object.rds"))
starv_sc_seurat <- readRDS(paste0(raw_object_dir, sample_ids_prev[2], "_sc_data_seurat_object.rds"))

# Add project name to meta.data
crtl_sc_seurat@meta.data$orig.ident <- "CRTL"
starv_sc_seurat@meta.data$orig.ident <- "STARV"

# Merge SeuratObjects
merged_sc_seurat <- merge(
  crtl_sc_seurat, y = starv_sc_seurat, add.cell.ids = c("CRTL", "STARV"), 
  project = "crtl_starv")
merged_sc_seurat

head(merged_sc_seurat@meta.data)
table(merged_sc_seurat$orig.ident)

# Save merged SeuratObject
saveRDS(merged_sc_seurat, paste0(raw_object_dir, sample_id, "_seurat_object.rds"))


# Plots for mapscvi #
PlotsMapSCVI <- function(sample_name, sc_data) {
    # Normal sized plots
    pdf(paste0(figure_dir, sample_name, "_plots_mapscvi_hypomap.pdf"), width = 15, height = 10)
    plot(plot_query_labels(query_seura_object=sc_data,reference_seurat=mapscvi::reference_hypoMap_downsample,overlay = FALSE,
                      label_col = "C66_named",labelonplot = FALSE))
    plot(plot_query_labels(query_seura_object=sc_data,reference_seurat=mapscvi::reference_hypoMap_downsample,overlay = FALSE,
                      label_col = "C66_named",labelonplot = TRUE))
    plot(plot_query_labels(query_seura_object=sc_data,reference_seurat=mapscvi::reference_hypoMap_downsample,label_col="C66_named",
                      overlay = TRUE,query_pt_size = 0.4,labelonplot = FALSE,label.size=1))
    plot(FeaturePlot(sc_data, features = "prediction_probability")+NoAxes())
    dev.off()

    # Big sized plots
    pdf(file = paste0(figure_dir, sample_name, "_plots_mapscvi_hypomap_big_size.pdf"), width = 25, height = 18)
    plot(plot_query_labels(query_seura_object=sc_data,reference_seurat=mapscvi::reference_hypoMap_downsample,overlay = FALSE,
                      label_col = "C66_named",labelonplot = TRUE))
    dev.off()
}

# Read in SeuratObject and run mapscvi
RunMapSCVI <- function(sample_name) {
    sc_data <- readRDS(paste0(raw_object_dir, sample_name, "_sc_data_seurat_object.rds"))
    print(sc_data) # should be 1 assay
    print(names(sc_data@reductions)) # should be NULL

    # Run MapSCVI
    sc_data <- mapscvi::map_new_seurat_hypoMap(sc_data, suffix = sample_name, label_col = "C66_named", max_epochs = 20)
    print(names(sc_data@reductions))
    print(head(sort(table(sc_data@meta.data$predicted),decreasing = TRUE),n = 10))

    # Save HypoMap SeuratObject
    saveRDS(sc_data, paste0(object_dir, sample_name, "_sc_data_seurat_object_hypomap_run.rds"))

    # Save predicted clusters CSV
    write.csv(sort(table(sc_data@meta.data$predicted),decreasing = TRUE), file = 
              paste0(predicted_cluster_dir, sample_name, "_file_names_clusters_mapscvi_hypomap.csv"))

    # Plot MapSCVI plots
    PlotsMapSCVI(sample_name = sample_name, sc_data = sc_data)
    
    return(sc_data)
}

# Run mapscvi
RunMapSCVI(sample_id)