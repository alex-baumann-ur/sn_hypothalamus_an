# Pathway enrichment analysis
# Set PROJECT_DIR, DATASET_ID, CONTROL_LABEL, and TREATMENT_LABEL before running.
# Group labels must match the input filenames and differential-expression tables.

# Configuration
project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/project")
dataset_id <- Sys.getenv("DATASET_ID", unset = "dataset")
control_label <- Sys.getenv("CONTROL_LABEL", unset = "control")
treatment_label <- Sys.getenv("TREATMENT_LABEL", unset = "treatment")
de_dir <- file.path(project_dir, "results", "differential_expression")
de_comparison_dir <- file.path(
  de_dir,
  paste0("de_markers_", tolower(control_label), "_vs_", tolower(treatment_label))
)
pathway_dir <- file.path(project_dir, "results", "pathway_analysis")
figure_dir <- file.path(project_dir, "results", "figures")

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
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichR)
library(stringr)
library(readxl)

de_markers_path <- file.path(de_comparison_dir, "filtered", "")
path_out <- file.path(pathway_dir, "")

marker_gene_dict <- list()
marker_gene_dict[[paste0("all_markers_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_all_markers_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("all_markers_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_all_markers_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("Oligodendrocytes_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Oligodendrocytes_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("Oligodendrocytes_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Oligodendrocytes_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("GABA_Neurons_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GABA_Neurons_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("GABA_Neurons_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GABA_Neurons_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("Astrocytes_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Astrocytes_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("Astrocytes_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Astrocytes_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("Microglia_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Microglia_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("Microglia_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Microglia_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("OPC_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_OPC_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("OPC_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_OPC_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("GLU_Neurons_", treatment_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GLU_Neurons_de_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("GLU_Neurons_", control_label)]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GLU_Neurons_de_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))

# specific questions with Agrp and Pomc
specific_markers_path <- file.path(de_comparison_dir, "other_explicit_de", "filtered", "")

marker_gene_dict[[paste0("Agrp_", treatment_label)]] <- read.csv(paste0(
  specific_markers_path,
  paste0("de_markers_agrp_neurons_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("Agrp_", control_label)]] <- read.csv(paste0(
  specific_markers_path,
  paste0("de_markers_agrp_neurons_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))
marker_gene_dict[[paste0("Pomc_", treatment_label)]] <- read.csv(paste0(
  specific_markers_path,
  paste0("de_markers_pomc_neurons_filtered_pos_log2fc0_5_pval005_pct01_", treatment_label, ".csv")
))
marker_gene_dict[[paste0("Pomc_", control_label)]] <- read.csv(paste0(
  specific_markers_path,
  paste0("de_markers_pomc_neurons_filtered_neg_log2fc0_5_pval005_pct01_", control_label, ".csv")
))

# with all filtered markers for up and down
marker_gene_dict[["all_markers_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_all_markers_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["Astrocytes_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Astrocytes_de_filtered_neg_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["GABA_Neurons_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GABA_Neurons_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["GLU_Neurons_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_GLU_Neurons_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["Microglia_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Microglia_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["Oligodendrocytes_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_Oligodendrocytes_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["OPC_ALL"]] <- read.csv(paste0(
  de_markers_path,
  paste0(dataset_id, "_OPC_de_filtered_pos_log2fc0_5_pval005_pct01_ALL.csv")
))
marker_gene_dict[["Agrp_ALL"]] <- read.csv(paste0(
  specific_markers_path,
  "de_markers_agrp_neurons_filtered_log2fc0_5_pval005_pct01_ALL.csv"
))
marker_gene_dict[["Pomc_ALL"]] <- read.csv(paste0(
  specific_markers_path,
  "de_markers_pomc_neurons_filtered_log2fc0_5_pval005_pct01_ALL.csv"
))


gene_set_lists <- list()
for (cell_type in names(marker_gene_dict)) {
  gene_set_lists[[cell_type]] <- marker_gene_dict[[cell_type]]$gene
}

head(gene_set_lists[["Agrp_ALL"]])

dbs <- listEnrichrDbs()
length(dbs)
head(dbs)

# find databases (column "libraryName") with mouse in the name
dbs[grep("Mouse", dbs$libraryName), ]


# Databases for mouse: Tabula_Muris, PerturbAtlas_MouseGenePerturbationSigs, WikiPathways_2024_Mouse, KOMP2_Mouse_Phenotypes_2022,
# HDSigDB_Mouse_2021, KEGG_2019_Mouse, WikiPathways_2019_Mouse, Mouse_Gene_Atlas

# gene_set_enrichtment_all_markers <- enrichr(gene_set_list_markers_all, databases = c(
#   "Tabula_Muris", "WikiPathways_2024_Mouse", "KOMP2_Mouse_Phenotypes_2022", #, "PerturbAtlas_MouseGenePerturbationSigs"
#   "HDSigDB_Mouse_2021", "KEGG_2019_Mouse", "Mouse_Gene_Atlas"))
# printEnrich(gene_set_enrichtment_all_markers, prefix = paste0(path_out, "All_Markers_enrichment"),
#               outFile = "excel")
# pdf(paste0(path_out, "All_Markers_enrichment", ".pdf"))
# for (db in names(gene_set_enrichtment_all_markers)) {
#   title <- paste0("All_Markers enrichment analysis ")
#   sub_title <- paste0(db, " database")
#   plot(plotEnrich(gene_set_enrichtment_all_markers[[db]])+ggtitle(label=title, subtitle=sub_title))
# }
# dev.off()

for (gene_set_list in c(paste0("Pomc_", treatment_label))) {
  # for (gene_set_list in c(paste0("Agrp_", control_label), paste0("Agrp_", treatment_label), paste0("Pomc_", control_label), paste0("Pomc_", treatment_label))) {
  # for (gene_set_list in names(gene_set_lists)) {
  gene_set_enrichment <- enrichr(gene_set_lists[[gene_set_list]], databases = c(
    "Tabula_Muris", "WikiPathways_2024_Mouse", "KOMP2_Mouse_Phenotypes_2022", # , "PerturbAtlas_MouseGenePerturbationSigs"
    "HDSigDB_Mouse_2021", "KEGG_2019_Mouse", "Mouse_Gene_Atlas", "JASPAR_PWM_Mouse_2025"
  ))
  printEnrich(gene_set_enrichment,
    prefix = paste0(path_out, gene_set_list, "_enrichment"),
    outFile = "excel"
  )
  # Name of the gene set list as a title for the plot
  pdf(paste0(path_out, gene_set_list, "_enrichment", ".pdf"))
  for (db in names(gene_set_enrichment)) {
    title <- paste0(gene_set_list, " enrichment analysis ")
    sub_title <- paste0(db, " database")
    plot(plotEnrich(gene_set_enrichment[[db]]) + ggtitle(label = title, subtitle = sub_title))
  }
  dev.off()
}


# ### Create VolcanoPlots for the cell types -> already done in 02_de_analysis.R
# raw_marker_gene_dict <- list()
# for (cell_type in names(marker_gene_dict)) {
#     raw_marker_gene_dict[[cell_type]] <- read.csv(paste0(de_markers_path,
#                                 paste0(paste0(dataset_id, "_"), cell_type, "_de_markers.csv")))
# }

# pdf(paste0(de_dirVolcanoPlots_cell_types.pdf"))
# for (cell_type in names(raw_marker_gene_dict)) {
#     top_genes <- raw_marker_gene_dict[[cell_type]][order(raw_marker_gene_dict[[cell_type]]$p_val_adj), ][1:10, ]
#     # Add a column with the significance
#     raw_marker_gene_dict[[cell_type]]$significance <- "NS"
#     raw_marker_gene_dict[[cell_type]]$significance[raw_marker_gene_dict[[cell_type]]$p_val_adj < 0.05 &
#         raw_marker_gene_dict[[cell_type]]$avg_log2FC > 0.5 &
#         (raw_marker_gene_dict[[cell_type]]$pct.1 > 0.1 | raw_marker_gene_dict[[cell_type]]$pct.2 > 0.1)] <- "Up"
#     raw_marker_gene_dict[[cell_type]]$significance[raw_marker_gene_dict[[cell_type]]$p_val_adj < 0.05 &
#         raw_marker_gene_dict[[cell_type]]$avg_log2FC < -0.5 &
#         (raw_marker_gene_dict[[cell_type]]$pct.1 > 0.1 | raw_marker_gene_dict[[cell_type]]$pct.2 > 0.1)] <- "Down"
#     plot(ggplot(raw_marker_gene_dict[[cell_type]], aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
#         geom_point(alpha = 0.6) +
#         geom_text_repel(data = top_genes, aes(label = gene), color = "black", max.overlaps = 100) +
#         scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
#         geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
#         geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#         labs(
#             x = "Average log2 Fold Change",
#             y = "-log10 P-value",
#             color = "Significance"
#         ) +
#         theme_minimal() +
#         ggtitle(paste0(cell_type, " markers"))+xlab("Average log2 Fold Change")+ylab("-log10 P-value"))
# }
# dev.off()


# genes_test <- gene_set_entrezid$converted_alias

# keytypes(org.Mm.eg.db) -> SYMBOL also possible as keytype
head(gene_set_lists[[paste0("Agrp_", treatment_label, "")]])

# ego <- enrichGO(gene          = gene_set_lists[[paste0("Agrp_", treatment_label, "")]],
#                 keyType       = "SYMBOL",
#                 OrgDb         = org.Mm.eg.db,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 readable      = TRUE)
# head(ego)
# write.csv(ego, paste0(path_out, paste0("Agrp_", treatment_label, "_GO_terms_BP_symbol.csv")))

### GO terms with SYMBOL as keytype
db_dict <- list("MF" = "Molecular Function", "BP" = "Biological Process", "CC" = "Cellular Component")

enrichgo_result <- function(gene_set_list, db) {
  ego_in <- enrichGO(
    gene = gene_set_lists[[gene_set_list]],
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db,
    ont = db,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  write.csv(ego_in, paste0(path_out, gene_set_list, "_GO_terms_", db, ".csv"))

  print(nrow(ego_in))

  # capitalize the first letter of the Description column
  if (nrow(ego_in) == 0) {
    return(NULL)
  }
  res <- ego_in@result
  res$Description <- sub("^([a-zA-Z])", "\\U\\1", res$Description, perl = TRUE)
  ego_in@result <- res

  return(ego_in)
}

enrichgo_figure <- function(ego_in, db, gene_set_list) {
  p <- dotplot(ego_in, showCategory = 15) +
    # scale_fill_gradientn(colors = c("#00838F", "#4DD0E1", "#E0F7FA")) +
    scale_fill_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
    ggtitle(
      label = paste0(gene_set_list, " enrichment analysis"),
      subtitle = paste0("GO ", db_dict[[db]], " database")
    )
  return(p)
}

for (gene_set_list in c(paste0("GABA_Neurons_", treatment_label, ""))) { # , "GLU_Neurons_ALL",
  # "Microglia_ALL", "Oligodendrocytes_ALL", "OPC_ALL")) {
  # for (gene_set_list in names(gene_set_lists)) {
  # for (gene_set_list in c(paste0("Oligodendrocytes_", control_label, ""), paste0("GABA_Neurons_", treatment_label, ""), paste0("GABA_Neurons_", control_label, ""), paste0("Astrocytes_", treatment_label, ""), paste0("Astrocytes_", control_label, ""),
  #             paste0("Microglia_", treatment_label, ""), paste0("Microglia_", control_label, ""), paste0("OPC_", treatment_label, ""), paste0("OPC_", control_label, ""), paste0("GLU_Neurons_", treatment_label, ""), paste0("GLU_Neurons_", control_label, ""), paste0("Agrp_", control_label, ""), paste0("Pomc_", treatment_label, ""), paste0("Pomc_", control_label, ""))) {
  in_list <- list()
  print(gene_set_list)
  for (db in c("BP")) {
    # for (db in c("MF", "BP", "CC")) {
    ego <- enrichgo_result(gene_set_list, db)
    if (!is.null(ego)) {
      print(db)
      in_list[[db]] <- enrichgo_figure(ego, db, gene_set_list)
    }
  }
  # pdf(paste0(path_out, gene_set_list, "_GO_terms_purple.pdf"), width = 10, height = 10)
  for (name in names(in_list)) {
    tiff(paste0(path_out, gene_set_list, "_GO_terms_", name, ".tiff"), units = "cm", width = 20, height = 20, res = 600)
    plot(in_list[[name]])
    dev.off()
  }
  # dev.off()
}

astro_ego <- enrichgo_result("Astrocytes_ALL", "MF")
head(astro_ego@result)

astro_ego <- enrichgo_result("Astrocytes_ALL", "BP")
head(astro_ego@result)

astro_ego <- enrichgo_result("Astrocytes_ALL", "CC")
head(astro_ego@result)


# Plots with GeneRatio vs. Term
gsea_path <- file.path(figure_dir, "gsea", "")

kegg_lists <- list()
for (gene_set_list in names(gene_set_lists)) {
  if (!(gene_set_list %in% c(
    paste0("Microglia_", control_label, ""), paste0("Pomc_", control_label, ""), "all_markers_ALL", "GABA_Neurons_ALL", "GLU_Neurons_ALL",
    "Microglia_ALL", "Oligodendrocytes_ALL", "OPC_ALL", "Agrp_ALL", "Pomc_ALL"
  ))) {
    kegg_lists[[gene_set_list]] <- read_xlsx(paste0(path_out, gene_set_list, "_enrichment.xlsx"), sheet = "KEGG_2019_Mouse")
  }
}

# for (gene_set_list in c(paste0("GABA_Neurons_", control_label, ""), paste0("Oligodendrocytes_", control_label, ""), paste0("all_markers_", control_label, ""), paste0("GLU_Neurons_", control_label, ""), paste0("OPC_", control_label, ""),
#   paste0("Astrocytes_", control_label, ""), paste0("Agrp_", treatment_label, ""), paste0("Agrp_", control_label, ""), paste0("Pomc_", treatment_label, ""))) {
#   kegg_lists[[gene_set_list]] <- read_xlsx(paste0(path_out, gene_set_list, "_enrichment.xlsx"), sheet = "KEGG_2019_Mouse")
# }
head(kegg_lists[[paste0("GABA_Neurons_", control_label, "")]])

names(kegg_lists)

plot_top_20_go_terms_kegg_padjust <- function(cell_type, database_used, title_between, path_output) {
  enrichment_result <- kegg_lists[[cell_type]]
  enrichment_result_top20 <- enrichment_result[1:20, ]
  # rename adjusted p-value to p.adjust
  enrichment_result_top20$Term <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Term, perl = TRUE)
  enrichment_result_top20$p.adjust <- enrichment_result_top20$Adjusted.P.value
  # filter for adjusted p-value < 0.05 and remove rows with NA
  enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$p.adjust < 0.05 & !is.na(enrichment_result_top20$Term), ]
  if (nrow(enrichment_result_top20) > 0) {
    enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[, 1])
    enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[, 2])
    enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
    # enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
    enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$GeneRatio)])
    if (nrow(enrichment_result_top20) > 20) {
      enrichment_result_top20 <- enrichment_result_top20[1:20, ]
    }
    p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Term, color = p.adjust, size = Count)) +
      geom_point() +
      theme_minimal() +
      theme(axis.title.y = element_blank()) +
      scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
      ggtitle(label = paste0(cell_type, " ", title_between, " enrichment analysis"), subtitle = paste0(str_replace_all(database_used, "_", " "), " database"))
    tiff(paste0(gsea_path, str_replace_all(cell_type, " ", "_"), "_", database_used, "padjust.tiff"), units = "cm", width = 20, height = 20, res = 600)
    plot(p)
    dev.off()
  }
}

# plot_top_20_go_terms_kegg_pval <- function(cell_type, database_used, title_between, path_output) {
#   enrichment_result <- kegg_lists[[cell_type]]
#   enrichment_result_top20 <- enrichment_result[1:20,]
#   # rename adjusted p-value to p.adjust
#   enrichment_result_top20$Term <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Term, perl = TRUE)
#   enrichment_result_top20$pvalue <- enrichment_result_top20$P.value
#   # filter for adjusted p-value < 0.05 and remove rows with NA
#   enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$pvalue < 0.05 & !is.na(enrichment_result_top20$Term),]
#   if (nrow(enrichment_result_top20) > 0) {
#     enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[,1])
#     enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[,2])
#     enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
#     # enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
#     enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$GeneRatio)])
#     if (nrow(enrichment_result_top20) > 20) {
#       enrichment_result_top20 <- enrichment_result_top20[1:20,]
#     }
#     p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Term, color = pvalue, size = Count)) +
#       geom_point() +
#       theme_minimal() +
#       theme(axis.title.y = element_blank()) +
#       scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
#       ggtitle(label = paste0(cell_type, " ", title_between, " enrichment analysis"), subtitle = paste0(str_replace_all(database_used, "_", " "), " database"))
#     tiff(paste0(gsea_path, str_replace_all(cell_type, " ", "_"), "_", database_used, "_pvalue.tiff"), units = "cm", width = 20, height = 20, res = 600)
#     plot(p)
#     dev.off()
#   }
# }

### Newest version of the plotting function for KEGG
# plot_top_20_go_terms_kegg_pval <- function(cell_type, cell_type_name, database_used, path_output) {
#   enrichment_result <- kegg_lists[[cell_type]]
#   enrichment_result_top20 <- enrichment_result[1:20,]
#   # rename adjusted p-value to p.adjust
#   enrichment_result_top20$Term <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Term, perl = TRUE)
#   enrichment_result_top20$pvalue <- enrichment_result_top20$P.value
#   # filter for adjusted p-value < 0.05 and remove rows with NA
#   enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$pvalue < 0.05 & !is.na(enrichment_result_top20$Term),]
#   if (nrow(enrichment_result_top20) > 0) {
#     enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[,1])
#     enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[,2])
#     enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
#     # enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
#     enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$GeneRatio)])
#     if (nrow(enrichment_result_top20) > 20) {
#       enrichment_result_top20 <- enrichment_result_top20[1:20,]
#     }
#     p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Term, color = pvalue, size = Count)) +
#       geom_point() + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
#       theme_minimal() +
#       theme(axis.title.y = element_blank(), axis.text.y = element_text(face = "bold")) +
#       scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
#       ggtitle(label = paste0(cell_type_name, ": ", str_replace_all(database_used, "_", " ")))
#     tiff(paste0(gsea_path, "Fig4_", str_replace_all(cell_type, " ", "_"), "_", database_used, "_pvalue.tiff"), units = "cm", width = 20, height = 20, res = 600)
#     plot(p)
#     dev.off()
#   }
# }


plot_top_20_go_terms_kegg_pval <- function(cell_type, cell_type_name, database_used, path_output) {
  enrichment_result <- kegg_lists[[cell_type]]
  enrichment_result_top20 <- enrichment_result[1:20, ]
  # rename adjusted p-value to p.adjust
  enrichment_result_top20$Term <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Term, perl = TRUE)
  enrichment_result_top20$pvalue <- enrichment_result_top20$P.value
  # filter for adjusted p-value < 0.05 and remove rows with NA
  enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$pvalue < 0.05 & !is.na(enrichment_result_top20$Term), ]
  if (nrow(enrichment_result_top20) > 0) {
    enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[, 1])
    enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$Overlap, "/", 2)[, 2])
    enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
    # enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
    enrichment_result_top20$Term <- factor(enrichment_result_top20$Term, levels = enrichment_result_top20$Term[order(enrichment_result_top20$GeneRatio)])
    if (nrow(enrichment_result_top20) > 20) {
      enrichment_result_top20 <- enrichment_result_top20[1:20, ]
    }
    p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Term, color = pvalue, size = Count)) +
      geom_point() +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      theme_minimal() +
      theme(axis.title.y = element_blank(), axis.text.y = element_text(face = "bold")) +
      scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
      ggtitle(label = paste0(cell_type_name, ": ", str_replace_all(database_used, "_", " "))) +
      guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))
    tiff(paste0(path_output, "Fig4_", str_replace_all(cell_type, " ", "_"), "_", database_used, "_pvalue.tiff"), units = "cm", width = 20, height = 23, res = 600)
    plot(p)
    dev.off()
  }
}


plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("GABA_Neurons_", treatment_label, ""), cell_type_name = "GABA neurons upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Oligodendrocytes_", treatment_label, ""), cell_type_name = "Oligodendrocytes upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("GLU_Neurons_", treatment_label, ""), cell_type_name = "GLU neurons upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("all_markers_", treatment_label, ""), cell_type_name = "All clusters upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("OPC_", treatment_label, ""), cell_type_name = "OPC upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Astrocytes_", treatment_label, ""), cell_type_name = "Astrocytes upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)

plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("GABA_Neurons_", control_label, ""), cell_type_name = "GABA neurons downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Oligodendrocytes_", control_label, ""), cell_type_name = "Oligodendrocytes downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("GLU_Neurons_", control_label, ""), cell_type_name = "GLU neurons downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("all_markers_", control_label, ""), cell_type_name = "All clusters downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("OPC_", control_label, ""), cell_type_name = "OPC downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Astrocytes_", control_label, ""), cell_type_name = "Astrocytes downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Agrp_", control_label, ""), cell_type_name = "Agrp neurons downregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Agrp_", treatment_label, ""), cell_type_name = "Agrp neurons upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)
plot_top_20_go_terms_kegg_pval(
  cell_type = paste0("Pomc_", treatment_label, ""), cell_type_name = "Pomc neurons upregulated", database_used = "KEGG_2019_Mouse",
  path_output = gsea_path
)

# for (cell_type in names(kegg_lists)) {
#   print(paste0("Plotting top 20 KEGG terms for ", cell_type))
#   plot_top_20_go_terms_kegg_pval(cell_type = cell_type, database_used = "KEGG_2019_Mouse", title_between = paste0("(vs. ", control_label, ")"), path_output = gsea_path)
#   plot_top_20_go_terms_kegg_padjust(cell_type = cell_type, database_used = "KEGG_2019_Mouse", title_between = paste0("(vs. ", control_label, ")"), path_output = gsea_path)
# }

### Plots with GeneRatio vs. Term for BP and CC
# bp_cc_path <- pathway_dirGO_terms/GO_tables/"
bp_cc_path <- file.path(figure_dir, "gsea", "")

names(gene_set_lists)

# bp_lists <- list()
# for (gene_set_list in c(paste0("GABA_Neurons_", treatment_label), paste0("Oligodendrocytes_", treatment_label), paste0("all_markers_", treatment_label), paste0("GLU_Neurons_", treatment_label))) {
#   bp_lists[[gene_set_list]] <- read.csv(paste0(bp_cc_path, gene_set_list, "_GO_terms_BP.csv"))
# }

bp_lists <- list()
for (gene_set_list in names(gene_set_lists)) {
  bp_lists[[gene_set_list]] <- read.csv(paste0(bp_cc_path, gene_set_list, "_GO_terms_BP.csv"))
}

names(bp_lists)

cc_lists <- list()
for (gene_set_list in names(gene_set_lists)) {
  cc_lists[[gene_set_list]] <- read.csv(paste0(bp_cc_path, gene_set_list, "_GO_terms_CC.csv"))
}

head(bp_lists[[paste0("GABA_Neurons_", control_label)]])

### This is the most recent version of the plotting function
# plot_top_20_go_terms_go_padjust <- function(cell_type, cell_type_name, database_used, path_output) {
#   if (database_used == "BP") {
#     enrichment_result <- bp_lists[[cell_type]]
#     db_full_name <- "GO Biological Process"
#   } else if (database_used == "CC") {
#     enrichment_result <- cc_lists[[cell_type]]
#     db_full_name <- "GO Cellular Component"
#   }
#   # enrichment_result_top20 <- enrichment_result[1:20,]
#   enrichment_result_top20 <- enrichment_result
#   enrichment_result_top20$Description <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Description, perl = TRUE)
#   # filter for adjusted p-value < 0.05 and remove rows with NA
#   enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$p.adjust < 0.05 & !is.na(enrichment_result_top20$Description),]
#   if (nrow(enrichment_result_top20) > 0) {
#     enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$GeneRatio, "/", 2)[,1])
#     enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$GeneRatio, "/", 2)[,2])
#     enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
#     # enrichment_result_top20$Description <- factor(enrichment_result_top20$Description, levels = enrichment_result_top20$Description[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
#     enrichment_result_top20$Description <- factor(enrichment_result_top20$Description, levels = enrichment_result_top20$Description[order(enrichment_result_top20$GeneRatio)])
#     if (nrow(enrichment_result_top20) > 20) {
#       enrichment_result_top20 <- enrichment_result_top20[1:20,]
#     }
#     p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
#       geom_point() + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
#       theme_minimal() +
#       theme(axis.title.y = element_blank(), axis.text.y = element_text(face = "bold")) +
#       scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
#       ggtitle(label = paste0(cell_type_name, ": ", db_full_name))
#     tiff(paste0(gsea_path, "Fig4_", str_replace_all(cell_type, " ", "_"), "_", database_used, "_padjust.tiff"), units = "cm", width = 20, height = 20, res = 600)
#     plot(p)
#     dev.off()
#   }
# }


plot_top_20_go_terms_go_padjust <- function(cell_type, cell_type_name, database_used, path_output) {
  if (database_used == "BP") {
    enrichment_result <- bp_lists[[cell_type]]
    db_full_name <- "GO Biological Process"
  } else if (database_used == "CC") {
    enrichment_result <- cc_lists[[cell_type]]
    db_full_name <- "GO Cellular Component"
  }
  # enrichment_result_top20 <- enrichment_result[1:20,]
  enrichment_result_top20 <- enrichment_result
  enrichment_result_top20$Description <- sub("^([a-zA-Z])", "\\U\\1", enrichment_result_top20$Description, perl = TRUE)
  # filter for adjusted p-value < 0.05 and remove rows with NA
  enrichment_result_top20 <- enrichment_result_top20[enrichment_result_top20$p.adjust < 0.05 & !is.na(enrichment_result_top20$Description), ]
  if (nrow(enrichment_result_top20) > 0) {
    enrichment_result_top20$Count <- as.numeric(str_split_fixed(enrichment_result_top20$GeneRatio, "/", 2)[, 1])
    enrichment_result_top20$All_Genes_Count <- as.numeric(str_split_fixed(enrichment_result_top20$GeneRatio, "/", 2)[, 2])
    enrichment_result_top20$GeneRatio <- enrichment_result_top20$Count / enrichment_result_top20$All_Genes_Count
    # enrichment_result_top20$Description <- factor(enrichment_result_top20$Description, levels = enrichment_result_top20$Description[order(enrichment_result_top20$p.adjust, decreasing = TRUE)])
    enrichment_result_top20$Description <- factor(enrichment_result_top20$Description, levels = enrichment_result_top20$Description[order(enrichment_result_top20$GeneRatio)])
    if (nrow(enrichment_result_top20) > 20) {
      enrichment_result_top20 <- enrichment_result_top20[1:20, ]
    }
    p <- ggplot(enrichment_result_top20, aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
      geom_point() +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      theme_minimal() +
      theme(axis.title.y = element_blank(), axis.text.y = element_text(face = "bold")) +
      scale_color_gradientn(colors = c("darkviolet", "orchid", "lightgrey")) +
      ggtitle(label = paste0(cell_type_name, ": ", db_full_name)) +
      guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))
    tiff(paste0(path_output, "Fig4_", str_replace_all(cell_type, " ", "_"), "_", database_used, "_padjust.tiff"), units = "cm", width = 20, height = 23, res = 600)
    plot(p)
    dev.off()
  }
}

plot_top_20_go_terms_go_padjust(cell_type = paste0("GABA_Neurons_", treatment_label), cell_type_name = "GABA neurons upregulated", database_used = "BP", path_output = bp_cc_path)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Oligodendrocytes_", treatment_label), cell_type_name = "Oligodendrocytes upregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("GLU_Neurons_", treatment_label), cell_type_name = "GLU neurons upregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("all_markers_", treatment_label), cell_type_name = "All clusters upregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("OPC_", treatment_label), cell_type_name = "OPC upregulated", database_used = "CC",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Astrocytes_", treatment_label), cell_type_name = "Astrocytes upregulated", database_used = "CC",
  path_output = bp_cc_path
)

plot_top_20_go_terms_go_padjust(cell_type = paste0("GABA_Neurons_", control_label), cell_type_name = "GABA neurons downregulated", database_used = "BP", path_output = bp_cc_path)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Oligodendrocytes_", control_label), cell_type_name = "Oligodendrocytes downregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("GLU_Neurons_", control_label), cell_type_name = "GLU neurons downregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("all_markers_", control_label), cell_type_name = "All clusters downregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Agrp_", control_label), cell_type_name = "Agrp neurons downregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Agrp_", treatment_label), cell_type_name = "Agrp neurons upregulated", database_used = "BP",
  path_output = bp_cc_path
)
plot_top_20_go_terms_go_padjust(
  cell_type = paste0("Pomc_", treatment_label), cell_type_name = "Pomc neurons upregulated", database_used = "BP",
  path_output = bp_cc_path
)

# for (cell_type in names(bp_lists)) {
#   print(paste0("Plotting top 20 GO terms for ", cell_type))
#   plot_top_20_go_terms_go_padjust(cell_type = cell_type, database_used = "BP", title_between = "(vs. Control)", path_output = bp_cc_path)
# }

# for (cell_type in names(cc_lists)) {
#   print(paste0("Plotting top 20 GO terms for ", cell_type))
#   plot_top_20_go_terms_go_padjust(cell_type = cell_type, database_used = "CC", title_between = "(vs. Control)", path_output = bp_cc_path)
# }


##############
library(msigdb)

de_markers_path <- file.path(de_comparison_dir, "unfiltered", "")
path_out <- file.path(pathway_dir, "enrichment_both", "")

de_genes_all_clusters <- read.csv(paste0(de_markers_path, paste0(dataset_id, "_all_markers_de_markers_", treatment_label, "_vs_", control_label, ".csv")))
head(de_genes_all_clusters)
ranks <- de_genes_all_clusters$avg_log2FC
names(ranks) <- de_genes_all_clusters$gene

head(ranks)

hallmarks_set <- msigdbr(
  species = "Mus musculus",
  collection = "MH",
  db_species = "MM"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)
