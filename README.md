# Single-cell RNA-seq analysis

Scripts for processing Cell Ranger outputs, mapping cells to the HypoMap reference, differential-expression analysis, pathway enrichment, and figure generation.

## Workflow

1. Run [`CellRanger/Snakefile_CellRanger`](CellRanger/Snakefile_CellRanger) with a tab-separated sample sheet and a Cell Ranger reference genome.
2. Run [`00_seurat_analysis_integration.R`](00_seurat_analysis_integration.R) to import the filtered matrices, create and merge Seurat objects, and run mapscvi.
3. Run [`01_hypomap_seurat_analysis.R`](01_hypomap_seurat_analysis.R) for Seurat integration and doublet detection.
4. Run [`02_differential_expression_analysis.R`](02_differential_expression_analysis.R), then [`03_filtering_files.py`](03_filtering_files.py) to filter marker tables.
5. Run [`04_pathway_analysis.R`](04_pathway_analysis.R) for enrichment analyses and [`05_figure_scripts.R`](05_figure_scripts.R) to generate figures.

[`03_filtering_files_files/03_filtering_files.html`](03_filtering_files_files/03_filtering_files.html) is a rendered report associated with the filtering workflow.

## Configuration

The scripts use project-relative `data/`, `metadata/`, and `results/` folders or configurable paths. Set the required project, dataset, condition, and sample values before running; see each script's configuration section. Input metadata, Cell Ranger outputs, and the Cell Ranger reference genome are not included here. Install the R and Python dependencies, including Snakemake, Cell Ranger, Seurat, and mapscvi, in the relevant environments.
