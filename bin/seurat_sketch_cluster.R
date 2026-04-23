required_packages <- c(
    "Matrix",
    "ggplot2", 
    "dplyr",
    "cowplot",
    "patchwork",
    "remotes",
    "devtools"
)

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

# Install Seurat 5
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}

# Verify installation
library(Seurat)
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")

library(dplyr)
library(ggplot2)
library(Matrix)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    cat("Usage: Rscript seurat_sketch_cluster.R <sample_id> <input_path> <output_dir> <prefix> <option>\n")
    cat("  sample_id:   Sample identifier\n")
    cat("  input_path:  Path to input directory (bin2cell output or spaceranger outs)\n")
    cat("  output_dir:  Output directory for results\n")
    cat("  prefix:      Prefix for cell IDs\n")
    cat("  option:      'Cellseg' for bin2cell output or 'Default' for spaceranger output\n")
    quit(status = 1)
}

sample_id <- args[1]
input_path <- args[2]
output_dir <- args[3]
prefix <- args[4]
option <- args[5]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Processing sample:", sample_id, "\n")
cat("Input path:", input_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("Cell ID prefix:", prefix, "\n")
cat("Pipeline option:", option, "\n")

load_bin2cell_data <- function(prefix, path, ident) {
    cat("Loading bin2cell segmented data from:", path, "\n")
    
    # List available files for debugging
    available_files <- list.files(path, pattern = "\\.(mtx|csv)$", full.names = FALSE)
    cat("Available files:", paste(available_files, collapse = ", "), "\n")
    
    count_files <- c("adata_count.mtx")
    genes_files <- c("adata_genes.csv")
    cells_files <- c("adata_cells.csv")
    meta_files <- c("adata_meta.csv")
    
    # Find the actual files
    count_file <- NULL
    for (f in count_files) {
        if (file.exists(file.path(path, f))) {
            count_file <- f
            break
        }
    }
    
    genes_file <- NULL
    for (f in genes_files) {
        if (file.exists(file.path(path, f))) {
            genes_file <- f
            break
        }
    }
    
    cells_file <- NULL
    for (f in cells_files) {
        if (file.exists(file.path(path, f))) {
            cells_file <- f
            break
        }
    }
    
    meta_file <- NULL
    for (f in meta_files) {
        if (file.exists(file.path(path, f))) {
            meta_file <- f
            break
        }
    }
    
    # Check if required files exist
    if (is.null(count_file)) {
        stop("Count matrix file not found. Looking for: ", paste(count_files, collapse = ", "))
    }
    if (is.null(genes_file)) {
        stop("Genes file not found. Looking for: ", paste(genes_files, collapse = ", "))
    }
    if (is.null(cells_file)) {
        stop("Cells file not found. Looking for: ", paste(cells_files, collapse = ", "))
    }
    
    # Read count matrix
    counts <- readMM(file.path(path, count_file))
    
    # Read genes
    genes_df <- read.csv(file.path(path, genes_file))
    if ("genes" %in% colnames(genes_df)) {
        genes <- genes_df$genes
    } else {
        genes <- genes_df[,1]
    }
    
    # Read cells
    cells_df <- read.csv(file.path(path, cells_file))
    if ("cells" %in% colnames(cells_df)) {
        cells <- cells_df$cells
    } else {
        cells <- cells_df[,1]
    } 
    
    # Transpose and convert to sparse matrix
    counts <- as(t(counts), "dgCMatrix")
    rownames(counts) <- genes
    colnames(counts) <- paste0(prefix, "_", cells)  # Add prefix to cell IDs
    
    cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")
    
    
    # Read metadata if available
    meta <- NULL
    if (!is.null(meta_file)) {
        meta_df <- read.csv(file.path(path, meta_file))
        rownames(meta_df) <- paste0(prefix, "_", meta_df[,1])  # First column as row names with prefix
        meta_df <- meta_df[,-1]  # Drop the first column
        meta_df$sample <- ident  # Assign identity label
        meta <- meta_df[colnames(counts), , drop = FALSE]  # Reorder to match cells
    }
    
    return(list(counts = counts, meta = meta))
}

load_spaceranger_data <- function(prefix, path, ident) {
    cat("Loading SpaceRanger data from:", path, "\n")
    
    # Look for the filtered feature matrix in different possible locations
    possible_paths <- c(
        file.path(path, "segmented_outputs", "filtered_feature_bc_matrix.h5"),
        file.path(path, "filtered_feature_bc_matrix.h5"),
        file.path(path, "binned_outputs", "square_008um", "filtered_feature_bc_matrix.h5")
    )
    
    # Find the first existing path
    data_path <- NULL
    for (p in possible_paths) {
        if (file.exists(p) || dir.exists(p)) {
            data_path <- p
            break
        }
    }
    
    if (is.null(data_path)) {
        stop("Could not find filtered_feature_bc_matrix in any of these locations:\n", 
             paste(possible_paths, collapse = "\n"))
    }
    
    cat("Found data at:", data_path, "\n")
    
    # Read the 10X data
    if (grepl("\\.h5$", data_path)) {
        counts <- Read10X_h5(data_path)
    } else {
        stop("Unsupported file format. Expected .h5 file, got: ", data_path)
    }
    
    cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")
    
    # Create minimal metadata
    meta <- data.frame(
        sample = ident,
        row.names = colnames(counts)
    )
        
    return(list(counts = counts, meta = meta))
}

seurat_sketch <- function(data, output_dir, sample_id) {
    cat("Creating Seurat object...\n")
    
    object <- CreateSeuratObject(
        counts = data$counts, 
        meta.data = data$meta, 
        min.cells = 0, 
        min.features = 0
    )
    
    cat("Initial object:", ncol(object), "cells,", nrow(object), "genes\n")
    
    # QC filtering
    hist_path <- file.path(output_dir, paste0(sample_id, "_hist.pdf"))
    pdf(hist_path)
    hist(object@meta.data$nFeature_RNA, breaks = 100, main = paste("Feature distribution -", sample_id))
    dev.off()

    filter_threshold <- quantile(object$nFeature_RNA,0.01)[[1]]
    cat("1% quantile of nFeature_RNA:", filter_threshold, "\n")

    object <- subset(object, subset = nFeature_RNA > filter_threshold)
    
    cat("After filtering (>", filter_threshold, " features):", ncol(object), "cells\n")

    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    
    # Sketch-based analysis
    cat("Performing sketch-based analysis...\n")
    sketch_size <- min(5000, ncol(object))  # Use fewer cells if dataset is small
    

    object <- SketchData(
        object = object,
        ncells = sketch_size,
        method = "LeverageScore",
        sketched.assay = "sketch",
        features = VariableFeatures(object)
    )
    
    DefaultAssay(object) <- "sketch"
    
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
    object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:30)
    object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 2)
    object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:30, min.dist = 0.2)
    
    # Project to full dataset
    object <- ProjectData(
        object = object,
        assay = "RNA",
        full.reduction = "full.pca.sketch",
        sketched.assay = "sketch",
        sketched.reduction = "pca.sketch",
        umap.model = "umap.sketch",
        dims = 1:30,
        refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
    )
    
    # Create plots
    DefaultAssay(object) <- "sketch"
    Idents(object) <- "seurat_cluster.sketched"
    p1 <- DimPlot(object, reduction = "umap.sketch", label = TRUE) + 
            ggtitle(paste("Sketched clustering (", sketch_size, " cells)", sep=""))
    
    DefaultAssay(object) <- "RNA"
    Idents(object) <- "seurat_cluster.projected"
    p2 <- DimPlot(object, reduction = "full.umap.sketch", label = TRUE) + 
            ggtitle("Projected clustering (full dataset)")
    
    # Save plots
    ggsave(file.path(output_dir, paste0(sample_id, "_umap_sketch.png")), 
            p1 + p2, width = 15, height = 5, dpi = 300)
    
    cat("Sketch-based analysis completed\n")

    return(object)
}

# Main execution
cat("Starting Seurat analysis for sample:", sample_id, "\n")

# Load data based on type
if (option == "Cellseg") {
    data <- load_bin2cell_data(prefix = prefix, path = input_path, ident = sample_id)
} else if (option == "Default") {
    data <- load_spaceranger_data(prefix = prefix, path = input_path, ident = sample_id)
} else {
    stop("Invalid option. Must be 'Cellseg' or 'Default', got: '", option, "'")
}

# Run sketch analysis
object <- seurat_sketch(data, output_dir = output_dir, sample_id = sample_id)

# Find marker genes
cat("Finding marker genes...\n")
DefaultAssay(object) <- 'RNA'

# Use the appropriate cluster identity
if ("seurat_cluster.projected" %in% colnames(object@meta.data)) {
    Idents(object) <- "seurat_cluster.projected"
} else {
    Idents(object) <- "seurat_clusters"
}

markers <- FindAllMarkers(object, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.1)
markers <- markers[markers$p_val_adj < 0.01, ]
markers <- markers[markers$avg_log2FC > 0.75, ]

cat("Found", nrow(markers), "significant marker genes\n")

# Save results
marker_file <- file.path(output_dir, paste0(sample_id, "_markers.csv"))
write.csv(markers, file = marker_file, row.names = FALSE)

meta_file <- file.path(output_dir, paste0(sample_id, "_seurat_meta.csv"))
write.csv(object@meta.data, meta_file)

# Save UMAP coordinates
umap_coords <- NULL
if ("full.umap.sketch" %in% names(object@reductions)) {
    umap_coords <- object@reductions$full.umap.sketch@cell.embeddings
}

if (!is.null(umap_coords)) {
    umap_file <- file.path(output_dir, paste0(sample_id, "_seurat_umap.csv"))
    write.csv(umap_coords, umap_file)
}

# Save Seurat object
rds_file <- file.path(output_dir, paste0(sample_id, "_seurat_object.rds"))
saveRDS(object, file = rds_file)

cat("Analysis completed successfully!\n")
cat("Output files:\n")
cat("  Markers:", marker_file, "\n")
cat("  Metadata:", meta_file, "\n")
if (!is.null(umap_coords)) {
    cat("  UMAP coords:", file.path(output_dir, paste0(sample_id, "_seurat_umap.csv")), "\n")
}
cat("  Seurat object:", rds_file, "\n")
cat("  Plots:", file.path(output_dir, paste0(sample_id, "_umap*.png")), "\n")
