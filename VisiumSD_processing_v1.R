#--- title: VisiumSD_processing
#--- author: Emilia Puig Lombardi
#--- date: 2024-07-12

#--- Environment setup
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(scater)
library(tibble)
library(SingleCellExperiment)
library(hdf5r) 
library(patchwork)
library(gtools)
library(clustree)
library(gridExtra)
library(tidyr)
library(magrittr)
library(STutility)
library(SpotClean)
library(S4Vectors)

# Set up directories and paths
setup_environment <- function(wd, sample_id) {
  datadir <- file.path(wd)
  datadir1 <- file.path(datadir, paste0(sample_id, "/"))
  outputdir <- file.path(datadir, paste0(sample_id, "_out"))
  
  if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }
  
  # File paths
  samples <- file.path(datadir1, "raw_feature_bc_matrix.h5")
  spotfiles <- file.path(datadir1, "spatial/tissue_positions.csv")
  imgs <- file.path(datadir1, "spatial/tissue_hires_image.png")
  json <- file.path(datadir1, "spatial/scalefactors_json.json")
  
  return(list(samples = samples, spotfiles = spotfiles, imgs = imgs, json = json, outputdir = outputdir))
}

# Main function to load and process data
process_sample <- function(wd, sample_id) {
  paths <- setup_environment(wd, sample_id)
  
  # Create info table for STutility input
  infoTable <- data.frame(samples = paths$samples, 
                          spotfiles = paths$spotfiles,
                          imgs = paths$imgs,
                          json = paths$json)
  
  # Load data
  obj <- InputFromTable(infoTable, platform = "Visium", disable.subset = FALSE)
  obj <- LoadImages(obj, time.resolve = FALSE)
  
  # Plot images
  ImagePlot(obj, method = "raster", annotate = FALSE)
  
  return(obj)
}

# Example usage
wd <- "HS_Visium_new"
sample_id <- "TIS4606_PP"
obj <- process_sample(wd, sample_id)

#--- Manual debris removal: opens ShinyApp locally
obj <- ManualAnnotation(obj)
# check labels
FeatureOverlay(obj, features = "labels")

#--- Metadata for later object creation 

# Calculate mitochondrial gene percentage
calculate_mito_percentage <- function(obj, assay = "RNA", pattern = "^MT-") {
  # Extract mitochondrial gene names based on the pattern
  mt.genes <- grep(pattern = pattern, x = rownames(obj), value = TRUE)
  
  # Calculate percent mitochondrial content
  mito_counts <- Matrix::colSums(LayerData(obj, assay = assay, layer = "counts")[mt.genes, ])
  total_counts <- Matrix::colSums(LayerData(obj, assay = assay, layer = "counts"))
  obj$percent.mito <- (mito_counts / total_counts) * 100
  
  return(obj)
}

# Clean and organize metadata in a df
clean_metadata <- function(obj) {
  # Retrieve and clean metadata
  meta.preqc <- obj@meta.data
  rownames(meta.preqc) <- gsub("_1", "", rownames(meta.preqc))
  meta.preqc$Barcode <- rownames(meta.preqc)
  
  # Select and reorder specific columns
  meta.preqc <- meta.preqc[, c(7, 5:6)]
  meta.preqc <- meta.preqc[order(meta.preqc$Barcode), ]
  
  return(meta.preqc)
}

# Example usage
obj <- calculate_mito_percentage(obj)
meta.preqc <- clean_metadata(obj)

#--- swapping correction / decontamination

# Read Visium data with SpotClean
load_visium_data <- function(count_file, tissue_csv, tissue_img, scale_factors) {
  obj_raw <- read10xRawH5(count_file)
  
  slide_info <- read10xSlide(
    tissue_csv_file = tissue_csv,
    tissue_img_file = tissue_img,
    scale_factor_file = scale_factors
  )
  
  obj <- createSlide(count_mat = obj_raw, slide_info = slide_info, gene_cutoff = 0)
  return(list(obj = obj, obj_raw = obj_raw))
}

# Vizu spot information
visualize_spot_data <- function(obj, obj_raw) {
  visualizeLabel(obj, "tissue")
  metadata(obj)$slide$total_counts <- Matrix::colSums(obj_raw)
  visualizeHeatmap(obj, "total_counts")
}

check_gene_heatmap <- function(obj, gene, viridis = FALSE) {
  visualizeHeatmap(obj, gene, viridis = viridis)
}

# IF NEEDED function to select highly variable or highly expressed genes
HiGenesKeep <- function(count_mat, top_high = 5000, mean_cutoff = 1, return_matrix = FALSE, verbose = TRUE) {
  mean_exp <- rowMeans(count_mat)
  top_genes <- rank(-mean_exp) <= top_high
  count_mat <- count_mat[top_genes, , drop = FALSE]
  mean_exp <- mean_exp[top_genes]
  high_exp_genes <- mean_exp >= mean_cutoff

  S_vf <- NormalizeData(CreateSeuratObject(count_mat), verbose = FALSE)
  S_vf <- FindVariableFeatures(S_vf, selection.method = "mvp", verbose = FALSE)

  if (as.integer(gsub("\\<(\\d+)\\.\\d+\\.\\d+", "\\1", S_vf@version)) >= 5) {
    high_variable_genes <- S_vf@assays$RNA@meta.data$vf_mvp_data_variable
  } else {
    high_variable_genes <- S_vf@assays$RNA@meta.features$mvp.variable
  }

  gene_tokeep <- high_variable_genes | high_exp_genes

  if (verbose) {
    message("Kept ", sum(gene_tokeep), " highly expressed or highly variable genes.")
  }

  if (return_matrix) {
    return(count_mat[gene_tokeep, , drop = FALSE])
  } else {
    return(names(which(gene_tokeep)))
  }
}

# Main decontamination function
decontaminate_spots <- function(obj, obj_raw, iterations = 10, candidate_radius = 20, seurat_v5 = TRUE, genes_keep = NULL) {
  if (seurat_v5) {
    decont_obj <- spotclean(obj, maxit = iterations, candidate_radius = candidate_radius)
  } else {
    decont_obj <- spotclean(obj, maxit = iterations, candidate_radius = candidate_radius, gene_keep = genes_keep)
  }
  return(decont_obj)
}

# Function to visualize contamination
visualize_contamination <- function(decont_obj) {
  visualizeHeatmap(decont_obj, metadata(decont_obj)$contamination_rate, logged = FALSE, 
                   legend_title = "contamination rate", legend_range = c(0, 1))
}

# Example usage
paths <- setup_environment(wd, sample_id)

# Load Visium data
visium_data <- load_visium_data(paths$samples, paths$spotfiles, paths$imgs, paths$json)
obj <- visium_data$obj
obj_raw <- visium_data$obj_raw

# Visualize spot data
visualize_spot_data(obj, obj_raw)

# Check specific gene heatmap
check_gene_heatmap(obj, "KRT14", viridis = FALSE)

# Calculate highly variable or highly expressed genes
genes_keep <- HiGenesKeep(obj_raw)

# Perform decontamination
decont_obj <- decontaminate_spots(obj, obj_raw, iterations = 10, candidate_radius = 20, seurat_v5 = TRUE)

# Visualize contamination rate
visualize_contamination(decont_obj)

# Visualize decontaminated gene
wrap_plots(check_gene_heatmap(obj, "KRT14", viridis = FALSE), check_gene_heatmap(decont_obj, "KRT14", viridis = FALSE))

# Convert SpotClean object to Seurat object and add metadata
convert_to_seurat <- function(decont_obj, spatial_path, sample_id, meta_data) {
  # Convert SpotClean object to Seurat object
  seurat_obj <- convertToSeurat(decont_obj, spatial_path)
  
  # Add contamination rate metadata
  seurat_obj$contamination_rate <- decont_obj@metadata$contamination_rate
  
  # Add original identifier metadata
  seurat_obj$orig.ident <- sample_id
  
  # Add custom metadata: labels and mitochondrial percentage
  labs <- meta_data$labels
  names(labs) <- meta_data$Barcode
  seurat_obj <- AddMetaData(seurat_obj, labs, col.name = "labels")
  
  pmt <- meta_data$percent.mito
  names(pmt) <- meta_data$Barcode
  seurat_obj <- AddMetaData(seurat_obj, pmt, col.name = "percent.mito")
  
  return(seurat_obj)
}

# Summarize label data
summarize_labels <- function(seurat_obj) {
  return(table(seurat_obj@meta.data$labels))
}

# Example usage
spatial_path <- file.path(wd, sample_id, "spatial/")
seurat_obj <- convert_to_seurat(decont_obj, spatial_path, sample_id, meta.preqc)

# Display the label summary
label_summary <- summarize_labels(seurat_obj)
print(label_summary)

#--- QC/filtering

# Function to plot QC metrics for a Seurat object
plot_qc_metrics <- function(seurat_obj, outdir = NULL, sampleID = NULL) {
  # Create spatial feature plots
  plot1 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial", pt.size.factor = 2) + theme(legend.position = "right")
  plot2 <- SpatialFeaturePlot(seurat_obj, features = "nFeature_Spatial", pt.size.factor = 2) + theme(legend.position = "right")
  plot3 <- SpatialFeaturePlot(seurat_obj, features = "percent.mito", pt.size.factor = 2) + theme(legend.position = "right")
  combined_plot <- wrap_plots(plot1, plot2, plot3)
  
  # Create violin plots
  plot4 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot5 <- VlnPlot(seurat_obj, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  plot6 <- VlnPlot(seurat_obj, features = "percent.mito", pt.size = 0.1) + NoLegend()
  combined_violin <- wrap_plots(plot4, plot5, plot6)

  if (!is.null(outdir) & !is.null(sampleID)) {
    # Save combined plots if outdir and sampleID provided
    ggsave(paste0(outdir, "/", sampleID, "_QC_spatial_plots.pdf"), combined_plot)
    ggsave(paste0(outdir, "/", sampleID, "_QC_violin_plots.pdf"), combined_violin)
  }
  
  return(list(spatial = combined_plot, violin = combined_violin))
}

# Function to filter data and generate QC distribution plots
# Generates distribution plots for each sample and dataset (features/counts) to observe 
# which values/how many values are outliers 
# Red   = lower threshold (-3*MAD)
# Blue  = upper threshold ( +3*MAD)
# Green = generic threshold (200 reads/genes)

filter_data_plots <- function(sr_ob, outdir, sampleID) {
  thresholds <- list(
    log_low_ncount_threshold = median(log(na.omit(sr_ob$nCount_Spatial))) - 2 * mad(log(na.omit(sr_ob$nCount_Spatial))),
    log_high_ncount_threshold = median(log(na.omit(sr_ob$nCount_Spatial))) + 2 * mad(log(na.omit(sr_ob$nCount_Spatial))),
    log_low_nfeature_threshold = median(log(na.omit(sr_ob$nFeature_Spatial))) - 2 * mad(log(na.omit(sr_ob$nFeature_Spatial))),
    log_high_nfeature_threshold = median(log(na.omit(sr_ob$nFeature_Spatial))) + 2 * mad(log(na.omit(sr_ob$nFeature_Spatial))),
    low_ncount_threshold = median(na.omit(sr_ob$nCount_Spatial)) - 2 * mad(na.omit(sr_ob$nCount_Spatial)),
    high_ncount_threshold = median(na.omit(sr_ob$nCount_Spatial)) + 2 * mad(na.omit(sr_ob$nCount_Spatial)),
    low_nfeature_threshold = median(na.omit(sr_ob$nFeature_Spatial)) - 2 * mad(na.omit(sr_ob$nFeature_Spatial)),
    high_nfeature_threshold = median(na.omit(sr_ob$nFeature_Spatial)) + 2 * mad(na.omit(sr_ob$nFeature_Spatial))
  )

  # Create QC distribution plots
  ncount_norm <- ggplot(sr_ob@meta.data, aes(x = nCount_Spatial)) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = thresholds$low_ncount_threshold), color = "red") +
    geom_vline(aes(xintercept = thresholds$high_ncount_threshold), color = "blue") +
    geom_vline(aes(xintercept = 200), color = "green")

  ncount_log <- ggplot(sr_ob@meta.data, aes(x = log(nCount_Spatial))) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = thresholds$log_low_ncount_threshold), color = "red") +
    geom_vline(aes(xintercept = thresholds$log_high_ncount_threshold), color = "blue") +
    geom_vline(aes(xintercept = log(200)), color = "green")

  nfeature_norm <- ggplot(sr_ob@meta.data, aes(x = nFeature_Spatial)) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = thresholds$low_nfeature_threshold), color = "red") +
    geom_vline(aes(xintercept = thresholds$high_nfeature_threshold), color = "blue") +
    geom_vline(aes(xintercept = 100), color = "green")

  nfeature_log <- ggplot(sr_ob@meta.data, aes(x = log(nFeature_Spatial))) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = thresholds$log_low_nfeature_threshold), color = "red") +
    geom_vline(aes(xintercept = thresholds$log_high_nfeature_threshold), color = "blue") +
    geom_vline(aes(xintercept = log(200)), color = "green")

  combined_qc_plots <- ncount_norm + ncount_log + nfeature_norm + nfeature_log +
    labs(caption = paste0("LRT: ", format(round(thresholds$low_ncount_threshold, 2), nsmall = 2), ", HRT: ", format(round(thresholds$high_ncount_threshold, 2), nsmall = 2), ", LFT: ", format(round(thresholds$low_nfeature_threshold, 2), nsmall = 2), ", HFT: ", format(round(thresholds$high_nfeature_threshold, 2), nsmall = 2), " / LOG -> LRT: ", format(round(thresholds$log_low_ncount_threshold, 2), nsmall = 2), ", HRT: ", format(round(thresholds$log_high_ncount_threshold, 2), nsmall = 2), ", LFT: ", format(round(thresholds$log_low_nfeature_threshold, 2), nsmall = 2), ", HFT: ", format(round(thresholds$log_high_nfeature_threshold, 2), nsmall = 2)))

  if (!is.null(outdir) & !is.null(sampleID)) {
    ggsave(paste0(outdir, "/", sampleID, "_QC_distribution_plots.pdf"), combined_qc_plots)
  }
  
  return(combined_qc_plots)
}

# Function to filter outlier data
filter_data <- function(sr_ob, outdir, sampleID) {
  # Select genes expressed in at least 3 spots
  selected_f <- rownames(sr_ob)[Matrix::rowSums(sr_ob) > 3]
  removed_f <- rownames(sr_ob)[Matrix::rowSums(sr_ob) <= 3]
  
  write.table(removed_f, file = paste0(outdir, "/", sampleID, "_removed_genes.txt"), sep = "\t", quote = FALSE, col.names = NA)
  
  # Subset data with selected genes
  data.filt <- subset(sr_ob, features = selected_f)
  
  # Identify outliers using MAD
  mad_reads.low <- isOutlier(data.filt$nCount_Spatial, nmads = 3, type = "lower")
  mad_reads.high <- isOutlier(data.filt$nCount_Spatial, nmads = 3, type = "higher")
  mad_gene.low <- isOutlier(data.filt$nFeature_Spatial, nmads = 3, type = "lower")
  mad_gene.high <- isOutlier(data.filt$nFeature_Spatial, nmads = 3, type = "higher")
  
  # Add metadata for outliers
  data.filt <- AddMetaData(data.filt, mad_reads.low | mad_gene.low, col.name = "low_outliers")
  data.filt <- AddMetaData(data.filt, mad_reads.high | mad_gene.high, col.name = "high_outliers")
  
  # Filter out poor quality spots based on feature counts
  filtered <- subset(data.filt, subset = nFeature_Spatial > 200 & percent.mito < 10)
  
  # Save QC plots
  plot2 <- SpatialDimPlot(filtered)
  plot3 <- SpatialDimPlot(data.filt)
  plot1 <- ggplot(data.frame(QCtype = c("RawTotalReads", "RawTotalGenes", "HighReadCounts", "LowReadCounts", "HighGene", "LowGene"), value = c(length(Cells(data.filt)), length(rownames(data.filt)), sum(mad_reads.high), sum(mad_reads.low), sum(mad_gene.high), sum(mad_gene.low))), aes(x = QCtype, y = value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = value), vjust = -0.3, size = 2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  if (!is.null(outdir)) {
    ggsave(paste0(outdir, "/", sampleID, "_QC_removed_values_nonlog.pdf"), plot1 | (plot2 / plot3))
  }
  
  return(filtered)
}

# Check outliers
check_outliers <- function(sr_ob) {
  p1 <- SpatialDimPlot(sr_ob, group.by = "low_outliers", pt.size.factor = 2)
  p2 <- SpatialDimPlot(sr_ob, group.by = "high_outliers", pt.size.factor = 2)
  p1 + p2
}

# Example execution for TIS4606_PP
outputdir <- "HS_Visium_new"
outputdir1 <- paste0(outputdir, "/Data_QC")
dir.create(outputdir1)

# Plot QC metrics
plot_qc_metrics(seurat_obj, outputdir, "TIS4606_PP")

#Remove debris spots
seurat_obj <- subset(seurat_obj, labels == "Selected")

# Filter data
obj_out <- filter_data(seurat_obj, outputdir1, "TIS4606_PP")

# Check outliers
check_outliers(obj_out)

# Visualize filtered data
VlnPlot(obj_out, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mito"), group.by = "orig.ident", pt.size = 0.1) & theme(axis.title.x = element_blank())
FeatureScatter(obj_out, "nCount_Spatial","nFeature_Spatial", pt.size = 0.5)

#--- Normalization and clustering

# Wrap ormalization, scaling, PCA, and UMAP
perform_normalization_and_clustering <- function(seurat_obj, outputdir, sampleID, npcs = 100) {
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, assay = "Spatial")

  # Find variable features, scale data, and run PCA
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = npcs)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:npcs)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:npcs)
  
  # Calculate variance explained by each PC
  total_variance <- seurat_obj@reductions$pca@misc$total.variance
  eigValues <- (seurat_obj[["pca"]]@stdev)^2
  varExplained <- eigValues / total_variance
  varExplained_cum <- cumsum(varExplained)
  
  # Variance information
  var_20pc <- sum(varExplained_cum <= 0.2)
  varpc_50PCA <- 100 * (varExplained_cum[50])
  message(paste0("The first 50 PCs explain ", round(varpc_50PCA), "% of the variance. 20% of the variance is explained by the first ", var_20pc, " PCs"))
  
  # Save variance plots
  variance_plots_dir <- file.path(outputdir, "Variance_Plots")
  dir.create(variance_plots_dir, showWarnings = FALSE)
  
  scree_plot <- ggplot(data.frame(PC = seq_along(varExplained), varExplained = varExplained), aes(x = PC, y = varExplained)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    ggtitle("Scree Plot") +
    ylab("Explained Variance")
  ggsave(filename = file.path(variance_plots_dir, paste0(sampleID, "_scree_plot.pdf")), plot = scree_plot)
  
  cumulative_variance_plot <- ggplot(data.frame(PC = seq_along(varExplained_cum), varExplained_cum = varExplained_cum), aes(x = PC, y = varExplained_cum)) +
    geom_point(size = 1) +
    theme_classic() +
    ggtitle("Cumulative Variance Explained") +
    xlab("PCs") +
    ylab("Cumulative Explained Variance") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = var_20pc, linetype = "dashed", color = "blue")
  ggsave(filename = file.path(variance_plots_dir, paste0(sampleID, "_cumulative_variance_plot.pdf")), plot = cumulative_variance_plot)
  
  ElbowPlot(seurat_obj, ndims = npcs, reduction = "pca") +
    theme_classic() +
    ggtitle("Elbow Plot") -> elbow_plot
  ggsave(filename = file.path(variance_plots_dir, paste0(sampleID, "_elbow_plot.pdf")), plot = elbow_plot)
  
  return(seurat_obj)
}

# Iterative clustering over a range of resolutions
perform_clustering <- function(seurat_obj, resolutions = seq(0, 1.6, by = 0.2), outputdir, sampleID) {
  # Directory for clustering results
  clustering_dir <- file.path(outputdir, "Clustering")
  dir.create(clustering_dir, showWarnings = FALSE)
  
  # Loop through resolutions and perform clustering
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  
  # Generate clustering tree plot to evaluate stability
  clustree_plot <- clustree(seurat_obj, prefix = "Spatial_snn_res.", node_colour = "sc3_stability")
  ggsave(filename = file.path(clustering_dir, paste0(sampleID, "_clustering_tree.pdf")), plot = clustree_plot)
  
  # Extract and save stability information
  stability <- clustree_plot$data[, c("Spatial_snn_res.", "sc3_stability")]
  write.table(stability, file = file.path(clustering_dir, paste0(sampleID, "_clustree_stability.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Calculate average stability for each resolution
  stability_ave <- aggregate(as.numeric(stability$sc3_stability), list(stability$Spatial_snn_res.), mean)
  colnames(stability_ave) <- c("Resolution", "Mean_Stability")
  stability_ave$Resolution <- as.numeric(as.character(stability_ave$Resolution))
  
  # Determine best resolution based on highest average stability
  best_res <- stability_ave$Resolution[which.max(stability_ave$Mean_Stability)]
  message(paste0("The best resolution based on stability is: ", best_res))
  
  # Save stability summary
  write.table(stability_ave, file = file.path(clustering_dir, paste0(sampleID, "_stability_summary.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  return(list(seurat_obj = seurat_obj, best_resolution = best_res))
}

# Example execution for TIS4606_PP
sampleID <- "TIS4606_PP"

# Perform normalization, PCA, and UMAP
obj_norm <- perform_normalization_and_clustering(obj_out, outputdir, sampleID)

# Perform clustering over a range of resolutions and determine best resolution
clustering_results <- perform_clustering(obj_norm, outputdir = outputdir, sampleID = sampleID)

# Visu
obj_norm <- FindClusters(obj_norm, resolution = clustering_results$best_resolution)
obj_norm <- RunUMAP(obj_norm, reduction = "pca", dims = 1:15)

SpatialDimPlot(obj_norm,pt.size.factor=3,image.alpha=0)

p1 <- DimPlot(obj_norm, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj_norm, label = TRUE, label.size = 3,pt.size.factor=2)
p1 + p2

#--- Identify spatially variable features

# Find spatially variable features using the 'moransi' method
obj_norm <- FindSpatiallyVariableFeatures(
  obj_norm, 
  assay = "Spatial", 
  features = VariableFeatures(obj_norm)[1:1000],
  selection.method = "moransi"
)

# Get the top spatially variable features
top_features <- head(SpatiallyVariableFeatures(obj_norm, selection.method = "moransi"), 6)

# Plot the top spatially variable features
SpatialFeaturePlot(
  obj_norm, 
  features = top_features, 
  ncol = 3, 
  alpha = c(0.1, 1), 
  pt.size = 1.5
)

# Plot specific features of interest
SpatialFeaturePlot(
  obj_norm, 
  features = c("IL17A", "IL17F", "TNF", "IL23A", "STAT4"), 
  ncol = 3, 
  pt.size = 2.5
)

SpatialFeaturePlot(
  obj_norm, 
  features = c("CXCL1", "FGR", "DOCK2"), 
  ncol = 3, 
  pt.size = 2
)

#--- OUTPUTS
output_dir <- "HS_Visium_new/Processed_Data"
dir.create(output_dir, showWarnings = FALSE)

rds_file <- file.path(output_dir, paste0(sampleID,"_preprocessed.rds"))
h5_file <- file.path(output_dir, paste0(sampleID,"_preprocessed.h5"))

# Save Seurat object to RDS file
saveRDS(obj_norm, file = rds_file)

# Save Seurat object to H5 file
# Seurat v4 and newer objects support saving as H5 files
SeuratDisk::SaveH5Seurat(obj_norm, filename = h5_file)

# Confirmation messages
message("Seurat object saved as RDS file: ", rds_file)
message("Seurat object saved as H5 file: ", h5_file)

