#----------------------------------------  NOTE -----------------------------------------#
# For basic VisiumHD data processing (post-SpaceRanger) using the 8um resolution as 
# recommended by 10x for first exploration.
# Exactly the same steps can be applied to the 16um resolution slot.
# Usage of the 2um data slot would be more interesting coupled with image segmentation to
# pick-up single-cell signal...
# ---------------------------------------------------------------------------------------#

library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
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
library(STutility)
library(magrittr)
library(SpotClean)
library(S4Vectors)

# Define paths and sample information
set_project_paths <- function(base_path, sample_id) {
  list(
    sample_h5 = file.path(base_path, sample_id, "binned_outputs/square_008um/raw_feature_bc_matrix.h5"),
    spotfile = file.path(base_path, sample_id, "binned_outputs/square_008um/tissue_positions.csv"),
    img = file.path(base_path, sample_id, "binned_outputs/square_008um/spatial/tissue_hires_image.png"),
    json = file.path(base_path, sample_id, "binned_outputs/square_008um/spatial/scalefactors_json.json")
  )
}

# Example usage on ADCA233 spheroid sample
project_base_path <- "cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PrecisionOncology/users/sandracurras/VisiumHD_data/240612_VH00207_47_AAFVTMCM5_spheroids/Output_spaceranger/VERSION2_json/"
sample_id <- "A1_ADCA233"
paths <- set_project_paths(project_base_path, sample_id)

#--------------------- Pre-treatment with STU - manual spot removal ---------------------#

# Create info table for input
infoTable <- data.frame(
  samples = paths$sample_h5,
  spotfiles = paths$spotfile,
  imgs = paths$img,
  json = paths$json
)

# Print the info table (for debugging purposes)
print(infoTable)

# Create the object
obj <- InputFromTable(infoTable, platform = "Visium", disable.subset = FALSE)

# Load images
obj <- LoadImages(obj, time.resolve = FALSE)

# Plot image
ImagePlot(obj, method = "raster", annotate = FALSE)

# Remove debris with manual annotation (opens Shiny app)
obj <- ManualAnnotation(obj)

# Feature overlay plot for labels
FeatureOverlay(obj, features = "labels")

# Calculate mitochondrial gene percentages
mt.genes <- grep("^MT-", rownames(obj), value = TRUE)
obj$percent.mito <- (Matrix::colSums(LayerData(obj, assay = "RNA", layer = "counts")[mt.genes, ]) /
                     Matrix::colSums(LayerData(obj, assay = "RNA", layer = "counts"))) * 100

# Pre-QC metadata processing
meta.preqc <- obj@meta.data
rownames(meta.preqc) <- gsub("_1", "", rownames(meta.preqc))
meta.preqc$Barcode <- rownames(meta.preqc)

# Subset and reorder metadata
meta.preqc <- meta.preqc[, c("Barcode", "percent.mito")]
meta.preqc <- meta.preqc[order(meta.preqc$Barcode), ]

# Display head of meta.preqc (for debugging)
head(meta.preqc)

#----------------------------- SpotClean bin decontamination ----------------------------#

# SpotClean decontamination workflow
run_spotclean <- function(base_path, sample_id, output_dir, max_iterations = 10, candidate_radius = 10) {
  # Generate paths for required files
  sample_paths <- set_project_paths(base_path, sample_id)
  
  # Load raw data
  obj_raw <- read10xRawH5(sample_paths$sample_h5)
  
  # Load slide information from tissue CSV, low-res image, and scale factors
  slide_info <- read10xSlide(
    tissue_csv_file   = sample_paths$spotfile,
    tissue_img_file   = file.path(base_path, sample_id, "binned_outputs/square_008um/spatial/tissue_lowres_image.png"),
    scale_factor_file = sample_paths$json
  )
  
  # Create SpotClean object
  obj <- createSlide(count_mat = obj_raw, slide_info = slide_info, gene_cutoff = 0)
  
  # Calculate total counts and visualize heatmap
  metadata(obj)$slide$total_counts <- Matrix::colSums(obj_raw)
  png(file.path(output_dir, paste0("SpotClean_", sample_id, "_totalCounts.png")))
  visualizeHeatmap(obj, "total_counts")
  dev.off()
  
  # Run SpotClean decontamination
  decont_obj <- spotclean(obj, maxit = max_iterations, candidate_radius = candidate_radius)
  
  # Visualize contamination rate
  png(file.path(output_dir, paste0("SpotClean_", sample_id, "_contaminationRate.png")))
  visualizeHeatmap(decont_obj, metadata(decont_obj)$contamination_rate, 
                   logged = FALSE, legend_title = "contamination rate", legend_range = c(0, 1))
  dev.off()
  
  # Visualize specific gene decontamination
  png(file.path(output_dir, paste0("SpotClean_", sample_id, "_decontKRT14.png")))
  wrap_plots(
    visualizeHeatmap(obj, "SERPINE1", viridis = FALSE),
    visualizeHeatmap(decont_obj, "KRT14", viridis = FALSE)
  )
  dev.off()
  
  # Convert to Seurat object for downstream analysis
  seurat_obj <- convertToSeurat(decont_obj, file.path(base_path, sample_id, "binned_outputs/square_008um/spatial/"))
  seurat_obj$contamination_rate <- decont_obj@metadata$contamination_rate
  seurat_obj$orig.ident <- paste0(sample_id, "_8um")
  
  return(seurat_obj)
}

# Example usage
# ! Added a tissue_positions.csv file (conversion from the parquet file in spatial/) for now as a workaround

output_dir <- "cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PrecisionOncology/users/Emilia/"
seurat_obj <- run_spotclean(project_base_path, sample_id, output_dir)

#--------------------------- IF skipping decontamination step ---------------------------#

# Load Seurat object and process without decontamination
process_spatial_data <- function(base_path, sample_id, meta_data, bin_size = 8) {
  # Load the spatial data
  seurat_obj <- Load10X_Spatial(data.dir = file.path(base_path, sample_id), bin.size = c(bin_size))
  
  # Add label metadata
  labs <- meta_data$labels
  names(labs) <- meta_data$Barcode
  seurat_obj <- AddMetaData(seurat_obj, labs, col.name = "labels")
  
  # Add mitochondrial percentage metadata
  pmt <- meta_data$percent.mito
  names(pmt) <- meta_data$Barcode
  seurat_obj <- AddMetaData(seurat_obj, pmt, col.name = "percent.mito")
  
  # Visualize metadata to verify
  head(seurat_obj@meta.data)
  
  # Set default assay to the appropriate spatial data
  DefaultAssay(seurat_obj) <- "Spatial.008um"
  
  # Plot spatial features
  plot1 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial.008um", pt.size.factor = 3, image.alpha = 0.5) + 
    theme(legend.position = "right")
  plot2 <- SpatialFeaturePlot(seurat_obj, features = "nFeature_Spatial.008um", pt.size.factor = 3, image.alpha = 0.5) + 
    theme(legend.position = "right")
  plot3 <- SpatialFeaturePlot(seurat_obj, features = "percent.mito", pt.size.factor = 3, image.alpha = 0.5) + 
    theme(legend.position = "right")
  
  wrap_plots(plot1, plot2, plot3)
  
  # Plot violin plots for the same features
  plot4 <- VlnPlot(seurat_obj, features = "nCount_Spatial.008um", pt.size = 0.1) + NoLegend()
  plot5 <- VlnPlot(seurat_obj, features = "nFeature_Spatial.008um", pt.size = 0.1) + NoLegend()
  plot6 <- VlnPlot(seurat_obj, features = "percent.mito", pt.size = 0.1) + NoLegend()
  
  wrap_plots(plot4, plot5, plot6)
  
  # Remove debris spots based on manual selection in metadata
  cleaned_obj <- subset(seurat_obj, labels == "Selected")
  
  # Print table of labels to check the result
  table(cleaned_obj@meta.data$labels)
  
  return(cleaned_obj)
}

# Example usage
cleaned_obj <- process_spatial_data(project_base_path, sample_id, meta.preqc)

#----------------------------------- QC and filtering -----------------------------------#

# Function to visualize QC plots and thresholds
filter_data_plots <- function(sr_ob, outdir, sampleID, spatial_col = "nCount_Spatial.008um", feature_col = "nFeature_Spatial.008um") {
  
  # Define thresholds using median and MAD
  log_low_ncount_threshold <- median(log(na.omit(sr_ob[[spatial_col]]))) - 2 * mad(log(na.omit(sr_ob[[spatial_col]])))
  log_high_ncount_threshold <- median(log(na.omit(sr_ob[[spatial_col]]))) + 2 * mad(log(na.omit(sr_ob[[spatial_col]])))
  log_low_nfeature_threshold <- median(log(na.omit(sr_ob[[feature_col]]))) - 2 * mad(log(na.omit(sr_ob[[feature_col]])))
  log_high_nfeature_threshold <- median(log(na.omit(sr_ob[[feature_col]]))) + 2 * mad(log(na.omit(sr_ob[[feature_col]])))
  
  # Thresholds without log transformation
  low_ncount_threshold <- median(na.omit(sr_ob[[spatial_col]])) - 2 * mad(na.omit(sr_ob[[spatial_col]]))
  high_ncount_threshold <- median(na.omit(sr_ob[[spatial_col]])) + 2 * mad(na.omit(sr_ob[[spatial_col]]))
  low_nfeature_threshold <- median(na.omit(sr_ob[[feature_col]])) - 2 * mad(na.omit(sr_ob[[feature_col]]))
  high_nfeature_threshold <- median(na.omit(sr_ob[[feature_col]])) + 2 * mad(na.omit(sr_ob[[feature_col]]))
  
  # Plot for normalized and log-transformed counts and features
  ncount_norm <- ggplot(sr_ob@meta.data, aes_string(x = spatial_col)) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = low_ncount_threshold), color = "red") +
    geom_vline(aes(xintercept = high_ncount_threshold), color = "blue") +
    geom_vline(aes(xintercept = 200), color = "green")
  
  ncount_log <- ggplot(sr_ob@meta.data, aes_string(x = paste0("log(", spatial_col, ")"))) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = log_low_ncount_threshold), color = "red") +
    geom_vline(aes(xintercept = log_high_ncount_threshold), color = "blue") +
    geom_vline(aes(xintercept = log(200)), color = "green")
  
  nfeature_norm <- ggplot(sr_ob@meta.data, aes_string(x = feature_col)) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = low_nfeature_threshold), color = "red") +
    geom_vline(aes(xintercept = high_nfeature_threshold), color = "blue") +
    geom_vline(aes(xintercept = 100), color = "green")
  
  nfeature_log <- ggplot(sr_ob@meta.data, aes_string(x = paste0("log(", feature_col, ")"))) +
    geom_histogram(aes(y = ..density..), bins = 200, fill = "grey") +
    geom_density() +
    geom_vline(aes(xintercept = log_low_nfeature_threshold), color = "red") +
    geom_vline(aes(xintercept = log_high_nfeature_threshold), color = "blue") +
    geom_vline(aes(xintercept = log(100)), color = "green")
  
  # Combine plots and add caption
  combined_plots <- ncount_norm + ncount_log + nfeature_norm + nfeature_log +
    labs(caption = paste0("LRT: ", format(round(low_ncount_threshold, 2), nsmall = 2), 
                          ", HRT: ", format(round(high_ncount_threshold, 2), nsmall = 2), 
                          ", LFT: ", format(round(low_nfeature_threshold, 2), nsmall = 2), 
                          ", HFT: ", format(round(high_nfeature_threshold, 2), nsmall = 2), 
                          "  /   LOG -> ", "LRT: ", format(round(log_low_ncount_threshold, 2), nsmall = 2), 
                          ", HRT: ", format(round(log_high_ncount_threshold, 2), nsmall = 2), 
                          ", LFT: ", format(round(log_low_nfeature_threshold, 2), nsmall = 2), 
                          ", HFT: ", format(round(log_high_nfeature_threshold, 2), nsmall = 2)))
  
  return(combined_plots)
}

# Function to filter data and identify outliers
filter_data <- function(sr_ob, outdir, sampleID, spatial_col = "nCount_Spatial.008um", feature_col = "nFeature_Spatial.008um") {
  
  # Filter genes expressed in at least 3 cells
  selected_f <- rownames(sr_ob)[Matrix::rowSums(sr_ob) > 3]
  removed_f <- rownames(sr_ob)[Matrix::rowSums(sr_ob) <= 3]
  
  write.table(removed_f, file = paste0(outdir, "/", sampleID, "_removed_genes.txt"), sep = "\t", quote = FALSE, col.names = NA)
  
  # Subset data for selected features
  data.filt <- subset(sr_ob, features = selected_f)
  
  # Identify outliers
  mad_reads.low <- isOutlier(data.filt[[spatial_col]], nmads = 3, type = "lower")
  mad_reads.high <- isOutlier(data.filt[[spatial_col]], nmads = 3, type = "higher")
  mad_gene.low <- isOutlier(data.filt[[feature_col]], nmads = 3, type = "lower")
  mad_gene.high <- isOutlier(data.filt[[feature_col]], nmads = 3, type = "higher")
  
  # Add outliers metadata
  data.filt <- AddMetaData(data.filt, mad_reads.low, col.name = "low_outliers")
  data.filt <- AddMetaData(data.filt, mad_reads.high, col.name = "high_outliers")
  data.filt <- AddMetaData(data.filt, mad_gene.low, col.name = "gene_low_outliers")
  data.filt <- AddMetaData(data.filt, mad_gene.high, col.name = "gene_high_outliers")
  
  # Create a QC summary
  qc_summary <- data.frame(QCtype = c("RawTotalReads", "RawTotalGenes", "HighReadCounts", "LowReadCounts", "HighGene", "LowGene", "Gene_In_<3_Cells"),
                           value = c(length(Cells(data.filt)), length(rownames(data.filt)), sum(mad_reads.high), sum(mad_reads.low), sum(mad_gene.high), sum(mad_gene.low), length(removed_f)))
  
  # Plot the QC summary
  plot1 <- ggplot(data = qc_summary, aes(x = QCtype, y = value)) + 
    geom_bar(stat = "identity", fill = "steelblue") + 
    geom_text(aes(label = value), vjust = -0.3, size = 2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  # Filter cells with more than 200 detected genes and less than 10% mitochondrial content
  filtered <- subset(data.filt, subset = nFeature_Spatial.008um > 200 & percent.mito < 10)
  
  plot2 <- SpatialDimPlot(filtered)
  plot3 <- SpatialDimPlot(data.filt)
  
  # Save QC plot
  ggsave(paste0(outdir, "/", sampleID, "_QC_removed_values_nonlog.pdf"))
  
  return(filtered)
}

# Example usage
outputdir <- sample_id
outputdir1 <- file.path(outputdir, "Data_QC")
dir.create(outputdir1)

obj_out <- filter_data(obj, outputdir1, sample_id)

# Function to check outliers visually
check_outliers <- function(sr_ob) {
  p1 <- SpatialDimPlot(sr_ob, group.by = "low_outliers", pt.size.factor = 1)
  p2 <- SpatialDimPlot(sr_ob, group.by = "high_outliers", pt.size.factor = 1)
  p3 <- SpatialDimPlot(sr_ob, group.by = "gene_low_outliers", pt.size.factor = 1)
  p4 <- SpatialDimPlot(sr_ob, group.by = "gene_high_outliers", pt.size.factor = 1)
  p1 + p2 + p3 + p4
}

# Check outliers for the filtered object
check_outliers(obj_out)

#Viz
VlnPlot(obj_out, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mito"), group.by = "orig.ident", pt.size = 0.1) & 
  theme(axis.title.x = element_blank())

FeatureScatter(obj_out, "nCount_Spatial.008um", "nFeature_Spatial.008um", pt.size = 0.5)

#----------------------------- Normalization and Clustering -----------------------------#

DefaultAssay(obj_out) <- "Spatial.008um"

# Normalize data
obj_norm <- NormalizeData(obj_out)

# Identify variable features and scale data
obj_norm <- FindVariableFeatures(obj_norm)
obj_norm <- ScaleData(obj_norm)

#--- Sketch Workflow for High-Dimensional Data

# Select 5,000 cells and create a new 'sketch' assay using Leverage Score method
obj_norm <- SketchData(
  object = obj_norm,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# Switch analysis to sketched cells
DefaultAssay(obj_norm) <- "sketch"

# Perform clustering on sketched cells
obj_norm <- FindVariableFeatures(obj_norm)
obj_norm <- ScaleData(obj_norm)

# Perform PCA and compute neighbors
obj_norm <- RunPCA(obj_norm, assay = "sketch", reduction.name = "pca.sketch")
obj_norm <- FindNeighbors(obj_norm, assay = "sketch", reduction = "pca.sketch", dims = 1:50)

# Calculate variance explained by each principal component (PC)
total_variance <- obj_norm@reductions$pca.sketch@misc$total.variance
eigValues <- (obj_norm[["pca.sketch"]]@stdev)^2
varExplained <- eigValues / total_variance
varExplained.cum <- cumsum(varExplained)

# Determine the number of PCs explaining 20% of variance and variance explained by first 50 PCs
var.20pc <- sum(varExplained.cum <= 0.2)
varpc.50PCA <- 100 * (varExplained.cum[50])

# Print the variance explained
print(paste0("The first 50 PCs explain ", round(varpc.50PCA), "% of the variance. 20% of the variance is explained by the first ", var.20pc, " PCs"))

# Visualization of variance explained

# Scree plot of explained variance by each PC
varExplained %>%
  enframe(name = "PC", value = "varExplained") %>%
  ggplot(aes(x = PC, y = varExplained)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ggtitle("Scree plot of explained variance") +
  ylab("Explained variance")

# Cumulative variance explained by PCs
ggplot(as.data.frame(varExplained.cum), aes(y = varExplained.cum, x = seq(1, length(varExplained.cum)))) +
  geom_point(size = 1) +
  theme_bw() +
  ggtitle("Cumulative variance explained by PCs") +
  xlab("PCs") +
  ylab("Cumulative explained variance") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = var.20pc, linetype = "dashed", color = "blue")

# Elbow plot to visualize the standard deviations of PCs
ElbowPlot(obj_norm, ndims = 50, reduction = "pca.sketch") +
  theme_bw() +
  ggtitle("Elbow plot of standard deviations of PCs")

# Create clustering output directory
int.outputdir <- file.path(outputdir, "Clustering")
dir.create(int.outputdir)

# Loop through resolutions to test clustering stability
resolutions <- seq(0.1, 3, by = 0.1)
df.2 <- FindClusters(obj_norm, resolution = 0) # Start with resolution = 0
df.2 <- FindClusters(df.2, resolution = resolutions)

# Clustering tree visualization using Clustree
clust <- clustree(na.omit(df.2@meta.data), prefix = "sketch_snn_res.", node_colour = "sc3_stability")
clust # Displays clustering at each resolution and how clusters evolve

# Extract and save stability data
stability <- clust$data[, c("sketch_snn_res.", "sc3_stability")]
write.table(stability, file = file.path(int.outputdir, "integrated_data_clustree_stability.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Calculate average stability for each resolution and determine the best resolution
stability.ave <- aggregate(as.numeric(stability$sc3_stability), list(stability$sketch_snn_res.), mean)
rownames(stability.ave) <- stability.ave$Group.1
stability.ave$Group.1 <- NULL
bestres <- as.numeric(rownames(stability.ave)[which.max(stability.ave$x)])
bestres # Best resolution

# Find clusters at the best resolution and run UMAP for visualization
obj_norm <- FindClusters(obj_norm, cluster.name = "seurat_cluster.sketched", resolution = bestres)
obj_norm <- RunUMAP(obj_norm, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:15)

# Visualize UMAP plot
SpatialDimPlot(obj_norm, pt.size.factor = 3, image.alpha = 0)

#--- Project sketched clustering onto the full dataset ---

obj_norm <- ProjectData(
  object = obj_norm,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# Switch back to full dataset
DefaultAssay(obj_norm) <- "Spatial"
Idents(obj_norm) <- "seurat_cluster.projected"

# Visual comparison of sketched clustering vs projected full dataset clustering
p1 <- DimPlot(obj_norm, reduction = "umap.sketch", label = FALSE) + 
  ggtitle("Sketched clustering (5,000 cells)") + 
  theme(legend.position = "bottom")

p2 <- DimPlot(obj_norm, reduction = "full.umap.sketch", label = FALSE) + 
  ggtitle("Projected clustering (full dataset)") + 
  theme(legend.position = "bottom")

p1 | p2 # Combine plots for visual comparison

# Additional detailed plots
p1 <- DimPlot(obj_norm, reduction = "full.umap.sketch", label = TRUE) + theme_bw()
p2 <- SpatialDimPlot(obj_norm, label = TRUE, label.size = 3, pt.size.factor = 3, image.alpha = 0.5)
p1 + p2 # Combined view of clustering and spatial distribution

#--------------------------------------- Markers ---------------------------------------#

# Create downsampled object to make visualization either
DefaultAssay(obj_norm) <- "Spatial"
Idents(obj_norm) <- "seurat_cluster.projected"
object_subset <- subset(obj_norm, cells = Cells(obj_norm[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

#Complete marker list
top_markers <- FindAllMarkers(obj_norm, only.pos = TRUE)
write.table(top_markers,"A1_ADCA233_8um_top_markers_unsupervised.txt",quote=F,sep="\t")

SCpubr::do_GroupwiseDEPlot(sample = obj_norm,
                           de_genes = tibble::tibble(top_markers),
                           top_genes = 5)

#-------------------------- Test label transfer for annotation --------------------------#

# This was tested for spheroid data (lung cancer) so that the scRNA-seq reference is a NSCLC dataset 
ref <- readRDS("NSCLC_GSE131907.rds")
colnames(ref@meta.data)[21] <- "Authors_cell_subtype"
colnames(ref@meta.data)[22] <- "Authors_cell_type_lvl1"
colnames(ref@meta.data)[23] <- "Authors_cell_type_lvl2"

Idents(ref) <- "Authors_cell_type_lvl1"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$Authors_cell_type_lvl1)
nUMI <- ref$nCount_RNA
cluster <- droplevels(cluster)

ref <- FindVariableFeatures(ref)

anchors <- FindTransferAnchors(reference = ref, query = obj_norm, 
							   normalization.method = "LogNormalize")

predictions.assay <- TransferData(anchorset = anchors, refdata = ref$Authors_cell_type_lvl1, 
								  prediction.assay = TRUE,
                                  weight.reduction = obj_norm[["full.pca.sketch"]], dims = 1:30)

obj_norm[["predictions"]] <- predictions.assay

DefaultAssay(obj_norm) <- "predictions"

SpatialFeaturePlot(obj_norm, features = c("Epithelial-cells", "Myeloid-cells","Fibroblasts", "Endothelial-cells", "Undetermined"),
				   pt.size.factor = 2, ncol = 3, 
				   crop = TRUE)

#---------------------------------------- Output ---------------------------------------#                        

saveRDS(obj_norm,paste0(sample_id,"_8um_preprocessed.rds"))



