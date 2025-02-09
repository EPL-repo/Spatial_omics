#Akoya data in Giotto

#--- Setup
library(Giotto)

results_folder = "~/Desktop/codex_eg/"

#python_path = NULL 
#if(is.null(python_path)) {
#  installGiottoEnvironment()
#}

#--- Global instructions
instrs = createGiottoInstructions(show_plot = FALSE,
                                  save_plot = TRUE,
                                  save_dir = results_folder)
                                
expr_path = paste0(results_folder, "PP/expr.csv")
loc_path  = paste0(results_folder, "PP/coords.csv")
meta_path = paste0(results_folder, "PP/meta.csv")

#tile size to adjust coordinates
xtilespan = 1344;
ytilespan = 1008;

#--- Create Giotto object and pre-process
codex_expression = readExprMatrix(expr_path, transpose = T)
codex_locations  = data.table::fread(loc_path)
codex_metadata   = data.table::fread(meta_path)

stitch_file = stitchTileCoordinates(location_file = codex_metadata, Xtilespan = xtilespan, Ytilespan = ytilespan);
codex_locations = stitch_file[,.(Xcoord, Ycoord)]

codex_test <- createGiottoObject(raw_exprs = codex_expression, 
                                 spatial_locs = codex_locations,
                                 instructions = instrs,
                                 cell_metadata = codex_metadata)

#An object of class giotto 
#>Active spat_unit:  cell 
#>Active feat_type:  rna 
#[SUBCELLULAR INFO]
#[AGGREGATE INFO]
#expression -----------------------
#  [cell][rna] raw
#spatial locations ----------------
#  [cell] raw


#DG_subset = subsetGiottoLocs(codex_test,
#                             x_max = 28497136, x_min = 10497136,
#                             y_max = 36563620, y_min = 20063620)

DG_subset <- codex_test

## visualize
spatPlot(gobject = DG_subset,point_size = 0.1, 
         coord_fix_ratio = NULL,point_shape = 'no_border',
         save_param = list(save_name = '1_c_spatPlot_SUBSET'))
                  
## filter
DG_subset <- filterGiotto(gobject = DG_subset,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 2,
                           expression_values = c('raw'),
                           verbose = T)
#Number of cells removed:  0  out of  88210 (stardist) / 207681 (deepcell) 
#Number of feats removed:  0  out of  33 

DG_subset <- normalizeGiotto(gobject = DG_subset, scalefactor = 6000, verbose = T,
                              log_norm = FALSE,library_size_norm = FALSE,
                              scale_feats = FALSE, scale_cells = TRUE)

## add gene & cell statistics
DG_subset <- addStatistics(gobject = DG_subset,expression_values = "normalized")

#--- Dimension reduction
# PCA
DG_subset <- runPCA(gobject = DG_subset,
                     expression_values = 'normalized',
                     scale_unit = T,
                     method = "factominer")
signPCA(DG_subset,
        scale_unit = T,
        scree_ylim = c(0, 3),
        save_param = list(save_name = '2_a_spatPlot'))

plotPCA(gobject = DG_subset,
        point_shape = 'no_border',
        point_size = 0.2,
        save_param = list(save_name = '2_b_PCA'))

# UMAP
DG_subset <- runUMAP(DG_subset,
                      dimensions_to_use = 1:14,
                      n_components = 2,
                      n_threads = 12)
plotUMAP(gobject = DG_subset,
         point_shape = 'no_border',
         point_size = 0.2,
         save_param = list(save_name = '2_c_UMAP'))

## sNN network (default)
DG_subset <- createNearestNetwork(gobject = DG_subset, 
								   dimensions_to_use = 1:14, k = 20)

#DG_subset <- doLeidenCluster(gobject = DG_subset, resolution = 0.4, 
#							  n_iterations = 100, name = 'leiden')

#codex_metadata = pDataDT(DG_subset)
#leiden_colors = Giotto:::getDistinctColors(length(unique(codex_metadata$leiden)))
#names(leiden_colors) = unique(codex_metadata$leiden)

#plotUMAP(gobject = DG_subset, 
#         cell_color = 'leiden', point_shape = 'no_border', point_size = 0.2, 
#         cell_color_code = leiden_colors,
#         save_param = list(save_name = '4_a_UMAP'))

#spatPlot(gobject = DG_subset, cell_color = 'leiden', point_shape = 'no_border', point_size = 0.2, 
#         cell_color_code = leiden_colors, coord_fix_ratio = 1,label_size =2,
#         legend_text = 5,legend_symbol_size = 2,
#         save_param = list(save_name = '4_b_spatplot'))

#spatDimPlot2D(gobject = DG_subset, cell_color = 'leiden', spat_point_shape = 'no_border', 
#              spat_point_size = 0.2, dim_point_shape = 'no_border', dim_point_size = 0.2, 
#              cell_color_code = leiden_colors,plot_alignment = c("horizontal"),
#              save_param = list(save_name = '5_a_spatdimplot'))


#--- Integrate Enable clustering
codex_metadata <- pDataDT(DG_subset)
clust_colors = Giotto:::getDistinctColors(length(unique(codex_metadata$Merged.clustering)))
names(clust_colors) = unique(codex_metadata$Merged.clustering)

plotUMAP(gobject = DG_subset, 
         cell_color = 'Merged.clustering', point_shape = 'no_border', point_size = 0.2,
         cell_color_code = clust_colors,
         save_param = list(save_name = '4_a_UMAP'))

spatPlot(gobject = DG_subset, cell_color = 'Merged.clustering', 
		 point_shape = 'no_border', point_size = 0.2, 
         cell_color_code = clust_colors, coord_fix_ratio = 1,label_size =2,
         legend_text = 5,legend_symbol_size = 2,
         save_param = list(save_name = '4_b_spatplot'))


            
#--- DE
cluster_column = 'Merged.clustering'
markers_scran = findMarkers_one_vs_all(gobject=DG_subset, method="scran",
                                       expression_values="normalized", 
                                       cluster_column=cluster_column, min_genes=3)
markergenes_scran = unique(markers_scran[, head(.SD, 5), by="cluster"][["feats"]])

plotMetaDataHeatmap(DG_subset, expression_values = "normalized", metadata_cols = c(cluster_column), 
                    selected_genes = markergenes_scran,
                    y_text_size = 8, show_values = 'zscores_rescaled',
                    save_param = list(save_name = '6_a_metaheatmap'))

# gini
markers_gini = findMarkers_one_vs_all(gobject=DG_subset, method="gini", expression_values="normalized",
                                      cluster_column=cluster_column, min_genes=5)
markergenes_gini = unique(markers_gini[, head(.SD, 5), by="cluster"][["feats"]])
plotMetaDataHeatmap(DG_subset, expression_values = "normalized", 
                    metadata_cols = c(cluster_column), selected_genes = markergenes_gini,
                    show_values = 'zscores_rescaled',
                    save_param = list(save_name = '6_c_metaheatmap'))
#---

plotMetaDataHeatmap(DG_subset, expression_values = 'scaled',
                    metadata_cols = c('Merged.clustering'),y_text_size = 6,
                    save_param = list(save_name = '7_a_metaheatmap'))
                    
spatDimGenePlot(DG_subset, 
                expression_values = 'scaled',
                genes = c("MPO","CD14","Pan-Cytokeratin", 
                		  "TCF-1","CD141","Podoplanin"),
                spat_point_shape = 'no_border',
                dim_point_shape = 'no_border',
                cell_color_gradient = c("darkblue", "white", "red"),
                save_param = list(save_name = '8_f_spatdimgeneplot'))
                
#--- Save object
saveGiotto(gobject = DG_subset,
           dir = results_folder,
           foldername = '19h1257-1-PP_Giotto')
