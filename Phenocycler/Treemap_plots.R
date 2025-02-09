setwd("Desktop/codex_eg/")

library(Giotto)
library(ggplot2)
library(treemapify)
library(dplyr)

giotto_obj = pp

cell_types <- meta$Merged.clustering
x <- meta$X.X
y <- meta$Y.Y

cell_type_counts <- table(cell_types)
cell_type_abundances <- cell_type_counts / sum(cell_type_counts)

cell_type_percentages <- cell_type_abundances * 100
cell_type_percentages <- round(cell_type_percentages)

# Convert abundances to a data frame
abundance_df <- data.frame(cell_type = names(cell_type_abundances),
						   abundance = cell_type_abundances)

abundance_df <- abundance_df %>%
arrange(desc(abundance.Freq))


color_mapping <- c(
     "Malignant A" = "#c67dc1",
     "Malignant B" = "#eeec9f",
     "Malignant C" = "#f90000",
     "Neutrophils" = "#199d77",
     "CD8+ T cells" = "#4dae4b",
     "Endothelial cells" = "#a65627",
     "Epithelial cells" = "#666666",
     "Exhausted T cells" = "#7fb1d3",
     "Plasma cells" = "#984ea3",
     "Possible B cells" = "#feb461",
     "Possible NK cells" = "#b3de68",
     "Helper T cells" = "#ffd82f",
     "Macrophages" = "#ff7f01",
     "mo-DCs" = "#e72a89",
     "Tregs" = "#377eb8"
)

treemap_plot <- ggplot(abundance_df, aes(area = abundance.Freq, fill = cell_type, label = cell_type)) +
     geom_treemap(size = 3.5, color = "white") +  # Increase line thickness and set line color to white
     geom_treemap_text(
         fontface = "italic",
         colour = "white",
         place = "centre",
         grow = TRUE
     ) +
     scale_fill_manual(values = color_mapping) +  # Set custom color mapping
     labs(title = "Cell Type Abundances") +
     theme(legend.position = "bottom")
 
# Show plot
print(treemap_plot)


