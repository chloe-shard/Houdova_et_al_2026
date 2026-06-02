# -----------------------------------------------------------------------------------------
#                           Neighbourhood composition score analysis
# -----------------------------------------------------------------------------------------

# load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

#### Neighborhood composition score - including positive spot (out of 7 neighbors) #####

# Seurat object names (Spatial T objects annotated with cell subtypes using RCTD analysis - see R file: Houdova_et_al_R_script_12_spatial_transcriptomics_RCTD_cell_type_annotation)
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Specify the thresholded cell type to include
thresholded_cell_type <- "FAP_PCs_thresholded" #swap for the other pericyte sub-populations e.g. "ECM_PCs_thresholded", "SMCs_thresholded" etc
neighbor_thresholded_cell_types <- c("M1_thresholded", "M2_thresholded", "Monocyte_thresholded", "Microglia_thresholded")

# Initialize a dataframe to store the results
results <- data.frame(SeuratObject = character(), NeighborCellType = character(), AverageFractionPositive = numeric(), stringsAsFactors = FALSE)

# Function to process each Seurat object
process_spatial_object <- function(spatial_obj, cell_type, neighbor_cell_type) {
  if (!cell_type %in% colnames(spatial_obj@meta.data)) {
    warning(paste("Annotation", cell_type, "not found in object:", deparse(substitute(spatial_obj))))
    return(NA)
  }
  if (!neighbor_cell_type %in% colnames(spatial_obj@meta.data)) {
    warning(paste("Neighbor annotation", neighbor_cell_type, "not found in object:", deparse(substitute(spatial_obj))))
    return(NA)
  }
  coords <- GetTissueCoordinates(spatial_obj)
  positive_spots <- WhichCells(spatial_obj, cells = rownames(spatial_obj@meta.data[spatial_obj@meta.data[[cell_type]] == "positive",]))
  if (length(positive_spots) == 0) {
    warning(paste("No positive spots found for cell type", cell_type, "in object:", deparse(substitute(spatial_obj))))
    return(NA)
  }
  positive_coords <- coords[positive_spots, ]
  distance_matrices <- list()
  calculate_distances <- function(spot_coords, all_coords) {
    spot_vector <- as.numeric(spot_coords)
    distances <- sqrt(rowSums((t(t(all_coords) - spot_vector)) ^ 2))
    return(distances)
  }
  for (spot in positive_spots) {
    spot_coords <- coords[spot, , drop = FALSE]
    distances <- calculate_distances(spot_coords, coords)
    distance_df <- data.frame(SpotID = rownames(coords), Distance = distances)
    distance_df <- distance_df[order(distance_df$Distance), ]
    distance_matrices[[spot]] <- distance_df
  }
  neighbor_matrix <- matrix(NA, nrow = length(positive_spots), ncol = 7)
  rownames(neighbor_matrix) <- positive_spots
  colnames(neighbor_matrix) <- paste0("Neighbor_", 1:7)
  for (i in seq_along(positive_spots)) {
    spot <- positive_spots[i]
    distance_df <- distance_matrices[[spot]]
    neighbor_matrix[i, 1:7] <- distance_df$SpotID[1:7]
  }
  neighbor_matrix_annotated <- neighbor_matrix
  for (i in 1:nrow(neighbor_matrix)) {
    for (j in 1:ncol(neighbor_matrix)) {
      spot_id <- neighbor_matrix[i, j]
      if (!is.na(spot_id)) {
        neighbor_matrix_annotated[i, j] <- spatial_obj@meta.data[spot_id, neighbor_cell_type]
      }
    }
  }
  positive_fractions <- apply(neighbor_matrix_annotated, 1, function(row) {
    positive_count <- sum(row == "positive", na.rm = TRUE)
    total_count <- sum(!is.na(row))
    fraction_positive <- positive_count / total_count
    return(fraction_positive)
  })
  average_fraction_positive <- mean(positive_fractions, na.rm = TRUE)
  return(average_fraction_positive)
}

# Process each Seurat object and store the results for each neighbor cell type
for (obj_name in names(seurat_objects)) {
  spatial_obj <- seurat_objects[[obj_name]]
  for (neighbor_cell_type in neighbor_thresholded_cell_types) {
    if (thresholded_cell_type %in% colnames(spatial_obj@meta.data)) {
      if (neighbor_cell_type %in% colnames(spatial_obj@meta.data)) {
        avg_fraction_positive <- process_spatial_object(spatial_obj, thresholded_cell_type, neighbor_cell_type)
        if (!is.na(avg_fraction_positive)) {
          results <- rbind(results, data.frame(SeuratObject = obj_name, NeighborCellType = neighbor_cell_type, AverageFractionPositive = avg_fraction_positive))
        }
      } else {
        warning(paste("Neighbor annotation", neighbor_cell_type, "not found in object:", obj_name))
      }
    } else {
      warning(paste("Annotation", thresholded_cell_type, "not found in object:", obj_name))
    }
  }
}

# Order the results by AverageFractionPositive
results <- results[order(-results$AverageFractionPositive), ]

# Print the results
print(results)

# clean up myeloid population names for plotting
results_plot <- results %>%
  mutate(
    MyeloidPopulation = gsub("_thresholded", "", NeighborCellType),
    MyeloidPopulation = factor(
      MyeloidPopulation,
      levels = c("M1", "M2", "Monocyte", "Microglia")
    )
  )

# Calculate mean neighbourhood composition score per myeloid population
mean_results <- results_plot %>%
  group_by(MyeloidPopulation) %>%
  summarise(
    Mean_NCS = mean(AverageFractionPositive, na.rm = TRUE),
    .groups = "drop"
  )

# Plot individual tissues + mean per myeloid population
ggplot(results_plot, aes(
  x = MyeloidPopulation,
  y = AverageFractionPositive
)) +
  geom_jitter(
    width = 0.15,
    size = 3,
    alpha = 0.75,
    aes(colour = SeuratObject)
  ) +
  geom_point(
    data = mean_results,
    aes(
      x = MyeloidPopulation,
      y = Mean_NCS
    ),
    inherit.aes = FALSE,
    size = 5,
    shape = 18,
    colour = "black"
  ) +
  geom_crossbar(
    data = mean_results,
    aes(
      x = MyeloidPopulation,
      y = Mean_NCS,
      ymin = Mean_NCS,
      ymax = Mean_NCS
    ),
    inherit.aes = FALSE,
    width = 0.45,
    colour = "black",
    linewidth = 0.6
  ) +
  labs(
    x = "Myeloid population",
    y = "Neighbourhood composition score",
    colour = "Tissue",
    title = "Myeloid neighbourhood composition around FAP+ pericyte spots"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Prepare results for statistics
results_stats <- results %>%
  mutate(
    MyeloidPopulation = gsub("_thresholded", "", NeighborCellType),
    MyeloidPopulation = factor(
      MyeloidPopulation,
      levels = c("M1", "M2", "Monocyte", "Microglia")
    )
  )

# One-way ANOVA
anova_model <- aov(
  AverageFractionPositive ~ MyeloidPopulation,
  data = results_stats
)

# Print ANOVA results
summary(anova_model)

# Tukey post-hoc test
tukey_results <- TukeyHSD(anova_model)

# Print Tukey results
print(tukey_results)

# Convert Tukey results into a dataframe
tukey_df <- as.data.frame(tukey_results$MyeloidPopulation)
tukey_df$comparison <- rownames(tukey_df)

# Reorder columns
tukey_df <- tukey_df %>%
  select(
    comparison,
    diff,
    lwr,
    upr,
    `p adj`
  )

print(tukey_df)
