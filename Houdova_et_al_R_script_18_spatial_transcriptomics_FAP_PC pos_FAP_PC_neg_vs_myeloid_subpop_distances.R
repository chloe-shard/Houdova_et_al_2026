############################################################################################################
#                   FAP_PC pos/ FAP_PC neg pericyte spots vs Myeloid subpop distances
############################################################################################################

# load libraries
library(Seurat)

#### Annotate FAP_PC negative pericyte spots (Pericyte positive / FAP_PC negative) ####

# List of Seurat object (vessel enriched)
spatial_object_names <- c(
  "T243", "T248", "T255", "T259", 
  "T260", "T262", "GBM_ZH881T1",
  "GBM_ZH916T1", "GBM_ZH1007nec", "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Loop over each Seurat object by name
for (obj_name in spatial_object_names) {
  obj <- get(obj_name)
  obj$FAP_PCs_neg_Pericyte_thresholded <- obj$Pericytes_broad_thresholded
  obj$FAP_PCs_neg_Pericyte_thresholded[obj$FAP_PCs_thresholded == "positive"] <- "negative"
  assign(obj_name, obj)
}

#### Myeloid cell subpop score vs distance to FAP_PCs or FAP_PCs_neg_Pericyte spots min max normalized ####

# Metadata coloumn of interest
metadata_column <- "FAP_PCs_thresholded"
# metadata_column <- "FAP_PCs_neg_Pericyte_thresholded"

# cell types of interest
cell_type_columns <- c("M1", "M2", "Monocyte", "Microglia")

# seurat objects (vessel enriched)
seurat_names <- c(
  "T243", "T248", "T255", "T259", 
  "T260", "T262", "GBM_ZH881T1",
  "GBM_ZH916T1", "GBM_ZH1007nec", "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Function to calculate distance-binned scores
calculate_distance_vs_score <- function(spatial_obj,
                                        metadata_column,
                                        cell_type_columns,
                                        positive_label = "positive",
                                        bin_width = 1) {
  if (!metadata_column %in% colnames(spatial_obj@meta.data)) {
    stop(paste("Metadata annotation column", metadata_column, "not found in the object"))
  }
  missing_cols <- setdiff(cell_type_columns, colnames(spatial_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing cell type score columns:", paste(missing_cols, collapse = ", ")))
  }
  coords <- GetTissueCoordinates(spatial_obj)
  common_spots <- intersect(rownames(coords), rownames(spatial_obj@meta.data))
  coords <- coords[common_spots, , drop = FALSE]
  meta <- spatial_obj@meta.data[common_spots, , drop = FALSE]
  positive_spots <- rownames(meta)[meta[[metadata_column]] == positive_label]
  if (length(positive_spots) == 0) {
    stop(paste("No", positive_label, "spots found in", metadata_column))
  }
  distance_score_list <- lapply(positive_spots, function(spot) {
    spot_coords <- coords[spot, , drop = FALSE]
    distances <- sqrt(rowSums((t(t(coords) - as.numeric(spot_coords)))^2))
    min_distance <- min(distances[distances > 0], na.rm = TRUE)
    normalized_distances <- distances / min_distance
    score_data <- data.frame(
      SpotID = rownames(coords),
      Distance = normalized_distances,
      meta[rownames(coords), cell_type_columns, drop = FALSE],
      stringsAsFactors = FALSE
    )
    return(score_data)
  })
  distance_score_df <- do.call(rbind, distance_score_list)
  max_distance <- ceiling(max(distance_score_df$Distance, na.rm = TRUE))
  distance_score_df$DistanceBin <- cut(
    distance_score_df$Distance,
    breaks = seq(1, max_distance + bin_width, by = bin_width),
    include.lowest = TRUE,
    right = FALSE
  )
  avg_scores <- aggregate(
    distance_score_df[, cell_type_columns, drop = FALSE],
    by = list(DistanceBin = distance_score_df$DistanceBin),
    FUN = mean,
    na.rm = TRUE
  )
  for (cell_type in cell_type_columns) {
    min_val <- min(avg_scores[[cell_type]], na.rm = TRUE)
    max_val <- max(avg_scores[[cell_type]], na.rm = TRUE)
    if (is.finite(min_val) && is.finite(max_val) && max_val != min_val) {
      avg_scores[[cell_type]] <- (avg_scores[[cell_type]] - min_val) / (max_val - min_val)
    } else {
      avg_scores[[cell_type]] <- NA
    }
  }
  return(avg_scores)
}


# Run across all tissues and produce final averaged output
run_distance_score_pipeline <- function(seurat_names,
                                        metadata_column,
                                        cell_type_columns,
                                        positive_label = "positive",
                                        bin_width = 1) {
  all_results <- list()
  for (obj_name in seurat_names) {
    message("Processing: ", obj_name)
    spatial_obj <- get(obj_name)
    result <- tryCatch({
      calculate_distance_vs_score(
        spatial_obj = spatial_obj,
        metadata_column = metadata_column,
        cell_type_columns = cell_type_columns,
        positive_label = positive_label,
        bin_width = bin_width
      )
    }, error = function(e) {
      message("Error processing ", obj_name, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(result)) {
      all_results[[obj_name]] <- result
    }
  }
  if (length(all_results) == 0) {
    stop("No Seurat objects were processed.")
  }
  cell_type_averages <- lapply(cell_type_columns, function(cell_type) {
    tissue_results <- lapply(names(all_results), function(obj_name) {
      object_data <- all_results[[obj_name]]
      if (!cell_type %in% colnames(object_data)) {
        return(NULL)
      }
      data.frame(
        Tissue = obj_name,
        DistanceBin = object_data$DistanceBin,
        Score = object_data[[cell_type]],
        stringsAsFactors = FALSE
      )
    })
    tissue_results <- do.call(rbind, tissue_results)
    averaged_data <- aggregate(
      Score ~ DistanceBin,
      data = tissue_results,
      FUN = mean,
      na.rm = TRUE
    )
    colnames(averaged_data)[2] <- paste0("Average_", cell_type)
    return(averaged_data)
  })
  names(cell_type_averages) <- cell_type_columns
  final_combined_averages <- Reduce(
    function(x, y) merge(x, y, by = "DistanceBin", all = TRUE),
    cell_type_averages
  )
  return(list(
    all_results = all_results,
    final_combined_averages = final_combined_averages
  ))
}

# Run pipeline
distance_score_output <- run_distance_score_pipeline(
  seurat_names = seurat_names,
  metadata_column = metadata_column,
  cell_type_columns = cell_type_columns,
  positive_label = "positive",
  bin_width = 1
)

# Extract final outputs
all_results <- distance_score_output$all_results
final_combined_averages <- distance_score_output$final_combined_averages

# print results
final_combined_averages

write.csv(
  final_combined_averages,
  file = "Myeloid_subpop_scores_vs_distance_FAP_PCs.csv",
  row.names = FALSE
)
