# -----------------------------------------------------------------------------------------
#                    Spatial Transcriptomics - Spot enrichment annotations
# -----------------------------------------------------------------------------------------

#Load libraries
library(Seurat)

#### Thresholding spot cell type scores to create binary enrichment annotations - Vessel enriched spots ####

# For this annotation, choose Option A in R file: Houdova_et_al_R_script_13_spatial_transcriptomics_RCTD_cell_type_annotation

# Define the cell types to combine into vascular score
vascular_cell_types <- c("Pericyte_broad", "Endothelial_broad")

# Function to create vascular score and threshold
process_seurat_object <- function(spatial_obj) {
  # Check that both required metadata columns exist
  missing_cols <- vascular_cell_types[!vascular_cell_types %in% colnames(spatial_obj@meta.data)]
  if (length(missing_cols) > 0) {
    warning(
      paste(
        "Missing columns:",
        paste(missing_cols, collapse = ", "),
        "- skipping this object."
      )
    )
    return(spatial_obj)
  }
  pericyte_values <- as.numeric(as.character(spatial_obj@meta.data$Pericyte_broad))
  endothelial_values <- as.numeric(as.character(spatial_obj@meta.data$Endothelial_broad))
  if (all(is.na(pericyte_values)) & all(is.na(endothelial_values))) {
    warning("Pericyte_broad and Endothelial_broad have no valid numeric values. Skipping this object.")
    return(spatial_obj)
  }
  pericyte_values[is.na(pericyte_values)] <- 0
  endothelial_values[is.na(endothelial_values)] <- 0
  spatial_obj@meta.data$Vascular_broad <- pericyte_values + endothelial_values
  vascular_values <- spatial_obj@meta.data$Vascular_broad
  mean_value <- mean(vascular_values, na.rm = TRUE)
  sd_value <- sd(vascular_values, na.rm = TRUE)
  threshold <- mean_value + sd_value
  spatial_obj@meta.data$Vascular_broad_thresholded <- ifelse(
    vascular_values > threshold,
    "positive",
    "negative"
  )
  return(spatial_obj)
}

# List of all Seurat objects
seurat_names <- c("T243", "T241", "T242", "T248", 
                  "T251", "T255", "T256", "T259", 
                  "T260", "T262", "T265", "T266", 
                  "T268", "T269", "T270", "T275", 
                  "T296", "T304", "T313", "T334",
                  "GBM_ZH881T1", "GBM_ZH916bulk", 
                  "GBM_ZH916T1", "GBM_ZH1007inf", 
                  "GBM_ZH1007nec", "GBM_ZH8812bulk",
                  "GBM_ZH1019inf", "GBM_ZH8811Abulk", 
                  "GBM_ZH8811Bbulk"
)

# Apply the function to each Seurat object
for (seurat_name in seurat_names) {
  message("Processing: ", seurat_name)
  spatial_obj <- get(seurat_name)
  processed_object <- process_seurat_object(spatial_obj)
  assign(seurat_name, processed_object)
}

#### Identify vessel enriched tissues - i.e. high Vascular_broad_thresholded aggregation score ####

# Specify the Seurat objects to analyse
seurat_names <- c("T243", "T241", "T242", "T248", 
                  "T251", "T255", "T256", "T259", 
                  "T260", "T262", "T265", "T266", 
                  "T268", "T269", "T270", "T275", 
                  "T296", "T304", "T313", "T334",
                  "GBM_ZH881T1", "GBM_ZH916bulk", 
                  "GBM_ZH916T1", "GBM_ZH1007inf", 
                  "GBM_ZH1007nec", "GBM_ZH8812bulk",
                  "GBM_ZH1019inf", "GBM_ZH8811Abulk", 
                  "GBM_ZH8811Bbulk"
)

# Specify the thresholded cell type / annotation to analyse
thresholded_cell_type <- "Vascular_broad_thresholded"

# Initialise dataframe to store results
results <- data.frame(
  SeuratObject = character(),
  AverageFractionPositive = numeric(),
  stringsAsFactors = FALSE
)

# Function to process each Seurat object
process_spatial_object <- function(spatial_obj, cell_type, obj_name = "unknown_object") {
  if (!cell_type %in% colnames(spatial_obj@meta.data)) {
    warning(paste("Annotation", cell_type, "not found in object:", obj_name))
    return(NA)
  }
  coords <- GetTissueCoordinates(spatial_obj)
  common_spots <- intersect(rownames(coords), rownames(spatial_obj@meta.data))
  if (length(common_spots) == 0) {
    warning(paste("No matching spot IDs between coordinates and metadata in object:", obj_name))
    return(NA)
  }
  coords <- coords[common_spots, , drop = FALSE]
  meta_data <- spatial_obj@meta.data[common_spots, , drop = FALSE]
  positive_spots <- rownames(meta_data)[meta_data[[cell_type]] == "positive"]
  if (length(positive_spots) == 0) {
    warning(paste("No positive spots found for", cell_type, "in object:", obj_name))
    return(NA)
  }
  if (nrow(coords) < 8) {
    warning(paste("Fewer than 8 spots available in object:", obj_name, "- cannot calculate 6-neighbour threshold."))
    return(NA)
  }
  calculate_distances <- function(spot_coords, all_coords) {
    spot_vector <- as.numeric(spot_coords)
    distances <- sqrt(rowSums((t(t(all_coords) - spot_vector)) ^ 2))
    return(distances)
  }
  distance_matrices <- list()
  for (spot in positive_spots) {
    spot_coords <- coords[spot, , drop = FALSE]
    distances <- calculate_distances(spot_coords, coords)
    distance_df <- data.frame(
      SpotID = rownames(coords),
      Distance = distances,
      stringsAsFactors = FALSE
    )
    distance_df <- distance_df[order(distance_df$Distance), ]
    distance_matrices[[spot]] <- distance_df
  }
  sixth_distances <- sapply(distance_matrices, function(df) df$Distance[7])
  seventh_distances <- sapply(distance_matrices, function(df) df$Distance[8])
  distance_threshold <- mean((sixth_distances + seventh_distances) / 2, na.rm = TRUE)
  if (is.na(distance_threshold)) {
    warning(paste("Distance threshold could not be calculated for object:", obj_name))
    return(NA)
  }
  neighbor_matrix <- matrix(
    NA,
    nrow = length(positive_spots),
    ncol = 6
  )
  rownames(neighbor_matrix) <- positive_spots
  colnames(neighbor_matrix) <- paste0("Neighbor_", 1:6)
  for (i in seq_along(positive_spots)) {
    spot <- positive_spots[i]
    distance_df <- distance_matrices[[spot]]
    neighbors <- distance_df$SpotID[
      distance_df$Distance <= distance_threshold &
        distance_df$SpotID != spot
    ]
    if (length(neighbors) > 0) {
      neighbor_matrix[i, 1:min(length(neighbors), 6)] <- neighbors[1:min(length(neighbors), 6)]
    }
  }
  neighbor_matrix_annotated <- neighbor_matrix
  for (i in 1:nrow(neighbor_matrix)) {
    for (j in 1:ncol(neighbor_matrix)) {
      spot_id <- neighbor_matrix[i, j]
      if (!is.na(spot_id)) {
        neighbor_matrix_annotated[i, j] <- meta_data[spot_id, cell_type]
      }
    }
  }
  positive_fractions <- apply(neighbor_matrix_annotated, 1, function(row) {
    positive_count <- sum(row == "positive", na.rm = TRUE)
    total_count <- sum(!is.na(row))
    if (total_count == 0) {
      return(NA)
    }
    fraction_positive <- positive_count / total_count
    return(fraction_positive)
  })
  average_fraction_positive <- mean(positive_fractions, na.rm = TRUE)
  return(average_fraction_positive)
}

# Process each Seurat object
for (seurat_name in seurat_names) {
  message("Processing: ", seurat_name)
  if (!exists(seurat_name)) {
    warning(paste("Object", seurat_name, "does not exist - skipping."))
    next
  }
  spatial_obj <- get(seurat_name)
  avg_fraction_positive <- process_spatial_object(
    spatial_obj = spatial_obj,
    cell_type = thresholded_cell_type,
    obj_name = seurat_name
  )
  if (!is.na(avg_fraction_positive)) {
    results <- rbind(
      results,
      data.frame(
        SeuratObject = seurat_name,
        AverageFractionPositive = avg_fraction_positive,
        stringsAsFactors = FALSE
      )
    )
  }
}

# print results
print(results)

#### Thresholding spot cell type scores to create binary enrichment annotations - Pericyte enriched spots ####

# For this annotation, choose Option A in R file: Houdova_et_al_R_script_13_spatial_transcriptomics_RCTD_cell_type_annotation

# Define the cell types
cell_types <- c("Pericyte_broad")

# Function to threshold cell type annotations
process_seurat_object <- function(spatial_obj) {
  for (cell_type in cell_types) {
    numeric_values <- as.numeric(as.character(spatial_obj@meta.data[[cell_type]]))
    if (all(is.na(numeric_values))) {
      warning(paste("Column", cell_type, "has no valid numeric values. Skipping."))
      next
    }
    mean_value <- mean(numeric_values, na.rm = TRUE)
    sd_value <- sd(numeric_values, na.rm = TRUE)
    threshold <- mean_value + sd_value
    new_annotation <- paste0(cell_type, "_thresholded")
    spatial_obj@meta.data[[new_annotation]] <- ifelse(numeric_values > threshold, "positive", "negative")
  }
  return(spatial_obj)
}

# List of all Seurat objects
seurat_names <- c("T243", "T241", "T242", "T248", 
                  "T251", "T255", "T256", "T259", 
                  "T260", "T262", "T265", "T266", 
                  "T268", "T269", "T270", "T275", 
                  "T296", "T304", "T313", "T334",
                  "GBM_ZH881T1", "GBM_ZH916bulk", 
                  "GBM_ZH916T1", "GBM_ZH1007inf", 
                  "GBM_ZH1007nec", "GBM_ZH8812bulk",
                  "GBM_ZH1019inf", "GBM_ZH8811Abulk", 
                  "GBM_ZH8811Bbulk"
)

# Apply the function to each Seurat object
for (seurat_name in seurat_names) {
  # Get the Seurat object by name
  spatial_obj <- get(seurat_name)
  processed_object <- process_seurat_object(spatial_obj)
  assign(seurat_name, processed_object)
}

#### Thresholding spot cell type scores to create binary enrichment annotations - Pericyte and Myeloid subpopulation enriched spots ####

# For this annotation, choose Option C in R file: Houdova_et_al_R_script_13_spatial_transcriptomics_RCTD_cell_type_annotation

# Define the cell types
cell_types <- c("Proliferative_PCs", "ECM_PCs", "FAP_PCs", "SMCs", "Interferon_PCs", "Transport_PCs", "M1", "M2", "Monocyte", "Microglia")

# Function to threshold cell type annotations
process_seurat_object <- function(spatial_obj) {
  for (cell_type in cell_types) {
    numeric_values <- as.numeric(as.character(spatial_obj@meta.data[[cell_type]]))
    if (all(is.na(numeric_values))) {
      warning(paste("Column", cell_type, "has no valid numeric values. Skipping."))
      next
    }
    mean_value <- mean(numeric_values, na.rm = TRUE)
    sd_value <- sd(numeric_values, na.rm = TRUE)
    threshold <- mean_value + sd_value
    new_annotation <- paste0(cell_type, "_thresholded")
    spatial_obj@meta.data[[new_annotation]] <- ifelse(numeric_values > threshold, "positive", "negative")
  }
  return(spatial_obj)
}

# List of all Seurat objects
seurat_names <- c("T243", "T241", "T242", "T248", 
                  "T251", "T255", "T256", "T259", 
                  "T260", "T262", "T265", "T266", 
                  "T268", "T269", "T270", "T275", 
                  "T296", "T304", "T313", "T334",
                  "GBM_ZH881T1", "GBM_ZH916bulk", 
                  "GBM_ZH916T1", "GBM_ZH1007inf", 
                  "GBM_ZH1007nec", "GBM_ZH8812bulk",
                  "GBM_ZH1019inf", "GBM_ZH8811Abulk", 
                  "GBM_ZH8811Bbulk"
)

# Apply the function to each Seurat object
for (seurat_name in seurat_names) {
  # Get the Seurat object by name
  spatial_obj <- get(seurat_name)
  processed_object <- process_seurat_object(spatial_obj)
  assign(seurat_name, processed_object)
}
