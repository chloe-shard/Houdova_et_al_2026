############### Gene Expression vs Distance to FAP_PC Positive vs Negative Spots #############

#load Libraries
library(Seurat)

# Define the metadata column
metadata_column <- "FAP_PC_enriched_neg_vessel_thresholded" # Or FAP_PC_enriched_pos_vessel_thresholded

# Macrophage markers
Genesofinterest <- c("CD68", "ITGAM", "CD14", "IL1B", "IL6", "IL12A", "IL12B", "IL23A", "TNF", "CD80", "CD86", "HLADR", "NOS2", "IRF5", "STAT1", "SOCS3", "MRC1", "IRF4", "CD163", "CD200R1", "MARCO", "MSR1", "STAT6", "SOCS2", "CCL2", "CCL7", "CSF1", "CSF2", "HGF", "MIF", "CXCL12", "CXCL1", "IL4", "IL34", "IFNG", "IL8", "IL10", "IL13", "TGFB1", "RARRES2", "SPHK1")

# List of Seurat object names (vessel enriched)
seurat_names <- c("T248", "T255", "T259",
                  "T262", "GBM_ZH881T1", "GBM_ZH916T1",
                  "GBM_ZH1007nec", "GBM_ZH8811Abulk", "GBM_ZH8811Bbulk",
                  "GBM_ZH8812bulk", "GBM_ZH1007inf", "GBM_ZH1019inf")

# Initialize a list to store results
all_results <- list()

# Calculate distances and gene expression for genes of interest
calculate_distance_vs_gene_expression <- function(seurat_obj, metadata_column, Genesofinterest) { 
  # Check if metadata annotation column is present
  if (!metadata_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Metadata annotation column", metadata_column, "not found in the object"))
  }
  tryCatch({
    gene_data <- FetchData(object = seurat_obj, vars = Genesofinterest, slot = "data")
  }, error = function(e) {
    stop(paste("Error fetching data for genes of interest:", e$message))
  })
  available_genes <- intersect(Genesofinterest, colnames(gene_data))
  if (length(available_genes) == 0) {
    stop("None of the genes of interest are found")
  } else if (length(available_genes) < length(Genesofinterest)) {
    warning("Some genes of interest are not found and will be skipped.")
  }
  coords <- GetTissueCoordinates(seurat_obj)
  positive_metadata_spots <- WhichCells(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[metadata_column]] == "positive", ]))
  if (length(positive_metadata_spots) == 0) {
    stop("No positive metadata spots found")
  }
  distance_expression_df <- data.frame()
  for (spot in positive_metadata_spots) {
    spot_coords <- coords[spot, , drop = FALSE]
    
    distances <- sqrt(rowSums((t(t(coords) - as.numeric(spot_coords))) ^ 2))
    min_distance <- min(distances[distances > 0], na.rm = TRUE)
    normalized_distances <- distances / min_distance
    spot_data <- data.frame(SpotID = rownames(coords), Distance = normalized_distances)
    for (gene in available_genes) {
      spot_data[[gene]] <- gene_data[spot_data$SpotID, gene]
    }
    distance_expression_df <- rbind(distance_expression_df, spot_data)
  }
  max_distance <- max(distance_expression_df$Distance, na.rm = TRUE)
  distance_expression_df$DistanceBin <- cut(distance_expression_df$Distance, breaks = seq(1, max_distance, by = 1), include.lowest = TRUE)
  avg_expression <- aggregate(. ~ DistanceBin, data = distance_expression_df[, c("DistanceBin", available_genes)], mean, na.rm = TRUE)
  for (gene in available_genes) {
    min_val <- min(avg_expression[[gene]], na.rm = TRUE)
    max_val <- max(avg_expression[[gene]], na.rm = TRUE)
    avg_expression[[gene]] <- (avg_expression[[gene]] - min_val) / (max_val - min_val)
  }
  
  return(avg_expression)
}

# Loop through each Seurat object and process for each gene
for (obj_name in seurat_names) {
  cat("Processing", obj_name, "\n")
  seurat_obj <- get(obj_name)
  
  tryCatch({
    result <- calculate_distance_vs_gene_expression(seurat_obj, metadata_column, Genesofinterest)
    all_results[[obj_name]] <- result
  }, error = function(e) {
    message(paste("Error processing object", obj_name, ":", e$message))
  })
}

# Print final results
print(all_results)

# Initialize an empty dataframe to store merged results
merged_results <- data.frame()

# Loop through each Seurat object's result in all_results and bind to the merged dataframe
for (obj_name in names(all_results)) {
  object_data <- all_results[[obj_name]]
  object_data$Object <- obj_name
  if (nrow(merged_results) == 0) {
    merged_results <- object_data
  } else {
    common_columns <- union(names(merged_results), names(object_data))
    for (col in setdiff(common_columns, names(object_data))) {
      object_data[[col]] <- NA
    }
    for (col in setdiff(common_columns, names(merged_results))) {
      merged_results[[col]] <- NA
    }
    object_data <- object_data[, names(merged_results)]
    merged_results <- rbind(merged_results, object_data)
  }
}

# Convert DistanceBin to a factor
merged_results$DistanceBin <- as.factor(merged_results$DistanceBin)

# Identify available genes in merged_results
available_genes <- intersect(Genesofinterest, names(merged_results))

# Filter out rows with all NA values
filtered_results <- merged_results[!apply(is.na(merged_results[available_genes]), 1, all), ]

# Ensure DistanceBin has valid values
filtered_results$DistanceBin <- droplevels(as.factor(filtered_results$DistanceBin))

# Initialize a list to store results for each gene
average_results_list <- list()

# Loop over each gene in available_genes to calculate its mean
for (gene in available_genes) {
  gene_data <- filtered_results[!is.na(filtered_results[[gene]]), ]
  if (nrow(gene_data) > 0) {
    gene_avg <- aggregate(gene_data[[gene]] ~ gene_data$DistanceBin, FUN = mean)
    colnames(gene_avg) <- c("DistanceBin", paste0("Average_", gene))
    average_results_list[[gene]] <- gene_avg
  } else {
    cat("No valid data for gene:", gene, "\n")
  }
}

# Print the results for each gene
for (gene in names(average_results_list)) {
  cat("Averaged results for", gene, "across distance bins:\n")
  print(average_results_list[[gene]])
}

# Initialize an empty dataframe to store merged results
merged_results <- data.frame()

# Loop through each Seurat object's result in all_results and bind to the merged dataframe
for (obj_name in names(all_results)) {
  object_data <- all_results[[obj_name]]
  object_data$Object <- obj_name
  if (nrow(merged_results) == 0) {
    merged_results <- object_data
  } else {
    common_columns <- union(names(merged_results), names(object_data))
    for (col in setdiff(common_columns, names(object_data))) {
      object_data[[col]] <- NA
    }
    for (col in setdiff(common_columns, names(merged_results))) {
      merged_results[[col]] <- NA
    }
    object_data <- object_data[, names(merged_results)]
    merged_results <- rbind(merged_results, object_data)
  }
}

# Convert DistanceBin to a factor to ensure it's treated as categorical
merged_results$DistanceBin <- as.factor(merged_results$DistanceBin)

# Identify available genes in merged_results
available_genes <- intersect(Genesofinterest, names(merged_results))

# Filter out rows with all NA values in the available_genes columns
filtered_results <- merged_results[!apply(is.na(merged_results[available_genes]), 1, all), ]

# Ensure DistanceBin has valid values and re-level it to drop any unused levels
filtered_results$DistanceBin <- droplevels(as.factor(filtered_results$DistanceBin))

# Initialize a list to store results for each gene
average_results_list <- list()

# Loop over each gene in available_genes to calculate its mean separately
for (gene in available_genes) {
  gene_data <- filtered_results[!is.na(filtered_results[[gene]]), ]
  if (nrow(gene_data) > 0) {
    gene_avg <- aggregate(gene_data[[gene]] ~ gene_data$DistanceBin, FUN = mean)
    colnames(gene_avg) <- c("DistanceBin", paste0("Average_", gene))
    average_results_list[[gene]] <- gene_avg
  } else {
    cat("No valid data for gene:", gene, "\n")
  }
}

# Merge all averaged data into one matrix
merged_averages <- Reduce(function(x, y) merge(x, y, by = "DistanceBin", all = TRUE), average_results_list)

# Save results
write.csv(merged_averages, "Path to file /Merged_averages_distance_to_FAP_PC_neg_vessel.csv", row.names=TRUE)
