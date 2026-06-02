# -----------------------------------------------------------------------------------------------------------
#                 Spatial Transcriptomics - FAP_PC enriched spots vs FAP_PC depleted vessel spots
# -----------------------------------------------------------------------------------------------------------

# load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(svglite)

##### annotate FAP+ PC depleted vascular spots #####

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
  obj$FAP_PCs_neg_vessel_thresholded <- obj$Vascular_broad_thresholded
  obj$FAP_PCs_neg_vessel_thresholded[obj$FAP_PCs_thresholded == "positive"] <- "negative"
  assign(obj_name, obj)
}

#### Plot both FAP_PC enriched spots and FAP_PCs depleted vessel spots ####

# List of Seurat object names as strings
spatial_object_names <- c(
  "T243", "T248", "T255", "T259", 
  "T260", "T262", "GBM_ZH881T1",
  "GBM_ZH916T1", "GBM_ZH1007nec", "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Define colors for each category
colors <- c(
  "Both_Positive" = "red", 
  "FAP_PCs_Positive" = "red", 
  "FAP_PCs_neg_vessel_Positive" = "blue", 
  "Both_Negative" = "grey"
)

# Loop over each Seurat object
for (obj_name in spatial_object_names) {
  message("Processing: ", obj_name)
  obj <- get(obj_name)
  required_cols <- c(
    "FAP_PCs_thresholded",
    "FAP_PCs_neg_vessel_thresholded"
  )
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    warning(
      "Skipping ", obj_name,
      " because missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
    next
  }
  combined_annotation <- with(
    obj@meta.data,
    ifelse(
      FAP_PCs_thresholded == "positive" & FAP_PCs_neg_vessel_thresholded == "positive",
      "Both_Positive",
      ifelse(
        FAP_PCs_thresholded == "positive" & FAP_PCs_neg_vessel_thresholded == "negative",
        "FAP_PCs_Positive",
        ifelse(
          FAP_PCs_thresholded == "negative" & FAP_PCs_neg_vessel_thresholded == "positive",
          "FAP_PCs_neg_vessel_Positive",
          "Both_Negative"
        )
      )
    )
  )
  # Add directly to metadata
  obj@meta.data$combined_annotation <- factor(
    combined_annotation,
    levels = c(
      "Both_Positive",
      "FAP_PCs_Positive",
      "FAP_PCs_neg_vessel_Positive",
      "Both_Negative"
    )
  )
  # Assign modified object back to global environment
  assign(obj_name, obj)
}

# Plot the SpatialDimPlot with the combined annotation and increased point size
SpatialDimPlot(GBM_ZH8811Abulk, group.by = "combined_annotation", cols = colors, pt.size.factor = 3.2) +
  ggtitle("SpatialDimPlot with FAP_PCs and FAP_PCs_neg_vessel Annotations")

#### Myeloid marker expression in FAP+ Per enriched vs depleted vascular neighborhoods (Log2 Fold Change) ####

# Define the two metadata columns and gene of interest
metadata_column1 <- "FAP_PCs_thresholded"      # First metadata column of interest
metadata_column2 <- "FAP_PCs_neg_vessel_thresholded" # Second metadata column of interest
genes_of_interest <- c("MRC1", "MSR1", "CD163")

# List of Seurat object (vessel enriched)
seurat_names <- c("T243", "T248", "T255",
                  "T259", "T260", "T262",
                  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
                  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
                  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Initialize a list to store results per Seurat object and gene
all_results <- list()

# Function to calculate gene expression in neighboring spots for two metadata columns
calculate_fold_change_expression <- function(seurat_obj, obj_name, metadata_column1, metadata_column2, gene) {
  if (!all(c(metadata_column1, metadata_column2) %in% colnames(seurat_obj@meta.data))) {
    stop(paste("One or both metadata annotation columns not found in the object"))
  }
  gene_data <- FetchData(object = seurat_obj, vars = gene, layer = "data")
  non_zero_spots <- gene_data[gene_data[[gene]] > 0, gene, drop = TRUE]
  gene_non_zero_mean <- mean(non_zero_spots, na.rm = TRUE)
  coords <- GetTissueCoordinates(seurat_obj)
  get_avg_expression <- function(metadata_column) {
    positive_spots <- WhichCells(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[metadata_column]] == "positive", ]))
    avg_expressions <- c()
    for (spot in positive_spots) {
      spot_coords <- coords[spot, , drop = FALSE]
      distances <- sqrt(rowSums((t(t(coords) - as.numeric(spot_coords))) ^ 2))
      closest_neighbors <- order(distances)[1:7]
      neighbor_gene_expression <- gene_data[closest_neighbors, gene, drop = TRUE]
      avg_expression <- mean(neighbor_gene_expression, na.rm = TRUE)
      avg_expressions <- c(avg_expressions, avg_expression)
    }
    mean(avg_expressions, na.rm = TRUE)
  }
  avg_expression1 <- get_avg_expression(metadata_column1)
  avg_expression2 <- get_avg_expression(metadata_column2)
  normalized_expression1 <- avg_expression1 / gene_non_zero_mean
  normalized_expression2 <- avg_expression2 / gene_non_zero_mean
  result_df <- data.frame(Object = obj_name, 
                          Gene = gene,
                          Avg_Expression_Column1 = avg_expression1, 
                          Avg_Expression_Column2 = avg_expression2,
                          Normalized_Expression_Column1 = normalized_expression1,
                          Normalized_Expression_Column2 = normalized_expression2)
  return(result_df)
}

# Loop through each Seurat object and each gene to calculate log2 fold change
for (obj_name in seurat_names) {
  cat("Processing", obj_name, "\n")
  seurat_obj <- get(obj_name)
  for (gene in genes_of_interest) {
    tryCatch({
      result <- calculate_fold_change_expression(seurat_obj, obj_name, metadata_column1, metadata_column2, gene)
      all_results[[paste(obj_name, gene, sep = "_")]] <- result
    }, error = function(e) {
      message(paste("Error processing object", obj_name, "and gene", gene, ":", e$message))
    })
  }
}

# Combine all results into a single dataframe
merged_results <- do.call(rbind, all_results)

# Calculate log2 fold change for each gene and perform statistical testing
final_results <- data.frame()
for (gene in genes_of_interest) {
  gene_data <- merged_results[merged_results$Gene == gene, ]
  gene_data$Log2_Fold_Change <- log2(gene_data$Avg_Expression_Column1 / gene_data$Avg_Expression_Column2)
  # Perform paired t-test or Wilcoxon test
  test_result <- tryCatch({
    t_test <- t.test(gene_data$Avg_Expression_Column1, gene_data$Avg_Expression_Column2, paired = TRUE)
    p_value <- t_test$p.value
  }, error = function(e) {
    message(paste("Error in test for gene", gene, ":", e$message))
    NA
  })
  gene_data$P_Value <- p_value
  final_results <- rbind(final_results, gene_data)
}

# Adjust p-values for multiple testing (Benjamini-Hochberg correction)
final_results$Adjusted_P_Value <- p.adjust(final_results$P_Value, method = "BH")

# Print the final results per tissue
print(final_results)

# Set plotting order by average log2 fold change
gene_order <- final_results %>%
  group_by(Gene) %>%
  summarise(Average_Log2_Fold_Change = mean(Log2_Fold_Change, na.rm = TRUE)) %>%
  arrange(desc(Average_Log2_Fold_Change)) %>%
  pull(Gene)

# Plot the log2 fold change for each gene
ggplot(final_results, aes(x = factor(Gene, levels = gene_order), y = Log2_Fold_Change)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without displaying outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Jittered points for each data point
  labs(title = "Log2 Fold Change of Gene Expression (Ordered by Average Fold Change)",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme_minimal()

# Export results
write.csv(final_results, "Path to file /Expression_of_interest_genes_of_interest_local_neighbourhoods_FAP_PCs_vs_FAP_PC_neg_vessel.csv", row.names=TRUE)

#### Ligand expression in FAP_PC enriched and FAP_PC depleted vessel spots (does not include neighbor spots) ####

# Define the two metadata columns and genes of interest
metadata_column1 <- "FAP_PCs_thresholded"             
metadata_column2 <- "FAP_PCs_neg_vessel_thresholded"
genes_of_interest <- c("WNT5A", "THBS2", "SERPINE1", "RARRES2", "CXCL12", "CXCL1", "CCL2", "CSF1", "SEMA3C")

# List of Seurat objects (vessel enriched)
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# List to store results
all_results <- list()

# Function to calculate gene expression only in positive spots for two metadata columns
calculate_fold_change_expression <- function(seurat_obj, obj_name, metadata_column1, metadata_column2, gene) {
  if (!all(c(metadata_column1, metadata_column2) %in% colnames(seurat_obj@meta.data))) {
    stop("One or both metadata annotation columns not found in the object")
  }
  if (!gene %in% rownames(seurat_obj)) {
    stop(paste("Gene", gene, "not found in the object"))
  }
  gene_data <- FetchData(object = seurat_obj, vars = gene, layer = "data")
  non_zero_spots <- gene_data[gene_data[[gene]] > 0, gene, drop = TRUE]
  gene_non_zero_mean <- mean(non_zero_spots, na.rm = TRUE)
  get_avg_expression_positive_spots <- function(metadata_column) {
    positive_spots <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[metadata_column]] == "positive"
    ]
    positive_spots <- intersect(positive_spots, rownames(gene_data))
    if (length(positive_spots) == 0) {
      warning(paste("No positive spots found for", metadata_column, "in", obj_name))
      return(NA)
    }
    positive_gene_expression <- gene_data[positive_spots, gene, drop = TRUE]
    avg_expression <- mean(positive_gene_expression, na.rm = TRUE)
    return(avg_expression)
  }
  avg_expression1 <- get_avg_expression_positive_spots(metadata_column1)
  avg_expression2 <- get_avg_expression_positive_spots(metadata_column2)
  normalized_expression1 <- avg_expression1 / gene_non_zero_mean
  normalized_expression2 <- avg_expression2 / gene_non_zero_mean
  result_df <- data.frame(
    Object = obj_name,
    Gene = gene,
    Avg_Expression_Column1 = avg_expression1,
    Avg_Expression_Column2 = avg_expression2,
    Normalized_Expression_Column1 = normalized_expression1,
    Normalized_Expression_Column2 = normalized_expression2
  )
  return(result_df)
}

# Loop through each Seurat object and each gene
for (obj_name in seurat_names) {
  cat("Processing", obj_name, "\n")
  seurat_obj <- get(obj_name)
  for (gene in genes_of_interest) {
    tryCatch({
      result <- calculate_fold_change_expression(
        seurat_obj,
        obj_name,
        metadata_column1,
        metadata_column2,
        gene
      )
      all_results[[paste(obj_name, gene, sep = "_")]] <- result
    }, error = function(e) {
      message(paste("Error processing object", obj_name, "and gene", gene, ":", e$message))
    })
  }
}

# Combine results
merged_results <- do.call(rbind, all_results)

# Calculate log2 fold change for each gene and perform statistical testing
final_results <- data.frame()
for (gene in genes_of_interest) {
  gene_data <- merged_results[merged_results$Gene == gene, ]
  gene_data$Log2_Fold_Change <- log2(
    gene_data$Avg_Expression_Column1 / gene_data$Avg_Expression_Column2
  )
  test_result <- tryCatch({
    t_test <- t.test(
      gene_data$Avg_Expression_Column1,
      gene_data$Avg_Expression_Column2,
      paired = TRUE
    )
    t_test$p.value
  }, error = function(e) {
    message(paste("Error in test for gene", gene, ":", e$message))
    NA
  })
  gene_data$P_Value <- test_result
  final_results <- rbind(final_results, gene_data)
}

# Adjust p-values for multiple testing
final_results$Adjusted_P_Value <- p.adjust(final_results$P_Value, method = "BH")

# Print final results
print(final_results)

# Set plotting order by average log2 fold change
gene_order <- final_results %>%
  group_by(Gene) %>%
  summarise(
    Average_Log2_Fold_Change = mean(Log2_Fold_Change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Average_Log2_Fold_Change)) %>%
  pull(Gene)

# Plot the log2 fold change for each gene
Ligand_boxplot <- ggplot(final_results, aes(x = factor(Gene, levels = gene_order), y = Log2_Fold_Change)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(
    title = "Log2 Fold Change of Gene Expression in Positive Spots Only",
    x = "Gene",
    y = "Log2 Fold Change"
  ) +
  theme_minimal()

# Export
ggsave(
  filename = "Ligand_boxplot.svg",
  plot = Ligand_boxplot,
  device = "svg",
  width = 12,
  height = 6,
  units = "in"
)

#### Myeloid marker gene expression vs distance to FAP_PC enriched or FAP_PC depleted vessel spots ####

# Define the metadata column and genes of interest
metadata_column <- "FAP_PCs_neg_vessel_thresholded" 
# metadata_column <- "FAP_PCs_thresholded"

Genesofinterest <- c(
  "CD68", "ITGAM", "CD14", "IL1B", "IL6", "IL12A", "IL12B", "IL23A",
  "TNF", "CD80", "CD86", "HLADR", "NOS2", "IRF5", "STAT1", "SOCS3",
  "MRC1", "IRF4", "CD163", "CD200R1", "MARCO", "MSR1", "STAT6",
  "SOCS2", "CCL2", "CCL7", "CSF1", "CSF2", "HGF", "MIF", "CXCL12",
  "CXCL1", "IL4", "IL34", "IFNG", "IL8", "IL10", "IL13", "TGFB1",
  "RARRES2", "SPHK1"
)

# List of Seurat object names (vessel enriched)
seurat_names <- c(
  "T248", "T255", "T259", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1",
  "GBM_ZH1007nec", "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Function to calculate distance-binned gene expression
calculate_distance_vs_gene_expression <- function(seurat_obj, metadata_column, Genesofinterest) {
  if (!metadata_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Metadata annotation column", metadata_column, "not found in the object"))
  }
  gene_data <- tryCatch(
    FetchData(object = seurat_obj, vars = Genesofinterest, layer = "data"),
    error = function(e) {
      stop(paste("Error fetching data for genes of interest:", e$message))
    }
  )
  available_genes <- intersect(Genesofinterest, colnames(gene_data))
  if (length(available_genes) == 0) {
    stop("None of the genes of interest are found")
  }
  if (length(available_genes) < length(Genesofinterest)) {
    warning("Some genes of interest are not found and will be skipped.")
  }
  coords <- GetTissueCoordinates(seurat_obj)
  positive_metadata_spots <- rownames(
    seurat_obj@meta.data[seurat_obj@meta.data[[metadata_column]] == "positive", , drop = FALSE]
  )
  positive_metadata_spots <- intersect(positive_metadata_spots, rownames(coords))
  if (length(positive_metadata_spots) == 0) {
    stop("No positive metadata spots found")
  }
  distance_expression_list <- list()
  for (spot in positive_metadata_spots) {
    spot_coords <- coords[spot, , drop = FALSE]
    distances <- sqrt(rowSums((t(t(coords) - as.numeric(spot_coords)))^2))
    min_distance <- min(distances[distances > 0], na.rm = TRUE)
    if (!is.finite(min_distance) || min_distance == 0) {
      next
    }
    normalized_distances <- distances / min_distance
    spot_data <- data.frame(
      SpotID = rownames(coords),
      Distance = normalized_distances
    )
    spot_data <- spot_data[spot_data$SpotID %in% rownames(gene_data), , drop = FALSE]
    for (gene in available_genes) {
      spot_data[[gene]] <- gene_data[spot_data$SpotID, gene]
    }
    distance_expression_list[[spot]] <- spot_data
  }
  distance_expression_df <- do.call(rbind, distance_expression_list)
  if (is.null(distance_expression_df) || nrow(distance_expression_df) == 0) {
    stop("No distance-expression data could be calculated.")
  }
  max_distance <- max(distance_expression_df$Distance, na.rm = TRUE)
  distance_expression_df$DistanceBin <- cut(
    distance_expression_df$Distance,
    breaks = seq(1, ceiling(max_distance), by = 1),
    include.lowest = TRUE
  )
  avg_expression <- aggregate(
    distance_expression_df[, available_genes, drop = FALSE],
    by = list(DistanceBin = distance_expression_df$DistanceBin),
    FUN = mean,
    na.rm = TRUE
  )
  for (gene in available_genes) {
    min_val <- min(avg_expression[[gene]], na.rm = TRUE)
    max_val <- max(avg_expression[[gene]], na.rm = TRUE)
    if (is.finite(min_val) && is.finite(max_val) && max_val > min_val) {
      avg_expression[[gene]] <- (avg_expression[[gene]] - min_val) / (max_val - min_val)
    } else {
      avg_expression[[gene]] <- NA
    }
  }
  return(avg_expression)
}

# Run analysis across all Seurat objects
all_results <- list()
for (obj_name in seurat_names) {
  cat("Processing", obj_name, "\n")
  seurat_obj <- get(obj_name)
  result <- tryCatch(
    calculate_distance_vs_gene_expression(
      seurat_obj = seurat_obj,
      metadata_column = metadata_column,
      Genesofinterest = Genesofinterest
    ),
    error = function(e) {
      message(paste("Error processing object", obj_name, ":", e$message))
      return(NULL)
    }
  )
  if (!is.null(result)) {
    result$Object <- obj_name
    all_results[[obj_name]] <- result
  }
}

# Merge all results into one dataframe
merged_results <- do.call(rbind, all_results)

# Clean DistanceBin
merged_results$DistanceBin <- droplevels(as.factor(merged_results$DistanceBin))

# Identify genes actually present in the merged output
available_genes <- intersect(Genesofinterest, colnames(merged_results))

# Remove rows where all gene values are NA
filtered_results <- merged_results[
  !apply(is.na(merged_results[, available_genes, drop = FALSE]), 1, all),
]

# Calculate average expression per DistanceBin across all Seurat objects
merged_averages <- aggregate(
  filtered_results[, available_genes, drop = FALSE],
  by = list(DistanceBin = filtered_results$DistanceBin),
  FUN = mean,
  na.rm = TRUE
)

# Rename gene columns
colnames(merged_averages)[-1] <- paste0("Average_", colnames(merged_averages)[-1])

# print result
print(merged_averages)

# Export
write.csv(merged_averages, "Path to file / Gene_of_interest_distance_to_FAP_PCs_neg_vessel.csv", row.names=TRUE)
