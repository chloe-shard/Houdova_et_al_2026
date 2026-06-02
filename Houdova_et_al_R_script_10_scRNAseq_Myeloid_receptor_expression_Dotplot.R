# -----------------------------------------------------------------------------------------------------------------------
#                                            Myeloid receptor expression - Dotplot
# -----------------------------------------------------------------------------------------------------------------------

# Load libraries
library(Seurat)

# Subset Myeloid clusters #
Idents(obj_Ebert) <- "Pooled_Cell_Type"
table(obj_Ebert$Pooled_Cell_Type)
obj_Ebert_Myeloid <- subset(obj_Ebert, idents = c("Myeloid"))

Idents(obj_Abdelfattah) <- "Pooled_Cell_Type"
table(obj_Abdelfattah$Pooled_Cell_Type)
obj_Abdelfattah_Myeloid <- subset(obj_Abdelfattah, idents = c("Myeloid"))

Idents(obj_LeBlanc) <- "Pooled_Cell_Type"
table(obj_LeBlanc$Pooled_Cell_Type)
obj_LeBlanc_Myeloid <- subset(obj_LeBlanc, idents = c("Myeloid"))

Idents(obj_Nomura) <- "Pooled_Cell_Type"
table(obj_Nomura$Pooled_Cell_Type)
obj_Nomura_Myeloid <- subset(obj_Nomura, idents = c("Myeloid"))

# Myeloid markers based on based on Pombo Antunes et al 2021
Monocyte <- c("FCN1", "VCAN", "CFP", "APOBEC3A")
Macrophage <- c("S100A6", "LYZ", "MARCO", "APOE", "CD14", "FCGR1A", "TGFBI", "CLEC12A", "FXYD5", "C1QA")
Microglia <- c("TMEM119", "P2RY12", "SALL1", "CX3CR1", "TMIGD3", "APOC2", "SCIN")

# Combine individual gene sets into a named list
signatures <- list(Monocyte = Monocyte, Macrophage = Macrophage, Microglia = Microglia) # Add more signatures as needed

seurat_obj <- obj_Abdelfattah_Myeloid # swap for the other datasets

# Normalize gene expression 
expr_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
gene_means <- rowMeans(expr_data)
normalized_expr_data <- sweep(expr_data, 1, gene_means, "/")

# Initialize a matrix to store the scores
num_cells <- ncol(seurat_obj)
num_signatures <- length(signatures)
scores_matrix <- matrix(NA, nrow = num_signatures, ncol = num_cells)

# Function to calculate state and non-state scores
calculate_scores <- function(normalized_expr_data, signature_genes) {
  valid_genes <- signature_genes[signature_genes %in% rownames(normalized_expr_data)]
  if (length(valid_genes) == 0) {
    return(rep(NA, ncol(normalized_expr_data)))
  }
  state_score <- colMeans(normalized_expr_data[valid_genes, , drop = FALSE], na.rm = TRUE)
  non_signature_genes <- setdiff(rownames(normalized_expr_data), valid_genes)
  non_state_score <- colMeans(normalized_expr_data[non_signature_genes, , drop = FALSE], na.rm = TRUE)
  ratio <- state_score / non_state_score
  return(ratio)
}

# Repeat 1-4 for each state signature
for (i in seq_along(signatures)) {
  signature_genes <- signatures[[i]]
  cat("Processing signature:", names(signatures)[i], "\n")
  cat("Valid genes:", signature_genes[signature_genes %in% rownames(normalized_expr_data)], "\n")
  scores_matrix[i, ] <- calculate_scores(normalized_expr_data, signature_genes)
}

# Assign state with maximum score to cell if state score/non-state score > 1
max_scores <- apply(scores_matrix, 2, max, na.rm = TRUE)
max_score_indices <- apply(scores_matrix, 2, function(x) ifelse(any(x > 1, na.rm = TRUE), which.max(x), NA))

# Create a vector for state assignment with signature names
signature_names <- names(signatures)
state_assignment <- sapply(seq_along(max_score_indices), function(idx) {
  if (!is.na(max_score_indices[idx]) && max_scores[idx] > 1) {
    return(signature_names[max_score_indices[idx]])
  } else {
    return("None")
  }
})


# Add state assignment to metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = state_assignment, col.name = "Myeloid_subpop")

# Check the table of state assignments
table(seurat_obj$Myeloid_subpop)

Idents(seurat_obj) <- "Myeloid_subpop"
DimPlot(seurat_obj)

obj_Abdelfattah_Myeloid  <- seurat_obj

# subset Macrophage cluster #
Idents(obj_Ebert_Myeloid) <- "Myeloid_subpop"
obj_Ebert_Mac <- subset(obj_Ebert_Myeloid, idents = c("Macrophage"))

Idents(obj_LeBlanc_Myeloid) <- "Myeloid_subpop"
obj_LeBlanc_Mac <- subset(obj_LeBlanc_Myeloid, idents = c("Macrophage"))

Idents(obj_Abdelfattah_Myeloid) <- "Myeloid_subpop"
obj_Abdelfattah_Mac <- subset(obj_Abdelfattah_Myeloid, idents = c("Macrophage"))

Idents(obj_Nomura_Myeloid) <- "Myeloid_subpop"
obj_Nomura_Mac <- subset(obj_Nomura_Myeloid, idents = c("Macrophage"))

# Score M1 and M2 
M1 <- c("IL1B", "IL6", "IL12A", "IL12B", "IL23A", "TNF", "CD80", "CD86", "HLA-DRA", "NOS2", "IRF5", "STAT1", "SOCS3")
M2 <- c("CD14", "MRC1", "IRF4", "CD163", "CD200R1", "MARCO", "MSR1", "STAT6", "SOCS2")

signatures <- list(M1 = M1, M2 = M2)

seurat_obj <- obj_Abdelfattah_Mac

# Normalize gene expression by dividing each value with the average expression of the gene across all cells
expr_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
gene_means <- rowMeans(expr_data)
normalized_expr_data <- sweep(expr_data, 1, gene_means, "/")

# Initialize a matrix to store the scores
num_cells <- ncol(seurat_obj)
num_signatures <- length(signatures)
scores_matrix <- matrix(NA, nrow = num_signatures, ncol = num_cells)

# Function to calculate state and non-state scores
calculate_scores <- function(normalized_expr_data, signature_genes) {
  valid_genes <- signature_genes[signature_genes %in% rownames(normalized_expr_data)]
  if (length(valid_genes) == 0) {
    return(rep(NA, ncol(normalized_expr_data)))
  }
  
  state_score <- colMeans(normalized_expr_data[valid_genes, , drop = FALSE], na.rm = TRUE)
  non_signature_genes <- setdiff(rownames(normalized_expr_data), valid_genes)
  non_state_score <- colMeans(normalized_expr_data[non_signature_genes, , drop = FALSE], na.rm = TRUE)
  
  ratio <- state_score / non_state_score
  return(ratio)
}

# Repeat 1-4 for each state signature
for (i in seq_along(signatures)) {
  signature_genes <- signatures[[i]]
  cat("Processing signature:", names(signatures)[i], "\n")
  cat("Valid genes:", signature_genes[signature_genes %in% rownames(normalized_expr_data)], "\n")
  scores_matrix[i, ] <- calculate_scores(normalized_expr_data, signature_genes)
}

# Assign state with maximum score to cell if state score/non-state score > 0.5
max_scores <- apply(scores_matrix, 2, max, na.rm = TRUE)

max_score_indices <- apply(scores_matrix, 2, function(x) ifelse(any(x > 0.5, na.rm = TRUE), which.max(x), NA))

# Create a vector for state assignment with signature names
signature_names <- names(signatures)
state_assignment <- sapply(seq_along(max_score_indices), function(idx) {
  if (!is.na(max_score_indices[idx]) && max_scores[idx] > 0.5) {
    return(signature_names[max_score_indices[idx]])
  } else {
    return("None")
  }
})

# Add state assignment to metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = state_assignment, col.name = "Macrophage_subpop")

# Check the table of state assignments
table(seurat_obj$Macrophage_subpop)

Idents(seurat_obj) <- "Macrophage_subpop"
DimPlot(seurat_obj)

obj_Abdelfattah_Mac <- seurat_obj

#### Add Macrophage annotations back to Myeloid metadata ####

# Define dataset prefixes
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

# Loop through each dataset
for (prefix in datasets) {
  myeloid_obj <- get(paste0("obj_", prefix, "_Myeloid"))
  mac_obj <- get(paste0("obj_", prefix, "_Mac"))
  meta_myeloid <- myeloid_obj@meta.data
  meta_mac <- mac_obj@meta.data
  macrophage_cells <- rownames(meta_myeloid)[meta_myeloid$Myeloid_subpop == "Macrophage"]
  common_cells <- intersect(macrophage_cells, rownames(meta_mac))
  meta_myeloid[common_cells, "Myeloid_subpop"] <- meta_mac[common_cells, "Macrophage_subpop"]
  myeloid_obj@meta.data <- meta_myeloid
  assign(paste0("obj_", prefix, "_Myeloid"), myeloid_obj)
}

#### Add Myeloid metadata to original object ####

# Dataset prefixes
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

# Loop through datasets
for (prefix in datasets) {
  full_obj <- get(paste0("obj_", prefix))
  myeloid_obj <- get(paste0("obj_", prefix, "_Myeloid"))
  new_annotation <- rep("None", ncol(full_obj))
  names(new_annotation) <- colnames(full_obj)
  myeloid_meta <- myeloid_obj@meta.data
  myeloid_annotations <- myeloid_meta$Myeloid_subpop
  names(myeloid_annotations) <- rownames(myeloid_meta)
  overlap_cells <- intersect(names(new_annotation), names(myeloid_annotations))
  new_annotation[overlap_cells] <- myeloid_annotations[overlap_cells]
  full_obj$Myeloid_subpop <- new_annotation
  assign(paste0("obj_", prefix), full_obj)
}

#### Receptor expression in myeloid subpopulations ####

genes_of_interest <- c("MCAM", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "RYK", "ROR2", "LRP5", "LRP6", "TLR4", "CD36", "CD47", "LRP1", "PLAU", "CMKLR1", "GPR1", "CCRL2", "CXCR4", "ACKR3", "CXCR1", "CXCR2", "CCR2", "CCR4", "CSF1R", "NRP1", "NRP2", "PLXNB1", "PLXND1")

# Function to calculate average expression (mean), z-score, and percentage expression per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    stop(paste("No genes of interest found in the dataset:", dataset_name))
  }
  assay_data <- GetAssayData(seurat_obj, layer = "data")  # Use 'data' for normalized counts
  if (!is.matrix(assay_data)) {
    assay_data <- as.matrix(assay_data)
  }
  expression_df <- as.data.frame(t(assay_data[genes_present, ]))  # Transpose for easier manipulation
  expression_df$cell_type <- seurat_obj@meta.data$Myeloid_subpop 
  avg_expression <- expression_df %>%
    group_by(cell_type) %>%
    summarise_at(vars(all_of(genes_present)), mean, na.rm = TRUE) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression")
  avg_expression_long <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = scale(avg_expression)) %>%
    ungroup()
  binary_expr <- assay_data[genes_present, ] > 0  # Convert to binary expression (expressed or not)
  binary_expr <- as.data.frame(t(binary_expr))  # Transpose for cell-wise rows
  binary_expr$cell_type <- seurat_obj@meta.data$Myeloid_subpop
  percentage_expr <- binary_expr %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))  # Clean up gene names
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "cell_type")) %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_"))  # Combine cell type and dataset
  return(plot_data)
}

# Process each dataset manually
plot_data_Nomura <- process_dataset(obj_Nomura_Myeloid, "Nomura")
plot_data_Ebert <- process_dataset(obj_Ebert_Myeloid, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc_Myeloid, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah_Myeloid, "Adelfattah")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
  plot_data_Nomura
)

# Specify the order of genes for plotting
custom_gene_order <- c("MCAM", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "RYK", "ROR2", "LRP5", "LRP6", "TLR4", "CD36", "CD47", "LRP1", "PLAU", "CMKLR1", "GPR1", "CCRL2", "CXCR4", "ACKR3", "CXCR1", "CXCR2", "CCR2", "CCR4", "CSF1R", "NRP1", "NRP2", "PLXNB1", "PLXND1")

# Specify the order of cell annotations for plotting
custom_cellstate_order <- c("Microglia", "M1", "M2", "Monocyte", "None")

# Ensure the 'gene' column respects the specified order
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)

# Ensure the 'Myeloid_subpop' in cell_type_dataset respects the specified order
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(Myeloid_subpop) {
      grep(paste0("^", Myeloid_subpop, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
    })
  )
)

custom_cellstate_order <- c("Microglia","M1","M2","Monocyte","None")
dataset_order <- c("Adelfattah","Ebert","LeBlanc","Nomura")

mat <- outer(custom_cellstate_order, dataset_order, paste, sep = "_")
wanted_levels <- as.vector(t(mat))  # <-- row-wise flatten

combined_plot_data$cell_type_dataset <- factor(
  as.character(combined_plot_data$cell_type_dataset),
  levels = wanted_levels
)

levels(combined_plot_data$cell_type_dataset)

combined_plot <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = percentage_expr, color = z_score)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
ggsave(
  filename = "Dotplot_Myeloid_cells_Receptors.svg",
  plot     = combined_plot,
  device   = "svg",
  width    = 8,
  height   = 5,
  units    = "in"
)
