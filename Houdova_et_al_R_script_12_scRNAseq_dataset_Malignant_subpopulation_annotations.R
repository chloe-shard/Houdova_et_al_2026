# -----------------------------------------------------------------------------------------
#                            Malignant sub-population annotations
# -----------------------------------------------------------------------------------------

# Load libraries
library(Seurat)
library(readr)

#### Subset Malignant cluster ####

Idents(obj_Ebert) <- "Pooled_Cell_Type"
table(obj_Ebert$Pooled_Cell_Type)
obj_Ebert_Malignant <- subset(obj_Ebert, idents = c("Malignant"))

Idents(obj_Abdelfattah) <- "Pooled_Cell_Type"
table(obj_Abdelfattah$Pooled_Cell_Type)
obj_Abdelfattah_Malignant <- subset(obj_Abdelfattah, idents = c("Malignant"))

Idents(obj_LeBlanc) <- "Pooled_Cell_Type"
table(obj_LeBlanc$Pooled_Cell_Type)
obj_LeBlanc_Malignant <- subset(obj_LeBlanc, idents = c("Malignant"))

Idents(obj_Nomura) <- "Pooled_Cell_Type"
table(obj_Nomura$Pooled_Cell_Type)
obj_Nomura_Malignant <- subset(obj_Nomura, idents = c("Malignant"))

#### Add Malignant annotations to dataset - cell scoring  ####

# Cell state markers obtained from Neftel et al., 2019, Cell
Neftel_cell_state_markers <- read_csv("Path to file/Neftel_cell_state_markers.csv")

# create gene marker lists
NPC1 <- na.omit(Neftel_cell_state_markers[["NPC1"]])
NPC2 <- na.omit(Neftel_cell_state_markers[["NPC2"]])
OPC  <- na.omit(Neftel_cell_state_markers[["OPC"]])
AC   <- na.omit(Neftel_cell_state_markers[["AC"]])
MES1 <- na.omit(Neftel_cell_state_markers[["MES1"]])
MES2 <- na.omit(Neftel_cell_state_markers[["MES2"]])

# Combine individual gene sets into a named list
signatures <- list(NPC1 = NPC1, NPC2 = NPC2, OPC = OPC, AC = AC, MES1 = MES1, MES2 = MES2) 

# Insert each dataset in turn
seurat_obj <- obj_LeBlanc_Malignant

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
seurat_obj <- AddMetaData(seurat_obj, metadata = state_assignment, col.name = "Malignant_subpop")

# Check the table of state assignments
table(seurat_obj$Malignant_subpop)

# Write back
obj_LeBlanc_Malignant  <- seurat_obj
