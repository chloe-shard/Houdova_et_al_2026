# -----------------------------------------------------------------------------------------
#                               Pericyte subpopulation annotations
# -----------------------------------------------------------------------------------------

# Load libraries
library(Seurat)
library(UCell)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)

########### Subset Pericyte cluster ####

Idents(obj_Ebert) <- "Pooled_Cell_Type"
table(obj_Ebert$Pooled_Cell_Type)
obj_Ebert_Pericyte <- subset(obj_Ebert, idents = c("Pericyte"))

Idents(obj_Abdelfattah) <- "Pooled_Cell_Type"
table(obj_Abdelfattah$Pooled_Cell_Type)
obj_Abdelfattah_Pericyte <- subset(obj_Abdelfattah, idents = c("Pericyte"))

Idents(obj_LeBlanc) <- "Pooled_Cell_Type"
table(obj_LeBlanc$Pooled_Cell_Type)
obj_LeBlanc_Pericyte <- subset(obj_LeBlanc, idents = c("Pericyte"))

Idents(obj_Nomura) <- "Pooled_Cell_Type"
table(obj_Nomura$Pooled_Cell_Type)
obj_Nomura_Pericyte <- subset(obj_Nomura, idents = c("Pericyte"))

########### reclustering of Pericyte cluster ####

#Ebert
DefaultAssay(obj_Ebert_Pericyte) <- 'RNA'
obj_Ebert_Pericyte  <- FindVariableFeatures(obj_Ebert_Pericyte, selection.method = "vst", nfeatures = 1000)

# Remove features with zero or negative values
variable_features <- VariableFeatures(obj_Ebert_Pericyte)
valid_features <- variable_features[variable_features > 0]

top20 <- head(VariableFeatures(obj_Ebert_Pericyte), 20)
plot1 <- VariableFeaturePlot(obj_Ebert_Pericyte)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2
obj_Ebert_Pericyte  <- ScaleData(obj_Ebert_Pericyte)
obj_Ebert_Pericyte  <- RunPCA(obj_Ebert_Pericyte , verbose = FALSE)

ElbowPlot(obj_Ebert_Pericyte)

obj_Ebert_Pericyte <- FindNeighbors(obj_Ebert_Pericyte, reduction = "pca", dims = 1:20)
obj_Ebert_Pericyte <- FindClusters(obj_Ebert_Pericyte, resolution = 0.2)
obj_Ebert_Pericyte <- RunUMAP(obj_Ebert_Pericyte, reduction = "pca", dims = 1:15)
UMAPPlot(obj_Ebert_Pericyte)

#Abdelfattah

DefaultAssay(obj_Abdelfattah_Pericyte) <- 'RNA'
obj_Abdelfattah_Pericyte  <- FindVariableFeatures(obj_Abdelfattah_Pericyte, selection.method = "vst", nfeatures = 1000)

# Remove features with zero or negative values
variable_features <- VariableFeatures(obj_Abdelfattah_Pericyte)
valid_features <- variable_features[variable_features > 0]

top20 <- head(VariableFeatures(obj_Abdelfattah_Pericyte), 20)
plot1 <- VariableFeaturePlot(obj_Abdelfattah_Pericyte)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2
obj_Abdelfattah_Pericyte  <- ScaleData(obj_Abdelfattah_Pericyte)
obj_Abdelfattah_Pericyte  <- RunPCA(obj_Abdelfattah_Pericyte , verbose = FALSE)

ElbowPlot(obj_Abdelfattah_Pericyte)

obj_Abdelfattah_Pericyte <- FindNeighbors(obj_Abdelfattah_Pericyte, reduction = "pca", dims = 1:20)
obj_Abdelfattah_Pericyte <- FindClusters(obj_Abdelfattah_Pericyte, resolution = 0.3)
obj_Abdelfattah_Pericyte <- RunUMAP(obj_Abdelfattah_Pericyte, reduction = "pca", dims = 1:35)
UMAPPlot(obj_Abdelfattah_Pericyte)

#LeBlanc

DefaultAssay(obj_LeBlanc_Pericyte) <- 'RNA'
obj_LeBlanc_Pericyte  <- FindVariableFeatures(obj_LeBlanc_Pericyte, selection.method = "vst", nfeatures = 1000)

# Remove features with zero or negative values
variable_features <- VariableFeatures(obj_LeBlanc_Pericyte)
valid_features <- variable_features[variable_features > 0]

top20 <- head(VariableFeatures(obj_LeBlanc_Pericyte), 20)
plot1 <- VariableFeaturePlot(obj_LeBlanc_Pericyte)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2
obj_LeBlanc_Pericyte  <- ScaleData(obj_LeBlanc_Pericyte)
obj_LeBlanc_Pericyte  <- RunPCA(obj_LeBlanc_Pericyte , verbose = FALSE)

ElbowPlot(obj_LeBlanc_Pericyte)

obj_LeBlanc_Pericyte <- FindNeighbors(obj_LeBlanc_Pericyte, reduction = "pca", dims = 1:20)
obj_LeBlanc_Pericyte <- FindClusters(obj_LeBlanc_Pericyte, resolution = 0.1)
obj_LeBlanc_Pericyte <- RunUMAP(obj_LeBlanc_Pericyte, reduction = "pca", dims = 1:10)
UMAPPlot(obj_LeBlanc_Pericyte)

#Nomura
DefaultAssay(obj_Nomura_Pericyte) <- 'RNA'
obj_Nomura_Pericyte  <- FindVariableFeatures(obj_Nomura_Pericyte, selection.method = "vst", nfeatures = 1000)

# Remove features with zero or negative values
variable_features <- VariableFeatures(obj_Nomura_Pericyte)
valid_features <- variable_features[variable_features > 0]

top20 <- head(VariableFeatures(obj_Nomura_Pericyte), 20)
plot1 <- VariableFeaturePlot(obj_Nomura_Pericyte)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2
obj_Nomura_Pericyte  <- ScaleData(obj_Nomura_Pericyte)
obj_Nomura_Pericyte  <- RunPCA(obj_Nomura_Pericyte , verbose = FALSE)

ElbowPlot(obj_Nomura_Pericyte)

obj_Nomura_Pericyte <- FindNeighbors(obj_Nomura_Pericyte, reduction = "pca", dims = 1:20)
obj_Nomura_Pericyte <- FindClusters(obj_Nomura_Pericyte, resolution = 0.2)
obj_Nomura_Pericyte <- RunUMAP(obj_Nomura_Pericyte, reduction = "pca", dims = 1:14)
UMAPPlot(obj_Nomura_Pericyte)

######### Add Pericyte annotations to dataset - cell scoring

Transport_PCs <- c("ITIH5", "P2RY14", "SLC6A13", "SLC6A12", "COL4A3", "ADARB2", "SLC30A10", "SLC19A1", "DIO3OS", "PTGDR2", "SHISA6")
ECM_PCs <- c("PRR16", "PIEZO2", "PRKG1", "THY1", "FAT1", "CD36", "IL1RAPL1", "NREP", "C2orf27A", "BMP1", "TGFBI", "TUSC3", "DOK6", "EBF2")
Interferon_PCs <- c("CSF2RA", "SAMSN1", "CD53", "DOCK2", "ADAM28", "PTPRC", "SLA", "RBM47", "HLA-DRB1", "VSIG4", "RNF144B", "ITGB2", "PLEK", "SLC11A1", "FYB", "HLA-DQB1", "PIKAP1", "CD86", "MSR1", "ITGAX", "SLC16A10", "MS4A4A", "SYK", "NLRP3", "FCGR3A", "AOAH", "FCER1G", "C1QC", "LCP1")
SMCs <- c("MYH11", "CSDC2", "DES", "AC097724.3", "PLN", "LMOD1", "CNN1", "NET1", "MYOCD")
FAP_PCs <- c("POSTN", "CA12", "CXCL6", "FBLN2", "CYP1B1", "FAP", "COL11A1", "PDGFRA", "SFRP2", "FLRT2", "THBS2", "SGCD", "FMOD", "LTBP2", "CLDN11", "GAS7", "PDPN", "ELN", "OGN", "DPSL3", "CLMP", "BICC1", "ROR2", "TMTM132C", "SFRP4", "GJA1", "OMD")
Proliferative_PCs <- c("ASPM", "PBK", "AURKB", "UBE2C", "FAM64A", "KIF2C", "SKA1", "SHCBP1", "CENPF", "RRM2", "KIF23", "CDC20", "KIF11", "TOP2A", "KIF15", "MAD2L1", "DLGAP5", "CENPM", "GTSE1", "CDCA8", "ANLN", "TPX2", "KIF14", "CCNA2", "NUF2", "KIF18B", "CEP55", "CCNB1", "CKAP2L", "DIAPH3", "TK1", "SGOL1", "NUSAP1")

# Combine individual gene sets into a named list
signatures <- list(Transport_PCs = Transport_PCs, ECM_PCs = ECM_PCs, Interferon_PCs = Interferon_PCs, SMCs = SMCs, FAP_PCs = FAP_PCs, Proliferative_PCs = Proliferative_PCs) # Add more signatures as needed

# Insert each dataset in turn
seurat_obj <- obj_LeBlanc_Pericyte

# 1. Normalize gene expression by dividing each value with the average expression of the gene across all cells
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

# 5. Repeat 1-4 for each state signature
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
seurat_obj <- AddMetaData(seurat_obj, metadata = state_assignment, col.name = "Pericyte_subpop")

# Check the table of state assignments
table(seurat_obj$Pericyte_subpop)

# Write back
obj_LeBlanc_Pericyte  <- seurat_obj

### UCell smoothing

# Assign proper row and column names to the UCell score matrix 
rownames(scores_matrix) <- names(signatures)
colnames(scores_matrix) <- colnames(seurat_obj)

# Transpose matrix so rows = cells and columns = signatures
score_df <- as.data.frame(t(scores_matrix))
rownames(score_df) <- colnames(seurat_obj)

# Add raw UCell signature scores to Seurat metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = score_df)

# Perform kNN smoothing of signature scores using PCA reduction
seurat_obj <- SmoothKNN(
  seurat_obj,
  signature.names = names(signatures),
  reduction = "pca"
)

# Extract smoothed signature score columns from metadata
knn_cols <- paste0(names(signatures), "_kNN")
knn_mat <- seurat_obj[[]] %>%
  dplyr::select(dplyr::all_of(knn_cols)) %>%
  as.matrix()

# Compute maximum smoothed signature score per cell
max_scores_knn <- apply(knn_mat, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

# Identify index of dominant signature per cell
max_idx_knn <- apply(knn_mat, 1, function(x) {
  if (all(is.na(x))) return(NA_integer_)
  which.max(replace(x, is.na(x), -Inf))
})
cutoff_knn <- 1  #(score ratio > 1)

# Assign each cell to dominant signature if above cutoff,
signature_names <- names(signatures)
state_assignment_knn <- ifelse(
  is.na(max_scores_knn) | max_scores_knn <= cutoff_knn,
  "None",
  signature_names[max_idx_knn]
)

# Save smoothed state assignment into Seurat metadata
seurat_obj$Pericyte_subpop_smooth <- state_assignment_knn

# Write back
obj_LeBlanc_Pericyte <- seurat_obj

######### Change the Pericyte subpop colour in UMAPs to be the same #####

# Define the desired order
state_order <- c(
  "Transport_PCs",
  "Proliferative_PCs",
  "ECM_PCs",
  "None",
  "SMCs",
  "Interferon_PCs",
  "FAP_PCs"
)

# Put all objects into list
seurat_list <- list(
  obj_Ebert_Pericyte,
  obj_Abdelfattah_Pericyte,
  obj_LeBlanc_Pericyte,
  obj_Nomura_Pericyte
)

# Enforce factor order for Pericyte_subpop_smooth
seurat_list <- lapply(seurat_list, function(obj) {
  obj$Pericyte_subpop_smooth <- factor(
    obj$Pericyte_subpop_smooth,
    levels = state_order
  )
  obj
})

# Write back to the original objects
obj_Ebert_Pericyte        <- seurat_list[[1]]
obj_Abdelfattah_Pericyte <- seurat_list[[2]]
obj_LeBlanc_Pericyte     <- seurat_list[[3]]
obj_Nomura_Pericyte     <- seurat_list[[4]]

# Build a color vector
present_states <- state_order[state_order %in% unique(unlist(
  lapply(seurat_list, function(obj) as.character(obj$Pericyte_subpop_smooth))
))]

default_colors <- hue_pal()(length(present_states))
names(default_colors) <- present_states

# Plotting function
create_umap_plot <- function(seurat_obj, colors, title) {
  DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "Pericyte_subpop_smooth",
    cols = colors
  ) +
    ggtitle(title) +
    theme_minimal()
}

# Create UMAPs for subpopulations
umap_plot1 <- create_umap_plot(obj_Ebert_Pericyte,        default_colors, "UMAP - Ebert (Pericyte)")
umap_plot2 <- create_umap_plot(obj_Abdelfattah_Pericyte, default_colors, "UMAP - Abdelfattah (Pericyte)")
umap_plot3 <- create_umap_plot(obj_LeBlanc_Pericyte,     default_colors, "UMAP - LeBlanc (Pericyte)")
umap_plot4 <- create_umap_plot(obj_Nomura_Pericyte,      default_colors, "UMAP - Nomura (Pericyte)")

# Display
umap_plot1 + umap_plot2 + umap_plot3 + umap_plot4

# Arrange plots in a 2x2 grid
combined_plot1 <- (umap_plot1 | umap_plot2) /
  (umap_plot3 | umap_plot4)

ggsave(
  filename = "UMAP_pericyte_subpop_all_datasets.svg",
  plot     = combined_plot1,
  device   = "svg",
  width    = 10,
  height   = 10,
  units    = "in"
)

# Create UMAPs for FAP in subpopulations
FAP_plot1 <- FeaturePlot(obj_Ebert_Pericyte, features = "FAP")
FAP_plot2 <- FeaturePlot(obj_Abdelfattah_Pericyte, features = "FAP")
FAP_plot3 <- FeaturePlot(obj_LeBlanc_Pericyte, features = "FAP")
FAP_plot4 <- FeaturePlot(obj_Nomura_Pericyte, features = "FAP")

# Display
FAP_plot1 + FAP_plot2 + FAP_plot3 + FAP_plot4

# Arrange plots in a 2x2 grid
combined_plot2 <- (FAP_plot1 | FAP_plot2) /
  (FAP_plot3 | FAP_plot4)

ggsave(
  filename = "UMAP_FAP_pericytes_all_datasets.svg",
  plot     = combined_plot2,
  device   = "svg",
  width    = 10,
  height   = 10,
  units    = "in"
)