# ------------------------------------------------------------------------------------------------------------------
#                                                   Ligand expression dotplots
# ------------------------------------------------------------------------------------------------------------------

#load libraries
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)


######### Ligand expression in FAP cells #############

# subset FAP positive cells
obj_Ebert_FAP_pos <- subset(obj_Ebert, subset = FAP > 0)
obj_LeBlanc_FAP_pos <- subset(obj_LeBlanc, subset = FAP > 0)
obj_Abdelfattah_FAP_pos <- subset(obj_Abdelfattah, subset = FAP > 0)
obj_Nomura_FAP_pos <- subset(obj_Nomura, subset = FAP > 0)

# specify ligands to plot
genes_of_interest <- c("WNT5A", "THBS2", "SERPINE1", "RARRES2", "CXCL12", "CXCL1", "CCL2", "CSF1", "SEMA3C")

# Function to calculate average expression (mean), z-score, and percentage expression per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    stop(paste("No genes of interest found in the dataset:", dataset_name))
  }
  assay_data <- GetAssayData(seurat_obj, layer = "data")
  if (!is.matrix(assay_data)) {
    assay_data <- as.matrix(assay_data)
  }
  expression_df <- as.data.frame(t(assay_data[genes_present, ]))
  expression_df$cell_type <- seurat_obj@meta.data$Pooled_Cell_Type 
  avg_expression <- expression_df %>%
    group_by(cell_type) %>%
    summarise_at(vars(all_of(genes_present)), mean, na.rm = TRUE) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression")
  avg_expression_long <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = scale(avg_expression)) %>%
    ungroup()
  binary_expr <- assay_data[genes_present, ] > 0
  binary_expr <- as.data.frame(t(binary_expr))  
  binary_expr$cell_type <- seurat_obj@meta.data$Pooled_Cell_Type
  percentage_expr <- binary_expr %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))  # Clean up gene names
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "cell_type")) %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_"))
  return(plot_data)
}

# Process each dataset
plot_data_Nomura <- process_dataset(obj_Nomura_FAP_pos, "Nomura")
plot_data_Ebert <- process_dataset(obj_Ebert_FAP_pos, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc_FAP_pos, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah_FAP_pos, "Adelfattah")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
  plot_data_Nomura
)

# Specify the order of cell types and genes for plotting
custom_gene_order <- c("WNT5A", "THBS2", "SERPINE1", "RARRES2", "CXCL12", "CXCL1", "CCL2", "CSF1", "SEMA3C")
custom_cellstate_order <- c("Pericyte", "Endothelial", "Malignant", "Myeloid", "Lymphocyte", "Oligodendrocyte", "OPC_non_malignant", "Astrocyte", "Neuron", "Other")
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(Pooled_Cell_Type) {
      grep(paste0("^", Pooled_Cell_Type, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
    })
  )
)

# Plot the combined dot plot with specified orders
combined_plot <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = percentage_expr, color = z_score)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Gene", y = "Cell Type (Dataset)", size = "Percentage Expressing", color = "Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(scales)

p <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = pmin(percentage_expr, 75),
                 color = z_score)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size_continuous(
    limits = c(0, 75),
    range  = c(0.5, 6),
    breaks = c(0, 25, 50, 75)
  ) +
  labs(
    x = "Gene", y = "Cell Type (Dataset)",
    size = "Percentage Expressing", color = "Z-score"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(
  filename = "Dotplot_FAP_positive_cells_ligands.svg",
  plot     = p,
  device   = "svg",
  width    = 8,
  height   = 6,
  units    = "in"
)

######### Ligand expression in FAP+ pericyte cells vs all cells #############

# Annotate FAP pos vs FAP neg pericytes
annotate_FAP_pericytes <- function(obj,
                                   pericyte_label = "Pericyte",
                                   celltype_col = "Pooled_Cell_Type",
                                   fap_var = "FAP",
                                   threshold = 0) {
  md <- obj@meta.data
  md$FAP_pericyte_status <- NA_character_
  is_pericyte <- as.character(md[[celltype_col]]) == pericyte_label
  fap_vec <- FetchData(obj, vars = fap_var)[, 1]
  md$FAP_pericyte_status[is_pericyte & !is.na(fap_vec) & fap_vec > threshold]  <- "FAP_pos_pericyte"
  md$FAP_pericyte_status[is_pericyte & !is.na(fap_vec) & fap_vec <= threshold] <- "FAP_neg_pericyte"
  obj@meta.data <- md
  obj
}

obj_Ebert        <- annotate_FAP_pericytes(obj_Ebert)
obj_LeBlanc      <- annotate_FAP_pericytes(obj_LeBlanc)
obj_Abdelfattah  <- annotate_FAP_pericytes(obj_Abdelfattah)
obj_Nomura       <- annotate_FAP_pericytes(obj_Nomura)

table(obj_Ebert$Pooled_Cell_Type, obj_Ebert$FAP_pericyte_status, useNA = "ifany")
table(obj_Ebert$FAP_pericyte_status, useNA = "ifany")

add_FAP_pericyte_celltype <- function(obj,
                                      celltype_col = "Pooled_Cell_Type",
                                      fap_status_col = "FAP_pericyte_status",
                                      new_col = "Pooled_Cell_Type_FAP") {
  md <- obj@meta.data
  md[[new_col]] <- as.character(md[[celltype_col]])
  is_pericyte <- md[[celltype_col]] == "Pericyte"
  md[[new_col]][is_pericyte] <- md[[fap_status_col]][is_pericyte]
  obj@meta.data <- md
  obj
}

obj_Ebert       <- add_FAP_pericyte_celltype(obj_Ebert)
obj_LeBlanc     <- add_FAP_pericyte_celltype(obj_LeBlanc)
obj_Abdelfattah <- add_FAP_pericyte_celltype(obj_Abdelfattah)
obj_Nomura      <- add_FAP_pericyte_celltype(obj_Nomura)

table(obj_Ebert$Pooled_Cell_Type_FAP, useNA = "ifany")

# specify ligands to plot
genes_of_interest <- c("WNT5A", "THBS2", "SERPINE1", "RARRES2", "CXCL12", "CXCL1", "CCL2", "CSF1", "SEMA3C")

# Function to calculate average expression (mean), z-score, and percentage expression per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  # Handle cases where no genes are present
  if (length(genes_present) == 0) {
    stop(paste("No genes of interest found in the dataset:", dataset_name))
  }
  assay_data <- GetAssayData(seurat_obj, layer = "data")
  if (!is.matrix(assay_data)) {
    assay_data <- as.matrix(assay_data)
  }
  expression_df <- as.data.frame(t(assay_data[genes_present, ]))
  expression_df$cell_type <- seurat_obj@meta.data$Pooled_Cell_Type_FAP 
  avg_expression <- expression_df %>%
    group_by(cell_type) %>%
    summarise_at(vars(all_of(genes_present)), mean, na.rm = TRUE) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression")
  avg_expression_long <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = scale(avg_expression)) %>%
    ungroup()
  binary_expr <- assay_data[genes_present, ] > 0
  binary_expr <- as.data.frame(t(binary_expr))  
  binary_expr$cell_type <- seurat_obj@meta.data$Pooled_Cell_Type_FAP
  percentage_expr <- binary_expr %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))  # Clean up gene names
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "cell_type")) %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_"))
  return(plot_data)
}

# Process each dataset
plot_data_Nomura <- process_dataset(obj_Nomura, "Nomura")
plot_data_Ebert <- process_dataset(obj_Ebert, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah, "Adelfattah")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
  plot_data_Nomura
)

# Specify the order of cell types and genes for plotting
custom_gene_order <- c("WNT5A", "THBS2", "SERPINE1", "RARRES2", "CXCL12", "CXCL1", "CCL2", "CSF1", "SEMA3C")
custom_cellstate_order <- c("FAP_pos_pericyte", "FAP_neg_pericyte", "Endothelial", "Malignant", "Myeloid", "Lymphocyte", "Oligodendrocyte", "OPC_non_malignant", "Astrocyte", "Neuron", "Other")
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(Pooled_Cell_Type) {
      grep(paste0("^", Pooled_Cell_Type, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
    })
  )
)

# Plot the combined dot plot with specified orders
combined_plot <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = percentage_expr, color = z_score)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Gene", y = "Cell Type (Dataset)", size = "Percentage Expressing", color = "Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(scales)

p <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = pmin(percentage_expr, 75),
                 color = z_score)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size_continuous(
    limits = c(0, 75),
    range  = c(0.5, 6),
    breaks = c(0, 25, 50, 75)
  ) +
  labs(
    x = "Gene", y = "Cell Type (Dataset)",
    size = "Percentage Expressing", color = "Z-score"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(
  filename = "Dotplot_FAP_positive_pericyte_vs_all_cells_ligands.svg",
  plot     = p,
  device   = "svg",
  width    = 8,
  height   = 6,
  units    = "in"
)
