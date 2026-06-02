# -----------------------------------------------------------------------------------------------------------------------
#                   Confirm author annotation of Pericyte population using classic pericyte markers - Dotplot
# -----------------------------------------------------------------------------------------------------------------------

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Pericyte marker genes (obtained from Ebert et al., 2020)
genes_of_interest <- c("PDGFRB", "THY1", "ACTA2", "CSPG4")

# calculate average expression (mean), z-score, and percentage expression per dataset
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
  expression_df$Pooled_Cell_Type <- seurat_obj@meta.data$Pooled_Cell_Type
  avg_expression <- expression_df %>%
    group_by(Pooled_Cell_Type) %>%
    summarise_at(vars(all_of(genes_present)), mean, na.rm = TRUE) %>%
    pivot_longer(-Pooled_Cell_Type, names_to = "gene", values_to = "avg_expression")
  avg_expression_long <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = scale(avg_expression)) %>%
    ungroup()
  binary_expr <- assay_data[genes_present, ] > 0 
  binary_expr <- as.data.frame(t(binary_expr))  
  binary_expr$Pooled_Cell_Type <- seurat_obj@meta.data$Pooled_Cell_Type
  percentage_expr <- binary_expr %>%
    group_by(Pooled_Cell_Type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-Pooled_Cell_Type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "Pooled_Cell_Type")) %>%
    mutate(cell_type_dataset = paste(Pooled_Cell_Type, dataset_name, sep = "_"))
  return(plot_data)
}

# Process each dataset
plot_data_Ebert <- process_dataset(obj_Ebert, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah, "Adelfattah")
plot_data_Nomura <- process_dataset(obj_Nomura, "Nomura")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Nomura,
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
)

# Specify the order of genes and cell types for plotting
custom_gene_order <- c("PDGFRB", "THY1", "ACTA2", "CSPG4")
custom_cellstate_order <- c("Pericyte", "Endothelial", "Myeloid", "Lymphocyte", "Oligodendrocyte", "OPC_non_malignant", "Astrocyte", "Neuron", "Malignant", "Other")
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(Pooled_Cell_Type) {
      grep(paste0("^", Pooled_Cell_Type, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
    })
  )
)

# Generate dot plot
Dotplot <- ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = percentage_expr, color = z_score)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Gene", y = "Cell Type (Dataset)", size = "Percentage Expressing", color = "Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# export
ggsave(
  filename = "Dotplot_Pericyte_cluster_markers_all_datasets.svg",
  plot     = Dotplot,
  device   = "svg",
  width    = 5,
  height   = 6,
  units    = "in"
)