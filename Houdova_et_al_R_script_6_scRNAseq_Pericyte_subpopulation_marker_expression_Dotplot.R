# ------------------------------------------------------------------------------------------------------------------
#                                                   Marker expression pericyte subpops
# ------------------------------------------------------------------------------------------------------------------

# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# OPTION 1 to plot

#Bejarano et al 2025 markers
PC_Transport_B <- c("SLC6A13", "SLC6A12", "ABCA9")
PC_ECM_B <- c("COL18A1", "OLFML2B", "PRSS23")
PC_interferon_B <- c("IFIT1", "IFIT2", "RSAD2")
PC_SMC_type1_B <- c("MYH11", "CASQ2", "ACTA2")
PC_SMC_type2_B <- c("PRKAA2", "NEBL", "ADIRF")
PC_Fibroblast_B <- c("CTHRC1", "ISLR", "PDGFRA", "FAP")
PC_Prolifrative_B <- c("UBE2C", "MKI67", "TOP2A")

genes_of_interest <- unique(c(
  PC_Transport_B,
  PC_ECM_B,
  PC_interferon_B,
  PC_SMC_type1_B,
  PC_SMC_type2_B,
  PC_Fibroblast_B,
  PC_Prolifrative_B
))

# OPTION 2 to plot

#New markers identified in this project
Transport_PCs <- c("ITIH5", "P2RY14", "SLC6A13", "SLC6A12", "COL4A3", "ADARB2", "SLC30A10", "SLC19A1", "DIO3OS", "PTGDR2", "SHISA6")
ECM_PCs <- c("PRR16", "PIEZO2", "PRKG1", "THY1", "FAT1", "CD36", "IL1RAPL1", "NREP", "C2orf27A", "BMP1", "TGFBI", "TUSC3", "DOK6", "EBF2")
Interferon_PCs <- c("CSF2RA", "SAMSN1", "CD53", "DOCK2", "ADAM28", "PTPRC", "SLA", "RBM47", "HLA-DRB1", "VSIG4", "RNF144B", "ITGB2", "PLEK", "SLC11A1", "FYB", "HLA-DQB1", "PIKAP1", "CD86", "MSR1", "ITGAX", "SLC16A10", "MS4A4A", "SYK", "NLRP3", "FCGR3A", "AOAH", "FCER1G", "C1QC", "LCP1")
SMCs <- c("MYH11", "CSDC2", "DES", "AC097724.3", "PLN", "LMOD1", "CNN1", "NET1", "MYOCD")
FAP_PCs <- c("POSTN", "CA12", "CXCL6", "FBLN2", "CYP1B1", "FAP", "COL11A1", "PDGFRA", "SFRP2", "FLRT2", "THBS2", "SGCD", "FMOD", "LTBP2", "CLDN11", "GAS7", "PDPN", "ELN", "OGN", "DPSL3", "CLMP", "BICC1", "ROR2", "TMTM132C", "SFRP4", "GJA1", "OMD")
Proliferative_PCs <- c("ASPM", "PBK", "AURKB", "UBE2C", "FAM64A", "KIF2C", "SKA1", "SHCBP1", "CENPF", "RRM2", "KIF23", "CDC20", "KIF11", "TOP2A", "KIF15", "MAD2L1", "DLGAP5", "CENPM", "GTSE1", "CDCA8", "ANLN", "TPX2", "KIF14", "CCNA2", "NUF2", "KIF18B", "CEP55", "CCNB1", "CKAP2L", "DIAPH3", "TK1", "SGOL1", "NUSAP1")

genes_of_interest <- unique(c(
  Transport_PCs,
  ECM_PCs,
  Interferon_PCs,
  SMCs,
  FAP_PCs,
  Proliferative_PCs
))

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
  expression_df$cell_type <- seurat_obj@meta.data$Pericyte_subpop_smooth
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
  binary_expr$cell_type <- seurat_obj@meta.data$Pericyte_subpop_smooth
  percentage_expr <- binary_expr %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "cell_type")) %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_"))
  return(plot_data)
}

# Process each dataset manually
plot_data_Nomura <- process_dataset(obj_Nomura_Pericyte, "Nomura")
plot_data_Ebert <- process_dataset(obj_Ebert_Pericyte, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc_Pericyte, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah_Pericyte, "Adelfattah")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
  plot_data_Nomura
)

# Specify the order of genes for plotting
custom_gene_order <- unique(c(
  PC_Transport_B,
  PC_ECM_B,
  PC_interferon_B,
  PC_SMC_type1_B,
  PC_SMC_type2_B,
  PC_Fibroblast_B,
  PC_Prolifrative_B
))

# Specify the order of pericyte annotations and genes for plotting
custom_cellstate_order <- c("Transport_PCs", "ECM_PCs", "Interferon_PCs", "SMCs", "FAP_PCs", "Proliferative_PCs", "None")
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(Pericyte_subpop_smooth) {
      grep(paste0("^", Pericyte_subpop_smooth, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
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

# Export
ggsave(
  filename = "Dotplot_Bejarano_2025_markers.svg",
  plot     = combined_plot,
  device   = "svg",
  width    = 10,
  height   = 8,
  units    = "in"
)
