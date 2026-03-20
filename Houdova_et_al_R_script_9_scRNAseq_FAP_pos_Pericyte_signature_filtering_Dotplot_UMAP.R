# -----------------------------------------------------------------------------------------------------------------------
#     Creating filtered FAP+ Pericyte signature for downstream Bulk seq analysis - marker Dotplot and signature UMAP
# -----------------------------------------------------------------------------------------------------------------------

#load libaries 
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)


################################## Annotate FAP pos vs FAP neg pericytes ##################################

annotate_FAP_pericytes <- function(obj,
                                   pericyte_label = "Pericyte",
                                   celltype_col = "Pooled_Cell_Type",
                                   fap_var = "FAP",
                                   threshold = 0) {
  md <- obj@meta.data
  md$FAP_pericyte_status <- NA_character_
  is_pericyte <- as.character(md[[celltype_col]]) == pericyte_label
  fap_vec <- FetchData(obj, vars = fap_var)[, 1]
  # (FetchData returns a 1-col data.frame; pull the vector)
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

######################### FAP+ PC signature expression across all cell populations ##############

target_cell_type <- "FAP_pos_pericyte"

# Exclude FAP_neg_pericyte population from from OFF-TARGET evaluation
exclude_cell_types <- c("FAP_neg_pericyte")

# Off-target limits (z-score + %)
max_other_z        <- 0
max_other_pct_expr <- 15

# "OR"  = Off target gene if either z or % exceed
# "AND" = off target gene only if both exceed
offtarget_rule <- "AND"

# Filter-out rule based on fraction of available datasets for that cell type
bad_fraction_threshold <- 0.50

# Require at least this many datasets for a given cell_type before using it to drop genes
min_celltype_datasets_required <- 1

# Only keep genes that exist in at least this many datasets overall
min_datasets_required <- 2

## Filtering code

# Off-target rows only (exclude target + excluded cell type)
offtarget_flags <- combined_plot_data %>%
  mutate(gene = as.character(gene)) %>%
  filter(!cell_type %in% c(target_cell_type, exclude_cell_types)) %>%
  group_by(gene, cell_type, dataset) %>%
  summarise(
    z   = z_score[1],
    pct = percentage_expr[1],
    .groups = "drop"
  ) %>%
  mutate(
    offtarget_bad = if (offtarget_rule == "OR") {
      (z > max_other_z) | (pct > max_other_pct_expr)
    } else {
      (z > max_other_z) & (pct > max_other_pct_expr)
    }
  )

# For each cell_type, compute fraction of datasets where each gene is bad
offtarget_bad_fraction <- offtarget_flags %>%
  group_by(gene, cell_type) %>%
  summarise(
    n_datasets_celltype = n(),
    n_bad = sum(offtarget_bad, na.rm = TRUE),
    frac_bad = n_bad / n_datasets_celltype,
    .groups = "drop"
  ) %>%
  filter(n_datasets_celltype >= min_celltype_datasets_required)

# 3) Mark genes to DROP if ANY off-target cell type has frac_bad >= threshold
genes_to_drop <- offtarget_bad_fraction %>%
  group_by(gene) %>%
  summarise(
    worst_frac_bad = max(frac_bad, na.rm = TRUE),
    worst_celltype_n = n_datasets_celltype[which.max(frac_bad)][1],
    worst_celltype = cell_type[which.max(frac_bad)][1],
    drop_gene = any(frac_bad >= bad_fraction_threshold),
    .groups = "drop"
  )

# 4) Enforce "gene must have at least X datasets overall"
gene_dataset_counts <- combined_plot_data %>%
  mutate(gene = as.character(gene)) %>%
  group_by(gene) %>%
  summarise(
    n_datasets_total = n_distinct(dataset),
    .groups = "drop"
  )

# 5) Final filtered gene list (KEPT genes)
filtered_genes <- gene_dataset_counts %>%
  left_join(genes_to_drop, by = "gene") %>%
  mutate(
    drop_gene = ifelse(is.na(drop_gene), FALSE, drop_gene),
    worst_frac_bad = ifelse(is.na(worst_frac_bad), 0, worst_frac_bad),
    worst_celltype_n = ifelse(is.na(worst_celltype_n), 0, worst_celltype_n),
    worst_celltype = ifelse(is.na(worst_celltype), NA_character_, worst_celltype)
  ) %>%
  filter(
    n_datasets_total >= min_datasets_required,
    drop_gene == FALSE
  ) %>%
  arrange(desc(n_datasets_total), desc(worst_frac_bad))

filtered_genes

# Optional: inspect drop reasons (fraction-based)
drop_reasons_fraction <- offtarget_bad_fraction %>%
  filter(frac_bad >= bad_fraction_threshold) %>%
  arrange(desc(frac_bad), desc(n_datasets_celltype), gene, cell_type)

drop_reasons_fraction

######################### FAP+ PC signature expression across all cell populations ##############

# specify FAP Pericyte signature genes to plot
genes_of_interest <- c("FAP", "LTBP2", "CYP1B1", "CLMP", "FBLN2", "OGN", "FMOD", "BICC1", "CXCL6", "SFRP2")

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
custom_gene_order <- c("FAP", "LTBP2", "CYP1B1", "CLMP", "FBLN2", "OGN", "FMOD", "BICC1", "CXCL6", "SFRP2")
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
  filename = "Dotplot_All_cells_FAP_Pericyte_signature.svg",
  plot     = p,
  device   = "svg",
  width    = 8,
  height   = 6,
  units    = "in"
)

################ Add signature to seurat objects 

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
})

# -----------------------------
# USER INPUT
# -----------------------------
FAP_PC_filtered <- c(
  "FAP", "LTBP2", "CYP1B1", "CLMP",
  "FBLN2", "OGN", "FMOD", "BICC1",
  "CXCL6", "SFRP2"
)

seurat_objects <- list(
  Ebert       = obj_Ebert,
  Abdelfattah = obj_Abdelfattah,
  LeBlanc     = obj_LeBlanc,
  Nomura      = obj_Nomura
)

assay_use <- "RNA"
score_name <- "FAP_PC_Score"

# ADD MODULE SCORE
seurat_objects <- imap(seurat_objects, function(seurat_obj, dataset_name) {
  DefaultAssay(seurat_obj) <- assay_use
  genes_present  <- intersect(FAP_PC_filtered, rownames(seurat_obj))
  genes_missing  <- setdiff(FAP_PC_filtered, genes_present)
  message("\nDataset: ", dataset_name)
  message("Genes present: ", length(genes_present))
  if (length(genes_missing) > 0) {
    message("Missing genes: ", paste(genes_missing, collapse = ", "))
  }
  seurat_obj <- AddModuleScore(
    object  = seurat_obj,
    features = list(genes_present),
    name    = score_name,
    assay   = assay_use
  )
  seurat_obj[[score_name]] <- seurat_obj[[paste0(score_name, "1")]]
  seurat_obj[[paste0(score_name, "1")]] <- NULL
  return(seurat_obj)
})

obj_Ebert       <- seurat_objects$Ebert
obj_Abdelfattah <- seurat_objects$Abdelfattah
obj_LeBlanc     <- seurat_objects$LeBlanc
obj_Nomura      <- seurat_objects$Nomura

# Add pericyte metadata to original object

# Dataset prefixes
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

for (prefix in datasets) {
  full_name <- paste0("obj_", prefix)
  peri_name <- paste0("obj_", prefix, "_Pericyte")
  full_obj     <- get(full_name)
  pericyte_obj <- get(peri_name)
  new_annotation <- setNames(rep("None", ncol(full_obj)), colnames(full_obj))
  pericyte_annotations <- as.character(pericyte_obj@meta.data$Pericyte_subpop_smooth)
  names(pericyte_annotations) <- rownames(pericyte_obj@meta.data)
  pericyte_annotations[pericyte_annotations %in% c("None", NA)] <- "Pericyte_None"
  overlap_cells <- intersect(names(new_annotation), names(pericyte_annotations))
  new_annotation[overlap_cells] <- pericyte_annotations[overlap_cells]
  full_obj$pericyte_subpop <- new_annotation
  assign(full_name, full_obj)
}

# Plot
p1 <- DimPlot(obj_Ebert, group.by = "pericyte_subpop")
p2 <- FeaturePlot(obj_Ebert, features = "FAP_PC_Score")

score_plot <- p1+p2

score_plot

# ----------------------------
# Save
# ----------------------------
ggsave(
  filename = "Pericyte_filtered_score_UMAPs.svg",
  plot     = score_plot,
  device   = "svg",
  width    = 12,
  height   = 5,
  units    = "in"
)
