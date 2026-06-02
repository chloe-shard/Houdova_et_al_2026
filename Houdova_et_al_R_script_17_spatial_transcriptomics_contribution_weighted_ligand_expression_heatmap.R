# -----------------------------------------------------------------------------------------
#                      Contribution weighted ligand expression in FAP+PC enriched niches
# -----------------------------------------------------------------------------------------

#load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(stringr)
library(patchwork)
library(purrr)
library(broom)
library(rstatix)

#### Step 1: calculate normalized expression of each ligand per cell type per dataset ####

# Ligands of interest
ligand_genes <- c("CCL2", "CSF1", "CXCL1", "CXCL12","RARRES2", "SEMA3C", "SERPINE1", "THBS2", "WNT5A")

# Cell type metadata column
celltype_col <- "Pooled_Cell_Type_FAP_updatedV3"

# list of Seurat objects (annotated with cell subtypes - see Houdova_et_al_R_script_11_scRNAseq_cell_subtype_annotations)
seurat_list <- list(
  Ebert = obj_Ebert,
  LeBlanc = obj_LeBlanc,
  Abdelfattah = obj_Abdelfattah,
  Nomura = obj_Nomura
)

# Cell types to exclude from plotting [cell types found only in one dataset, with low cell numbers or not belonging to a defined subtype - i.e. "None"]
exclude_celltypes_plot <- c( "OPC_non_malignant", 
                             "Astrocyte", 
                             "Neuron", 
                             "Other",
                             "None_Malignant",
                             "None_Myeloid",
                             "Other",
                             "Neuron")

# Function to calculate normalised expression
calculate_ligand_norm <- function(seurat_obj, dataset_name, ligand_genes, celltype_col) {
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    warning(paste("Skipping", dataset_name, "- metadata column not found:", celltype_col))
    return(NULL)
  }
  genes_present <- ligand_genes[ligand_genes %in% rownames(seurat_obj)]
  genes_missing <- setdiff(ligand_genes, genes_present)
  
  if (length(genes_missing) > 0) {
    warning(
      paste(
        "In", dataset_name, "genes were not found and will be skipped:",
        paste(genes_missing, collapse = ", ")
      )
    )
  }
  
  if (length(genes_present) == 0) {
    warning(paste("Skipping", dataset_name, "- no ligand genes present."))
    return(NULL)
  }
  
  expression_df <- FetchData(
    seurat_obj,
    vars = c(genes_present, celltype_col)
  )
  
  expression_long <- expression_df %>%
    rename(cell_type = all_of(celltype_col)) %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "gene",
      values_to = "expression"
    )
  
  global_means <- expression_long %>%
    group_by(gene) %>%
    summarise(
      global_mean_expression = mean(expression, na.rm = TRUE),
      .groups = "drop"
    )
  
  norm_by_celltype <- expression_long %>%
    group_by(cell_type, gene) %>%
    summarise(
      avg_expression = mean(expression, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    left_join(global_means, by = "gene") %>%
    mutate(
      dataset = dataset_name,
      normalised_expression = avg_expression / global_mean_expression
    ) %>%
    select(
      dataset,
      cell_type,
      gene,
      n_cells,
      avg_expression,
      global_mean_expression,
      normalised_expression
    )
  
  return(norm_by_celltype)
}

# Run across all datasets and merge into one table
ligand_norm_all_datasets <- bind_rows(
  lapply(
    names(seurat_list),
    function(dataset_name) {
      calculate_ligand_norm(
        seurat_obj = seurat_list[[dataset_name]],
        dataset_name = dataset_name,
        ligand_genes = ligand_genes,
        celltype_col = celltype_col
      )
    }
  )
)

# View full merged table with all cell types retained
ligand_norm_all_datasets

# filter cells for the heatmap
ligand_norm_for_heatmap <- ligand_norm_all_datasets %>%
  filter(!cell_type %in% exclude_celltypes_plot) %>%
  mutate(
    dataset_celltype = paste(dataset, cell_type, sep = " | ")
  )

# Set the dataset order for plotting
dataset_order <- c("Nomura", "Abdelfattah", "LeBlanc", "Ebert")

# Set ligand/gene order for plotting
gene_order <- c("CCL2", "CSF1", "CXCL1", "CXCL12","RARRES2", "SEMA3C", "SERPINE1", "THBS2", "WNT5A")

# Set cell type order for plotting
celltype_order <- c(
  "Lymphocyte",
  "Microglia",
  "Monocyte",
  "M1",
  "NPC",
  "Oligodendrocyte",
  "OPC",
  "FAP_neg_pericyte",
  "M2",
  "Endothelial",
  "FAP_pos_pericyte",
  "AC",
  "MES"
)

# Prepare heatmap dataframe
ligand_norm_for_heatmap_grouped <- ligand_norm_for_heatmap %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    cell_type = factor(cell_type, levels = celltype_order),
    
    celltype_dataset = paste(cell_type, dataset, sep = " | ")
  ) %>%
  arrange(cell_type, dataset) %>%
  mutate(
    celltype_dataset = factor(
      celltype_dataset,
      levels = unique(celltype_dataset)
    ),
    gene = factor(
      gene,
      levels = rev(gene_order)
    )
  )

# plot heatmap
heatmap_ligand_normalised <- ligand_norm_for_heatmap_grouped %>%
  ggplot(
    aes(
      x = celltype_dataset,
      y = gene,
      fill = normalised_expression
    )
  ) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint =5,
    na.value = "grey90"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Ligand expression normalised within each dataset",
    fill = "Normalised\nexpression"
  )

#print heatmap
heatmap_ligand_normalised

# export 
ggsave(
  filename = "heatmap_ligand_normalised.svg",
  plot = heatmap_ligand_normalised,
  device = svglite::svglite,
  width = 16,
  height = 8,
  units = "in"
)

#### Step 2: calculate mean normalized expression of each ligand per cell type across all 4 datasets ####

library(dplyr)
library(ggplot2)

# Cell types to exclude from plotting [cell types found only in one dataset, with low cell numbers or not belonging to a defined subtype - i.e. "None"]
exclude_celltypes_plot <- c("OPC_non_malignant", 
                            "Astrocyte", 
                            "Neuron", 
                            "Other",
                            "None_Malignant",
                            "None_Myeloid",
                            "Other",
                            "Neuron")

# Average normalised expression across datasets
Ligand_average_expression <- ligand_norm_all_datasets %>%
  filter(!cell_type %in% exclude_celltypes_plot) %>%
  group_by(cell_type, gene) %>%
  summarise(
    mean_normalised_expression = mean(normalised_expression, na.rm = TRUE),
    sd_normalised_expression = sd(normalised_expression, na.rm = TRUE),
    n_datasets = n_distinct(dataset[!is.na(normalised_expression)]),
    .groups = "drop"
  )

# Set gene order
gene_order <- c("CCL2", "CSF1", "CXCL1", "CXCL12","RARRES2", "SEMA3C", "SERPINE1", "THBS2", "WNT5A")

# Set cell type order
celltype_order <- c(
  "Lymphocyte",
  "Microglia",
  "Monocyte",
  "M1",
  "NPC",
  "Oligodendrocyte",
  "M2",
  "FAP_pos_pericyte",
  "OPC",
  "Endothelial",
  "FAP_neg_pericyte",
  "MES",
  "AC"
)

# plot heatmap of mean normalized expression
heatmap_ligand_average <- Ligand_average_expression %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    gene = factor(gene, levels = rev(gene_order))
  ) %>%
  ggplot(
    aes(
      x = cell_type,
      y = gene,
      fill = mean_normalised_expression
    )
  ) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 3,
    na.value = "grey90"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Mean ligand expression normalised across datasets",
    fill = "Mean normalised\nexpression"
  )

# print heatmap
heatmap_ligand_average

# export heatmap
ggsave(
  filename = "Ligand_heatmap_mean_normalised_exp.pdf",
  plot = heatmap_ligand_average,
  width = 8,
  height = 8
)

#### Step 3: Calculate cell proportions in FAP_PCs enriched spots (does not include neighbour spots) ####

# Seurat object names (Spatial T objects annotated with cell subtypes using RCTD analysis - see R file: Houdova_et_al_R_script_13_spatial_transcriptomics_RCTD_cell_type_annotation)
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Metadata columns containing cell proportion scores
cell_prop_cols <- c(
  "MES_Houdova",
  "NPC_Houdova",
  "Monocyte_Houdova",
  "AC_Houdova",
  "Endothelial_Houdova",
  "FAP_neg_pericyte_Houdova",
  "FAP_pos_pericyte_Houdova",
  "Lymphocyte_Houdova",
  "M1_Houdova",
  "M2_Houdova",
  "Microglia_Houdova",
  "Oligodendrocyte_Houdova",
  "OPC_Houdova"
)

# Thresholded spots to analyse (niche of interest)
threshold_cols <- c(
  "FAP_PCs_thresholded"
)

# Empty list to store all results
all_results_list <- list()

# Function to calculate proportions
for (threshold_col in threshold_cols) {
  message("Running analysis for: ", threshold_col)
  results_list <- list()
  for (obj_name in seurat_names) {
    message("Processing: ", obj_name)
    seurat_obj <- get(obj_name)
    if (!threshold_col %in% colnames(seurat_obj@meta.data)) {
      warning(threshold_col, " not found in ", obj_name, ". Skipping.")
      next
    }
    missing_cols <- setdiff(cell_prop_cols, colnames(seurat_obj@meta.data))
    if (length(missing_cols) > 0) {
      warning(
        "The following cell proportion columns are missing in ",
        obj_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
    present_cell_prop_cols <- intersect(cell_prop_cols, colnames(seurat_obj@meta.data))
    
    if (length(present_cell_prop_cols) == 0) {
      warning("No cell proportion columns found in ", obj_name, ". Skipping.")
      next
    }
    meta <- seurat_obj@meta.data
    
    # Identify positive spots for this threshold column
    positive_spots <- rownames(meta)[meta[[threshold_col]] == "positive"]
    if (length(positive_spots) == 0) {
      warning("No positive spots found for ", threshold_col, " in ", obj_name, ". Skipping.")
      next
    }
    avg_props <- colMeans(
      meta[positive_spots, present_cell_prop_cols, drop = FALSE],
      na.rm = TRUE
    )
    results_list[[obj_name]] <- data.frame(
      SeuratObject = obj_name,
      ThresholdGroup = threshold_col,
      n_positive_spots = length(positive_spots),
      t(avg_props),
      row.names = NULL
    )
  }
  all_results_list[[threshold_col]] <- bind_rows(results_list)
}

# Combine both analyses into one dataframe
positive_avg_props_all <- bind_rows(all_results_list)

# View combined result
positive_avg_props_all

# prepare datafram for plotting
plot_data_all <- positive_avg_props_all %>%
  pivot_longer(
    cols = -c(SeuratObject, ThresholdGroup, n_positive_spots),
    names_to = "CellType",
    values_to = "AverageProportion"
  )

# clean labels
plot_data_all <- plot_data_all %>%
  mutate(
    ThresholdGroup = case_when(
      ThresholdGroup == "FAP_PCs_thresholded" ~ "FAP_PCs enriched spots",
      TRUE ~ ThresholdGroup
    )
  )

# Order cell types by median average proportion
plot_data_all <- plot_data_all %>%
  mutate(
    CellType = reorder(CellType, AverageProportion, FUN = median, na.rm = TRUE)
  )

# Plot boxplots for FAP_PCs enriched spots only
FAP_PCs_enriched_boxplot <- ggplot(
  plot_data_all,
  aes(x = CellType, y = AverageProportion)
) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey80",
    colour = "black"
  ) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    title = "Cell proportions in FAP_PCs enriched spots",
    x = "Cell type",
    y = "Average proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# print boxplot
FAP_PCs_enriched_boxplot

# export
ggsave(
  filename = "FAP_PCs_enriched_boxplot.svg",
  plot = FAP_PCs_enriched_boxplot,
  device = "svg",
  width = 12,
  height = 6,
  units = "in"
)

# Output mean proportion
mean_celltype_props_by_group <- plot_data_all %>%
  group_by(ThresholdGroup, CellType) %>%
  summarise(
    mean_proportion = mean(AverageProportion, na.rm = TRUE),
    median_proportion = median(AverageProportion, na.rm = TRUE),
    sd_proportion = sd(AverageProportion, na.rm = TRUE),
    sem_proportion = sd(AverageProportion, na.rm = TRUE) / sqrt(sum(!is.na(AverageProportion))),
    n_tissues = sum(!is.na(AverageProportion)),
    .groups = "drop"
  ) %>%
  arrange(ThresholdGroup, desc(mean_proportion))

# export mean
write.csv(
  mean_celltype_props_by_group,
  "Cell_type_proportions_mean.csv",
  row.names = FALSE
)

# export all tissues
write.csv(
  plot_data_all,
  "Cell_type_proportions_per_tissue.csv",
  row.names = FALSE
)

#### Step 4: Calculate contribution weighted expression score ####

# load in cell proportions csv files
Cell_type_proportions_mean <- read_csv("Path to file/Cell_type_proportions_mean.csv")
Cell_type_proportions_per_tissue <- read_csv("Path to file/Cell_type_proportions_per_tissue.csv")

# 1. Check exact ThresholdGroup labels
Cell_type_proportions_mean %>%
  distinct(ThresholdGroup)

# Set spots of interest (FAP_PCs_pos)
threshold_group_of_interest <- "FAP_PCs positive spots"

# 2. Create cell type mapping to merge the nomenclature from single cell and spatial data
celltype_map <- tibble(
  cell_type = c(
    "AC",
    "Endothelial",
    "FAP_neg_pericyte",
    "FAP_pos_pericyte",
    "Lymphocyte",
    "M1",
    "M2",
    "MES",
    "Microglia",
    "Monocyte",
    "NPC",
    "OPC",
    "Oligodendrocyte"
  ),
  CellType = c(
    "AC_Houdova",
    "Endothelial_Houdova",
    "FAP_neg_pericyte_Houdova",
    "FAP_pos_pericyte_Houdova",
    "Lymphocyte_Houdova",
    "M1_Houdova",
    "M2_Houdova",
    "MES_Houdova",
    "Microglia_Houdova",
    "Monocyte_Houdova",
    "NPC_Houdova",
    "OPC_Houdova",
    "Oligodendrocyte_Houdova"
  )
)

# 3. Extract mean proportions from FAP_PCs_pos spots
FAP_PCs_positive_proportions <- Cell_type_proportions_mean %>%
  filter(ThresholdGroup == threshold_group_of_interest) %>%
  select(
    CellType,
    mean_proportion,
    median_proportion,
    sd_proportion,
    sem_proportion,
    n_tissues
  )

# 4. Join proportions onto ligand table and weight expression
Ligand_average_expression_normalised <- Ligand_average_expression %>%
  mutate(
    cell_type = as.character(cell_type)
  ) %>%
  left_join(
    celltype_map,
    by = "cell_type"
  ) %>%
  left_join(
    FAP_PCs_positive_proportions,
    by = "CellType"
  ) %>%
  mutate(
    mean_normalised_expression_weighted_by_mean_proportion =
      mean_normalised_expression * mean_proportion
  )

# 5. Check output
Ligand_average_expression_normalised %>%
  select(
    cell_type,
    CellType,
    gene,
    mean_normalised_expression,
    mean_proportion,
    mean_normalised_expression_weighted_by_mean_proportion,
    sd_normalised_expression,
    n_datasets
  ) %>%
  head()

# Check whether any cell types failed to match
Ligand_average_expression_normalised %>%
  filter(is.na(mean_proportion)) %>%
  distinct(cell_type, CellType)

# Prepare heatmap dataframe
heatmap_df <- Ligand_average_expression_normalised %>%
  select(
    cell_type,
    gene,
    mean_normalised_expression,
    mean_normalised_expression_weighted_by_mean_proportion
  ) %>%
  pivot_longer(
    cols = c(
      mean_normalised_expression,
      mean_normalised_expression_weighted_by_mean_proportion
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mean_normalised_expression =
        "Original mean normalised expression",
      mean_normalised_expression_weighted_by_mean_proportion =
        "Contribution-weighted mean normalised expression"
    )
  )

# Set cell type order for plotting
celltype_order <- c(
  "Lymphocyte",
  "Microglia",
  "Monocyte",
  "M1",
  "NPC",
  "Oligodendrocyte",
  "M2",
  "FAP_pos_pericyte",
  "OPC",
  "Endothelial",
  "FAP_neg_pericyte",
  "MES",
  "AC"
)

# set gene order for plotting
gene_order <- Ligand_average_expression_normalised %>%
  distinct(gene) %>%
  pull(gene)

# plot contribution-weighted heatmap
heatmap_weighted_only <- Ligand_average_expression_normalised %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    gene = factor(gene, levels = rev(gene_order))
  ) %>%
  ggplot(
    aes(
      x = cell_type,
      y = gene,
      fill = mean_normalised_expression_weighted_by_mean_proportion
    )
  ) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0.3,
    na.value = "grey90"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Contribution-weighted ligand normalised expression",
    fill = "Weighted\nnormalised expression"
  )

#print heatmap
heatmap_weighted_only

# export
ggsave(
  filename = "Ligand_heatmap_average_normalised_and_contibution_weighted_exp.svg",
  plot = heatmap_ligand_average + heatmap_weighted_only,
  device = svglite::svglite,
  width = 16,
  height = 8,
  units = "in"
)


#### Step 5: Statistics on contribution weighted expression - ANOVA ####

# 1. Define threshold group of interest
threshold_group_of_interest <- "FAP_PCss positive spots"

# 2. Create helper functions
get_significance <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

is_valid_anova_df <- function(df) {
  n_distinct(df$cell_type) >= 2 &&
    nrow(df) > n_distinct(df$cell_type)
}

# 3. Create cell type mapping to merge the nomenclature from single cell and spatial data
celltype_map <- tibble(
  cell_type = c(
    "AC",
    "Endothelial",
    "FAP_neg_pericyte",
    "FAP_pos_pericyte",
    "Lymphocyte",
    "M1",
    "M2",
    "MES",
    "Microglia",
    "Monocyte",
    "NPC",
    "OPC",
    "Oligodendrocyte"
  ),
  CellType = c(
    "AC_Houdova",
    "Endothelial_Houdova",
    "FAP_neg_pericyte_Houdova",
    "FAP_pos_pericyte_Houdova",
    "Lymphocyte_Houdova",
    "M1_Houdova",
    "M2_Houdova",
    "MES_Houdova",
    "Microglia_Houdova",
    "Monocyte_Houdova",
    "NPC_Houdova",
    "OPC_Houdova",
    "Oligodendrocyte_Houdova"
  )
)

# 4. Prepare per-tissue cell type proportions
FAP_PCs_positive_proportions_per_tissue <- Cell_type_proportions_per_tissue %>%
  filter(ThresholdGroup == threshold_group_of_interest) %>%
  select(
    SeuratObject,
    ThresholdGroup,
    n_positive_spots,
    CellType,
    AverageProportion
  )

# Check number of tissues per cell type
FAP_PCs_positive_proportions_per_tissue %>%
  count(CellType, name = "n_tissues") %>%
  arrange(CellType)

# 5. Join ligand expression to per-tissue proportions
Ligand_weighted_expression_per_tissue <- Ligand_average_expression %>%
  mutate(
    cell_type = as.character(cell_type)
  ) %>%
  left_join(
    celltype_map,
    by = "cell_type"
  ) %>%
  left_join(
    FAP_PCs_positive_proportions_per_tissue,
    by = "CellType"
  ) %>%
  mutate(
    weighted_expression =
      mean_normalised_expression * AverageProportion
  ) %>%
  filter(
    !is.na(SeuratObject),
    !is.na(AverageProportion),
    !is.na(weighted_expression)
  )

# Check output
Ligand_weighted_expression_per_tissue %>%
  select(
    SeuratObject,
    gene,
    cell_type,
    CellType,
    mean_normalised_expression,
    AverageProportion,
    weighted_expression
  ) %>%
  head()

# Cell types from ligand table that did not map to spatial CellType names
Ligand_average_expression %>%
  mutate(cell_type = as.character(cell_type)) %>%
  left_join(celltype_map, by = "cell_type") %>%
  filter(is.na(CellType)) %>%
  distinct(cell_type)

# Ligand/cell types that did not get per-tissue proportions
Ligand_average_expression %>%
  mutate(cell_type = as.character(cell_type)) %>%
  left_join(celltype_map, by = "cell_type") %>%
  left_join(FAP_PCs_positive_proportions_per_tissue, by = "CellType") %>%
  filter(is.na(AverageProportion)) %>%
  distinct(cell_type, CellType)

# 6. Create clean weighted expression object for statistics
Ligand_weighted_expression_clean <- Ligand_weighted_expression_per_tissue %>%
  filter(
    !is.na(weighted_expression),
    !is.na(cell_type)
  )

# 7. ANOVA per ligand
weighted_anova_by_ligand <- Ligand_weighted_expression_clean %>%
  group_by(gene) %>%
  group_modify(~ {
    
    df <- .x
        if (!is_valid_anova_df(df)) {
      return(tibble(
        term = NA_character_,
        df = NA_real_,
        sumsq = NA_real_,
        meansq = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_
      ))
    }
    
    aov(weighted_expression ~ cell_type, data = df) %>%
      broom::tidy()
  }) %>%
  ungroup() %>%
  filter(term == "cell_type") %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    significance = get_significance(p_adj_BH)
  ) %>%
  arrange(p_adj_BH)

# Print ANOVA results
weighted_anova_by_ligand

# 8. Tukey post hoc test per ligand
weighted_tukey_by_ligand <- Ligand_weighted_expression_clean %>%
  group_by(gene) %>%
  group_modify(~ {
    
    df <- .x
    
    if (!is_valid_anova_df(df)) {
      return(tibble())
    }
    
    fit <- aov(weighted_expression ~ cell_type, data = df)
    tukey <- TukeyHSD(fit, "cell_type")
    
    as.data.frame(tukey$cell_type) %>%
      tibble::rownames_to_column("comparison") %>%
      separate(
        comparison,
        into = c("group1", "group2"),
        sep = "-"
      ) %>%
      rename(
        mean_difference = diff,
        conf_low = lwr,
        conf_high = upr,
        p_adj_Tukey = `p adj`
      )
  }) %>%
  ungroup() %>%
  mutate(
    significance = get_significance(p_adj_Tukey)
  ) %>%
  arrange(gene, p_adj_Tukey)

# Print Tukey results
weighted_tukey_by_ligand

# 9. Prepare Tukey results for ligand-cell type significance summary
weighted_tukey_sig_input <- weighted_tukey_by_ligand %>%
  filter(!is.na(p_adj_Tukey))

# Create one row per ligand-celltype-comparison direction
tukey_long <- bind_rows(
  weighted_tukey_sig_input %>%
    transmute(
      gene,
      cell_type = group1,
      compared_cell_type = group2,
      p_adj_Tukey,
      mean_difference,
      significant_p001 = p_adj_Tukey < 0.001
    ),
  
  weighted_tukey_sig_input %>%
    transmute(
      gene,
      cell_type = group2,
      compared_cell_type = group1,
      p_adj_Tukey,
      mean_difference = -mean_difference,
      significant_p001 = p_adj_Tukey < 0.001
    )
)

# 10. Annotate significant ligand-cell type combinations
tukey_sig_summary_by_ligand_celltype <- tukey_long %>%
  group_by(gene, cell_type) %>%
  summarise(
    n_comparisons_total = n_distinct(compared_cell_type),
    n_comparisons_p001 = sum(significant_p001, na.rm = TRUE),
    
    significant_vs_all_other_celltypes_p001 =
      n_comparisons_p001 == n_comparisons_total,
    
    significant_vs_at_least_10_celltypes_p001 =
      n_comparisons_p001 >= 10,
    
    significant_vs_at_least_3_celltypes_p001 =
      n_comparisons_p001 >= 8,
    
    significant_comparisons_p001 = paste(
      unique(compared_cell_type[significant_p001]),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    significance_category = case_when(
      significant_vs_all_other_celltypes_p001 ~
        "Significant compared to all other cell types, p_adj_Tukey < 0.001",
      
      significant_vs_at_least_10_celltypes_p001 ~
        "Significant compared to at least 10 cell types, p_adj_Tukey < 0.001",
      
      significant_vs_at_least_8_celltypes_p001 ~
        "Significant compared to at least 8 cell types, p_adj_Tukey < 0.001",
      
      TRUE ~
        "Not significant compared to at least 8 cell types at p_adj_Tukey < 0.001"
    )
  ) %>%
  arrange(
    gene,
    desc(significant_vs_all_other_celltypes_p001),
    desc(significant_vs_at_least_10_celltypes_p001),
    desc(significant_vs_at_least_8_celltypes_p001),
    desc(n_comparisons_p001),
    cell_type
  )

# print results
tukey_sig_summary_by_ligand_celltype

# export
write.csv(
  tukey_sig_summary_by_ligand_celltype,
  "tukey_sig_summary_by_ligand_celltype.csv",
  row.names = FALSE
)