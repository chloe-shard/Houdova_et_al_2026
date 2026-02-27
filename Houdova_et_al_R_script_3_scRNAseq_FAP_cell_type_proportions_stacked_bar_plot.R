# ------------------------------------------------------------------------------------------------------------------
#             Quantify cell type proportions of FAP expressing cells per dataset - stacked bar plot
# ------------------------------------------------------------------------------------------------------------------


#Load libraries 
library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)

# Seurat objects
seurat_objects <- list(
  Ebert = obj_Ebert,
  LeBlanc = obj_LeBlanc,
  Abdelfattah = obj_Abdelfattah,
  Nomura = obj_Nomura
)

# Compute per-patient cell-level contributions
calculate_FAP_proportions <- function(seurat_obj, dataset_name) {
  
  seurat_obj <- subset(seurat_obj, subset = FAP > 0)
  
  fap_expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")["FAP", ]
  
  meta <- seurat_obj@meta.data %>%
    mutate(FAP_expr = as.numeric(fap_expr)) %>%
    filter(FAP_expr > 0)
  
  proportions <- meta %>%
    group_by(Patient_ID, Pooled_Cell_Type) %>%
    summarise(
      n_cells = n(),
      total_expr = sum(FAP_expr),
      .groups = "drop"
    ) %>%
    group_by(Patient_ID) %>%
    mutate(
      prop_cells = n_cells / sum(n_cells),
      prop_expr  = total_expr / sum(total_expr),
      dataset    = dataset_name
    ) %>%
    ungroup()
  
  proportions
}

results_list <- imap_dfr(seurat_objects, calculate_FAP_proportions)

# View per-patient result
head(results_list)

# Average per dataset
dataset_summary <- results_list %>%
  group_by(dataset, Pooled_Cell_Type) %>%
  summarise(
    avg_prop_cells = mean(prop_cells),
    avg_prop_expr = mean(prop_expr),
    .groups = "drop"
  )

# View final dataset-level summary
head(dataset_summary)

# Normalize proportions so each dataset sums to 1
dataset_summary_norm <- dataset_summary %>%
  group_by(dataset) %>%
  mutate(
    norm_prop_cells = avg_prop_cells / sum(avg_prop_cells),
    norm_prop_expr = avg_prop_expr / sum(avg_prop_expr)
  ) %>%
  ungroup()

custom_colors <- c(
  "Endothelial"       = "#F8766D",
  "Pericyte"          = "#00A9FF",
  "Oligodendrocyte"   = "#0CB702",
  "Malignant"         = "#C77CFF",
  "Myeloid"           = "#E68613",
  "Lymphocyte"        = "#ABA300",
  "OPC_non_malignant" = "#00BFC4",
  "Astrocyte"         = "#FF61CC",
  "Other"             = "#999999",
  "Neuron"            = "#66C2A5"
)

ggplot(dataset_summary_norm, aes(x = dataset, y = norm_prop_cells, fill = Pooled_Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = custom_colors,
    breaks = names(custom_colors),
    drop = FALSE                     
  ) +
  labs(
    title = "Normalised proportion of FAP+ cells by cell type",
    x = "Dataset",
    y = "Proportion of FAP+ Cells (Normalised)",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
