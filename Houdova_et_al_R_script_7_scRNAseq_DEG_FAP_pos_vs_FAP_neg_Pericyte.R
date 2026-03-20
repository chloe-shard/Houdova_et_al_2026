# ------------------------------------------------------------------------------------------------------------------
#                                                   DEG FAP pos vs FAP neg Pericyte
# ------------------------------------------------------------------------------------------------------------------

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(VennDiagram)
  library(grid)
})

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

############### Differential gene expression FAP+ vs FAP- #######

# USER INPUTS
seurat_objects <- list(
  Ebert = obj_Ebert_Pericyte,
  Abdelfattah = obj_Abdelfattah_Pericyte,
  LeBlanc = obj_LeBlanc_Pericyte,
  Nomura = obj_Nomura_Pericyte
)

assay_use   <- "RNA"
fap_gene    <- "FAP"
min_pct     <- 0.25
fdr_cutoff  <- 0.05
lfc_cutoff  <- 0.5     # logFC threshold for calling up/down
k_of_n      <- 3       # "3 of 4" logic (>= 3 datasets)

# HELPERS
get_expression <- function(seurat_object, gene, assay = "RNA") {
  assay_obj <- seurat_object[[assay]]
  
  # Case 1: Old Seurat versions where @data exists
  if ("data" %in% slotNames(assay_obj)) {
    if (gene %in% rownames(assay_obj@data)) {
      return(assay_obj@data[gene, ])
    } else {
      stop(paste("Gene", gene, "not found in RNA @data"))
    }
  }
  
  # Case 2: Seurat v5 (Assay5) using layers
  layers_available <- SeuratObject::Layers(assay_obj)
  if ("data" %in% layers_available) {
    expr <- GetAssayData(seurat_object, assay = assay, layer = "data")
    if (gene %in% rownames(expr)) {
      return(expr[gene, ])
    } else {
      stop(paste("Gene", gene, "not found in RNA/data layer"))
    }
  }
  
  stop("No usable expression layer (@data or RNA/data layer) found.")
}

get_fc_col <- function(markers_df) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  fc_col <- intersect(fc_candidates, colnames(markers_df))[1]
  if (is.na(fc_col) || length(fc_col) == 0) {
    stop("Could not find a logFC column in FindMarkers output.")
  }
  fc_col
}

run_de_one_dataset <- function(seurat_object,
                               dataset_name,
                               assay = "RNA",
                               fap_gene = "FAP",
                               min_pct = 0.25,
                               fdr_cutoff = 0.05,
                               lfc_cutoff = 0.5) {
  
  # FAP status
  fap_expr <- get_expression(seurat_object, fap_gene, assay = assay)
  seurat_object$FAP_status <- ifelse(fap_expr > 0, "Positive", "Negative")
  Idents(seurat_object) <- "FAP_status"
  
  # DE
  de <- FindMarkers(
    object  = seurat_object,
    ident.1 = "Positive",
    ident.2 = "Negative",
    min.pct = min_pct
  )
  
  fc_col <- get_fc_col(de)
  
  stats_tbl <- de %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    transmute(
      gene = gene,
      logFC = .data[[fc_col]],
      pval  = p_val,
      FDR   = p_val_adj
    )
  
  up_genes <- stats_tbl %>%
    filter(FDR < fdr_cutoff, logFC >  lfc_cutoff) %>%
    pull(gene)
  
  down_genes <- stats_tbl %>%
    filter(FDR < fdr_cutoff, logFC < -lfc_cutoff) %>%
    pull(gene)
  
  list(
    dataset = dataset_name,
    stats   = stats_tbl,
    up      = up_genes,
    down    = down_genes
  )
}

# RUN ACROSS DATASETS
dge_results <- imap(seurat_objects, ~ run_de_one_dataset(
  seurat_object = .x,
  dataset_name  = .y,
  assay         = assay_use,
  fap_gene      = fap_gene,
  min_pct       = min_pct,
  fdr_cutoff    = fdr_cutoff,
  lfc_cutoff    = lfc_cutoff
))

# UP/DOWN LISTS + COMMON GENES
up_lists   <- lapply(dge_results, `[[`, "up")
down_lists <- lapply(dge_results, `[[`, "down")

# In all datasets
common_up_all   <- Reduce(intersect, up_lists)
common_down_all <- Reduce(intersect, down_lists)

# In >= k_of_n datasets (default 3 of 4)
up_counts   <- table(unlist(up_lists))
down_counts <- table(unlist(down_lists))

common_up_kofn   <- names(up_counts[up_counts >= k_of_n])
common_down_kofn <- names(down_counts[down_counts >= k_of_n])

# WIDE TABLE: logFC/pval/FDR PER DATASET FOR GENES IN >= k_of_n
wide_tbl <- reduce(dge_results, function(acc, res) {
  ds <- res$dataset
  tmp <- res$stats %>%
    rename(
      !!paste0(ds, "_logFC") := logFC,
      !!paste0(ds, "_pval")  := pval,
      !!paste0(ds, "_FDR")   := FDR
    )
  if (is.null(acc)) tmp else full_join(acc, tmp, by = "gene")
}, .init = NULL)

genes_of_interest <- union(common_up_kofn, common_down_kofn)
dataset_names <- names(seurat_objects)

final_tbl <- wide_tbl %>%
  filter(gene %in% genes_of_interest) %>%
  rowwise() %>%
  mutate(
    n_up = sum(sapply(dataset_names, function(ds) {
      fc  <- get(paste0(ds, "_logFC"))
      fdr <- get(paste0(ds, "_FDR"))
      isTRUE(!is.na(fc) && !is.na(fdr) && fdr < fdr_cutoff && fc >  lfc_cutoff)
    })),
    n_down = sum(sapply(dataset_names, function(ds) {
      fc  <- get(paste0(ds, "_logFC"))
      fdr <- get(paste0(ds, "_FDR"))
      isTRUE(!is.na(fc) && !is.na(fdr) && fdr < fdr_cutoff && fc < -lfc_cutoff)
    })),
    direction_kofn = case_when(
      n_up   >= k_of_n ~ paste0("Up_in_>=", k_of_n, "of", length(dataset_names)),
      n_down >= k_of_n ~ paste0("Down_in_>=", k_of_n, "of", length(dataset_names)),
      TRUE             ~ "Other"
    )
  ) %>%
  ungroup() %>%
  arrange(desc(n_up), desc(n_down), gene)

# VENN DIAGRAMS (UP + DOWN)
venn_up <- venn.diagram(
  x = up_lists,
  category.names = names(up_lists),
  filename = NULL
)

venn_down <- venn.diagram(
  x = down_lists,
  category.names = names(down_lists),
  filename = NULL
)

# Draw (in RStudio plot pane)
grid.newpage(); grid.draw(venn_up)
grid.newpage(); grid.draw(venn_down)

# OUTPUT OBJECTS
# Gene vectors:
# common_up_all, common_down_all
# common_up_kofn, common_down_kofn
#
# Main per-gene per-dataset stats table for genes in >=k_of_n:
final_tbl

#write to disk
write.csv(final_tbl, "Path_to_file/FAP_DE_genes_3of4_per_dataset_stats.csv", row.names = FALSE)

# Z score plot

# Universal gene expression accessor (handles both old and new Seurat versions)
get_expression <- function(seurat_object, gene, assay = "RNA") {
  assay_obj <- seurat_object[[assay]]
  
  # Old Seurat: @data exists
  if ("data" %in% slotNames(assay_obj)) {
    if (gene %in% rownames(assay_obj@data)) {
      return(assay_obj@data[gene, ])
    } else {
      warning(paste("Gene", gene, "not found in RNA @data"))
      return(rep(NA, ncol(seurat_object)))
    }
  }
  
  # Seurat v5: uses layers
  layers_available <- SeuratObject::Layers(assay_obj)
  if ("data" %in% layers_available) {
    expr <- GetAssayData(seurat_object, assay = assay, layer = "data")
    if (gene %in% rownames(expr)) {
      return(expr[gene, ])
    } else {
      warning(paste("Gene", gene, "not found in RNA/data layer"))
      return(rep(NA, ncol(seurat_object)))
    }
  }
  
  stop("No usable expression layer (@data or layer='data') found.")
}

# Function to calculate average expression and z-score
calculate_avg_zscores <- function(genes, seurat_obj) {
  zscores_positive <- numeric(length(genes))
  zscores_negative <- numeric(length(genes))
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    
    # Get expression for all cells
    all_expr <- get_expression(seurat_obj, gene)
    
    # Skip gene if expression retrieval failed
    if (all(is.na(all_expr))) {
      zscores_positive[i] <- NA
      zscores_negative[i] <- NA
      next
    }
    
    overall_mean <- mean(all_expr)
    overall_sd <- sd(all_expr)
    
    positive_cells <- colnames(seurat_obj)[seurat_obj$FAP_Status == 'Positive']
    negative_cells <- colnames(seurat_obj)[seurat_obj$FAP_Status == 'Negative']
    
    avg_positive <- mean(all_expr[positive_cells])
    avg_negative <- mean(all_expr[negative_cells])
    
    zscores_positive[i] <- (avg_positive - overall_mean) / overall_sd
    zscores_negative[i] <- (avg_negative - overall_mean) / overall_sd
  }
  
  return(list(Positive = zscores_positive, Negative = zscores_negative))
}

# Seurat objects
seurat_objects <- list(
  Dataset1 = obj_Ebert_Pericyte,
  Dataset2 = obj_Abdelfattah_Pericyte,
  Dataset3 = obj_LeBlanc_Pericyte,
  Dataset4 = obj_Nomura_Pericyte
)

# Combine up + down gene list
genes_of_interest <- c(common_upregulated_3of4, common_downregulated_3of4)

# Initialise output dataframe
avg_zscores_df <- data.frame(Gene = genes_of_interest)

# Loop over datasets
for (dataset_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Create FAP_Status metadata using version-safe accessor
  fap_expr <- get_expression(seurat_obj, "FAP")
  seurat_obj$FAP_Status <- ifelse(fap_expr > 0, "Positive", "Negative")
  
  # Compute Z-scores
  zscores <- calculate_avg_zscores(genes_of_interest, seurat_obj)
  
  # Append to output dataframe
  avg_zscores_df[[paste0(dataset_name, "_Positive")]] <- zscores$Positive
  avg_zscores_df[[paste0(dataset_name, "_Negative")]] <- zscores$Negative
}

# Output the Z-scores data frame
print(avg_zscores_df)

write.csv(avg_zscores_df, "avg_zscores_FAP_pos_vs_neg_df.csv")
