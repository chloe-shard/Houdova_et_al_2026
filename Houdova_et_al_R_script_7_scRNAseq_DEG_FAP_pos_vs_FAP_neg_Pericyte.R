# ------------------------------------------------------------------------------------------------------------------
#                                                   DEG FAP pos vs FAP neg Pericyte
# ------------------------------------------------------------------------------------------------------------------

# Load libraries
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(VennDiagram)
  library(grid)

#### Subset Pericyte cluster ####

# Subset Ebert dataset to Pericyte cells
Idents(obj_Ebert) <- "Pooled_Cell_Type"
table(obj_Ebert$Pooled_Cell_Type)
obj_Ebert_Pericyte <- subset(obj_Ebert, idents = c("Pericyte"))

# Subset Abdelfattah dataset to Pericyte cells
Idents(obj_Abdelfattah) <- "Pooled_Cell_Type"
table(obj_Abdelfattah$Pooled_Cell_Type)
obj_Abdelfattah_Pericyte <- subset(obj_Abdelfattah, idents = c("Pericyte"))

# Subset LeBlanc dataset to Pericyte cells
Idents(obj_LeBlanc) <- "Pooled_Cell_Type"
table(obj_LeBlanc$Pooled_Cell_Type)
obj_LeBlanc_Pericyte <- subset(obj_LeBlanc, idents = c("Pericyte"))

# Subset Nomura dataset to Pericyte cells
Idents(obj_Nomura) <- "Pooled_Cell_Type"
table(obj_Nomura$Pooled_Cell_Type)
obj_Nomura_Pericyte <- subset(obj_Nomura, idents = c("Pericyte"))


############### Differential gene expression FAP+ vs FAP- #######

# Set input datasets and analysis parameters
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
lfc_cutoff  <- 0.5
k_of_n      <- 3       # Number of datasets required for shared DEG call, e.g. 3 of 4


# Extract gene expression from old or new Seurat object structures
get_expression <- function(seurat_object, gene, assay = "RNA") {
  assay_obj <- seurat_object[[assay]]
  
  # Use @data slot for older Seurat objects
  if ("data" %in% slotNames(assay_obj)) {
    if (gene %in% rownames(assay_obj@data)) {
      return(assay_obj@data[gene, ])
    } else {
      stop(paste("Gene", gene, "not found in RNA @data"))
    }
  }
  
  # Use data layer for Seurat v5 objects
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


# Identify the log fold-change column name from FindMarkers output
get_fc_col <- function(markers_df) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  fc_col <- intersect(fc_candidates, colnames(markers_df))[1]
  
  if (is.na(fc_col) || length(fc_col) == 0) {
    stop("Could not find a logFC column in FindMarkers output.")
  }
  fc_col
}


# Run FAP-positive vs FAP-negative DEG analysis for one dataset
run_de_one_dataset <- function(seurat_object,
                               dataset_name,
                               assay = "RNA",
                               fap_gene = "FAP",
                               min_pct = 0.25,
                               fdr_cutoff = 0.05,
                               lfc_cutoff = 0.5) {
    fap_expr <- get_expression(seurat_object, fap_gene, assay = assay)
  seurat_object$FAP_status <- ifelse(fap_expr > 0, "Positive", "Negative")
  Idents(seurat_object) <- "FAP_status"
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
    filter(FDR < fdr_cutoff, logFC > lfc_cutoff) %>%
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

# Run DEG code across datasets
dge_results <- imap(seurat_objects, ~ run_de_one_dataset(
  seurat_object = .x,
  dataset_name  = .y,
  assay         = assay_use,
  fap_gene      = fap_gene,
  min_pct       = min_pct,
  fdr_cutoff    = fdr_cutoff,
  lfc_cutoff    = lfc_cutoff
))

# Extract upregulated and downregulated gene lists from each dataset
up_lists   <- lapply(dge_results, `[[`, "up")
down_lists <- lapply(dge_results, `[[`, "down")

# Identify genes shared across all datasets
common_up_all   <- Reduce(intersect, up_lists)
common_down_all <- Reduce(intersect, down_lists)

# Count how many datasets each DEG appears in
up_counts   <- table(unlist(up_lists))
down_counts <- table(unlist(down_lists))

# Identify genes significant in at least k datasets, e.g. 3 of 4
common_up_kofn   <- names(up_counts[up_counts >= k_of_n])
common_down_kofn <- names(down_counts[down_counts >= k_of_n])

# Create wide table with logFC, p-value and FDR for each dataset
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

# Keep genes that are consistently upregulated or downregulated across k datasets
genes_of_interest <- union(common_up_kofn, common_down_kofn)
dataset_names <- names(seurat_objects)

# Add summary columns showing how many datasets support each DEG direction
final_tbl <- wide_tbl %>%
  filter(gene %in% genes_of_interest) %>%
  rowwise() %>%
  mutate(
    n_up = sum(sapply(dataset_names, function(ds) {
      fc  <- get(paste0(ds, "_logFC"))
      fdr <- get(paste0(ds, "_FDR"))
      isTRUE(!is.na(fc) && !is.na(fdr) && fdr < fdr_cutoff && fc > lfc_cutoff)
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

# Generate Venn diagram for upregulated genes
venn_up <- venn.diagram(
  x = up_lists,
  category.names = names(up_lists),
  filename = NULL
)

# Generate Venn diagram for downregulated genes
venn_down <- venn.diagram(
  x = down_lists,
  category.names = names(down_lists),
  filename = NULL
)

# Draw Venn diagrams in RStudio plot pane
grid.newpage(); grid.draw(venn_up)
grid.newpage(); grid.draw(venn_down)

# View final DEG summary table
final_tbl

# Save file
write.csv(final_tbl, "Path_to_file/FAP_DE_genes_3of4_per_dataset_stats.csv", row.names = FALSE)

#### Z-score plot #####

# Extract gene expression from old or new Seurat object structures
get_expression <- function(seurat_object, gene, assay = "RNA") {
  assay_obj <- seurat_object[[assay]]
 
   # Use @data slot for older Seurat objects
  if ("data" %in% slotNames(assay_obj)) {
    if (gene %in% rownames(assay_obj@data)) {
      return(assay_obj@data[gene, ])
    } else {
      warning(paste("Gene", gene, "not found in RNA @data"))
      return(rep(NA, ncol(seurat_object)))
    }
  }
  
  # Use data layer for Seurat v5 objects
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


# Calculate average expression Z-scores for FAP-positive and FAP-negative cells
calculate_avg_zscores <- function(genes, seurat_obj) {
  zscores_positive <- numeric(length(genes))
  zscores_negative <- numeric(length(genes))
  for (i in seq_along(genes)) {
    gene <- genes[i]
    # Extract expression values for the current gene
    all_expr <- get_expression(seurat_obj, gene)
    # Skip genes not found in the dataset
    if (all(is.na(all_expr))) {
      zscores_positive[i] <- NA
      zscores_negative[i] <- NA
      next
    }
    # Calculate dataset-wide mean and standard deviation
    overall_mean <- mean(all_expr)
    overall_sd <- sd(all_expr)
    # Identify FAP-positive and FAP-negative cells
    positive_cells <- colnames(seurat_obj)[seurat_obj$FAP_Status == "Positive"]
    negative_cells <- colnames(seurat_obj)[seurat_obj$FAP_Status == "Negative"]
    # Calculate average expression in FAP-positive and FAP-negative cells
    avg_positive <- mean(all_expr[positive_cells])
    avg_negative <- mean(all_expr[negative_cells])
    # Convert average expression values to Z-scores
    zscores_positive[i] <- (avg_positive - overall_mean) / overall_sd
    zscores_negative[i] <- (avg_negative - overall_mean) / overall_sd
  }
  return(list(Positive = zscores_positive, Negative = zscores_negative))
}

# Store Pericyte datasets for Z-score analysis
seurat_objects <- list(
  Ebert = obj_Ebert_Pericyte,
  Abdelfattah = obj_Abdelfattah_Pericyte,
  LeBlanc = obj_LeBlanc_Pericyte,
  Nomura = obj_Nomura_Pericyte
)

# Combine shared upregulated and downregulated DEG lists
genes_of_interest <- c(common_up_kofn, common_down_kofn)

# Initialise output dataframe for Z-scores
avg_zscores_df <- data.frame(Gene = genes_of_interest)

# Calculate FAP-positive and FAP-negative Z-scores for each dataset
for (dataset_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[dataset_name]]
  # Create FAP status metadata
  fap_expr <- get_expression(seurat_obj, "FAP")
  seurat_obj$FAP_Status <- ifelse(fap_expr > 0, "Positive", "Negative")
  # Calculate Z-scores for selected genes
  zscores <- calculate_avg_zscores(genes_of_interest, seurat_obj)
  # Add Z-score results to output dataframe
  avg_zscores_df[[paste0(dataset_name, "_Positive")]] <- zscores$Positive
  avg_zscores_df[[paste0(dataset_name, "_Negative")]] <- zscores$Negative
}

# View Z-score dataframe
print(avg_zscores_df)

# Save Z-score dataframe to disk
write.csv(avg_zscores_df, "avg_zscores_FAP_pos_vs_neg_df.csv")