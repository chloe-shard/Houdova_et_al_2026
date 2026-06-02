# -------------------------------------------------------------------------------------------------------------
#                                      RCTD annotation of cell types in spatial T data
# -------------------------------------------------------------------------------------------------------------

#Load libraries
library(Seurat)
library(spacexr)
library(progress)
library(dplyr)

#### OPTION A: Abdelfattah scRNAseq object annotated with "broad" cell type annotations (e.g. Malignant, pericyte, Endothelial etc) ####

# load Abdelfattah dataset
obj_Abdelfattah <- readRDS("G:/TAOF/Raw_seq_data/Processed_rds_files/Single_cell/GBM_patients/Abdelfattah.rds")
table(obj_Abdelfattah$Patient)
Idents(obj_Abdelfattah) <- "Patient"
DimPlot(obj_Abdelfattah)

#remove recurrent samples
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("ndGBM-01", "ndGBM-02", "ndGBM-04", "ndGBM-05", "ndGBM-06", "ndGBM-07", "ndGBM-08", "ndGBM-09", "ndGBM-10", "ndGBM-11"))
Idents(obj_Abdelfattah) <- "Assignment"
table(obj_Abdelfattah$Assignment)

# List of Seurat objects and their corresponding cell type column names
seurat_list <- list(obj_Abdelfattah)
cell_type_columns <- c("Assignment")

# Define cell type mappings for each Seurat object
cell_type_mappings <- list(
  # Mapping for obj_Abdelfattah
  c(
    "BCells" = "Lymphocyte",
    "Endo" = "Endothelial",
    "Glioma" = "Malignant",
    "Oligo" = "Oligodendrocyte",
    "Pericytes" = "Pericyte",
    "TCells" = "Lymphocyte",
    "Other" = "Other"
  )
)

# Loop through Seurat object to add the new annotations
for (i in 1:length(seurat_list)) {
  seurat_obj <- seurat_list[[i]]
  cell_type_column <- cell_type_columns[i]
  cell_type_mapping <- cell_type_mappings[[i]]
  seurat_obj@meta.data$Pooled_Cell_Type <- recode(
    seurat_obj@meta.data[[cell_type_column]], !!!cell_type_mapping
  )
  if (i == 1) {
    obj_Abdelfattah <- seurat_obj
  } 
}

# updated metadata
table(obj_Abdelfattah$Pooled_Cell_Type, useNA = "ifany")

# subset to remove "Other" annotated cells
obj_Abdelfattah <- subset(
  obj_Abdelfattah,
  subset = !Pooled_Cell_Type %in% c("Other")
)

# remove unused cell type annotations
obj_Abdelfattah$Pooled_Cell_Type <- droplevels(
  obj_Abdelfattah$Pooled_Cell_Type
)

# check metadata
table(obj_Abdelfattah$Pooled_Cell_Type)

#### OPTION B: Abdelfattah scRNAseq object annotated with all cell subtype annotations (e.g. FAP_pos_pericyte, FAP_neg_pericyte, AC, MES, OPC, NPC, M1, M2 etc) ####

# To annotate all cell subtype annotateions see R File: Houdova_et_al_R_script_11_scRNAseq_dataset_cell_subtype_annotations

# subset to remove Other, and "None" cells
obj_Abdelfattah <- subset(
  obj_Abdelfattah,
  subset = !Pooled_Cell_Type_FAP_updatedV3 %in% c("Other", "None_Malignant", "None_Myeloid", "None")
)

# remove unused cell type annotations
obj_Abdelfattah$Pooled_Cell_Type_FAP_updatedV3 <- droplevels(
  obj_Abdelfattah$Pooled_Cell_Type_FAP_updatedV3
)

# check metadata
table(obj_Abdelfattah$Pooled_Cell_Type_FAP_updatedV3)

#### OPTION C: Abdelfattah scRNAseq object annotated with all pericyte, myeloid and malignant subtype annotations (e.g. Transport_PCs, ECM_PCs, Interferon_PCs, SMCs, FAP_PCs, Proliferative_PCs, M1, M2, Monocyte, Microglia, NPC, OPC, MES, AC) ####

# To annotate Pericyte subpopulations see R File: Houdova_et_al_R_script_5_scRNAseq_Pericyte_subpopulation_annotations_UMAP
# To annotate Myeloid subpopulations see R File: Houdova_et_al_R_script_10_scRNAseq_Myeloid_receptor_expression_Dotplot
# To annotate Malignant subpopulations see R File: Houdova_et_al_R_script_12_scRNAseq_Malignant_subpopulation_annotations

# Merge Myeloid annotations with "broad" cell annotations

# list datasets
datasets <- c("Abdelfattah")

# function
for (prefix in datasets) {
  message("Processing: ", prefix)
  full_obj_name    <- paste0("obj_", prefix)
  myeloid_obj_name <- paste0("obj_", prefix, "_Myeloid")
  
  # Check objects exist
  if (!exists(full_obj_name) || !exists(myeloid_obj_name)) {
    warning("Missing object(s) for: ", prefix)
    next
  }
  
  full_obj    <- get(full_obj_name)
  myeloid_obj <- get(myeloid_obj_name)
  
  # Check required metadata columns exist
  required_full_cols <- "Pooled_Cell_Type"
  required_myeloid_cols <- "Myeloid_subpop"
  
  if (!required_full_cols %in% colnames(full_obj@meta.data)) {
    warning(required_full_cols, " not found in: ", full_obj_name)
    next
  }
  
  if (!required_myeloid_cols %in% colnames(myeloid_obj@meta.data)) {
    warning(required_myeloid_cols, " not found in: ", myeloid_obj_name)
    next
  }
  
  # Initialise all cells as "None"
  full_obj$Myeloid_subpop <- "None"
  
  # Find overlapping cells
  overlap_cells <- intersect(
    colnames(full_obj),
    colnames(myeloid_obj)
  )
  
  # Add myeloid subtype annotations to matching cells
  full_obj$Myeloid_subpop[overlap_cells] <- as.character(
    myeloid_obj$Myeloid_subpop[overlap_cells]
  )
  
  # Start from existing pooled annotation
  pooled <- as.character(full_obj$Pooled_Cell_Type_FAP)
  
  # Replace only cells labelled "Myeloid"
  myeloid_idx <- pooled == "Myeloid"
  pooled[myeloid_idx] <- as.character(full_obj$Myeloid_subpop[myeloid_idx])
  
  # If a cell was labelled Myeloid but has no subtype annotation
  pooled[myeloid_idx & full_obj$Myeloid_subpop == "None"] <- "None_Myeloid"
  
  # Save updated annotation
  full_obj$Pooled_Cell_Type_Myeloid <- factor(pooled)
  
  # Put object back into environment
  assign(full_obj_name, full_obj)
  
  # Checks
  message("Finished: ", prefix)
  message("Matched Myeloid cells: ", length(overlap_cells))
  print(table(full_obj$Myeloid_subpop, useNA = "ifany"))
  print(table(full_obj$Pooled_Cell_Type_Myeloid, useNA = "ifany"))
}

# check metadata
table(obj_Ebert$Pooled_Cell_Type_Myeloid)

# Merge Pericyte annotations with "broad" cell annotations

# list datasets
datasets <- c("Abdelfattah")

# function
for (prefix in datasets) {
  message("Processing: ", prefix)
  full_obj_name    <- paste0("obj_", prefix)
  Pericyte_obj_name <- paste0("obj_", prefix, "_Pericyte")
  
  # Check objects exist
  if (!exists(full_obj_name) || !exists(Pericyte_obj_name)) {
    warning("Missing object(s) for: ", prefix)
    next
  }
  
  full_obj    <- get(full_obj_name)
  Pericyte_obj <- get(Pericyte_obj_name)
  
  # Check required metadata columns exist
  required_full_cols <- "Pooled_Cell_Type_Myeloid"
  required_Pericyte_cols <- "Pericyte_subpop_smooth"
  
  if (!required_full_cols %in% colnames(full_obj@meta.data)) {
    warning(required_full_cols, " not found in: ", full_obj_name)
    next
  }
  
  if (!required_Pericyte_cols %in% colnames(Pericyte_obj@meta.data)) {
    warning(required_Pericyte_cols, " not found in: ", Pericyte_obj_name)
    next
  }
  
  # Initialise all cells as "None"
  full_obj$Pericyte_subpop <- "None"
  
  # Find overlapping cells
  overlap_cells <- intersect(
    colnames(full_obj),
    colnames(Pericyte_obj)
  )
  
  # Add Pericyte subtype annotations to matching cells
  full_obj$Pericyte_subpop[overlap_cells] <- as.character(
    Pericyte_obj$Pericyte_subpop[overlap_cells]
  )
  
  # Start from existing pooled annotation
  pooled <- as.character(full_obj$Pooled_Cell_Type_FAP)
  
  # Replace only cells labelled "Pericyte"
  Pericyte_idx <- pooled == "Pericyte"
  pooled[Pericyte_idx] <- as.character(full_obj$Pericyte_subpop[Pericyte_idx])
  
  # If a cell was labelled Pericyte but has no subtype annotation
  pooled[Pericyte_idx & full_obj$Pericyte_subpop == "None"] <- "None_Pericyte"
  
  # Save updated annotation
  full_obj$Pooled_Cell_Type_Pericyte_Pericyte <- factor(pooled)
  
  # Put object back into environment
  assign(full_obj_name, full_obj)
  
  # Checks
  message("Finished: ", prefix)
  message("Matched Pericyte cells: ", length(overlap_cells))
  print(table(full_obj$Pericyte_subpop, useNA = "ifany"))
  print(table(full_obj$Pooled_Cell_Type_Pericyte_Pericyte, useNA = "ifany"))
}

# check metadata
table(obj_Ebert$Pooled_Cell_Type_Myeloid_Pericyte)

# Merge Malignant annotations with "broad" cell annotations

# list datasets
datasets <- c("Abdelfattah")

# Rename map: merge MES and NPC subpopulations
rename_map <- c(
  "MES1" = "MES",
  "MES2" = "MES",
  "NPC1" = "NPC",
  "NPC2" = "NPC"
)

# Rename Malignant_subpop values
obj_Abdelfattah$Malignant_subpop_merged <- obj_Abdelfattah$Malignant_subpop

obj_Abdelfattah$Malignant_subpop_merged[
  obj_Abdelfattah$Malignant_subpop_merged %in% names(rename_map)
] <- rename_map[
  obj_Abdelfattah$Malignant_subpop_merged[
    obj_Abdelfattah$Malignant_subpop_merged %in% names(rename_map)
  ]
]

# Check output
table(obj_Abdelfattah$Malignant_subpop_merged, useNA = "ifany")

# function
for (prefix in datasets) {
  message("Processing: ", prefix)
  full_obj_name    <- paste0("obj_", prefix)
  Malignant_obj_name <- paste0("obj_", prefix, "_Malignant")
  
  # Check objects exist
  if (!exists(full_obj_name) || !exists(Malignant_obj_name)) {
    warning("Missing object(s) for: ", prefix)
    next
  }
  
  full_obj    <- get(full_obj_name)
  Malignant_obj <- get(Malignant_obj_name)
  
  # Check required metadata columns exist
  required_full_cols <- "Pooled_Cell_Type_Myeloid_Malignant"
  required_Malignant_cols <- "Malignant_subpop_merged"
  
  if (!required_full_cols %in% colnames(full_obj@meta.data)) {
    warning(required_full_cols, " not found in: ", full_obj_name)
    next
  }
  
  if (!required_Malignant_cols %in% colnames(Malignant_obj@meta.data)) {
    warning(required_Malignant_cols, " not found in: ", Malignant_obj_name)
    next
  }
  
  # Initialise all cells as "None"
  full_obj$Malignant_subpop <- "None"
  
  # Find overlapping cells
  overlap_cells <- intersect(
    colnames(full_obj),
    colnames(Malignant_obj)
  )
  
  # Add Malignant subtype annotations to matching cells
  full_obj$Malignant_subpop[overlap_cells] <- as.character(
    Malignant_obj$Malignant_subpop[overlap_cells]
  )
  
  # Start from existing pooled annotation
  pooled <- as.character(full_obj$Pooled_Cell_Type_FAP)
  
  # Replace only cells labelled "Malignant"
  Malignant_idx <- pooled == "Malignant"
  pooled[Malignant_idx] <- as.character(full_obj$Malignant_subpop[Malignant_idx])
  
  # If a cell was labelled Malignant but has no subtype annotation
  pooled[Malignant_idx & full_obj$Malignant_subpop == "None"] <- "None_Malignant"
  
  # Save updated annotation
  full_obj$Pooled_Cell_Type_Myeloid_Pericyte_Malignant <- factor(pooled)
  
  # Put object back into environment
  assign(full_obj_name, full_obj)
  
  # Checks
  message("Finished: ", prefix)
  message("Matched Malignant cells: ", length(overlap_cells))
  print(table(full_obj$Malignant_subpop, useNA = "ifany"))
  print(table(full_obj$Pooled_Cell_Type_Myeloid_Pericyte_Malignant, useNA = "ifany"))
}

# check metadata
table(obj_Ebert$Pooled_Cell_Type_Myeloid_Pericyte_Malignant)

# subset to remove Other, and "None" cells
obj_Abdelfattah <- subset(
  obj_Abdelfattah,
  subset = !Pooled_Cell_Type_Myeloid_Pericyte_Malignant %in% c("Other", "None_Pericyte", "None_Myeloid", "None_Malignant", "None")
)

# remove unused cell type annotations
obj_Abdelfattah$Pooled_Cell_Type_Myeloid_Pericyte_Malignant <- droplevels(
  obj_Abdelfattah$Pooled_Cell_Type_Myeloid_Pericyte_Malignant
)

# check metadata
table(obj_Abdelfattah$Pooled_Cell_Type_Myeloid_Pericyte_Malignant)


#### Downsample so each group has a in a chosen metadata column has at most 800 cells ####

set.seed(123)

# Specify metadata column and max cells per group
metadata_col <- "Pooled_Cell_Type_FAP_updatedV3" # this is for OPTION B - change for other options
max_cells <- 800

# safety check
if (!metadata_col %in% colnames(obj_Abdelfattah_subset@meta.data)) {
  stop(paste(metadata_col, "is not a metadata column in the Seurat object."))
}

# Count cells per metadata group
celltype_counts <- table(obj_Abdelfattah_subset@meta.data[[metadata_col]])

# Downsample within each metadata group
cells_to_keep <- unlist(lapply(names(celltype_counts), function(ct) {
  
  ct_cells <- colnames(obj_Abdelfattah_subset)[
    obj_Abdelfattah_subset@meta.data[[metadata_col]] == ct
  ]
  
  if (length(ct_cells) > max_cells) {
    sample(ct_cells, max_cells)
  } else {
    ct_cells
  }
}))

# Subset object
obj_Abdelfattah_downsampled <- subset(
  obj_Abdelfattah_subset,
  cells = cells_to_keep
)

# Check counts
table(obj_Abdelfattah_downsampled@meta.data[[metadata_col]])


#### Run RCTD analysis ####

# function to run RCTD annotation
run_RCTD <- function(
    sc_obj,
    spatial_obj,
    cell_type_col = "Pooled_Cell_Type_FAP_updatedV3", # this is for OPTION B - change for other options
    suffix = "_Houdova_OptionB" # this is for OPTION B - change for other options
) {
  if (!cell_type_col %in% colnames(sc_obj@meta.data)) {
    stop(paste0(
      "The metadata column '", cell_type_col,
      "' was not found in sc_obj@meta.data."
    ))
  }
  reference_counts <- GetAssayData(
    sc_obj,
    assay = "RNA",
    layer = "counts"
  )
  cell_type_metadata <- sc_obj@meta.data[[cell_type_col]]
  valid_cells <- !is.na(cell_type_metadata)
  reference_counts <- reference_counts[, valid_cells, drop = FALSE]
  cell_types <- as.factor(cell_type_metadata[valid_cells])
  names(cell_types) <- colnames(reference_counts)
  spatial_counts <- GetAssayData(
    spatial_obj,
    assay = "Spatial",
    layer = "counts"
  )
  spatial_coords <- GetTissueCoordinates(spatial_obj)
  shared_genes <- intersect(
    rownames(reference_counts),
    rownames(spatial_counts)
  )
  reference_counts <- reference_counts[shared_genes, , drop = FALSE]
  spatial_counts <- spatial_counts[shared_genes, , drop = FALSE]
  spatial_nUMI <- Matrix::colSums(spatial_counts)
  shared_barcodes <- Reduce(intersect, list(
    colnames(spatial_counts),
    names(spatial_nUMI),
    rownames(spatial_coords)
  ))
  spatial_counts <- spatial_counts[, shared_barcodes, drop = FALSE]
  spatial_nUMI <- spatial_nUMI[shared_barcodes]
  spatial_coords <- spatial_coords[shared_barcodes, , drop = FALSE]
  reference <- Reference(reference_counts, cell_types)
  puck <- SpatialRNA(
    coords = spatial_coords,
    counts = spatial_counts,
    nUMI = spatial_nUMI
  )
  rctd <- create.RCTD(puck, reference)
  myRCTD <- run.RCTD(rctd, doublet_mode = "full")
  cell_type_proportions <- myRCTD@results$weights
  common_barcodes <- intersect(
    rownames(cell_type_proportions),
    colnames(spatial_obj)
  )
  cell_type_proportions <- cell_type_proportions[
    common_barcodes,
    ,
    drop = FALSE
  ]
  for (cell_type in colnames(cell_type_proportions)) {
    new_col <- rep(NA_real_, ncol(spatial_obj))
    names(new_col) <- colnames(spatial_obj)
    
    new_col[common_barcodes] <- cell_type_proportions[, cell_type]
    
    new_name <- paste0(cell_type, suffix)
    
    spatial_obj[[new_name]] <- new_col
  }
  return(spatial_obj)
}

# spatial object names
spatial_object_names <- c(
  "T243", 
  "T248",
  "T255",
  "T259", 
  "T260",
  "T262",
  "GBM_ZH881T1",
  "GBM_ZH916T1",
  "GBM_ZH1007nec", 
  "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", 
  "GBM_ZH1007inf", 
  "GBM_ZH1019inf"
)

# Create a progress bar
pb <- progress_bar$new(
  format = "  Running [:bar] :percent in :elapsed ETA: :eta",
  total = length(spatial_object_names),
  width = 60
)

# Run RCTD and overwrite objects in the global environment
for (name in spatial_object_names) {
  message("Running RCTD for: ", name)
  spatial_obj <- get(name, envir = .GlobalEnv)
  updated_obj <- run_RCTD(
    sc_obj = obj_Abdelfattah_downsampled,
    spatial_obj = spatial_obj,
    suffix = "_Houdova"
  )
  assign(name, updated_obj, envir = .GlobalEnv)
  pb$tick()
}

#### Save objects after checking ####

spatial_object_names <- c(
  "T243", 
  "T248",
  "T255",
  "T259", 
  "T260",
  "T262",
  "GBM_ZH881T1",
  "GBM_ZH916T1",
  "GBM_ZH1007nec", 
  "GBM_ZH8811Abulk",
  "GBM_ZH8812bulk", 
  "GBM_ZH1007inf", 
  "GBM_ZH1019inf"
)

output_base_path <- "Path to file"

#save annotated objects
for (obj_name in spatial_object_names) {
  if (!exists(obj_name)) {
    warning("Object not found: ", obj_name)
    next
  }
  message("Saving: ", obj_name)
  saveRDS(
    get(obj_name),
    file = file.path(output_base_path, paste0(obj_name, ".rds"))
  )
}