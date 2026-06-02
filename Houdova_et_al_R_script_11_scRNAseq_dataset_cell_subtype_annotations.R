# -------------------------------------------------------------------------------------------------------------
#                               scRNAseq dataset annotations of cell subtype populations
# -------------------------------------------------------------------------------------------------------------

# load libraries
library(dplyr)
library(Seurat)

#### Annotate FAP pos vs FAP neg pericytes and add to Seurat objects ####

# function to annotate FAP pos vs FAP neg pericytes
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

# apply across datasets
obj_Ebert        <- annotate_FAP_pericytes(obj_Ebert)
obj_LeBlanc      <- annotate_FAP_pericytes(obj_LeBlanc)
obj_Abdelfattah  <- annotate_FAP_pericytes(obj_Abdelfattah)
obj_Nomura       <- annotate_FAP_pericytes(obj_Nomura)


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

# Check metatdata column
table(obj_Ebert$Pooled_Cell_Type_FAP, useNA = "ifany")

#### Add myeloid annotations from myeloid-subset objects back into full objects ####

# For adding Myeloid_subpop metadata column see R file : Houdova_et_al_R_script_10_scRNAseq_Myeloid_receptor_expression_Dotplot

# list datasets
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

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
  required_full_cols <- "Pooled_Cell_Type_FAP"
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
  full_obj$Pooled_Cell_Type_FAP_updated <- factor(pooled)
  
  # Put object back into environment
  assign(full_obj_name, full_obj)
  
  # Checks
  message("Finished: ", prefix)
  message("Matched Myeloid cells: ", length(overlap_cells))
  print(table(full_obj$Myeloid_subpop, useNA = "ifany"))
  print(table(full_obj$Pooled_Cell_Type_FAP_updated, useNA = "ifany"))
}

# check metadata
table(obj_Ebert$Pooled_Cell_Type_FAP_updated)

#### Add malignant annotations from malignant-subset objects back into full objects - annotated as "Malignant_subpop" ####

# For adding Malignant_subpop metadata column see R file : Houdova_et_al_R_script_12_scRNAseq_dataset_Malignant_subpopulation_annotations

# Datasets
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

# malignant Malignant_subpop column to use for each dataset
Malignant_subpop_column_map <- c(
  Ebert = "Malignant_subpop",
  LeBlanc = "Malignant_subpop",
  Abdelfattah = "Malignant_subpop",
  Nomura = "Malignant_subpop"
)

# function
for (prefix in datasets) {
  message("Processing: ", prefix)
  obj_name <- paste0("obj_", prefix)
  if (!exists(obj_name)) {
    warning("Object not found: ", obj_name)
    next
  }
  obj <- get(obj_name)
  if (!"Pooled_Cell_Type_FAP_updated" %in% colnames(obj@meta.data)) {
    warning("Pooled_Cell_Type_FAP_updated column not found in: ", obj_name)
    print(colnames(obj@meta.data))
    next
  }
  Malignant_subpop_col <- Malignant_subpop_column_map[[prefix]]
  
  if (!Malignant_subpop_col %in% colnames(obj@meta.data)) {
    warning(Malignant_subpop_col, " column not found in: ", obj_name)
    print(colnames(obj@meta.data))
    next
  }
  pooled <- as.character(obj@meta.data$Pooled_Cell_Type_FAP_updated)
  Malignant_subpop <- as.character(obj@meta.data[[Malignant_subpop_col]])
  malignant_idx <- pooled == "Malignant"
  pooled[malignant_idx] <- Malignant_subpop[malignant_idx]
  pooled[malignant_idx & Malignant_subpop == "None"] <- "None_Malignant"
  obj$Pooled_Cell_Type_FAP_updatedV2 <- factor(pooled)
  assign(obj_name, obj)
  message("Finished: ", prefix)
  message("Used malignant annotation column: ", Malignant_subpop_col)
  print(table(obj$Pooled_Cell_Type_FAP_updatedV2))
}

# rename cell states
# Dataset prefixes
datasets <- c("Ebert", "LeBlanc", "Abdelfattah", "Nomura")

# rename map (merge MES and NPC subpopulations)
rename_map <- list(
  
  Ebert = c(
    "MES1"   = "MES",
    "MES2"   = "MES",
    "NPC1"   = "NPC",
    "NPC2"   = "NPC",
  ),
  
  LeBlanc = c(
    "MES1"   = "MES",
    "MES2"   = "MES",
    "NPC1"   = "NPC",
    "NPC2"   = "NPC",
  ),
  
  Abdelfattah = c(
    "MES1"   = "MES",
    "MES2"   = "MES",
    "NPC1"   = "NPC",
    "NPC2"   = "NPC",
  ),
  
  Nomura = c(
    "MES1"   = "MES",
    "MES2"   = "MES",
    "NPC1"   = "NPC",
    "NPC2"   = "NPC",
  )
)

# Loop through datasets
for (prefix in datasets) {
  message("Processing: ", prefix)
  obj_name <- paste0("obj_", prefix)
  if (!exists(obj_name)) {
    warning("Object not found: ", obj_name)
    next
  }
  obj <- get(obj_name)
  if (!"Pooled_Cell_Type_FAP_updatedV2" %in% colnames(obj@meta.data)) {
    warning("Pooled_Cell_Type_FAP_updatedV2 column not found in: ", obj_name)
    print(colnames(obj@meta.data))
    next
  }
  pooled <- as.character(obj@meta.data$Pooled_Cell_Type_FAP_updatedV2)
  current_map <- rename_map[[prefix]]
  matched_labels <- intersect(names(current_map), unique(pooled))
  if (length(matched_labels) == 0) {
    warning("No matching labels to rename for: ", obj_name)
  } else {
    for (old_label in matched_labels) {
      pooled[pooled == old_label] <- current_map[[old_label]]
    }
  }
  obj$Pooled_Cell_Type_FAP_updatedV3 <- factor(pooled)
  assign(obj_name, obj)
  message("Finished: ", prefix)
  print(table(obj$Pooled_Cell_Type_FAP_updatedV3))
}

# Check output
table(obj_Ebert$Pooled_Cell_Type_FAP_updatedV3, useNA = "ifany")
