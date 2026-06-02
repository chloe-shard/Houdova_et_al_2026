# -------------------------------------------------------------------------------------------------------------
#                               scRNAseq dataset import and initial annotations for downstream analysis
# -------------------------------------------------------------------------------------------------------------

# load libraries
library(dplyr)
library(Seurat)

##### load GBM datasets and subset for primary GBM tumours ####

#Ebert
obj_Ebert <- readRDS("file_path/Ebert.rds")
Idents(obj_Ebert) <- "combined_clusters_annot"
table(obj_Ebert$combined_clusters_annot)

#Abdelfattah
obj_Abdelfattah <- readRDS("file_path/Abdelfattah.rds")
table(obj_Abdelfattah$Patient)
Idents(obj_Abdelfattah) <- "Patient"
DimPlot(obj_Abdelfattah)
#remove recurrent samples
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("ndGBM-01", "ndGBM-02", "ndGBM-04", "ndGBM-05", "ndGBM-06", "ndGBM-07", "ndGBM-08", "ndGBM-09", "ndGBM-10", "ndGBM-11"))
Idents(obj_Abdelfattah) <- "Assignment"
table(obj_Abdelfattah$Assignment)

#LeBlanc
obj_LeBlanc <- readRDS("file_path/LeBlanc.rds")
DefaultAssay(obj_LeBlanc) <- 'RNA'
#remove recurrent samples
table(obj_LeBlanc$sample)
Idents(obj_LeBlanc) <- "sample"
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1"))
Idents(obj_LeBlanc) <- "cell_type"
table(obj_LeBlanc$cell_type)

#Normura
obj_Nomura <- readRDS("file_path/Nomura.rds")
#remove recurrent samples
Idents(obj_Nomura) <- "PvsR"
obj_Nomura <- subset(obj_Nomura, idents = grep("Primary$", levels(Idents(obj_Nomura)), value = TRUE))
table(obj_Nomura$PvsR)
table(obj_Nomura$CellType)

########## Align author cell annotations by adding new metadata column "Pooled_Cell_Type" to harmonize nomenclature ####

# List of Seurat objects and their corresponding cell type column names
seurat_list <- list(obj_Ebert, obj_Abdelfattah, obj_LeBlanc, obj_Nomura)
cell_type_columns <- c("combined_clusters_annot", "Assignment", "cell_type", "CellType")

# Define cell type mappings for each Seurat object
cell_type_mappings <- list(
  # Mapping for obj_Ebert
  c(
    "Endothelial cells_2" = "Endothelial",
    "Macrophages_5" = "Myeloid",
    "Pericyte cells_4" = "Pericyte",
    "T cells_1" = "Lymphocyte",
    "Tumor cells_0" = "Malignant",
    "Tumor/prolif. cells_3" = "Myeloid"
  ),
  # Mapping for obj_Abdelfattah
  c(
    "BCells" = "Lymphocyte",
    "Endo" = "Endothelial",
    "Glioma" = "Malignant",
    "Oligo" = "Oligodendrocyte",
    "Pericytes" = "Pericyte",
    "TCells" = "Lymphocyte",
    "Other" = "Other"
  ),
  # Mapping for obj_LeBlanc
  c(
    "endothelial" = "Endothelial",
    "fibroblast" = "Pericyte",
    "immune" = "Myeloid",
    "malignant" = "Malignant",
    "oligodendrocyte" = "Oligodendrocyte",
    "neuron" = "Neuron"
  ),
  # Mapping for obj_Nomura
  c(
    "Endothel" = "Endothelial",
    "Excitatory neuron" = "Neuron",
    "Inhibitory neuron" = "Neuron",
    "Lymphocyte" = "Lymphocyte",
    "Oligodendrocyte" = "Oligodendrocyte",
    "OPC" = "OPC_non_malignant",
    "Pericyte" = "Pericyte",
    "Astrocyte" = "Astrocyte",
    "TAM" = "Myeloid",
    "Malignant" = "Malignant"
  )
)

# Loop through each Seurat object to add the new annotations
for (i in 1:length(seurat_list)) {
  # Get the current Seurat object, cell type column name, and mapping
  seurat_obj <- seurat_list[[i]]
  cell_type_column <- cell_type_columns[i]
  cell_type_mapping <- cell_type_mappings[[i]]
  
  # Recode the cell type annotations to the new pooled annotations
  seurat_obj@meta.data$Pooled_Cell_Type <- recode(
    seurat_obj@meta.data[[cell_type_column]], !!!cell_type_mapping
  )
  
  # Update the original Seurat object explicitly
  if (i == 1) {
    obj_Ebert <- seurat_obj
  } else if (i == 2) {
    obj_Abdelfattah <- seurat_obj
  } else if (i == 3) {
    obj_LeBlanc <- seurat_obj
  } else if (i == 4) {
    obj_Nomura <- seurat_obj
  }
}

# Example: View the updated metadata for obj_Ebert
table(obj_Ebert$Pooled_Cell_Type, useNA = "ifany")

#### Align Patient column in metadata across datasets ####

# Define dataset names and corresponding patient ID column names
dataset_names <- c("Ebert", "Abdelfattah", "LeBlanc", "Nomura")
patient_id_columns <- c("orig.ident", "Patient", "patient", "Patient")

# Loop through datasets and assign Patient_ID metadata
for (i in seq_along(dataset_names)) {
  seurat_obj <- get(paste0("obj_", dataset_names[i]))
  patient_col <- patient_id_columns[i]
  
  # Create the new Patient_ID column
  seurat_obj$Patient_ID <- seurat_obj@meta.data[[patient_col]]
  
  # Assign it back to the original object name
  assign(paste0("obj_", dataset_names[i]), seurat_obj)
}

#Check
table(obj_Ebert$Patient_ID)
table(obj_Abdelfattah$Patient_ID)
table(obj_LeBlanc$Patient_ID)
table(obj_Nomura$Patient_ID)
