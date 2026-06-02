# ------------------------------------------------------------------------------------------------------------------
#                                                   FAP+ UMAPs
# ------------------------------------------------------------------------------------------------------------------

# load libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Inputs for analysis
seurat_list <- list(
  Ebert            = obj_Ebert,
  Abdelfattah      = obj_Abdelfattah,
  LeBlanc          = obj_LeBlanc,
  Nomura           = obj_Nomura
)

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

# FAP cutoff limits
fap_min <- 0
fap_max <- 2

# check palette coverage
all_annotations <- unique(unlist(lapply(seurat_list, function(obj) unique(obj@meta.data$Pooled_Cell_Type))))
missing_colors <- setdiff(all_annotations, names(custom_colors))
if (length(missing_colors) > 0) {
  warning("The following annotations are missing colours: ", paste(missing_colors, collapse = ", "))
}

get_umap_limits <- function(seurat_obj, reduction = "umap") {
  emb <- Embeddings(seurat_obj, reduction = reduction)
  list(
    x = range(emb[, 1], na.rm = TRUE),
    y = range(emb[, 2], na.rm = TRUE)
  )
}

create_umap_plot <- function(seurat_obj, colors, title, lims) {
  DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by  = "Pooled_Cell_Type",
    cols      = colors,
    pt.size   = pt.size
  ) +
    ggtitle(title) +
    coord_cartesian(xlim = lims$x, ylim = lims$y) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10)
    )
}

create_fap_plot <- function(seurat_obj, title, lims, fap_min, fap_max, raster = FALSE) {
  FeaturePlot(
    seurat_obj,
    features   = "FAP",
    min.cutoff = fap_min,
    max.cutoff = fap_max,
    raster     = raster
  ) +
    ggtitle(title) +
    coord_cartesian(xlim = lims$x, ylim = lims$y) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )
}

# Build per-dataset paired panels
rows <- lapply(names(seurat_list), function(ds) {
  
  obj  <- seurat_list[[ds]]
  lims <- get_umap_limits(obj, reduction = "umap")
  
  p_celltype <- create_umap_plot(
    seurat_obj = obj,
    colors     = custom_colors,
    title      = paste0(ds, " - Cell type"),
    lims       = lims,
  )
  
  p_fap <- create_fap_plot(
    seurat_obj = obj,
    title      = paste0(ds, " - FAP"),
    lims       = lims,
    fap_min    = fap_min,
    fap_max    = fap_max,
    raster     = FALSE
  )
  
  # One row: Cell type | FAP
  (p_celltype | p_fap) + plot_layout(guides = "collect")
})

final_plot <- wrap_plots(rows, ncol = 1) &
  theme(legend.position = "right")

# export
ggsave(
  filename = "final_plot.svg",
  plot     = final_plot,
  device   = "svg",
  width    = 10,
  height   = 10,
  units    = "in"
)
