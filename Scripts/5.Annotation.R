############################################
# Annotation Prep Script for Parsef Integrated Data
# Purpose: Explore marker expression before annotation
############################################

# ---------------------------- #
# Load Required Libraries
# ---------------------------- #
library(Seurat)
library(SeuratExtend)   # For DimPlot2, VlnPlot2
library(ggplot2)
library(dplyr)

# ---------------------------- #
# Paths
# ---------------------------- #
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
annotation_dir <- file.path(base_dir, "Annotation")
dir.create(annotation_dir, showWarnings = FALSE)

featureplot_dir <- file.path(annotation_dir, "FeaturePlots")
vlnplot_dir <- file.path(annotation_dir, "VlnPlots")

dir.create(featureplot_dir, showWarnings = FALSE)
dir.create(vlnplot_dir, showWarnings = FALSE)

# ---------------------------- #
# Load Integrated Object
# ---------------------------- #
seu_parsef_integrated <- readRDS(file.path(base_dir, "saved_R_data", "seu_parsef_integrated.rds"))


# ---------------------------- #
# Set Idents & Remove Small Clusters
# ---------------------------- #
Idents(seu_parsef_integrated) <- "mnn.snn.louvianmlr_1.5"  # Adjust clustering column if needed

# Remove clusters with <100 cells
cluster_sizes <- table(Idents(seu_parsef_integrated))
large_clusters <- names(cluster_sizes[cluster_sizes >= 100])
seu_parsef_integrated <- subset(seu_parsef_integrated, idents = large_clusters)

# ---------------------------- #
# Marker List
# ---------------------------- #
rna.features <- c(
  'ASCL2','BATF','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD200','CD27',
  'CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B','CTLA4','CTSW',
  'CXCL8','CXCR3','CXCR5','EBI3','FABP5','FCGR2B','FCGR3A','FOXP3','GNLY','GZMK','HAVCR2',
  'HIST1H4C','HLA-DPA1','HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHM',
  'IKZF2','IL10','IL18BP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','KLRD1','KLRC1',
  'LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67','MT-ND1','MT-ND5','MS4A1',
  'NELL2','NCAM1','NKG7','NR4A1','PDCD1','PRDM1','PRF1','SELL','SERPINA1','SERPING1','SH2D1A',
  'TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18','TNFRSF4','TNFRSF9','TOX','TRBC1','TRDV1','TRDV2',
  'TRGC2','TRGV9','XCL1','XCL2'
)

# ---------------------------- #
# Generate Plots for Each Marker
# ---------------------------- #



# Feature Plots
for (gene in rna.features) {
  if (gene %in% rownames(seu_parsef_integrated)) {
    p_feat <- DimPlot2(
      seu_parsef_integrated,
      features = gene,
      reduction = "umap.mnn.rna",   # adjust if different reduction
      cols = 'A'
    )
    
    ggsave(
      filename = file.path(featureplot_dir, paste0("Feature_", gene, ".png")),
      plot = p_feat,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

# Violin Plots
for (gene in rna.features) {
  if (gene %in% rownames(seu_parsef_integrated)) {
    p_vln <- VlnPlot2(
      seu_parsef_integrated,
      features = gene,
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) 
    
    ggsave(
      filename = file.path(vlnplot_dir, paste0("Vln_", gene, ".png")),
      plot = p_vln,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

# ---------------------------- #
# Save Post-Subset Object
# ---------------------------- #
saveRDS(seu_parsef_integrated, file.path(annotation_dir, "seu_parsef_integrated_annotationReady.rds"))
