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
  'CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
  'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
  'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
  'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
  'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
  'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
  'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
  'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67',
  'CXCR5','PDCD1','BCL6','IL21','BATF','MAF','ASCL2','SH2D1A','CD200','TOX','CCR2','PRDM1','IL10'
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
