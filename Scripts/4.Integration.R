# ============================
# Integration for Parse Dataset
# ============================

# Load required libraries
library(Seurat)
library(harmony)         # Optional for batch correction
library(patchwork)
library(scCustomize)
library(SeuratWrappers)
library(SeuratExtend)

# Define base directory
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"

# ----------------------------
# Load merged parse dataset
# ----------------------------
seu_parsef_clean <- readRDS(
  file.path(base_dir, "saved_R_data", "seu_parsef_filtered_CellCycle_DoubletClean.rds")
)
# Create Integration folder
integration_dir <- file.path(base_dir, "Integration_Plots")
dir.create(integration_dir, showWarnings = FALSE)

# ----------------------------
# Split object by sample for integration
# ----------------------------
seu_parsef_clean[["RNA"]] <- JoinLayers(seu_parsef_clean[["RNA"]])
seu_parsef_clean[["RNA"]] <- split(seu_parsef_clean[["RNA"]], f = seu_parsef_clean$sample)

# ----------------------------
# Standard preprocessing per sample
# ----------------------------
seu_parsef_clean <- NormalizeData(seu_parsef_clean)
seu_parsef_clean <- FindVariableFeatures(seu_parsef_clean)
seu_parsef_clean <- ScaleData(seu_parsef_clean)
seu_parsef_clean <- RunPCA(seu_parsef_clean)
seu_parsef_clean <- FindNeighbors(seu_parsef_clean, dims = 1:30, reduction = "pca")
seu_parsef_clean <- FindClusters(seu_parsef_clean, resolution = 2, cluster.name = "unintegrated_clusters")
seu_parsef_clean <- RunUMAP(seu_parsef_clean, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


# ----------------------------
# Integrate data
# ----------------------------

### Fast CCA Integration

seu_parsef_clean <- IntegrateLayers(
  seu_parsef_clean,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca.rna"
)

seu_parsef_clean <- IntegrateLayers(
  object = seu_parsef_clean, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony"
)

### Fast MNN Integration

seu_parsef_clean <- IntegrateLayers(
  seu_parsef_clean,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn.rna"
)


# ----------------------------
# Clustering and Neighbours
# ----------------------------

# CCA

seu_parsef_clean <- FindNeighbors(seu_parsef_clean, reduction = "integrated.cca.rna", dims = 1:30)
seu_parsef_clean <- RunUMAP(seu_parsef_clean, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")



seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 0.5, cluster.name = 'cca.snn.louvianmlr_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1, cluster.name = 'cca.snn.louvianmlr_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1.5, cluster.name = 'cca.snn.louvianmlr_1.5')

seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 0.5, cluster.name = 'cca.snn.slm_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1, cluster.name = 'cca.snn.slm_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1.5, cluster.name = 'cca.snn.slm_1.5')


DimPlot2(seu_parsef_clean, label = TRUE, box = TRUE, label.color = "black", repel = TRUE,  reduction = "umap.cca.rna",
         group.by = "cca.snn.louvianmlr_0.5",cols = 'default',pt.size=1) 


# ---------------------------- #
#   CCA Integration Plots      #
# ---------------------------- #

# Row 1: Azimuth Predictions
p1 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "predicted.celltype.l1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "predicted.celltype.l3",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L3")

# Row 2: LouvainMLR
p4 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.louvianmlr_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.louvianmlr_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.louvianmlr_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.louvianmlr_1.5")

# Row 3: SLM
p7 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.slm_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.slm_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.cca.rna",
  group.by = "cca.snn.slm_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("cca.snn.slm_1.5")

# Combine into a 3x3 grid
combined_plot_cca <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)


# Save to Integration folder
ggsave(
  filename = file.path(integration_dir, "CCA_Clustering_Grid.png"),
  plot = combined_plot_cca,
  width = 24,   # Adjust width as needed
  height = 17,  # Adjust height as needed
  dpi = 300
)

# ============================ #
#       FastMNN Clustering    #
# ============================ #

# Neighbors and UMAP
seu_parsef_clean <- FindNeighbors(seu_parsef_clean, reduction = "integrated.mnn.rna", dims = 1:30)
seu_parsef_clean <- RunUMAP(seu_parsef_clean, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna")

# LouvainMLR Clustering
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 0.5, cluster.name = 'mnn.snn.louvianmlr_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1, cluster.name = 'mnn.snn.louvianmlr_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1.5, cluster.name = 'mnn.snn.louvianmlr_1.5')

# SLM Clustering
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 0.5, cluster.name = 'mnn.snn.slm_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1, cluster.name = 'mnn.snn.slm_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1.5, cluster.name = 'mnn.snn.slm_1.5')

# ---------------------------- #
#   FastMNN Integration Plots  #
# ---------------------------- #

# Row 1: Azimuth Predictions
p1 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "predicted.celltype.l1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "predicted.celltype.l3",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L3")

# Row 2: LouvainMLR
p4 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.louvianmlr_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.louvianmlr_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.louvianmlr_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.louvianmlr_1.5")

# Row 3: SLM
p7 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.slm_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.slm_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.mnn.rna",
  group.by = "mnn.snn.slm_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("mnn.snn.slm_1.5")

# Combine into 3x3 grid
combined_plot_mnn <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)

# Save to Integration folder
ggsave(
  filename = file.path(integration_dir, "MNN_Clustering_Grid.png"),
  plot = combined_plot_mnn,
  width = 24,
  height = 17,
  dpi = 300
)

# ============================ #
#       Harmony Integration    #
# ============================ #

# Neighbors and UMAP
seu_parsef_clean <- FindNeighbors(seu_parsef_clean, reduction = "harmony", dims = 1:30)
seu_parsef_clean <- RunUMAP(seu_parsef_clean, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# LouvainMLR Clustering
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 0.5, cluster.name = 'harmony.snn.louvianmlr_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1, cluster.name = 'harmony.snn.louvianmlr_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 2, resolution = 1.5, cluster.name = 'harmony.snn.louvianmlr_1.5')

# SLM Clustering
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 0.5, cluster.name = 'harmony.snn.slm_0.5')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1, cluster.name = 'harmony.snn.slm_1')
seu_parsef_clean <- FindClusters(seu_parsef_clean, algorithm = 3, resolution = 1.5, cluster.name = 'harmony.snn.slm_1.5')

# ---------------------------- #
#   Harmony Integration Plots  #
# ---------------------------- #

# Row 1: Azimuth Predictions
p1 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l3",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Azimuth L3")

# Row 2: LouvainMLR
p4 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.louvianmlr_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.louvianmlr_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.louvianmlr_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.louvianmlr_1.5")

# Row 3: SLM
p7 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.slm_0.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.slm_1",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction = "umap.harmony",
  group.by = "harmony.snn.slm_1.5",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("harmony.snn.slm_1.5")

# Combine into 3x3 grid
combined_plot_harmony <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)

# Save to Integration folder
ggsave(
  filename = file.path(integration_dir, "Harmony_Clustering_Grid.png"),
  plot = combined_plot_harmony,
  width = 24,
  height = 17,
  dpi = 300
)

seu_parsef_clean[["RNA"]] <- JoinLayers(seu_parsef_clean[["RNA"]])

# ----------------------------
# Save integrated object
# ----------------------------
saveRDS(
  seu_parsef_clean,
  file.path(base_dir, "saved_R_data", "seu_parsef_integrated.rds")
)