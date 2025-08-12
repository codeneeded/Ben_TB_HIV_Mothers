library(Seurat)
library(dplyr)
library(SeuratExtend)
library(cowplot)
library(ggplot2)
library(EnhancedVolcano)

# Paths
base_dir       <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
saved_dir      <- file.path(base_dir, "saved_R_data")
annotation_dir <- file.path(base_dir, "Annotation")

rds_in  <- file.path(saved_dir, "seu_parsef_integrated_annotationReady.rds")
csv_in  <- file.path(annotation_dir, "IGRA_annotation.csv")

# Load Seurat object
seu <- readRDS(rds_in)

# Load annotation table
ann <- read.csv(csv_in, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Remove columns with empty names
ann <- ann[, colnames(ann) != ""]

# Now rename properly
ann <- ann %>%
  rename(cluster_id = `Cluster ID`,
         IGRA_Annotation = Annotation,
         IGRA_Notes = Notes) %>%
  mutate(cluster_id = as.character(cluster_id))


# Create mapping vectors
annot_map <- setNames(ann$IGRA_Annotation, ann$cluster_id)
notes_map <- setNames(ann$IGRA_Notes, ann$cluster_id)

# Map annotations into metadata
cluster_vec <- as.character(seu$mnn.snn.louvianmlr_1.5)
seu$IGRA_Annotation <- unname(annot_map[cluster_vec])
seu$IGRA_Notes      <- unname(notes_map[cluster_vec])

# Clean annotation text
seu$IGRA_Annotation <- trimws(seu$IGRA_Annotation)  # remove leading/trailing spaces
seu$IGRA_Annotation <- as.character(seu$IGRA_Annotation)  # force to character
seu$IGRA_Annotation <- factor(seu$IGRA_Annotation)        # back to factor

# Check merge
table(seu$IGRA_Annotation, useNA = "ifany")

# --- Settings ---
genes_for_split <- c("CXCR5", "ICOS")
new_label <- "CD4+ Tfh"

# Identify cells in cluster 1
is_cluster1 <- seu$IGRA_Annotation == "CD4+ Activated"

# Get expression (log-normalized data slot)
expr <- GetAssayData(seu, layer = "data", assay = DefaultAssay(seu))  # updated for Seurat v5

# Co-expression: both genes > 0
coexp <- Matrix::colSums(expr[genes_for_split, , drop = FALSE] > 0) == length(genes_for_split)

# Convert to character to allow new label
tmp_ann <- as.character(seu$IGRA_Annotation)
tmp_ann[is_cluster1 & coexp] <- new_label

# Write back as factor
seu$IGRA_Annotation <- factor(tmp_ann)

# Check result
table(seu$IGRA_Annotation)

### Set Idents for Plotting
Idents(seu) <- "IGRA_Annotation"

# Path to save
out_path <- file.path(annotation_dir, "Annotated_UMAP.png")

# Create plot
p <- DimPlot2(
  seu,
  reduction = "umap.mnn.rna",
  group.by = "IGRA_Annotation",
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
) + ggtitle("Annotated Plot")

# Save
ggsave(
  filename = out_path,
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)


# Proportion of cells in each IGRA_Annotation
p1 <- ClusterDistrBar(
  origin  = seu$IGRA_status,
  cluster = seu$IGRA_Annotation
)

ggsave(
  filename = file.path(annotation_dir, "ClusterDistribution_Proportion.png"),
  plot = p1, width = 10, height = 6, dpi = 300
)

# Absolute cell counts
p2 <- ClusterDistrBar(
  origin  = seu$IGRA_status,
  cluster = seu$IGRA_Annotation,
  percent = FALSE
)

ggsave(
  filename = file.path(annotation_dir, "ClusterDistribution_Absolute.png"),
  plot = p2, width = 10, height = 6, dpi = 300
)

# IGRA pos vs neg, reversed, normalized
p3 <- ClusterDistrBar(
  origin     = seu$IGRA_status,
  cluster    = seu$IGRA_Annotation,
  rev        = TRUE,
  normalize  = TRUE
)

ggsave(
  filename = file.path(annotation_dir, "ClusterDistribution_IGRA_PosVsNeg.png"),
  plot = p3, width = 10, height = 6, dpi = 300
)

################################ DGE ##############################################################################

# Paths
dge_dir <- file.path(base_dir, "DGE")
dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)

# Make sure Idents are set to IGRA_Annotation
Idents(seu) <- "IGRA_Annotation"

# Loop through clusters
for (clust in levels(seu$IGRA_Annotation)) {
  
  message("Running DGE for cluster: ", clust)
  
  # Subset to this cluster
  seu_sub <- subset(seu, idents = clust)
  
  # Ensure IGRA_status is a factor
  seu_sub$IGRA_status <- factor(seu_sub$IGRA_status)
  
  # Run MAST test controlling for library size
  markers <- FindMarkers(
    seu_sub,
    ident.1 = "Positive",           # IGRA positive group
    ident.2 = "Negative",           # IGRA negative group
    group.by = "IGRA_status",
    test.use = "MAST",
    latent.vars = c("nCount_RNA"),  # control for library size
    logfc.threshold = 0,            # keep all genes
    min.pct = 0
  )
  
  # Save results
  out_file <- file.path(dge_dir, paste0("DGE_", gsub(" ", "_", clust), "_IGRApos_vs_neg.csv"))
  write.csv(markers, out_file)
}

############################################ Volcano Plots ####################################################################


# Paths
dge_dir <- file.path(base_dir, "DGE")
plot_dir <- file.path(dge_dir, "Volcano_Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# List all cluster DGE CSV files
dge_files <- list.files(dge_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in dge_files) {
  # Read results
  deg <- read.csv(f, stringsAsFactors = FALSE)
  
  # Ensure gene column exists
  # If the Seurat FindMarkers CSV has rownames as gene names, fix:
  if (!"hgnc_symbol" %in% colnames(deg)) {
    deg$hgnc_symbol <- rownames(deg)
  }
  
  # Remove genes without symbols
  deg <- deg %>% filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
  
  # Set coloring rules
  keyvals <- ifelse(
    abs(deg$avg_log2FC) > 1.5 & deg$p_val_adj < 0.01, "#CD0BBC",
    ifelse(deg$p_val_adj < 0.01, "#28E2E5", "gray30")
  )
  keyvals[is.na(keyvals)] <- "gray30"
  names(keyvals)[keyvals == "gray30"] <- "NS"
  names(keyvals)[keyvals == "#28E2E5"] <- "adj(p-value) < 0.01"
  names(keyvals)[keyvals == "#CD0BBC"] <- "FC > 1.5"
  
  # Extract cluster name from filename
  clust_name <- gsub("^DGE_|_IGRApos_vs_neg\\.csv$", "", basename(f))
  
  # Make volcano plot
  vp <- EnhancedVolcano(
    deg,
    lab = deg$hgnc_symbol,
    x = "avg_log2FC",         # Seurat FindMarkers column
    y = "p_val",              # raw p-value column
    pCutoffCol = "p_val_adj", # adjusted p-value column
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 0.01,
    FCcutoff = 1.5,
    pointSize = 4.0,
    labSize = 3.0,
    labCol = "black",
    labFace = "bold",
    colAlpha = 4/5,
    legendPosition = "right",
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = "black",
    title = paste0(clust_name, ": IGRA+ vs IGRA-"),
    subtitle = "Differential Gene Expression",
    colCustom = keyvals
  )
  
  # Save
  ggsave(
    filename = file.path(plot_dir, paste0("Volcano_", clust_name, ".png")),
    plot = vp + guides(color = guide_legend(reverse = TRUE)),
    dpi = 500,
    width = 10,
    height = 7
  )
}
