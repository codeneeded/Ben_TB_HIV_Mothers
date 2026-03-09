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

# Check mergebs download project -i 491055580 -o /media/akshay-iyer/Elements/10x_TB/BONE-19321 --concurrency high


table(seu$IGRA_Annotation, useNA = "ifany")

# --- Settings ---
genes_for_split <- c("CXCR5", "ICOS")
new_label <- "CD4+ Tfh"


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

#saveRDS(seu, file = "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/saved_R_data/seu_annotated.rds")

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

p
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
  normalize  = TRUE,
  cols=c('#A3F8A9','#FF918A')
)
p3
ggsave(
  filename = file.path(annotation_dir, "ClusterDistribution_IGRA_PosVsNeg.png"),
  plot = p3, width = 10, height = 6, dpi = 300
)


p4 <- ClusterDistrPlot(
  origin = seu$orig.ident,
  cluster = seu$IGRA_Annotation,
  condition = seu$IGRA_status,
  cols = c('#A3F8A9','#FF918A')
)
ggsave(
  filename = file.path(annotation_dir, "ClusterFrequency_Sample_IGRA_PosVsNeg.png"),
  plot = p4, width = 10, height = 6, dpi = 300
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
  names(keyvals)[keyvals == "#28E2E5"] <- "adj(p-value) < 0.05"
  names(keyvals)[keyvals == "#CD0BBC"] <- "FC > 1.5"
  
  # Extract cluster name from filename
  clust_name <- gsub("^DGE_|_IGRApos_vs_neg\\.csv$", "", basename(f))
  
  # Make volcano plot
  vp <- EnhancedVolcano(
    deg,
    lab = deg$X,
    x = "avg_log2FC",         # Seurat FindMarkers column
    y = "p_val",              # raw p-value column
    pCutoffCol = "p_val_adj", # adjusted p-value column
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 0.05,
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
# ─── Custom-label volcano plots for two specific clusters ───────────────────

custom_labels <- list(
  "CD4+_TCM" = c(
    "RELA","JAK3","RPBJ","DNMT1","EZH1","SIN3A","HIF1A","DDIT4",
    "TSC2","SQSTM1","FYCO1","BCL2","ZAP70","AKT3","DNAJB9","GPR15"
  ),
  "CD14+_Monocytes" = c(
    "SMAD3","JARID2","PCH2","PCGF3","STAB1","THBD","PPARD","VEGFA",
    "VCAN","SDC2","COL23A1","LUCAT1","STAT1","STAT2","TRIM25",
    "RNF213","CYBB","NCF2","FGR","CTSS","LYZ","HLA-B"
  )
)
for (clust_name in names(custom_labels)) {
  
  # Find matching file (pattern matches the cluster name inside filename)
  matched_file <- dge_files[grepl(clust_name, dge_files, fixed = TRUE)]
  
  if (length(matched_file) == 0) {
    message("No DGE file found for cluster: ", clust_name)
    next
  }
  if (length(matched_file) > 1) {
    message("Multiple files matched for: ", clust_name, " — using first: ", matched_file[1])
    matched_file <- matched_file[1]
  }
  
  # Read & prep (same as main loop)
  deg <- read.csv(matched_file, stringsAsFactors = FALSE)
  if (!"hgnc_symbol" %in% colnames(deg)) deg$hgnc_symbol <- rownames(deg)
  deg <- deg %>% filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
  
  # Colour logic
  keyvals <- ifelse(
    abs(deg$avg_log2FC) > 0.3 & deg$p_val_adj < 0.05, "#CD0BBC",
    ifelse(deg$p_val_adj < 0.05, "#28E2E5", "gray30")
  )
  keyvals[is.na(keyvals)] <- "gray30"
  names(keyvals)[keyvals == "gray30"]    <- "NS"
  names(keyvals)[keyvals == "#28E2E5"]   <- "adj(p-value) < 0.05"
  names(keyvals)[keyvals == "#CD0BBC"]   <- "FC > 0.3"
  
  # Only label the genes of interest; all others get ""
  genes_to_label <- custom_labels[[clust_name]]
  labels_vec <- ifelse(deg$X %in% genes_to_label, deg$X, "")
  
  vp <- EnhancedVolcano(
    deg,
    lab              = labels_vec,
    selectLab        = genes_to_label,
    x                = "avg_log2FC",
    y                = "p_val",
    pCutoffCol       = "p_val_adj",
    xlab             = bquote(~Log[2]~ 'fold change'),
    pCutoff          = 0.05,
    FCcutoff         = 0.3,
    pointSize        = 4.0,
    labSize          = 3.5,
    labCol           = "black",
    labFace          = "bold",
    colAlpha         = 4/5,
    legendPosition   = "right",
    legendLabSize    = 14,
    legendIconSize   = 4.0,
    drawConnectors   = TRUE,
    widthConnectors  = 1.0,
    colConnectors    = "black",
    typeConnectors   = "open",
    endsConnectors   = "last",
    lengthConnectors = unit(0.01, "npc"),
    boxedLabels      = TRUE,
    max.overlaps     = Inf,
    title            = paste0(clust_name, ": IGRA+ vs IGRA-"),
    subtitle         = "Differential Gene Expression — selected genes labelled",
    colCustom        = keyvals
  )
  
  # Safe filename (spaces → underscores, + → plus)
  safe_name <- gsub("[+ ]+", "_", clust_name)
  safe_name <- gsub("_+$", "", safe_name)   # trim trailing underscores
  
  ggsave(
    filename = file.path(plot_dir, paste0("Volcano_", safe_name, "_customLabels.png")),
    plot     = vp + guides(color = guide_legend(reverse = TRUE)),
    dpi      = 500,
    width    = 10,
    height   = 7
  )
  
  message("Saved custom-label volcano for: ", clust_name)
}
##################### Pseudobulking and DEG #######################################################################
## ---- Setup ----
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

## ---- 1) Pseudobulk by sample (orig.ident) ----
pseudobulk_by_sample <- function(
    seu,
    sample_col = "orig.ident",
    group_col  = "IGRA_status",
    assay      = "RNA",
    min_cells_per_sample = 20,   # tweak as needed
    keep_groups = c("Negative","Positive") # we compare these two levels
){
  stopifnot(sample_col %in% colnames(seu@meta.data))
  stopifnot(group_col  %in% colnames(seu@meta.data))
  
  X <- GetAssayData(seu, assay = assay, slot = "counts")  # genes x cells (dgCMatrix)
  md <- seu@meta.data[, c(sample_col, group_col), drop = FALSE]
  colnames(md) <- c("sample", "group")
  md$group <- droplevels(factor(md$group))
  
  ## keep only the two groups of interest
  keep_cells <- md$group %in% keep_groups & !is.na(md$group) & !is.na(md$sample)
  X <- X[, keep_cells, drop = FALSE]
  md <- md[keep_cells, , drop = FALSE]
  
  ## filter samples with very few cells
  n_cells_by_sample <- table(md$sample)
  good_samples <- names(n_cells_by_sample)[n_cells_by_sample >= min_cells_per_sample]
  keep_cells2 <- md$sample %in% good_samples
  X <- X[, keep_cells2, drop = FALSE]
  md <- md[keep_cells2, , drop = FALSE]
  
  ## ensure each sample maps to a single group (warn & drop ambiguous)
  grp_by_sample <- tapply(as.character(md$group), md$sample, function(v) unique(v))
  bad <- vapply(grp_by_sample, function(v) length(v) != 1, logical(1))
  if (any(bad)) {
    warning("Dropping samples with mixed IGRA_status: ",
            paste(names(bad)[bad], collapse = ", "))
  }
  keep_samples <- names(grp_by_sample)[!bad]
  keep_cells3 <- md$sample %in% keep_samples
  X <- X[, keep_cells3, drop = FALSE]
  md <- md[keep_cells3, , drop = FALSE]
  
  ## sparse aggregation via sample model matrix (genes x cells %*% cells x samples)
  smp_f <- factor(md$sample, levels = unique(md$sample))
  M <- Matrix::sparse.model.matrix(~ 0 + smp_f)   # columns are samples
  colnames(M) <- levels(smp_f)
  PB <- X %*% M                                   # genes x samples (dgCMatrix)
  
  ## sample metadata
  sample_md <- data.frame(
    sample = colnames(PB),
    group  = vapply(colnames(PB), function(s) as.character(unique(md$group[md$sample == s]))[1], character(1)),
    n_cells = as.integer(n_cells_by_sample[colnames(PB)]),
    row.names = colnames(PB),
    stringsAsFactors = FALSE
  )
  ## keep only samples from the two groups and both groups present
  sample_md$group <- factor(sample_md$group, levels = keep_groups)
  PB <- PB[, rownames(sample_md), drop = FALSE]
  
  if (nlevels(droplevels(sample_md$group)) < 2)
    stop("Need at least one sample in each of: ", paste(keep_groups, collapse = " & "))
  
  list(counts = PB, sample_md = sample_md)
}

## ---- 2) edgeR DGE (Positive vs Negative) ----
run_edger_two_group <- function(pb_counts, sample_md, ref = "Negative"){
  stopifnot(all(colnames(pb_counts) == rownames(sample_md)))
  group <- droplevels(factor(sample_md$group))
  ## set ref to Negative so logFC is Positive vs Negative
  group <- relevel(group, ref = ref)
  
  y <- DGEList(counts = pb_counts, samples = sample_md, group = group)
  keep <- filterByExpr(y, group = group)   # sensible filtering by group
  y <- y[keep,, keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group)          # Intercept + groupPositive
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "groupPositive")
  
  top <- topTags(qlf, n = Inf)$table
  ## tidy column names and add FDR
  top$gene <- rownames(top)
  top <- top[, c("gene","logFC","logCPM","F","PValue","FDR")]
  list(table = top, qlf = qlf, fit = fit, y = y, design = design)
}

## ---- 3) Put it together & run ----
## ---- Pseudobulk-level DGE (save to same folder) ----

# 1) Create pseudobulk counts by sample
pb <- pseudobulk_by_sample(
  seu,
  sample_col = "orig.ident",
  group_col  = "IGRA_status",
  assay      = "RNA",
  min_cells_per_sample = 20,
  keep_groups = c("Negative","Positive")
)

# 2) Run DGE (Positive vs Negative)
dge <- run_edger_two_group(pb$counts, pb$sample_md, ref = "Negative")

# 3) Save results to the same DGE folder
out_file_pb <- file.path(dge_dir, "Pseudobulk_DGE_IGRApos_vs_neg.csv")
write.csv(dge$table, out_file_pb, row.names = FALSE)

message("Saved pseudobulk DGE results to: ", out_file_pb)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(EnhancedVolcano)
})

## ---- Paths ----
dge_dir  <- file.path(base_dir, "DGE")
plot_dir <- file.path(dge_dir, "Volcano_Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Load pseudobulk DGE (edgeR) ----
# Default filename produced earlier; change if you used a different name.
pb_csv <- file.path(dge_dir, "Pseudobulk_DGE_IGRApos_vs_neg.csv")
if (!file.exists(pb_csv)) {
  # fall back: pick the first CSV that looks like pseudobulk
  cand <- list.files(dge_dir, pattern = "Pseudobulk.*\\.csv$", full.names = TRUE)
  if (length(cand) == 0) stop("No pseudobulk CSV found in: ", dge_dir)
  pb_csv <- cand[1]
}

deg <- read.csv(pb_csv, stringsAsFactors = FALSE)

## Expecting edgeR columns: gene, logFC, logCPM, F, PValue, FDR
req <- c("gene", "logFC", "PValue", "FDR")
missing <- setdiff(req, colnames(deg))
if (length(missing)) stop("Missing columns in pseudobulk CSV: ", paste(missing, collapse = ", "))

deg <- deg %>%
  filter(!is.na(gene) & gene != "") %>%
  mutate(
    # Safety: coerce numerics
    logFC  = as.numeric(logFC),
    PValue = as.numeric(PValue),
    FDR    = as.numeric(FDR)
  )

## ---- Coloring rules (same style as before, adapted to edgeR cols) ----
keyvals <- ifelse(
  abs(deg$logFC) > 1.5 & deg$FDR < 0.01, "#CD0BBC",
  ifelse(deg$FDR < 0.01, "#28E2E5", "gray30")
)
keyvals[is.na(keyvals)] <- "gray30"
names(keyvals)[keyvals == "gray30"]  <- "NS"
names(keyvals)[keyvals == "#28E2E5"] <- "adj(p-value) < 0.01"
names(keyvals)[keyvals == "#CD0BBC"] <- "FC > 1.5"

## ---- Volcano plot (edgeR: x=logFC, y=PValue, padj=FDR) ----
vp <- EnhancedVolcano(
  deg,
  lab           = deg$gene,
  x             = "logFC",
  y             = "PValue",
  pCutoffCol    = "FDR",
  pCutoff       = 0.01,
  FCcutoff      = 1.5,
  xlab          = bquote(~Log[2]~ 'fold change (IGRA+ vs IGRA-)'),
  pointSize     = 4.0,
  labSize       = 3.0,
  labCol        = "black",
  labFace       = "bold",
  colAlpha      = 4/5,
  legendPosition= "right",
  legendLabSize = 14,
  legendIconSize= 4.0,
  drawConnectors= TRUE,
  widthConnectors = 1.0,
  colConnectors   = "black",
  title         = "Pseudobulk: IGRA+ vs IGRA-",
  subtitle      = "edgeR (sample-level pseudobulk)",
  colCustom     = keyvals
)

## ---- Save ----
out_png <- file.path(plot_dir, "Volcano_Pseudobulk_IGRApos_vs_neg.png")
ggsave(
  filename = out_png,
  plot     = vp + guides(color = guide_legend(reverse = TRUE)),
  dpi      = 500,
  width    = 10,
  height   = 7
)

message("Saved pseudobulk volcano plot: ", out_png)
table(pb$sample_md$group)
summary(dge$table$PValue < 0.05)
summary(dge$table$FDR < 0.05)
