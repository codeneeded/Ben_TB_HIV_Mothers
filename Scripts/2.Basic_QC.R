# Quality Control Visualizations
##Ref -> https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html
##Ref -> http://bioconductor.org/books/3.15/OSCA.intro/getting-scrna-seq-datasets.html
#Ref -> https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html
# Ref -> https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html

##############################################################
# QC COMMENTS (Applicable to Both RAW and PARSEF Pipelines)
#
# The following QC steps are applied to all Seurat objects 
# (RAW and PARSEF). Each metric is assessed against expected 
# biological and technical thresholds:
#
# 1. nGenes (nFeature_RNA):
#    - High-quality data should have a single main peak.
#    - Shoulders or bimodal shapes can indicate failed cells 
#      or heterogeneous biology (e.g. quiescent populations).
#    - Cutoff example: > 250 features retained.
#
# 2. UMI Counts (nCount_RNA):
#    - Low UMI counts (< 500) indicate low sequencing depth.
#    - Between 500–1000 is borderline but may still be usable.
#
# 3. Complexity Score (log10GenesPerUMI):
#    - Ratio of genes detected per UMI (novelty).
#    - Low complexity (< 0.8) indicates over-sequencing a 
#      narrow set of transcripts (e.g. poor cell diversity).
#
# 4. Percent Mitochondrial Reads (percent_mito):
#    - High % mitochondrial reads (> 15–20%) can indicate 
#      stressed/dying cells.
#    - May vary by tissue type; context-dependent.
#
# 5. Percent Ribosomal Reads (percent_ribo):
#    - Very low ribosomal RNA (< 5%) can indicate poor quality.
#    - Very high ribosomal RNA may indicate unremoved debris 
#      or technical bias.
# Ribosomal content (percent_ribo) is expected to be low (~1%) in Parse data due to rRNA depletion.
# No hard cutoff applied. Cells with percent_ribo > 5% will be reviewed manually,
# but filtering will only occur if high-ribo clusters are identified as artifacts.
# 10x Genomics scRNA-seq ribosomal content:
# Typical cells show ~5–15% ribosomal gene expression.
# High percent_ribo (>20%) may indicate stressed, damaged, or highly proliferative cells.
# Unlike Parse (where ~1% is normal), 10x baseline ribo content is higher.
# 6. Percent Hemoglobin Genes (percent_hb):
#    - High levels (> 20%) indicate strong RBC contamination.
#
# 7. Percent Platelet Genes (percent_plat):
#    - High levels (> 2%) indicate platelet contamination.
#
# All cutoffs should be adjusted based on dataset-specific
# distributions, protocol differences, and biological context.
##############################################################
########################################
# QC for RAW MERGED SEURAT OBJECT
########################################
# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)

# Base directories
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
qc_dir_raw <- file.path(base_dir, "QC_Plots", "RAW_QC")
dir.create(qc_dir_raw, recursive = TRUE, showWarnings = FALSE)

# Load RAW Seurat object
seu_raw <- readRDS(file.path(base_dir, "saved_R_data", "seu_raw_merged.rds"))
Idents(seu_raw) <- "sample"

########################################
# 1. Exploration
########################################

exploratory_dir_raw <- file.path(qc_dir_raw, "1_Exploration")
dir.create(exploratory_dir_raw, showWarnings = FALSE)

# Cells per sample
png(file.path(exploratory_dir_raw, "Cells_per_Sample.png"), width = 1800, height = 1200)
seu_raw@meta.data %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Number of Cells per Sample")
dev.off()

# Cells per IGRA group
png(file.path(exploratory_dir_raw, "Cells_per_IGRA.png"), width = 1200, height = 1000)
seu_raw@meta.data %>%
  ggplot(aes(x = IGRA_status, fill = IGRA_status)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of Cells per IGRA Group")
dev.off()

########################################
# 2. FeatureScatter QC
########################################

fs_dir_raw <- file.path(qc_dir_raw, "2_FeatureScatter_QC")
dir.create(fs_dir_raw, showWarnings = FALSE)

# UMI vs Gene colored by %mito
png(file.path(fs_dir_raw, "UMI_vs_Gene_colored_by_Mito.png"), width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object = seu_raw,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  low_cutoff_UMI = 500,
  high_cutoff_UMI = 40000,
  meta_gradient_low_cutoff = 20,
  combination = TRUE
)
dev.off()

# Gene vs %mito
png(file.path(fs_dir_raw, "Gene_vs_Mito.png"), width = 1800, height = 1200)
QC_Plot_GenevsFeature(
  seurat_object = seu_raw,
  feature1 = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  high_cutoff_feature = 20
)
dev.off()

########################################
# 3. Histogram QC
########################################

hist_dir_raw <- file.path(qc_dir_raw, "3_Histogram_Plots")
dir.create(hist_dir_raw, showWarnings = FALSE)

metadata_raw <- seu_raw@meta.data
metadata_raw$log10GenesPerUMI <- log10(metadata_raw$nFeature_RNA / metadata_raw$nCount_RNA)

## UMI Count
png(file = file.path(hist_dir_raw, "UMI_Count.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 500) + ggtitle("UMI Count Distribution")
dev.off()

## nGenes
png(file = file.path(hist_dir_raw, "nGenes.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 250) + ggtitle("Gene Count Distribution")
dev.off()

## Complexity Score
png(file = file.path(hist_dir_raw, "Complexity_Score.png"), width = 1800, height = 1200)
metadata_raw %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic() +
  geom_vline(xintercept = 0.8) + ggtitle("Complexity Score (Genes per UMI)")
dev.off()

## Mito Ratio
png(file = file.path(hist_dir_raw, "Mito_Ratio.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = percent_mito, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 15) + ggtitle("Mitochondrial Ratio Distribution")
dev.off()

## Ribo Ratio
png(file = file.path(hist_dir_raw, "Ribo_Ratio.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = percent_ribo, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 5) + ggtitle("Ribosomal Ratio Distribution")
dev.off()

## Heme Ratio
png(file = file.path(hist_dir_raw, "Heme_Ratio.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = percent_hb, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 20) + ggtitle("Hemoglobin Ratio Distribution")
dev.off()

## Platelet Ratio
png(file = file.path(hist_dir_raw, "Platelet_Ratio.png"), width = 1800, height = 1200)
metadata_raw %>% 
  ggplot(aes(color = orig.ident, x = percent_plat, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 2) + ggtitle("Platelet Ratio Distribution")
dev.off()


########################################
# QC for PARSEF MERGED SEURAT OBJECT
########################################

qc_dir_parsef <- file.path(base_dir, "QC_Plots", "PARSEF_QC")
dir.create(qc_dir_parsef, recursive = TRUE, showWarnings = FALSE)

# Load PARSEF Seurat object
seu_parsef <- readRDS(file.path(base_dir, "saved_R_data", "seu_parsef_merged.rds"))
Idents(seu_parsef) <- "sample"

########################################
# 1. Exploration
########################################

exploratory_dir_parsef <- file.path(qc_dir_parsef, "1_Exploration")
dir.create(exploratory_dir_parsef, showWarnings = FALSE)

png(file.path(exploratory_dir_parsef, "Cells_per_Sample.png"), width = 1800, height = 1200)
seu_parsef@meta.data %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Number of Cells per Sample")
dev.off()

png(file.path(exploratory_dir_parsef, "Cells_per_IGRA.png"), width = 1200, height = 1000)
seu_parsef@meta.data %>%
  ggplot(aes(x = IGRA_status, fill = IGRA_status)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of Cells per IGRA Group")
dev.off()

########################################
# 2. FeatureScatter QC
########################################

fs_dir_parsef <- file.path(qc_dir_parsef, "2_FeatureScatter_QC")
dir.create(fs_dir_parsef, showWarnings = FALSE)

png(file.path(fs_dir_parsef, "UMI_vs_Gene_colored_by_Mito.png"), width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object = seu_parsef,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  low_cutoff_UMI = 500,
  high_cutoff_UMI = 40000,
  meta_gradient_low_cutoff = 20,
  combination = TRUE
)
dev.off()

png(file.path(fs_dir_parsef, "Gene_vs_Mito.png"), width = 1800, height = 1200)
QC_Plot_GenevsFeature(
  seurat_object = seu_parsef,
  feature1 = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  high_cutoff_feature = 20
)
dev.off()

########################################
# 3. Histogram QC
########################################

hist_dir_parsef <- file.path(qc_dir_parsef, "3_Histogram_Plots")
dir.create(hist_dir_parsef, showWarnings = FALSE)

metadata_parsef <- seu_parsef@meta.data
metadata_parsef$log10GenesPerUMI <- log10(metadata_parsef$nFeature_RNA / metadata_parsef$nCount_RNA)

png(file = file.path(hist_dir_parsef, "UMI_Count.png"), width = 1800, height = 1200)
metadata_parsef %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 500) + ggtitle("UMI Count Distribution")
dev.off()

png(file = file.path(hist_dir_parsef, "nGenes.png"), width = 1800, height = 1200)
metadata_parsef %>% png(file.path(exploratory_dir_parsef, "Cells_per_Sample.png"), width = 1800, height = 1200)
seu_parsef@meta.data %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Number of Cells per Sample")
dev.off()

png(file.path(exploratory_dir_parsef, "Cells_per_IGRA.png"), width = 1200, height = 1000)
seu_parsef@meta.data %>%
  ggplot(aes(x = IGRA_status, fill = IGRA_status)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of Cells per IGRA Group")
dev.off()

  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 250) + ggtitle("Gene Count Distribution")
dev.off()

png(file = file.path(hist_dir_parsef, "Complexity_Score.png"), width = 1800, height = 1200)
metadata_parsef %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic() +
  geom_vline(xintercept = 0.8) + ggtitle("Complexity Score (Genes per UMI)")
dev.off()

png(file = file.path(hist_dir_parsef, "Mito_Ratio.png"), width = 1800, height = 1200)
metadata_parsef %>% 
  ggplot(aes(color = orig.ident, x = percent_mito, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 15) + ggtitle("Mitochondrial Ratio Distribution")
dev.off()

png(file = file.path(hist_dir_parsef, "Ribo_Ratio.png"), width = 1800, height = 1200)
metadata_parsef %>% 
  ggplot(aes(color = orig.ident, x = percent_ribo, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 5) + ggtitle("Ribosomal Ratio Distribution")
dev.off()

png(file = file.path(hist_dir_parsef, "Heme_Ratio.png"), width = 1800, height = 1200)
metadata_parsef %>% 
  ggplot(aes(color = orig.ident, x = percent_hb, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 20) + ggtitle("Hemoglobin Ratio Distribution")
dev.off()

png(file = file.path(hist_dir_parsef, "Platelet_Ratio.png"), width = 1800, height = 1200)
metadata_parsef %>% 
  ggplot(aes(color = orig.ident, x = percent_plat, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 2) + ggtitle("Platelet Ratio Distribution")
dev.off()

############################### FINAL QC ##########################################################
# ------------------- #
# 1. Setup QC folders #
# ------------------- #
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/QC_Plots"
qc_final_dir <- file.path(base_dir, "FINAL_QC")
dir.create(qc_final_dir, showWarnings = FALSE)
# ---------------------- #
# 2. Pre-QC Visualization
# ---------------------- #

feats_parse <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_hb", "percent_plat")

png(file.path(qc_final_dir, "PreQC_Features_Grouped.png"), width = 1800, height = 1200)
VlnPlot(seu_parsef, group.by = "orig.ident", features = feats_parse, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

# ---------------------------------------- #
# 3. QC Filtering - Parse-specific cutoffs #
# ---------------------------------------- #
seu_parsef <- JoinLayers(seu_parsef)

# 📝 Parse-specific cutoff notes:
# nUMI >= 500 
# nGene >= 600 
# percent_mito < 15
# percent_hb < 20
# percent_plat < 2 

# Apply filters
filtered_parse <- subset(
  seu_parsef,
  subset = (nCount_RNA >= 500) &
    (nFeature_RNA >= 600) &
    (percent_mito < 15) &
    (percent_hb < 20) &
    (percent_plat < 2)
)


# ----------------------- #
# 4. Post-QC Visualization
# ----------------------- #
# Cells per sample after filtering
png(file.path(qc_final_dir, "Cells_per_Sample_PostQC.png"), width = 1800, height = 1200)
filtered_parse@meta.data %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Number of Cells per Sample (Post-QC)")
dev.off()

# Cells per IGRA group after filtering
png(file.path(qc_final_dir, "Cells_per_IGRA_PostQC.png"), width = 1000, height = 700)
filtered_parse@meta.data %>%
  ggplot(aes(x = IGRA_status, fill = IGRA_status)) +
  geom_bar(width = 0.6, color = "black", alpha = 0.9) +
  scale_fill_manual(
    values = c("Negative" = "#A3F8A9", "Positive" = "#FF918A"),
    labels = c("Negative" = "IGRA Negative", "Positive" = "IGRA Positive"),
    name = "IGRA Status"
  ) +
  geom_text(
    stat = "count",
    aes(label = scales::comma(..count..)),
    vjust = -0.5,
    size = 6
  ) +
  labs(
    title = "Number of Cells per IGRA Group (Post-QC)",
    x = "IGRA Status",
    y = "Cell Count"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold")
  ) +
  ylim(0, NA)  # ensures text labels fit above bars
dev.off()

##############


##############
png(file.path(qc_final_dir, "PostQC_Features_Grouped.png"), width = 1800, height = 1200)
VlnPlot(filtered_parse, group.by = "orig.ident", features = feats_parse, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

# ---------------------- #
# 5. Save Post-QC Object #
# ---------------------- #
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
saveRDS(filtered_parse, file.path(base_dir, "saved_R_data", "seu_parsef_filtered.rds"))
