suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readr)
})

# -----------------------------
# Paths
# -----------------------------
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
saved_dir <- file.path(base_dir, "saved_R_data")
rds_in <- file.path(saved_dir, "seu_annotated.rds")

out_dir <- file.path(base_dir, "Module_Scoring")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(rds_in)
DefaultAssay(seu) <- ifelse("RNA" %in% Assays(seu), "RNA", DefaultAssay(seu))
assay_use <- DefaultAssay(seu)
cluster_col <- "IGRA_Annotation"
group_col   <- "IGRA_status"

# -----------------------------
# Module definitions (line by line, grouped in families)
# -----------------------------

modules <- list(
  
  # ---- TGF-beta ----
  TGFb = list(
    ReceptorsLigands = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3","ENG","ACVRL1"),
    SMAD_Core        = c("SMAD2","SMAD3","SMAD4","SMAD7","SKI","SKIL","PMEPA1"),
    Downstream       = c("FOXP3","IKZF2","CTLA4","LRRC32","ITGAE","CCR7","SELL","SERPINE1","SERPINE1","PTPN14","COL1A1","COL3A1"),
    Noncanonical     = c("MAPK14","MAPK8","MAPK9","MAP3K7","FOS","JUN","PIK3CD","AKT1","MTOR","RPTOR","NFKBIA","RELB")
  ),
  
  # ---- IL-10 / STAT3 ----
  IL10_STAT3 = list(
    ReceptorsTransducers = c("IL10RA","IL10RB","STAT3"),
    TFH_Adjacent_TFs     = c("BCL6","MAF","BATF","PIM1"),
    Treg_Overlap         = c("FOXP3","CTLA4","IKZF2"),
    Myeloid_M2_Hue       = c("MRC1","MSR1","CD163")
  ),
  
  # ---- Reservoir CD4 ----
  Reservoir_CD4 = list(
    Tfh_Core         = c("CXCR5","PDCD1","ICOS","BCL6","MAF","BATF","SH2D1A","IL21"),
    Tph_Discriminator= c("PRDM1","CXCL13","CCR2","PDCD1","ICOS"),
    Reservoir_Memory = c("TCF7","CCR7","SELL","BCL2","IL7R")
  ),
  
  # ---- ISG ----
  ISG = list(
    TypeI  = c("ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","MX2","OAS1","OAS2","OAS3","RSAD2",
               "IFITM1","IFITM2","IFITM3","BST2","TRIM5","APOBEC3G","SAMHD1"),
    TypeII = c("STAT1","IRF1","GBP1","GBP5","CXCL9","CXCL10")
  ),
  
  # ---- Myeloid ----
  Myeloid = list(
    Inflammatory  = c("IL1B","TNF","NFKBIA","NFKBIZ","CCL2","CXCL8","SERPINA1","S100A8","S100A9"),
    Regulatory    = c("IL10RA","IL10RB","STAT3","SOCS3","TGFB1","TGFBR1","TGFBR2","MRC1","MSR1","CD163"),
    DC_Maturation = c("CCR7","LAMP3","FSCN1","CCL19","CD40","RELB")
  ),
  
  # ---- NK ----
  NK = list(
    Inhibition  = c("KLRC1","KLRD1"),
    Effector    = c("PRF1","GNLY","NKG7","GZMK"),
    Activation  = c("IFNG","XCL1","XCL2","CCL5")
  ),
  
  # ---- CD8 ----
  CD8 = list(
    Effector    = c("PRF1","GNLY","NKG7","CCL5","XCL1","XCL2","IFNG"),
    Memory      = c("TCF7","CCR7","SELL","BCL2","IL7R"),
    Exhaustion  = c("PDCD1","TIGIT","LAG3","HAVCR2","TOX","EOMES")
  ),
  
  # ---- B cells ----
  Bcell = list(
    GC_like     = c("BCL6","CXCR5","AICDA","S1PR2","RGS13","MEF2B","BACH2"),
    Atypical_ABC= c("TBX21","ITGAX","FCRL5","ZEB2"),
    Plasmablast = c("PRDM1","XBP1","JCHAIN","MZB1","SDC1","TNFRSF17")
  ),
  
  # ---- Autophagy ----
  Autophagy = list(
    Initiation   = c("ULK1","RB1CC1","ATG13","ATG101","BECN1","PIK3C3","PIK3R4","ATG14","WIPI1","WIPI2","ATG2A","ATG2B"),
    Conjugation  = c("ATG5","ATG12","ATG16L1","ATG7","ATG3","MAP1LC3A","MAP1LC3B","MAP1LC3B2","GABARAP","GABARAPL1","GABARAPL2"),
    Fusion       = c("RAB7A","STX17","SNAP29","VAMP8","EPG5"),
    Lysosome     = c("LAMP1","LAMP2","CTSD","CTSB","ATP6V1A","ATP6V1B2","ATP6V0D1","PSAP","MCOLN1"),
    CLEAR_TFs    = c("TFEB","TFE3","MITF"),
    Cargo        = c("SQSTM1","NBR1","OPTN","CALCOCO2","TAX1BP1","TOLLIP"),
    Mitophagy    = c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13"),
    Xenophagy    = c("NOD2","RIPK2","TBK1","OPTN","CALCOCO2","SQSTM1","TAX1BP1","LRSAM1","ATG16L1","ATG5","ATG7"),
    cGAS_STING   = c("STING1","CGAS"),
    mTORC        = c("MTOR","RPTOR","MLST8","RHEB","TSC1","TSC2","DEPDC5"),
    AMPK         = c("PRKAA1","PRKAA2","PRKAB1","PRKAG1","STK11")
  )
)

# -----------------------------
# Scoring + plotting
# -----------------------------
score_module <- function(obj, genes, name) {
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0) {
    obj[[paste0("MS_", name)]] <- NA
    return(obj)
  }
  obj <- AddModuleScore(obj, features = list(present), name = paste0("MS_", name), assay = assay_use)
  col <- paste0("MS_", name, "1")
  obj[[paste0("MS_", name)]] <- obj[[col]]
  obj[[col]] <- NULL
  return(obj)
}

plot_module <- function(obj, module_field, folder, title) {
  md <- obj@meta.data
  md <- md[!is.na(md[[group_col]]) & !is.na(md[[cluster_col]]), ]
  gg <- ggplot(md, aes_string(x = group_col, y = module_field, fill = group_col)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.3) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(as.formula(paste0("~", cluster_col)), scales = "free_y") +
    theme_minimal(base_size = 14) +
    xlab("") + ylab("Module Score") +
    ggtitle(title)
  ggsave(file.path(folder, paste0(module_field, ".png")), gg, width = 20, height = 18, dpi = 300, bg="white")
}

# Helper to choose reference group (pos vs neg) else first level
pick_ref_levels <- function(vec) {
  lv <- unique(as.character(vec))
  ref_idx <- grep("(pos|positive|igra\\+)", lv, ignore.case = TRUE)
  if (length(ref_idx) > 0) {
    ref <- lv[ref_idx[1]]
    other <- setdiff(lv, ref)
    if (length(other) == 1) return(c(ref, other))
  }
  # Fallback: first two unique levels in order encountered
  if (length(lv) >= 2) return(lv[1:2])
  return(lv)
}

# Compute Wilcoxon p and HL shift per module x cluster
module_stats <- list()

coverage <- data.frame()

for (family in names(modules)) {
  family_dir <- file.path(out_dir, family)
  dir.create(family_dir, showWarnings = FALSE)
  for (sub in names(modules[[family]])) {
    mod_name <- paste(family, sub, sep="_")
    field <- paste0("MS_", mod_name)
    seu <- score_module(seu, modules[[family]][[sub]], mod_name)
    plot_module(seu, field, family_dir, paste(family, "-", sub))
    
    present <- intersect(modules[[family]][[sub]], rownames(seu))
    missing <- setdiff(modules[[family]][[sub]], present)
    coverage <- rbind(coverage,
                      data.frame(Module=mod_name, Present=length(present), Missing=length(missing),
                                 Present_Genes=paste(present, collapse=";"),
                                 Missing_Genes=paste(missing, collapse=";"))
    )
    
    # Collect stats per cluster
    md <- seu@meta.data
    if (!field %in% colnames(md)) next
    md <- md[!is.na(md[[group_col]]) & !is.na(md[[cluster_col]]) & !is.na(md[[field]]), 
             c(group_col, cluster_col, field)]
    if (nrow(md) == 0) next
    
    for (cl in sort(unique(md[[cluster_col]]))) {
      subdf <- md[md[[cluster_col]] == cl, ]
      grps <- unique(as.character(subdf[[group_col]]))
      if (length(grps) != 2) {
        module_stats[[length(module_stats)+1]] <- data.frame(
          Module = mod_name, Cluster = cl, PValue = NA_real_, HL_Shift = NA_real_
        )
        next
      }
      
      # Pick consistent order: ref (pos) vs other; HL_Shift is ref - other
      ord <- pick_ref_levels(subdf[[group_col]])
      if (length(ord) != 2) {
        module_stats[[length(module_stats)+1]] <- data.frame(
          Module = mod_name, Cluster = cl, PValue = NA_real_, HL_Shift = NA_real_
        )
        next
      }
      
      x <- subdf[subdf[[group_col]] == ord[1], field]
      y <- subdf[subdf[[group_col]] == ord[2], field]
      
      # Require minimal n per group
      if (length(x) < 3 || length(y) < 3) {
        module_stats[[length(module_stats)+1]] <- data.frame(
          Module = mod_name, Cluster = cl, PValue = NA_real_, HL_Shift = NA_real_
        )
        next
      }
      
      wt <- tryCatch(
        wilcox.test(x, y, alternative = "two.sided", exact = FALSE, conf.int = TRUE),
        error = function(e) NULL
      )
      if (is.null(wt)) {
        module_stats[[length(module_stats)+1]] <- data.frame(
          Module = mod_name, Cluster = cl, PValue = NA_real_, HL_Shift = NA_real_
        )
      } else {
        # wt$estimate is the Hodges–Lehmann location shift (x - y) when conf.int=TRUE
        hl <- unname(if (!is.null(wt$estimate)) wt$estimate else median(x) - median(y))
        module_stats[[length(module_stats)+1]] <- data.frame(
          Module = mod_name,
          Cluster = cl,
          PValue = wt$p.value,
          HL_Shift = hl
        )
      }
    }
  }
}

# -----------------------------
# Outputs
# -----------------------------
# Coverage table
write_tsv(coverage, file.path(out_dir, "Module_Gene_Coverage.tsv"))

# Stats table (Module x Cluster with p-value and HL shift)
stats_df <- dplyr::bind_rows(module_stats)
stats_csv <- file.path(out_dir, "Module_Stats_byCluster.csv")
readr::write_csv(stats_df, stats_csv)

