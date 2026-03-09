
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(viridis)
library(stringr)
library(scales)
library(tidyr)
library(SeuratExtend)

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
    "TGFB–TGFBR module"             = c("ACVRL1","ENG","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3"),
    "SMAD2/3–SMAD4 module"          = c("PMEPA1","SKI","SKIL","SMAD2","SMAD3","SMAD4","SMAD7"),
    "TGFβ transcriptional targets"  = c("CCR7","COL1A1","COL3A1","CTLA4","FOXP3","IKZF2","ITGAE","LRRC32","PTPN14","SELL","SERPINE1"),
    "TGFβ noncanonical (MAPK/PI3K)" = c("AKT1","FOS","JUN","MAP3K7","MAPK14","MAPK8","MAPK9","MTOR","NFKBIA","PIK3CD","RELB","RPTOR")
  ),
  
  # ---- IL-10 / STAT3 ----
  IL10_STAT3 = list(
    "IL10RA/IL10RB–JAK1/TYK2–STAT3 module" = c("IL10RA","IL10RB","STAT3"),
    "STAT3-Tfh-linked targets"             = c("BATF","BCL6","MAF","PIM1"),
    "STAT3-regulatory/Treg-overlap"        = c("CTLA4","FOXP3","IKZF2"),
    "Alternative-activation (M2-like) signature" = c("CD163","MRC1","MSR1")
  ),
  
  # ---- Reservoir CD4 ----
  CD4_Reservoir = list(
    Tfh_Core          = c("BCL6","BATF","CXCR5","ICOS","IL21","MAF","PDCD1","SH2D1A"),
    Tph_Discriminator = c("CCR2","CXCL13","ICOS","PDCD1","PRDM1"),
    Reservoir_Memory  = c("BCL2","CCR7","IL7R","SELL","TCF7")
  ),
  
  # ---- ISG ----
  ISG = list(
    TypeI  = c("APOBEC3G","BST2","IFI44L","IFI6","IFIT1","IFIT3","IFITM1","IFITM2","IFITM3",
               "ISG15","MX1","MX2","OAS1","OAS2","OAS3","RSAD2","SAMHD1","TRIM5"),
    TypeII = c("CIITA","CXCL10","CXCL11","GBP1","GBP2","GBP3","GBP4","GBP5","GBP6","HLA-DPA1",
               "HLA-DPB1","HLA-DRA","HLA-DRB1","IRF1","IRF8","PSMB8","PSMB9","SOCS1","STAT1","TAP1","TAP2")
  ),
  
  # ---- Myeloid ----
  Myeloid = list(
    Inflammatory  = c("CCL2","CXCL8","IL1B","NFKBIA","NFKBIZ","S100A8","S100A9","SERPINA1","TNF"),
    Regulatory    = c("CD163","IL10RA","IL10RB","MRC1","MSR1","SOCS3","STAT3","TGFB1","TGFBR1","TGFBR2"),
    DC_Maturation = c("CCL19","CD40","CCR7","FSCN1","LAMP3","RELB")
  ),
  
  # ---- NK ----
  NK = list(
    Inhibition  = c("KLRC1","KLRD1"),
    Effector    = c("GNLY","GZMK","NKG7","PRF1"),
    Activation  = c("CCL5","IFNG","XCL1","XCL2")
  ),
  
  # ---- CD8 ----
  CD8 = list(
    Effector    = c("CCL5","GNLY","IFNG","NKG7","PRF1","XCL1","XCL2"),
    Memory      = c("BCL2","CCR7","IL7R","SELL","TCF7"),
    Exhaustion  = c("EOMES","HAVCR2","LAG3","PDCD1","TIGIT","TOX")
  ),
  
  # ---- B cells ----
  Bcell = list(
    GC_like      = c("AICDA","BACH2","BCL6","CXCR5","MEF2B","RGS13","S1PR2"),
    Atypical_ABC = c("FCRL5","ITGAX","TBX21","ZEB2"),
    Plasmablast  = c("JCHAIN","MZB1","PRDM1","SDC1","TNFRSF17","XBP1")
  ),
  
  # ---- Autophagy ----
  Autophagy = list(
    Initiation   = c("ATG101","ATG13","ATG14","ATG2A","ATG2B","BECN1","PIK3C3","PIK3R4","RB1CC1","ULK1","WIPI1","WIPI2"),
    Conjugation  = c("ATG12","ATG16L1","ATG3","ATG5","ATG7","GABARAP","GABARAPL1","GABARAPL2",
                     "MAP1LC3A","MAP1LC3B","MAP1LC3B2"),
    Fusion       = c("EPG5","RAB7A","SNAP29","STX17","VAMP8"),
    Lysosome     = c("ATP6V0D1","ATP6V1A","ATP6V1B2","CTSB","CTSD","LAMP1","LAMP2","MCOLN1","PSAP"),
    CLEAR_TFs    = c("MITF","TFE3","TFEB"),
    Cargo        = c("CALCOCO2","NBR1","OPTN","SQSTM1","TAX1BP1","TOLLIP"),
    Mitophagy    = c("BCL2L13","BNIP3","BNIP3L","FKBP8","FUNDC1","PHB2","PINK1","PRKN"),
    Xenophagy    = c("ATG16L1","ATG5","ATG7","CALCOCO2","LRSAM1","NOD2","OPTN","RIPK2","SQSTM1","TAX1BP1","TBK1"),
    cGAS_STING   = c("CGAS","STING1"),
    mTORC        = c("DEPDC5","MLST8","MTOR","RHEB","RPTOR","TSC1","TSC2"),
    AMPK         = c("PRKAA1","PRKAA2","PRKAB1","PRKAG1","STK11")
  ),
  
  # ---- Reservoir Favoring (all former "#1" sets) ----
  Reservoir_Favoring = list(
    NFkB_activity = sort(c(
      "BCL2A1","BIRC2","BIRC3","CCL2","CCL3","CCL4","CCL5","CXCL8","ICAM1",
      "IL1B","IL6","NFKB1","NFKBIA","PTGS2","RELA","RELB","SELE","TNF","TNFAIP2",
      "TNFAIP3","TRAF1","TRAF2","VCAM1","NOD2"
    )),
    STAT_reactivation_axis = sort(c(
      "CD40LG","CISH","CXCL9","CXCL10","IFNGR1","IFNGR2","IL6R","IL6ST","IRF1",
      "ICOS","JAK1","JAK2","OSMR","PIM1","SOCS1","SOCS3","STAT1","STAT3","STAT5A","STAT5B","TYK2"
    )),
    MYC_targets = sort(c(
      "CDK4","E2F1","EIF4A1","EIF4E","EIF5A","HSPD1","HSP90AA1","LDHA","MAX","MYC",
      "NCL","NPM1","PABPC1","PCNA","RPL11","RPLP0","RPS6","SLC2A1","SRSF1","ODC1"
    )),
    Cell_cycle_E2F_G2M = sort(c(
      "AURKA","AURKB","BUB1","BUB1B","CCNA2","CCNB1","CCNB2","CDC20","CDC25A","CDC25C",
      "CDK1","CDK2","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MKI67","PCNA","PLK1","TK1","TOP2A"
    )),
    Glycolysis = sort(c(
      "ALDOA","ENO1","GAPDH","GPI","HK2","LDHA","LDHB","PDK1","PDK3","PFKFB3","PFKP","PGAM1",
      "PGK1","PKM","SLC2A1","SLC2A3","TPI1"
    )),
    OXPHOS = sort(c(
      "ATP5F1A","ATP5F1B","ATP5MC1","ATP5PF","COX4I1","COX5B","COX6C","CYC1","NDUFA5","NDUFA9",
      "NDUFB8","NDUFS1","SDHB","SDHC","UQCRC1","UQCRC2","UQCRQ"
    )),
    Apoptosis = sort(c(
      "BCL2","BCL2A1","BCL2L1","BCL2L2","BIRC2","BIRC3","CFLAR","MDM2","MCL1","TRAF1","TRAF2","XIAP"
    )),
    Chromatin_repression = sort(c(
      "CHD4","DNMT1","DNMT3A","DNMT3B","EED","EZH2","GATAD2A","HDAC1","HDAC2",
      "KDM1A","MBD2","MTA1","MTA2","NCOR1","NCOR2","RCOR1","SETDB1","SIN3A","SIN3B","SUZ12"
    ))
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

# -------- Safe IDs (ASCII-only) for column names --------
safe_id <- function(x) {
  y <- iconv(x, to = "ASCII//TRANSLIT")           # turn “β” → "b", etc.
  y <- gsub("[^A-Za-z0-9]+", "_", y, perl = TRUE) # non-alnum → "_"
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}

# Build a safe Seurat meta field name from family + submodule (pretty)
module_field_name <- function(family, sub_pretty) {
  paste0("MS_", paste(family, safe_id(sub_pretty), sep = "_"))
}

# -------- Pretty display names (unchanged) --------
pretty_module_name <- function(mod_name) {
  x <- mod_name
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\bNFkB\\b", "NF-κB", x)
  x <- gsub("\\bOXPHOS\\b", "Oxidative phosphorylation", x)
  x <- gsub("\\bSTAT reactivation axis\\b", "STAT reactivation axis", x)
  x <- gsub("\\bCell cycle E2F G2M\\b", "Cell-cycle/E2F/G2M", x)
  x
}

# -------- Fixed: safe filenames (no “invalid range” error) --------
safe_filename <- function(label) {
  x <- iconv(label, to = "ASCII//TRANSLIT")
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[^A-Za-z0-9_. -]", "", x, perl = TRUE)  # put "-" at end of [] or use perl
  x <- gsub("\\s+", "_", x)
  x
}


plot_module <- function(obj, module_field, folder, title) {
  md <- obj@meta.data
  
  # filter valid rows
  keep <- !is.na(md[[group_col]]) & !is.na(md[[cluster_col]]) & !is.na(md[[module_field]])
  md <- md[keep, , drop = FALSE]
  if (nrow(md) == 0) return(invisible(NULL))
  
  # build a safe DF with standard column names to avoid parse() on weird names
  df <- data.frame(
    .group   = md[[group_col]],
    .cluster = md[[cluster_col]],
    .score   = md[[module_field]],
    stringsAsFactors = FALSE
  )
  
  # force order if exactly these levels exist
  if (all(c("Negative", "Positive") %in% unique(df$.group))) {
    df$.group <- factor(df$.group, levels = c("Negative", "Positive"))
  }
  
  pretty_title <- pretty_module_name(title)
  file_stub    <- safe_filename(pretty_module_name(module_field))
  
  gg <- ggplot(df, aes(x = .group, y = .score, fill = .group)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.18, size = 0.4, alpha = 0.25) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.x.npc = "center",
      label.y.npc = 0.9,
      size = 7,
      color = "red"
    ) +
    facet_wrap(~ .cluster, scales = "free_y") +
    scale_fill_manual(
      values = c("Positive" = "#FF918A", "Negative" = "#A3F8A9"),
      labels = c("Negative" = "IGRA Negative", "Positive" = "IGRA Positive"),
      name   = "IGRA Status"
    ) +
    labs(
      title = pretty_title,
      x = "IGRA Status",
      y = "Module Score"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 16) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 14, face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(folder, paste0(file_stub, ".png")),
    plot = gg, width = 20, height = 18, dpi = 300, bg = "white"
  )
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

# -----------------------------
# Compute Wilcoxon p and HL shift per module x cluster
# -----------------------------
module_stats <- list()
coverage <- data.frame()

for (family in names(modules)) {
  family_dir <- file.path(out_dir, safe_id(family))
  dir.create(family_dir, showWarnings = FALSE, recursive = TRUE)
  
  el <- modules[[family]]
  
  # ---------- CASE A: family is a single flat module (character vector) ----------
  if (is.character(el)) {
    genes <- el
    # Use a safe field based only on the family name (avoid duplication family_family)
    field     <- paste0("MS_", safe_id(family))
    mod_id    <- sub("^MS_", "", field)
    mod_pretty<- pretty_module_name(family)
    # Title for flat modules = just the module pretty name
    mod_label <- mod_pretty
    
    # score + plot
    seu <- score_module(seu, genes, mod_id)
    plot_module(seu, field, family_dir, mod_label)
    
    # coverage
    present <- intersect(genes, rownames(seu))
    missing <- setdiff(genes, rownames(seu))
    coverage <- rbind(
      coverage,
      data.frame(
        Module         = mod_id,
        Module_Pretty  = mod_pretty,
        Family         = family,
        Submodule      = NA_character_,
        Present        = length(present),
        Missing        = length(missing),
        Present_Genes  = paste(present, collapse = ";"),
        Missing_Genes  = paste(missing, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
    
    # stats per cluster
    md <- seu@meta.data
    if (field %in% colnames(md)) {
      keep <- !is.na(md[[group_col]]) & !is.na(md[[cluster_col]]) & !is.na(md[[field]])
      if (any(keep)) {
        md2 <- md[keep, c(group_col, cluster_col, field), drop = FALSE]
        for (cl in sort(unique(md2[[cluster_col]]))) {
          subdf <- md2[md2[[cluster_col]] == cl, , drop = FALSE]
          grps  <- unique(as.character(subdf[[group_col]]))
          
          if (length(grps) != 2 || nrow(subdf) < 6) {
            module_stats[[length(module_stats) + 1]] <- data.frame(
              Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
              PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
            )
            next
          }
          
          ord <- pick_ref_levels(subdf[[group_col]])
          if (length(ord) != 2) {
            module_stats[[length(module_stats) + 1]] <- data.frame(
              Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
              PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
            )
            next
          }
          
          x <- subdf[subdf[[group_col]] == ord[1], field]
          y <- subdf[subdf[[group_col]] == ord[2], field]
          if (length(x) < 3 || length(y) < 3) {
            module_stats[[length(module_stats) + 1]] <- data.frame(
              Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
              PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
            )
            next
          }
          
          wt <- tryCatch(wilcox.test(x, y, alternative = "two.sided", exact = FALSE, conf.int = TRUE),
                         error = function(e) NULL)
          if (is.null(wt)) {
            module_stats[[length(module_stats) + 1]] <- data.frame(
              Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
              PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
            )
          } else {
            hl <- unname(if (!is.null(wt$estimate)) wt$estimate else median(x) - median(y))
            module_stats[[length(module_stats) + 1]] <- data.frame(
              Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
              PValue = wt$p.value, HL_Shift = hl, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    next  # done with flat family; move to next family
  }
  
  # ---------- CASE B: family is a list of submodules ----------
  for (sub_pretty in names(el)) {
    genes    <- el[[sub_pretty]]
    
    # Safe meta-data field name
    field     <- module_field_name(family, sub_pretty)
    mod_id    <- sub("^MS_", "", field)
    mod_pretty<- paste(family, sub_pretty, sep = " : ")
    
    # Titles:
    # show only submodule for TGFb, IL10_STAT3, Reservoir_CD4; otherwise "Family - Sub"
    if (family %in% c("IL10_STAT3", "TGFb", "Reservoir_CD4")) {
      mod_label <- sub_pretty
    } else {
      mod_label <- paste(family, "-", sub_pretty)
    }
    
    # score + plot
    seu <- score_module(seu, genes, mod_id)
    plot_module(seu, field, family_dir, mod_label)
    
    # coverage
    present <- intersect(genes, rownames(seu))
    missing <- setdiff(genes, rownames(seu))
    coverage <- rbind(
      coverage,
      data.frame(
        Module         = mod_id,
        Module_Pretty  = mod_pretty,
        Family         = family,
        Submodule      = sub_pretty,
        Present        = length(present),
        Missing        = length(missing),
        Present_Genes  = paste(present, collapse = ";"),
        Missing_Genes  = paste(missing, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
    
    # stats per cluster
    md <- seu@meta.data
    if (!field %in% colnames(md)) next
    keep <- !is.na(md[[group_col]]) & !is.na(md[[cluster_col]]) & !is.na(md[[field]])
    if (!any(keep)) next
    md2 <- md[keep, c(group_col, cluster_col, field), drop = FALSE]
    
    for (cl in sort(unique(md2[[cluster_col]]))) {
      subdf <- md2[md2[[cluster_col]] == cl, , drop = FALSE]
      grps  <- unique(as.character(subdf[[group_col]]))
      
      if (length(grps) != 2) {
        module_stats[[length(module_stats) + 1]] <- data.frame(
          Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
          PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
        )
        next
      }
      
      ord <- pick_ref_levels(subdf[[group_col]])
      if (length(ord) != 2) {
        module_stats[[length(module_stats) + 1]] <- data.frame(
          Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
          PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
        )
        next
      }
      
      x <- subdf[subdf[[group_col]] == ord[1], field]
      y <- subdf[subdf[[group_col]] == ord[2], field]
      if (length(x) < 3 || length(y) < 3) {
        module_stats[[length(module_stats) + 1]] <- data.frame(
          Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
          PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
        )
        next
      }
      
      wt <- tryCatch(wilcox.test(x, y, alternative = "two.sided", exact = FALSE, conf.int = TRUE),
                     error = function(e) NULL)
      if (is.null(wt)) {
        module_stats[[length(module_stats) + 1]] <- data.frame(
          Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
          PValue = NA_real_, HL_Shift = NA_real_, stringsAsFactors = FALSE
        )
      } else {
        hl <- unname(if (!is.null(wt$estimate)) wt$estimate else median(x) - median(y))
        module_stats[[length(module_stats) + 1]] <- data.frame(
          Module = mod_id, Module_Pretty = mod_pretty, Cluster = cl,
          PValue = wt$p.value, HL_Shift = hl, stringsAsFactors = FALSE
        )
      }
    }
  }
}
# -----------------------------
# Outputs
# -----------------------------

# -----------------------------
# Composite module scores (per-cell) + plots
# -----------------------------
md <- seu@meta.data

# Find submodule fields (per-cell)
find_fields <- function(prefix) grep(paste0("^MS_", prefix, "_"), colnames(md), value = TRUE)

# Composite IL-10 and TGFb per cell (average all submodules inside each family)
il10_fields <- find_fields("IL10_STAT3")
tgfb_fields <- find_fields("TGFb")

if (length(il10_fields) > 0) {
  md$MS_IL10_composite <- rowMeans(md[, il10_fields, drop = FALSE], na.rm = TRUE)
}
if (length(tgfb_fields) > 0) {
  md$MS_TGFb_composite <- rowMeans(md[, tgfb_fields, drop = FALSE], na.rm = TRUE)
}

# Save back
seu@meta.data <- md

# Plot the two per-cell composites (boxplot IGRA +/- faceted by cluster)
composite_dir <- file.path(out_dir, "Composites")
dir.create(composite_dir, showWarnings = FALSE, recursive = TRUE)
if ("MS_IL10_composite" %in% colnames(seu@meta.data)) {
  plot_module(seu, "MS_IL10_composite", composite_dir, "IL-10 composite")
}
if ("MS_TGFb_composite" %in% colnames(seu@meta.data)) {
  plot_module(seu, "MS_TGFb_composite", composite_dir, "TGFβ composite")
}

# -----------------------------
# Per-cell RPI and mTOR–MAPK (MAI) (z-scored within cluster) + plots
# -----------------------------
# -----------------------------
# Per-cell RPI (updated for Reservoir_Favoring submodules)
# -----------------------------

# Helper to build SAFE module-score column names you already use
rf  <- function(sub) module_field_name("Reservoir_Favoring", sub)
mf  <- function(fam, sub_pretty) module_field_name(fam, sub_pretty)  # keep existing

# RPI positives (Reservoir_Favoring submodules)
pos_components <- c(
  rf("NFkB_activity"),
  rf("STAT_reactivation_axis"),
  rf("MYC_targets"),
  rf("Cell_cycle_E2F_G2M")
)

# Glycolysis / OXPHOS are also under Reservoir_Favoring now
gly_ox <- c(
  rf("Glycolysis"),
  rf("OXPHOS")
)

# RPI negatives
neg_tgfb_target <- mf("TGFb", "TGFβ transcriptional targets")
neg_il10_regs   <- c(
  mf("IL10_STAT3", "STAT3-regulatory/Treg-overlap"),
  mf("IL10_STAT3", "IL10RA/IL10RB–JAK1/TYK2–STAT3 module")
)
neg_other <- c(
  rf("Apoptosis"),
  rf("Chromatin_repression")
)

# Build per-cell derived pieces
md <- seu@meta.data
stopifnot(all(c(group_col, cluster_col) %in% colnames(md)))

# Per-cell means for GlyOx and IL-10-regulatory block
has_glyox <- gly_ox[gly_ox %in% colnames(md)]
if (length(has_glyox) > 0) {
  md$MS_GlyOx_mean <- rowMeans(md[, has_glyox, drop = FALSE], na.rm = TRUE)
}

has_il10_regs <- neg_il10_regs[neg_il10_regs %in% colnames(md)]
if (length(has_il10_regs) > 0) {
  md$MS_IL10reg_mean <- rowMeans(md[, has_il10_regs, drop = FALSE], na.rm = TRUE)
}

# z-scale within cluster for all needed columns
need_cols_cell <- unique(c(pos_components, "MS_GlyOx_mean", neg_tgfb_target, "MS_IL10reg_mean", neg_other))
need_cols_cell <- need_cols_cell[need_cols_cell %in% colnames(md)]
md <- z_within_cluster_cells(md, need_cols_cell, cluster_col)

# Compute per-cell RPI = z(pos + GlyOx) - z(neg)
z_pos_cols <- paste0("z_", c(pos_components, "MS_GlyOx_mean"))
z_neg_cols <- paste0("z_", c(neg_tgfb_target, "MS_IL10reg_mean", neg_other))
z_pos_cols <- z_pos_cols[z_pos_cols %in% colnames(md)]
z_neg_cols <- z_neg_cols[z_neg_cols %in% colnames(md)]

md$MS_RPI <- NA_real_
if (length(z_pos_cols) > 0 && length(z_neg_cols) > 0) {
  zpos <- if (length(z_pos_cols) == 1) md[[z_pos_cols]] else rowMeans(md[, z_pos_cols, drop = FALSE], na.rm = TRUE)
  zneg <- if (length(z_neg_cols) == 1) md[[z_neg_cols]] else rowMeans(md[, z_neg_cols, drop = FALSE], na.rm = TRUE)
  md$MS_RPI <- zpos - zneg
}


# Per-cell MAI = z(pro-autophagy + MAPK noncanonical) - z(mTORC)
mapk_noncanon <- mf("TGFb", "TGFβ noncanonical (MAPK/PI3K)")
auto_pro_cols <- c("MS_Autophagy_Initiation","MS_Autophagy_Conjugation","MS_Autophagy_Fusion",
                   "MS_Autophagy_Lysosome","MS_Autophagy_CLEAR_TFs","MS_Autophagy_Cargo",
                   "MS_Autophagy_Mitophagy","MS_Autophagy_Xenophagy","MS_Autophagy_AMPK",
                   mapk_noncanon)
auto_anti_col <- "MS_Autophagy_mTORC"

# z-scale inputs for MAI (if present)
need_mai <- unique(c(auto_pro_cols, auto_anti_col))
md <- z_within_cluster_cells(md, need_mai, cluster_col)

z_auto_pro <- paste0("z_", auto_pro_cols)
z_auto_pro <- z_auto_pro[z_auto_pro %in% colnames(md)]
z_auto_anti <- paste0("z_", auto_anti_col)

md$MS_MAI <- NA_real_
if (length(z_auto_pro) > 0 && z_auto_anti %in% colnames(md)) {
  zpro <- if (length(z_auto_pro) == 1) md[[z_auto_pro]] else rowMeans(md[, z_auto_pro, drop = FALSE], na.rm = TRUE)
  zanti <- md[[z_auto_anti]]
  md$MS_MAI <- zpro - zanti
}

# Save back and plot per-cell RPI/MAI
seu@meta.data <- md
plot_module(seu, "MS_RPI", composite_dir, "Reactivation Propensity Index")
plot_module(seu, "MS_MAI", composite_dir, "mTOR–MAPK Autophagy Index")

# -----------------------------
# SAMPLE LEVEL (safe version)
# -----------------------------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

sample_col  <- "orig.ident"
group_col   <- if (exists("group_col")) group_col else "IGRA_status"
cluster_col <- if (exists("cluster_col")) cluster_col else "seurat_clusters"

# Output
sample_out_dir <- file.path(composite_dir, "Sample_Level")
dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Pull meta.data and sanity checks
# -----------------------------
md <- seu@meta.data
need_cols <- c(sample_col, group_col, cluster_col, "MS_RPI", "MS_MAI")
stopifnot(all(need_cols %in% colnames(md)))

# Ensure types are friendly
md[[sample_col]]  <- as.character(md[[sample_col]])
md[[group_col]]   <- as.character(md[[group_col]])
md[[cluster_col]] <- as.character(md[[cluster_col]])

# Helper: majority Group per sample (ties broken deterministically)
sample_group <- md %>%
  filter(!is.na(.data[[group_col]]), .data[[group_col]] != "") %>%
  count(.data[[sample_col]], .data[[group_col]], name = "n") %>%
  group_by(.data[[sample_col]]) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(all_of(sample_col), all_of(group_col))

# -----------------------------
# Summary helpers
# -----------------------------
summarise_sample_metric <- function(df, metric) {
  stopifnot(metric %in% colnames(df))
  out <- df %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      Mean   = mean(.data[[metric]], na.rm = TRUE),
      SD     = sd(.data[[metric]],   na.rm = TRUE),
      NCells = dplyr::n(),
      .groups = "drop"
    ) %>%
    left_join(sample_group, by = setNames(sample_col, sample_col)) %>%
    rename(Sample = !!sample_col, Group = !!group_col) %>%
    mutate(
      Group = ifelse(is.na(Group) | Group == "", "Unknown", Group)
    ) %>%
    # drop rows with non-finite means (boxplot can't handle them)
    filter(is.finite(Mean))
  out
}

summarise_sample_cluster_metric <- function(df, metric) {
  stopifnot(metric %in% colnames(df))
  out <- df %>%
    group_by(.data[[sample_col]], .data[[cluster_col]]) %>%
    summarise(
      Mean   = mean(.data[[metric]], na.rm = TRUE),
      SD     = sd(.data[[metric]],   na.rm = TRUE),
      NCells = dplyr::n(),
      .groups = "drop"
    ) %>%
    left_join(sample_group, by = setNames(sample_col, sample_col)) %>%
    rename(Sample = !!sample_col, Group = !!group_col, Cluster = !!cluster_col) %>%
    mutate(
      Group = ifelse(is.na(Group) | Group == "", "Unknown", Group)
    ) %>%
    filter(is.finite(Mean))
  out
}

# -----------------------------
# Compute summaries
# -----------------------------
rpi_sample_overall <- summarise_sample_metric(md, "MS_RPI")
mai_sample_overall <- summarise_sample_metric(md, "MS_MAI")
rpi_sample_by_cluster <- summarise_sample_cluster_metric(md, "MS_RPI")
mai_sample_by_cluster <- summarise_sample_cluster_metric(md, "MS_MAI")

# Save CSVs
write.csv(rpi_sample_overall,  file.path(sample_out_dir, "RPI_per_sample_overall.csv"),  row.names = FALSE)
write.csv(mai_sample_overall,  file.path(sample_out_dir, "MAI_per_sample_overall.csv"),  row.names = FALSE)
write.csv(rpi_sample_by_cluster, file.path(sample_out_dir, "RPI_per_sample_by_cluster.csv"), row.names = FALSE)
write.csv(mai_sample_by_cluster, file.path(sample_out_dir, "MAI_per_sample_by_cluster.csv"), row.names = FALSE)

# -----------------------------
# Plot (safe) helpers
# -----------------------------
plot_sample_overall <- function(tab, title_y, fname) {
  dat <- tab %>% filter(is.finite(Mean))
  if (nrow(dat) == 0) {
    message("[", title_y, "] No finite data for overall plot; skipping.")
    return(invisible(NULL))
  }
  # If only one group, avoid empty boxplot stats
  have_box <- length(unique(dat$Group)) > 0 && any(!is.na(dat$Mean))
  p <- ggplot(dat, aes(x = Group, y = Mean)) +
    { if (have_box) geom_boxplot(outlier.shape = NA, width = 0.6) } +
    geom_jitter(aes(size = NCells, color = Sample),
                width = 0.15, alpha = 0.85, show.legend = TRUE) +
    scale_size_continuous(name = "Cells per sample") +
    labs(title = paste0(title_y, " — Per Sample (overall)"),
         x = NULL, y = paste0("Mean ", title_y, " per sample")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  ggsave(file.path(sample_out_dir, fname), p, width = 9, height = 6.5, dpi = 350)
}

plot_sample_by_cluster <- function(tab, title_y, fname, min_pts_per_cluster = 3) {
  dat0 <- tab %>% filter(is.finite(Mean))
  if (nrow(dat0) == 0) {
    message("[", title_y, "] No finite data for by-cluster plot; skipping.")
    return(invisible(NULL))
  }
  # Keep clusters with enough points to form a boxplot/jitter panel
  keep_clusters <- dat0 %>%
    count(Cluster, name = "n_pts") %>%
    filter(n_pts >= min_pts_per_cluster) %>%
    pull(Cluster)
  dat <- dat0 %>% filter(Cluster %in% keep_clusters)
  if (nrow(dat) == 0) {
    message("[", title_y, "] All clusters filtered due to low n; skipping.")
    return(invisible(NULL))
  }
  p <- ggplot(dat, aes(x = Group, y = Mean)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(size = NCells, color = Sample),
                width = 0.15, alpha = 0.85, show.legend = FALSE) +
    facet_wrap(~ Cluster, scales = "free_y") +
    labs(title = paste0(title_y, " — Per Sample by Cluster"),
         x = NULL, y = paste0("Mean ", title_y, " per sample (cluster)")) +
    theme_bw(base_size = 11)
  ggsave(file.path(sample_out_dir, fname), p, width = 12, height = 8.5, dpi = 350)
}

# -----------------------------
# Make & save plots (safe)
# -----------------------------
plot_sample_overall(rpi_sample_overall, "RPI", "RPI_per_sample_overall.png")
plot_sample_overall(mai_sample_overall, "MAI", "MAI_per_sample_overall.png")
plot_sample_by_cluster(rpi_sample_by_cluster, "RPI", "RPI_per_sample_by_cluster.png")
plot_sample_by_cluster(mai_sample_by_cluster, "MAI", "MAI_per_sample_by_cluster.png")

message("Saved sample-level CSVs and plots to: ", sample_out_dir)

# -----------------------------
# Sample-level IL-10 vs TGFβ correlation (one dot per sample), colored by IGRA
# -----------------------------
sample_df <- seu@meta.data %>%
  dplyr::filter(!is.na(orig.ident), !is.na(.data[[group_col]])) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    IL10_comp = mean(MS_IL10_composite, na.rm = TRUE),
    TGFb_comp = mean(MS_TGFb_composite,  na.rm = TRUE),
    IGRA      = dplyr::first(.data[[group_col]]),
    .groups   = "drop"
  ) %>%
  tidyr::drop_na(IL10_comp, TGFb_comp)

# scatter (sample-level)
p_scatter <- ggplot(sample_df, aes(x = TGFb_comp, y = IL10_comp, color = IGRA)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  scale_color_manual(
    values = c("Positive" = "#FF918A", "Negative" = "#A3F8A9"),
    labels = c("Negative" = "IGRA Negative", "Positive" = "IGRA Positive"),
    name   = "IGRA Status"
  ) +
  labs(
    title = "IL-10 composite vs TGFβ composite (sample means)",
    x = "TGFβ composite (sample mean)",
    y = "IL-10 composite (sample mean)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(composite_dir, "IL10_vs_TGFb_composites_sample_scatter.png"),
       p_scatter, width = 10, height = 7, dpi = 300, bg = "white")
# -----------------------------
# Sample-level IL-10 vs TGFβ scatter split by cell type (facet by cluster)
# -----------------------------
sample_by_cluster_df <- seu@meta.data %>%
  dplyr::filter(!is.na(orig.ident), !is.na(.data[[group_col]]), !is.na(.data[[cluster_col]])) %>%
  dplyr::group_by(orig.ident, .data[[cluster_col]]) %>%
  dplyr::summarise(
    IL10_comp = mean(MS_IL10_composite, na.rm = TRUE),
    TGFb_comp = mean(MS_TGFb_composite,  na.rm = TRUE),
    IGRA      = dplyr::first(.data[[group_col]]),
    .groups   = "drop"
  ) %>%
  tidyr::drop_na(IL10_comp, TGFb_comp)

p_scatter_bycell <- ggplot(sample_by_cluster_df,
                           aes(x = TGFb_comp, y = IL10_comp, color = IGRA)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "black") +
  scale_color_manual(
    values = c("Positive" = "#FF918A", "Negative" = "#A3F8A9"),
    labels = c("Negative" = "IGRA Negative", "Positive" = "IGRA Positive"),
    name   = "IGRA Status"
  ) +
  facet_wrap(stats::as.formula(paste0("~", cluster_col)), scales = "free") +
  labs(
    title = "IL-10 composite vs TGFβ composite (sample means), split by cell type",
    x = "TGFβ composite (sample mean)",
    y = "IL-10 composite (sample mean)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(composite_dir, "IL10_vs_TGFb_composites_sample_scatter_byCellType.png"),
  p_scatter_bycell, width = 22, height = 16, dpi = 300, bg = "white"
)

# -----------------------------
# Write module coverage & stats
# -----------------------------
# Coverage table
write_tsv(coverage, file.path(out_dir, "Module_Gene_Coverage.tsv"))

# Stats table (Module x Cluster with p-value and HL shift)
stats_df <- dplyr::bind_rows(module_stats)
stats_csv <- file.path(out_dir, "Module_Stats_byCluster.csv")
readr::write_csv(stats_df, stats_csv)

# -----------------------------
# Reviewer-friendly outputs
# -----------------------------
coverage %>%
  dplyr::mutate(Module_Reviewer = pretty_module_name(Module)) %>%
  readr::write_csv(file.path(out_dir, "Module_Gene_Coverage_pretty.csv"))

stats_df %>%
  dplyr::mutate(Module_Reviewer = pretty_module_name(Module)) %>%
  readr::write_csv(file.path(out_dir, "Module_Stats_byCluster_pretty.csv"))

# ===========================
# Overall Summary (per family)
# ===========================
stats_df$Module
plot_module_summary <- function(stats_df, out_dir) {
  if (is.null(stats_df) || nrow(stats_df) == 0) return(invisible(NULL))
  
  df <- stats_df %>%
    dplyr::filter(!is.na(PValue)) %>%
    dplyr::mutate(
      # Use the pretty label if available; else build one dynamically
      Module_Pretty = dplyr::coalesce(
        if ("Module_Pretty" %in% names(.)) .data$Module_Pretty else NULL,
        pretty_module_name(.data$Module)
      ),
      # Gracefully handle absence of Family column
      Family = if ("Family" %in% names(.)) .data$Family else sub("^([^_]+).*", "\\1", .data$Module),
      Direction = ifelse(HL_Shift >= 0, "Positive", "Negative"),
      SizeMag   = pmin(abs(HL_Shift), 2),
      LogP      = -log10(PValue),
      Signif    = dplyr::case_when(
        PValue <= 1e-4 ~ "****",
        PValue <= 1e-3 ~ "***",
        PValue <= 1e-2 ~ "**",
        PValue <= 5e-2 ~ "*",library(dplyr)
        library(ggplot2)
        library(ggpubr)
        
        # Load saved CSV
        rpi_sample_overall <- read.csv(file.path(sample_out_dir, "RPI_per_sample_overall.csv"),
                                       stringsAsFactors = FALSE)
        
        # Group colours
        group_colours <- c("pos" = "#FF918A", "neg" = "#A3F8A9")  # adjust to match your exact Group labels
        
        # Check what your group labels actually are:
        # unique(rpi_sample_overall$Group)
        # Then update the names above accordingly e.g. c("IGRA+" = "#FF918A", "IGRA-" = "#A3F8A9")
        
        p_rpi_overall <- ggplot(rpi_sample_overall, aes(x = Group, y = Mean, fill = Group)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
          geom_jitter(width = 0.15, size = 2.5, colour = "black") +  # all dots black
          scale_fill_manual(values = group_colours) +
          stat_compare_means(method = "wilcox.test", label = "p.format", 
                             label.x = 1.5, vjust = -0.5) +
          labs(
            title = "RPI per Sample (Overall)",
            x     = NULL,
            y     = "Mean MS_RPI"
          ) +
          theme_classic(base_size = 14) +
          theme(legend.position = "none")   # remove legend
        
        ggsave(
          filename = file.path(sample_out_dir, "RPI_per_sample_overall_replot.png"),
          plot     = p_rpi_overall,
          dpi      = 500,
          width    = 5,
          height   = 5
        )
        TRUE ~ ""
      )
    )
  
  
  # keep a stable module order within each family (by pretty label)
  df <- df %>%
    dplyr::group_by(Family) %>%
    dplyr::mutate(Module_Pretty = factor(Module_Pretty, levels = sort(unique(Module_Pretty), decreasing = TRUE))) %>%
    dplyr::ungroup()
  
  df$Cluster <- factor(df$Cluster, levels = sort(unique(df$Cluster)))
  fill_vals  <- c("Positive" = "#FF918A", "Negative" = "#A3F8A9")
  
  families   <- unique(df$Family)
  for (fam in families) {
    subdf  <- df[df$Family == fam, , drop = FALSE]
    n_rows <- length(unique(subdf$Module_Pretty))
    
    # --- Dynamic sizing (inches) ---
    row_height_in <- 0.45
    baseline_in   <- 2.5
    min_h_in      <- 6
    max_h_in      <- 26
    fig_h         <- max(min_h_in, min(max_h_in, baseline_in + row_height_in * n_rows))
    fig_w         <- 14
    
    # --- Font/marker tuning by density ---
    base_size <- max(12, min(16, 16 - 0.08 * (n_rows - 15)))
    size_min  <- 3.0
    size_max  <- max(6.5, min(9.5, 9.5 - 0.07 * (n_rows - 15)))
    star_size <- 3.8
    star_nudge   <- 0.16
    top_expand   <- 0.22
    bottom_expand<- 0.12
    
    if (n_rows <= 3) {
      size_max     <- 7.2; star_nudge <- 0.18; top_expand <- 0.30; bottom_expand <- 0.22
    } else if (n_rows == 4) {
      size_max     <- 7.8; star_nudge <- 0.18; top_expand <- 0.28; bottom_expand <- 0.20
    } else if (n_rows == 11) {
      size_max     <- min(size_max, 8.8); star_nudge <- 0.27; star_size <- 3.6
      top_expand   <- 0.16; bottom_expand <- 0.10
    } else if (n_rows <= 6) {
      star_nudge   <- 0.16; top_expand <- 0.30; bottom_expand <- 0.20
    } else if (n_rows <= 10) {
      star_nudge   <- 0.18; top_expand <- 0.22; bottom_expand <- 0.12
    } else if (n_rows <= 18) {
      star_nudge   <- 0.22; top_expand <- 0.15; bottom_expand <- 0.08
    } else {
      star_nudge   <- 0.25; top_expand <- 0.12; bottom_expand <- 0.06
    }
    
    p <- ggplot2::ggplot(subdf, ggplot2::aes(x = Cluster, y = Module_Pretty)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = Direction, size = SizeMag),
        shape = 21, color = "black", stroke = 0.6, alpha = 1
      ) +
      ggplot2::geom_text(
        data = subset(subdf, Signif != ""),
        ggplot2::aes(label = Signif),
        color = "red", fontface = "bold",
        size = star_size, vjust = 0,
        position = ggplot2::position_nudge(y = star_nudge)
      ) +
      ggplot2::scale_fill_manual(values = fill_vals, name = "Direction") +
      ggplot2::scale_size_continuous(range = c(size_min, size_max), name = "|HL Shift|") +
      ggplot2::labs(
        title = paste0(fam, " Module Summary: IGRA+ vs IGRA−"),
        x = "Cluster", y = "Module"
      ) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(bottom_expand, top_expand))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0.5),
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y  = ggplot2::element_text(color = "black"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.position    = "top",
        legend.box         = "horizontal",
        legend.title       = ggplot2::element_text(face = "bold"),
        legend.key.height  = grid::unit(10, "pt"),
        legend.key.width   = grid::unit(24, "pt"),
        legend.box.margin  = ggplot2::margin(t = 4, r = 4, b = 4, l = 4),
        legend.background  = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.4),
        plot.margin        = ggplot2::margin(10, 16, 14, 10)
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          override.aes = list(size = 6, stroke = 0.6),
          order = 1, title.position = "top", label.position = "right"
        ),
        size = ggplot2::guide_legend(
          order = 2, title.position = "top", label.position = "right"
        )
      )
    
    ggplot2::ggsave(
      filename  = file.path(out_dir, paste0("Module_Summary_", safe_id(fam), ".png")),
      plot      = p, width = fig_w, height = fig_h, units = "in", dpi = 300, bg = "white",
      limitsize = FALSE
    )
  }
}

# Run it
plot_module_summary(stats_df, out_dir)


############################# Post Annotation VLN Plots ###################
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers"
annotation_dir <- file.path(base_dir, "Annotation")
vlnplot_dir <- file.path(annotation_dir, "Post-Annotation/VlnPlots")
vlnplot_dir_2 <- file.path(annotation_dir, "Post-Annotation/VlnPlots_Extra")

dir.create(vlnplot_dir, showWarnings = FALSE)
dir.create(vlnplot_dir_2, showWarnings = FALSE)

igra_order <- c(
  "CD4+ naïve","CD4+ Activated","CD4+ Tfh","CD4+ Treg",
  "CD8+ TCM","CD8+ EM","CD8+ TEMRA",
  "dnT",
  "NKCD56bright","NK Cytotoxic",
  "cDC2","moDC precursor","moDC","DC-like APCs",
  "CD14+ Monocytes","CD16+ NCM",
  "Naïve B","Memory B","Intermediate B","GC-like B"
)
seu$IGRA_Annotation <- factor(seu$IGRA_Annotation, levels = igra_order)



# ---------------------------- #
# Marker List
# ---------------------------- #
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

rna.features.2 <- c(
  'TRAC','TRBC2','IL7R',
  'FOS','FOSB','JUN','JUNB','JUND','EGR1','EGR2','EGR3','NR4A2','NR4A3',
  'DUSP1','DUSP2','DUSP4','DUSP5','IER2','IER3','TNFAIP3',
  'HSP90AA1','HSPA1A','HSPA1B',
  'TNFRSF4','TNFRSF9','TNFRSF18',
  'CCR7','LEF1','KLF2','KLF3',
  'CCL4','CXCR4','ITGA4','ITGB1','IL12RB2','CCR6','IL23R','TOX2','IL21',
  'GZMB','CTSW','FGFBP2',
  'TOP2A','HMGB2','TYMS','PCNA','STMN1',
  'TUBB','UBE2C','CENPF',
  'S100A4'
)
for (gene in rna.features) {
  vln_file <- file.path(vlnplot_dir, paste0("Vln_", gene, ".png"))
  
  # Skip if plot already exists
  if (!file.exists(vln_file) && gene %in% rownames(seu)) {
    p_vln <- VlnPlot2(
      seu,
      features = gene,
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) +
      theme(plot.margin = unit(c(1, 1, 1, 2), "cm")) # Increase left margin (4th value)
    
    ggsave(
      filename = vln_file,
      plot = p_vln,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

for (gene in rna.features.2) {
  vln_file <- file.path(vlnplot_dir_2, paste0("Vln_", gene, ".png"))
  
  # Skip if plot already exists
  if (!file.exists(vln_file) && gene %in% rownames(seu)) {
    p_vln <- VlnPlot2(
      seu,
      features = gene,
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) +
      theme(plot.margin = unit(c(1, 1, 1, 2), "cm")) # Increase left margin (4th value)
    
    ggsave(
      filename = vln_file,
      plot = p_vln,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

###############
library(dplyr)
library(ggplot2)
library(ggpubr)

sample_out_dir <- '/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/Module_Scoring/Composites/Sample_Level'
# Load saved CSV
rpi_sample_overall <- read.csv(file.path(sample_out_dir, "RPI_per_sample_overall.csv"),
                               stringsAsFactors = FALSE)

# Group colours
group_colours <- c("Positive" = "#FF918A", "Negative" = "#A3F8A9")  # adjust to match your exact Group labels

# Check what your group labels actually are:
# unique(rpi_sample_overall$Group)
# Then update the names above accordingly e.g. c("IGRA+" = "#FF918A", "IGRA-" = "#A3F8A9")

p_rpi_overall <- ggplot(rpi_sample_overall, aes(x = Group, y = Mean, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.5, colour = "black") +  # all dots black
  scale_fill_manual(values = group_colours) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.x = 1.5, vjust = -0.5) +
  labs(
    title = "RPI per Sample (Overall)",
    x     = NULL,
    y     = "Mean MS_RPI"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")   # remove legend

ggsave(
  filename = file.path(sample_out_dir, "RPI_per_sample_overall_replot.png"),
  plot     = p_rpi_overall,
  dpi      = 500,
  width    = 5,
  height   = 5
)
