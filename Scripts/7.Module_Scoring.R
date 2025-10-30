
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(viridis)
library(stringr)
library(scales)

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
  
  if (all(c("Negative", "Positive") %in% unique(md[[group_col]]))) {
    md[[group_col]] <- factor(md[[group_col]], levels = c("Negative", "Positive"))
  }
  
  gg <- ggplot(md, aes_string(x = group_col, y = module_field, fill = group_col)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.18, size = 0.4, alpha = 0.25) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.x.npc = "center",
      label.y.npc = 0.9,      # lowered from the facet strip
      size = 7,
      color = "red"           # red significance stars
    ) +
    facet_wrap(stats::as.formula(paste0("~", cluster_col)), scales = "free_y") +
    scale_fill_manual(
      values = c("Positive" = "#FF918A", "Negative" = "#A3F8A9"),
      labels = c("Negative" = "IGRA Negative", "Positive" = "IGRA Positive"),
      name   = "IGRA Status"
    ) +
    labs(
      title = title,
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
    filename = file.path(folder, paste0(module_field, ".png")),
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


############# Overall Summary

plot_module_summary <- function(stats_df, out_dir) {
  if (nrow(stats_df) == 0) return(invisible(NULL))
  
  df <- stats_df %>%
    dplyr::filter(!is.na(PValue)) %>%
    dplyr::mutate(
      Family    = sub("^([^_]+).*", "\\1", Module),
      Direction = ifelse(HL_Shift >= 0, "Positive", "Negative"),
      SizeMag   = pmin(abs(HL_Shift), 2),
      LogP      = -log10(PValue),
      Signif    = dplyr::case_when(
        PValue <= 1e-4 ~ "****",
        PValue <= 1e-3 ~ "***",
        PValue <= 1e-2 ~ "**",
        PValue <= 5e-2 ~ "*",
        TRUE ~ ""
      )
    )
  
  df$Cluster <- factor(df$Cluster, levels = sort(unique(df$Cluster)))
  df$Module  <- factor(df$Module,  levels = rev(unique(df$Module)))
  fill_vals  <- c("Positive" = "#FF918A", "Negative" = "#A3F8A9")
  families   <- unique(df$Family)
  
  for (fam in families) {
    subdf  <- df[df$Family == fam, , drop = FALSE]
    n_rows <- length(unique(subdf$Module))
    
    # --- Dynamic sizing (inches) ---
    row_height_in <- 0.45
    baseline_in   <- 2.5
    min_h_in      <- 6
    max_h_in      <- 26
    fig_h         <- max(min_h_in, min(max_h_in, baseline_in + row_height_in * n_rows))
    fig_w         <- 14
    
    # --- Defaults (good for 5–10 rows) ---
    base_size <- max(12, min(16, 16 - 0.08 * (n_rows - 15)))
    size_min  <- 3.0
    size_max  <- max(6.5, min(9.5, 9.5 - 0.07 * (n_rows - 15)))
    star_size <- 3.8
    star_nudge   <- 0.16
    top_expand   <- 0.22
    bottom_expand<- 0.12
    
    # --- Special handling for your row counts ---
    if (n_rows <= 3) {
      # Very few rows → big dots; slightly limit max size & raise stars a touch
      size_max     <- 7.2
      star_nudge   <- 0.18
      top_expand   <- 0.30
      bottom_expand<- 0.22
    } else if (n_rows == 4) {
      size_max     <- 7.8
      star_nudge   <- 0.18
      top_expand   <- 0.28
      bottom_expand<- 0.20
    } else if (n_rows == 11) {
      # Dense but not extreme → lift stars slightly & give more top space
      size_max     <- min(size_max, 8.8)  # keep dots substantial but not overwhelming
      star_nudge   <- 0.27                # a bit higher to clear big dots
      star_size    <- 3.6                 # tiny reduction for extra clearance
      top_expand   <- 0.16
      bottom_expand<- 0.10
    } else if (n_rows <= 6) {
      star_nudge   <- 0.16
      top_expand   <- 0.30
      bottom_expand<- 0.20
    } else if (n_rows <= 10) {
      star_nudge   <- 0.18
      top_expand   <- 0.22
      bottom_expand<- 0.12
    } else if (n_rows <= 18) {
      star_nudge   <- 0.22
      top_expand   <- 0.15
      bottom_expand<- 0.08
    } else {
      star_nudge   <- 0.25
      top_expand   <- 0.12
      bottom_expand<- 0.06
    }
    
    p <- ggplot2::ggplot(subdf, ggplot2::aes(x = Cluster, y = Module)) +
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
      filename  = file.path(out_dir, paste0("Module_Summary_", fam, ".png")),
      plot      = p,
      width     = fig_w,
      height    = fig_h,
      units     = "in",
      dpi       = 300,
      bg        = "white",
      limitsize = FALSE
    )
    
    message(
      "✅ Saved: ", file.path(out_dir, paste0("Module_Summary_", fam, ".png")),
      " (", sprintf("%.1f", fig_w), "×", sprintf("%.1f", fig_h), " in; rows=", n_rows,
      "; star_nudge=", star_nudge, "; star_size=", star_size,
      "; size_max=", sprintf("%.1f", size_max),
      "; expandTop=", top_expand, "; expandBottom=", bottom_expand, ")"
    )
  }
}

plot_module_summary(stats_df, out_dir)

###########################3
## ---- Prep a test subset ----
if (!"Family" %in% names(stats_df)) {
  stats_df$Family <- sub("^([^_]+).*", "\\1", stats_df$Module)
}
test_fam <- unique(stats_df$Family)[1]

subdf <- subset(stats_df, Family == test_fam & !is.na(PValue))
subdf$Direction <- ifelse(subdf$HL_Shift >= 0, "Positive", "Negative")
subdf$SizeMag   <- pmin(abs(subdf$HL_Shift), 2)   # cap for nicer scaling

# Significance stars (common cutoffs); blank if not significant
subdf$Signif <- ifelse(subdf$PValue <= 1e-4, "****",
                       ifelse(subdf$PValue <= 1e-3, "***",
                              ifelse(subdf$PValue <= 1e-2, "**",
                                     ifelse(subdf$PValue <= 5e-2, "*", "")
                              )))

# Factor ordering for tidy axes
subdf$Cluster <- factor(subdf$Cluster, levels = sort(unique(subdf$Cluster)))
subdf$Module  <- factor(subdf$Module,  levels = rev(unique(subdf$Module)))

fill_vals <- c("Positive" = "#FF918A", "Negative" = "#A3F8A9")

n_mods <- length(unique(subdf$Module))
fig_h  <- max(6, min(20, n_mods * 0.5))

## ---- Build test plot ----
p <- ggplot(subdf, aes(x = Cluster, y = Module)) +
  # Dots: fully opaque, outlined for contrast
  geom_point(
    aes(fill = Direction, size = SizeMag),
    shape = 21, color = "black", stroke = 0.5, alpha = 1
  ) +
  # Stars: moved clearly above the dot using a Y nudge on discrete axis
  geom_text(
    data = subset(subdf, Signif != ""),
    aes(label = Signif),
    color = "red", fontface = "bold",
    size = 4, vjust = 0,                    # anchor bottom of text
    position = position_nudge(y = 0.18)       # push above the dot
  ) +
  scale_fill_manual(values = fill_vals, name = "Direction") +
  scale_size_continuous(
    range = c(3.2, 9.5),
    name  = "|HL Shift|"
  ) +
  labs(
    title = paste0(test_fam, " Module Summary: IGRA+ vs IGRA−"),
    x = "Cluster", y = "Module"
  ) +
  # Add top/bottom breathing room so nudged stars aren’t clipped
  scale_y_discrete(expand = expansion(mult = c(0.04, 0.12))) +
  coord_cartesian(clip = "off") +
  # Publication theme
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    
    # Legend polish
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.height = unit(10, "pt"),
    legend.key.width  = unit(24, "pt"),
    legend.box.margin = margin(t = 4, r = 4, b = 4, l = 4),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 6, stroke = 0.6),
      order = 1, title.position = "top", label.position = "right"
    ),
    size = guide_legend(
      order = 2, title.position = "top", label.position = "right"
    )
  )

p

## ---- Save test figure ----
ggsave(
  filename = file.path(out_dir, paste0("Module_Summary_TEST_", test_fam, ".png")),
  plot = p, width = 14, height = fig_h, dpi = 300, bg = "white"
)

