library(dplyr)
library(enrichR)
library(openxlsx)
library(ggplot2)


# -----------------------------
# Enrichr databases
# -----------------------------
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)

tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

# -----------------------------
# Input directories for DGE CSVs
# -----------------------------
input_dirs <- list(
  IGRAPosvsNeg = "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/DGE"
)

# -----------------------------
# Output base directory
# -----------------------------
base_output <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/Pathway_Analysis_EnrichR"

# -----------------------------
# Run Enrichr on significant genes (p_adj < 0.05)
# -----------------------------
run_enrichment <- function(gene_df, label, base_output) {
  gene_list <- rownames(gene_df)
  if (length(gene_list) == 0) {
    message("No genes found for ", label)
    return(NULL)
  }
  
  enrichment <- enrichr(gene_list, databases)
  
  # Create output subfolders
  dir.create(file.path(base_output, "CSVs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_output, "Plots"), recursive = TRUE, showWarnings = FALSE)
  
  # Save Excel results
  excel_out <- file.path(base_output, "CSVs", paste0(label, "_Enrichment.xlsx"))
  wb <- createWorkbook()
  for (db in names(enrichment)) {
    addWorksheet(wb, substr(db, 1, 31))
    writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  saveWorkbook(wb, excel_out, overwrite = TRUE)
  
  # Process top TF results
  top_tf_list <- list()
  top_pathway_list <- list()
  
  for (db_name in names(enrichment)) {
    db_results <- enrichment[[db_name]]
    
    # Rename Combined Score if needed
    if ("Combined Score" %in% colnames(db_results)) {
      db_results <- db_results %>% rename(Combined.Score = `Combined Score`)
    }
    
    if (!"Combined.Score" %in% colnames(db_results)) {
      message("Skipping ", db_name, " for ", label, " — no Combined.Score.")
      next
    }
    
    sig_results <- db_results %>% filter(Adjusted.P.value < 0.05)
    
    if (nrow(sig_results) > 0) {
      top_terms <- sig_results %>%
        arrange(desc(Combined.Score)) %>%
        slice_head(n = 10) %>%
        mutate(Database = db_name)
      
      if (db_name %in% tf_databases) {
        top_tf_list[[db_name]] <- top_terms
      } else if (db_name %in% pathway_databases) {
        top_pathway_list[[db_name]] <- top_terms
      }
    }
  }
  
  # Plot TFs
  tf_df <- bind_rows(top_tf_list)
  if ("Combined.Score" %in% colnames(tf_df) && nrow(tf_df) > 0) {
    tf_df <- tf_df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    
    p_tf <- ggplot(tf_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      scale_y_log10() +  # <-- log scale here
      coord_flip() +
      labs(
        title = paste("Top Transcription Factors -", label),
        x = "TF Term", y = "log10(Combined Score)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, "Plots", paste0(label, "_Transcription_Factors.png")),
      plot = p_tf, width = 12, height = 10, dpi = 300, bg = "white"
    )
  }
  
  # Plot Pathways
  pathway_df <- bind_rows(top_pathway_list)
  if ("Combined.Score" %in% colnames(pathway_df) && nrow(pathway_df) > 0) {
    pathway_df <- pathway_df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    
    p_path <- ggplot(pathway_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      scale_y_log10() +  # <-- log scale here
      coord_flip() +
      labs(
        title = paste("Top Pathways -", label),
        x = "Pathway Term", y = "log10(Combined Score)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, "Plots", paste0(label, "_Pathways.png")),
      plot = p_path, width = 12, height = 10, dpi = 300, bg = "white"
    )
  }
}

# -----------------------------
# Loop over all comparisons and clusters
# -----------------------------
# Loop over all comparisons and clusters
for (comparison in names(input_dirs)) {
  input_dir <- input_dirs[[comparison]]
  comparison_output <- file.path(base_output, comparison)
  
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    message("No CSV files found in: ", input_dir)
    next
  }
  
  for (csv in csv_files) {
    dge <- tryCatch(read.csv(csv, row.names = 1, check.names = FALSE),
                    error = function(e) { message("Read failed: ", csv, " — ", e$message); return(NULL) })
    if (is.null(dge)) next
    
    if (!"p_val_adj" %in% colnames(dge)) {
      message("Skipping (no p_val_adj): ", basename(csv), " | cols: ", paste(colnames(dge), collapse = ", "))
      next
    }
    
    sig_genes <- dge %>% dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05)
    if (nrow(sig_genes) == 0) {
      message("No significant genes in ", basename(csv))
      next
    }
    
    label <- tools::file_path_sans_ext(basename(csv))
    message("Running Enrichr for: ", label)
    run_enrichment(sig_genes, label, comparison_output)
  }
}

######## Specific Plots #######################
library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)


read_tf_enrichr <- function(xlsx_path, tf_databases) {
  shs <- getSheetNames(xlsx_path)
  tf_sheets <- intersect(shs, substr(tf_databases, 1, 31))
  if (!length(tf_sheets)) stop("No TF sheets found in: ", xlsx_path)
  
  bind_rows(lapply(tf_sheets, function(s) {
    df <- readWorkbook(xlsx_path, sheet = s)
    df %>%
      select(Term, Adjusted.P.value, Combined.Score, Genes) %>%
      mutate(
        Database = s,
        Term = as.character(Term),
        Adjusted.P.value = as.numeric(Adjusted.P.value),
        Combined.Score = as.numeric(Combined.Score)
      )
  }))
}

# Escape for literal substring OR allow regex by setting use_regex = TRUE
.build_pattern <- function(terms, use_regex = FALSE) {
  terms <- unique(str_trim(terms))
  if (!use_regex) {
    esc <- function(x) str_replace_all(x, "([\\W])", "\\\\\\1")
    terms <- vapply(terms, esc, "", USE.NAMES = FALSE)
  }
  paste0("(", paste(terms, collapse = "|"), ")")
}

plot_curated_tfs <- function(label,
                             comparison = "IGRAPosvsNeg",
                             curated_terms,
                             base_output,
                             tf_databases,
                             match_mode = c("pattern","exact"),
                             use_regex = FALSE) {
  
  match_mode <- match.arg(match_mode)
  xlsx_path <- file.path(base_output, comparison, "CSVs", paste0(label, "_Enrichment.xlsx"))
  if (!file.exists(xlsx_path)) stop("Excel not found: ", xlsx_path)
  
  tf_all <- read_tf_enrichr(xlsx_path, tf_databases) %>%
    mutate(Term_norm = str_trim(Term))
  
  if (match_mode == "exact") {
    curated_norm <- str_trim(curated_terms)
    tf_cur <- tf_all %>%
      filter(tolower(Term_norm) %in% tolower(curated_norm),
             !is.na(Adjusted.P.value), Adjusted.P.value < 0.05)
    missing <- setdiff(curated_norm, unique(tf_all$Term_norm))
  } else {
    pat <- .build_pattern(curated_terms, use_regex = use_regex)
    tf_cur <- tf_all %>%
      filter(str_detect(Term_norm, regex(pat, ignore_case = TRUE)),
             !is.na(Adjusted.P.value), Adjusted.P.value < 0.05)
    # Terms never matched by substring:
    matched_set <- unique(tf_cur$Term_norm)
    missing <- curated_terms[!vapply(curated_terms, function(t)
      any(str_detect(matched_set, regex(.build_pattern(t, use_regex), ignore_case = TRUE))), logical(1))]
  }
  
  tf_cur <- tf_cur %>%
    distinct(Term, Database, .keep_all = TRUE) %>%
    arrange(desc(Combined.Score))
  
  if (!nrow(tf_cur)) {
    message("No curated TF terms (q<0.05) found for: ", label)
    if (length(missing)) message("Not present/matched: ", paste(missing, collapse = " | "))
    return(invisible(NULL))
  }
  
  tf_cur <- tf_cur %>%
    mutate(Term = factor(Term, levels = unique(Term[order(Combined.Score)])))
  
  p <- ggplot(tf_cur, aes(x = Term, y = Combined.Score, fill = Database)) +
    geom_bar(stat = "identity") +
    scale_y_log10() +
    coord_flip() +
    labs(
      title = paste("Curated Transcription Factors -", label),
      x = "TF Term", y = "log10(Combined Score)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  out_dir <- file.path(base_output, comparison, "Plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outfile_png <- file.path(out_dir, paste0(label, "_TFs_Curated.png"))
  ggsave(outfile_png, plot = p, width = 12, height = 10, dpi = 300, bg = "white")
  
  message("Saved: ", outfile_png)
  if (length(missing)) message("Not present/matched: ", paste(missing, collapse = " | "))
  
  invisible(p)
}

# -----------------------------
# Define your curated term lists
# -----------------------------
curated_CD4_activated <- c(
  "STAT3 mouse",
  "STAT3 20064451 ChIP-Seq CD4+T Mouse",
  "YY1 33199912 ChIP-Seq 293T Human KidneyEmbryo",
  "YY1 21170310 ChIP-Seq MESCs Mouse",
  "CREM 20920259 ChIP-Seq GC1-SPG Mouse",
  "FOXO1 32281255 ChIP-Seq Chondrocytes Human Osteoarthritis",
  "FOXO1 25302145 ChIP-Seq T-LYMPHOCYTE Mouse",
  "BCL3 23251550 ChIP-Seq MUSCLE Mouse",
  "RUNX1 21571218 ChIP-Seq MEGAKARYOCYTES Human",
  "RUNX1 30185409 ChIP-Seq HPC Mouse BoneMarrow Leukemia",
  "RUNX1 22343733 ChIP-Seq Kasumi1 Human PeriphBlood",
  "CHD4 35695185 ChIP-Seq nicBasalRootGanglia Mouse Embryo",
  "SIN3A 21632747 ChIP-Seq MESCs Mouse",
  "SIN3B 21632747 ChIP-Seq MESCs Mouse",
  "SETDB1 19884255 ChIP-Seq MESCs Mouse"
)

curated_CD14_mono <- c(
  "STAT3 20064451 ChIP-Seq CD4+T Mouse",
  "STAT3 23295773 ChIP-Seq U87 Human",
  "STAT3 18555785 Chip-Seq ESCs Mouse",
  "NR3C1 27076634 ChIP-Seq BEAS2B Human Lung Inflammation",
  "VDR 24787735 ChIP-Seq THP-1 Human",
  "SMAD2 18955504 ChIP-ChIP HaCaT Human",
  "SMAD2/3 21741376 ChIP-Seq EPCs Human",
  "SMAD3 18955504 ChIP-ChIP HaCaT Human",
  "SMAD3 30307970 ChIP-Seq HCASMC Human CoronaryArtery",
  "SMAD3 21741376 ChIP-Seq ESCs Human",
  "SMAD3 21741376 ChIP-Seq HESCs Human",
  "SMAD4 21799915 ChIP-Seq A2780 Human",
  "SMAD4 21741376 ChIP-Seq HESCs Human",
  "SMAD7 mouse",
  "NCOR 22465074 ChIP-Seq MACROPHAGES Mouse",
  "NCOR1 26117541 ChIP-Seq K562 Human",
  "SMRT 22465074 ChIP-Seq MACROPHAGES Mouse",
  "SMRT 27268052 Chip-Seq Bcells Human"
)

# -----------------------------
# Make curated plots for the two labels you mentioned
# -----------------------------
label_CD4  <- "DGE_CD4+_Activated_IGRApos_vs_neg"
label_CD14 <- "DGE_CD14+_Monocytes_IGRApos_vs_neg"
comparison_name <- "IGRAPosvsNeg"

# Pattern (substring) matching recommended to catch TRRUST-style entries like "STAT3"
plot_curated_tfs(label_CD4,  comparison_name, curated_CD4_activated, base_output, tf_databases, match_mode = "pattern")
plot_curated_tfs(label_CD14, comparison_name, curated_CD14_mono,     base_output, tf_databases, match_mode = "pattern")

