library(Seurat)


############################ Read in Unfiltered Parse Objects and Create Merged Seurat ##################################################

# Base path containing sample folders like "sample_1273088", "sample_1273505", etc.
base_dir <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/Ben_PARSE_TB_Data/GEX_Output/combined"

# List all sample folders
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^sample_", basename(sample_dirs))]

# Initialize list to hold Seurat objects
seurat_list <- list()

# Loop through all sample folders
for (sample_path in sample_dirs) {
  sample_name <- basename(sample_path)
  sample_short <- sub("^sample_", "", sample_name)
  dge_path <- file.path(sample_path, "DGE_unfiltered")
  
  # Skip if required files are missing
  if (!file.exists(file.path(dge_path, "count_matrix.mtx"))) {
    warning(paste("Skipping", sample_name, "- missing DGE_unfiltered data"))
    next
  }
  
  # Read matrix
  mat <- ReadParseBio(dge_path)
  
  # Fix empty gene names
  if (any(rownames(mat) == "")) {
    rownames(mat)[rownames(mat) == ""] <- "unknown"
  }
  
  # Read cell metadata
  cell_meta <- read.csv(file.path(dge_path, "cell_metadata.csv"), row.names = 1)
  
  # Create Seurat object
  seu <- CreateSeuratObject(
    counts = mat,
    meta.data = cell_meta,
    project = sample_name
  )
  
  # Assign sample information
  seu$orig.ident <- sample_name
  seu$sample <- sample_short
  
  # Store
  seurat_list[[sample_short]] <- seu
  
  cat("✅ Loaded", sample_name, "\n")
}

# Merge all Seurat objects
seu_merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = "Parse_Merged"
)

#################################################### Fix Metadata ########################################################################

# Create cleaned IGRA metadata
igra_df <- data.frame(
  sample = c("1273088", "1273531", "6072090", "6056780", "6072210",
             "1273601", "6060100", "6065400", "1273505", "6060150", "6060180", "6060280"),
  IGRA_status = c("Negative", "Negative", "Negative", "Negative", "Negative", "Negative",
                  "Positive", "Positive", "Positive", "Positive", "Positive", "Positive"),
  stringsAsFactors = FALSE
)

# Map IGRA status to cells in merged object
seu_merged$IGRA_status <- igra_df$IGRA_status[match(seu_merged$sample, igra_df$sample)]

# Define output path
save_path <- "/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/saved_R_data/seu_merged.rds"

# Save the object
saveRDS(seu_merged, file = save_path)
