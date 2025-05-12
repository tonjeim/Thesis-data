# Load required libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(glue)
library(rtracklayer)
options(contrasts = c("contr.treatment", "contr.poly"))

# Load metadata
load("C:/Users/tonje/OneDrive/Data/tidy_meta.RData", verbose = TRUE)

# Define species and their corresponding paths
species_info <- tribble(
  ~species_code, ~species_name, ~annotation_file,
  "bdis", "Brachypodium distachyon", "Bdistachyon_556_v3.2.gene_exons.gtf",
  "bsyl", "Brachypodium sylvaticum", "BsylvaticumAin_1_721_v2.1.gene_exons.gtf",
  "hvul", "Hordeum vulgare", "HvulgareMorex_702_V3.gene.gtf",
  "hjub", "Hordeum jubatum", "Hordeum_jubatum_BCC2055.pgsb.r1.Mar2024.gtf"
)

# Function to process each species and return both tidy data and dds objects
process_species <- function(species_code, species_name, annotation_file) {
  # Get salmon files
  salmon_path <- glue::glue("C:/Users/tonje/OneDrive/Data/quantification/salmon-{species_code}/")
  
  sfs <- list.files(
    path = salmon_path,
    pattern = "quant.sf",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Set names based on directory structure
  names(sfs) <- basename(dirname(sfs))
  
  sfs <- sfs[!names(sfs) %in% excluded_samples]
  
  # Import annotation and create transcript map
  anno_path <- glue::glue("C:/Users/tonje/OneDrive/Data/annotations/{annotation_file}")
  anno <- rtracklayer::import(anno_path)
  
  tx_map <- anno |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "transcript") |>
    select(transcript_id, gene_id)
  
  # Import with tximport
  txi <- tximport(files = sfs, type = "salmon", tx2gene = tx_map)
  
  # Prepare metadata
  meta <- all.meta |>
    filter(species == species_name) |>
    column_to_rownames("library") |>
    as.data.frame(stringsAsFactors = TRUE) |>
    (\(x) x[match(colnames(txi$counts), rownames(x)), ])()
  
  # Create DESeq object
  dds <- DESeqDataSetFromTximport(txi, meta, ~ time.point)
  
  # Function to process tissue data and return both dds and tidy data
  process_tissue <- function(tissue_type) {
    n <- 10  # min number of counts
    k <- 15  # min number of samples
    
    # Filter and process tissue-specific data
    dds_tissue <- dds[, meta$tissue == tissue_type] |>
      (\(x) x[rowSums(counts(x) >= n) >= k, ])()
    
    # Run DESeq
    dds_tissue <- DESeq(dds_tissue, fitType = "glmGamPoi", test = "LRT", reduced = ~ 1)
    
    # Get normalized counts and format for tidy data
    tidy_data <- vst(dds_tissue, blind = FALSE) |>
      assay() |>
      as_tibble(rownames = "GeneID") |>
      pivot_longer(-GeneID, names_to = "library", values_to = "vst.count") |>
      left_join(meta |> as_tibble(rownames = "library"), by = join_by("library"))
    
    # Return both the DESeq object and the tidy data
    return(list(
      dds = dds_tissue,
      tidy_data = tidy_data
    ))
  }
  
  # Process both tissues
  leaf_results <- process_tissue("leaf")
  root_results <- process_tissue("root")
  
  # Combine tidy data
  tidy_combined <- bind_rows(
    leaf = leaf_results$tidy_data,
    root = root_results$tidy_data,
    .id = "tissue"
  )
  
  # Return both dds objects and combined tidy data
  return(list(
    leaf_dds = leaf_results$dds,
    root_dds = root_results$dds,
    tidy_data = tidy_combined
  ))
}

# Process all species and store results
results <- species_info |>
  pmap(process_species) |>
  set_names(species_info$species_code)
dds_obj <- results$bsyl$root_dds

# Extract PCA data
pca_data <- plotPCA(vst(dds_obj, blind = FALSE), intgroup = c("time.point"), returnData = TRUE, ntop = 1500)

# Create the PCA plot with ggplot2
ggplot(pca_data, aes(PC1, PC2, color = 'red')) +
  geom_point(size = 3) +  # Add points
  geom_text(aes(label = rownames(pca_data)), vjust = -1, size = 3) +  # Add sample ID labels
  labs(title = "PCA of vst-normalized data", x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")) +
  theme_minimal()

ggplot(pca_data, aes(PC1, PC2, color = time.point)) +
  geom_point(size = 3) +  # Add points
  geom_text(aes(label = rownames(pca_data)), vjust = -1, size = 3) +  # Add sample ID labels
  labs(title = "PCA of vst-normalized data", x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")) +
  theme_minimal()

