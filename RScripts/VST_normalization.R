# Load required libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(glue)
library(rtracklayer)
options(contrasts = c("contr.treatment", "contr.poly"))

outliers <- c('BDIS05_T1_L', 'BDIS17_T2_L', 'BDIS58_T4_L', 'BDIS75_T5_L',
              'BDIS17_T2_R', 'BDIS37_T3_R', 'BDIS58_T4_R', 'BDIS75_T5_R',
              'BSYL75_T5_R', 'BSYL77_T5_R', 'BSYL80_T5_R',
              'HVUL37_T3_R')

# Load metadata
load("Data/tidy_meta.RData", verbose = TRUE)

# Define species and their corresponding paths
species_info <- tribble(
  ~species_code, ~species_name, ~annotation_file,
  "bdis", "Brachypodium distachyon", "Bdistachyon_556_v3.2.gene_exons.gtf",
  "bsyl", "Brachypodium sylvaticum", "BsylvaticumAin_1_721_v2.1.gene_exons.gtf",
  "hvul", "Hordeum vulgare", "HvulgareMorex_702_V3.gene.gtf",
  "hjub", "Hordeum jubatum", "Hordeum_jubatum_BCC2055.pgsb.r1.Mar2024.gtf"
)

# Function to process each species
process_species <- function(species_code, species_name, annotation_file) {
  # Get salmon files
  salmon_path <- glue::glue("Data/quantification/salmon-{species_code}/")
  
  sfs <- list.files(
    path = salmon_path,
    pattern = "quant.sf",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Set names based on directory structure
  names(sfs) <- basename(dirname(sfs))
  
  # Remove outliers from salmon files    
  sfs <- sfs[!names(sfs) %in% outliers]
  
  # Import annotation and create transcript map
  anno_path <- glue::glue("Data//annotations/{annotation_file}")
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
  
  # Create DESeq object and process by tissue
  dds <- DESeqDataSetFromTximport(txi, meta, ~ time.point)
  
  
  # Function to process tissue data
  process_tissue <- \(tissue_type) {
    n <- 10  # min number of counts
    k <- 4  # min number of samples
  
    
    # Filter and process tissue-specific data
    dds_tissue <- dds[, meta$tissue == tissue_type] |>
      (\(x) x[rowSums(counts(x) >= n) >= k, ])()
    
    # Calculate size factors specifically for this tissue
    dds_tissue <- estimateSizeFactors(dds_tissue, type = "poscount")
    
    # Run DESeq 
    res_tissue <- DESeq(dds_tissue, fitType = "glmGamPoi", test = "LRT", reduced = ~ 1) 
    
    
    # Get normalized counts and format
    vst(res_tissue, blind = FALSE) |>
      assay() |>
      as_tibble(rownames = "GeneID") |>
      pivot_longer(-GeneID, names_to = "library", values_to = "vst.count") |>
      left_join(meta |> as_tibble(rownames = "library"), by = join_by("library"))
  }
  
  # Process both tissues and combine results
  bind_rows(
    leaf = process_tissue("leaf"),
    root = process_tissue("root"),
    .id = "tissue"
  )
}

# Process all species
results <- species_info |>
  pmap(process_species) |>
  set_names(species_info$species_code)

# Re-convert to wide format:
# BDIS
results$bdis |>
  filter(tissue == "leaf") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'BDIS_L_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

results$bdis |>
  filter(tissue == "root") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'BDIS_R_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)


# BSYL
results$bsyl |>
  filter(tissue == "leaf") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'BSYL_L_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

results$bsyl |>
  filter(tissue == "root") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'BSYL_R_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# HVUL
results$hvul |>
  filter(tissue == "leaf") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'HVUL_L_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

results$hvul |>
  filter(tissue == "root") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'HVUL_R_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# HJUB
results$hjub |>
  filter(tissue == "leaf") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'HJUB_L_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

results$hjub |>
  filter(tissue == "root") |>
  select(GeneID, library, vst.count) |>
  pivot_wider(names_from = library, values_from = vst.count) |>
  pivot_longer(cols = -GeneID, names_to = "TimePoint", values_to = "Value") %>%
  mutate(Time = str_extract(TimePoint, "T[1-5]")) %>% 
  group_by(GeneID, Time) %>%
  summarize(MeanValue = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = MeanValue) |>
  write.table(file = 'HJUB_R_vst-normalized.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
