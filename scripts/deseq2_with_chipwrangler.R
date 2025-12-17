#!/usr/bin/env Rscript
### CHANGED: explicit shebang for subprocess safety

deseq_with_chipwrangler <- function(
    counts_file,
    sample_metadata,
    conditions,             # list of 2 strings, e.g. list("TRP_4hr","TRP_0hr")
    condition_col = "condition",
    spike_genomes
) {

  if (length(conditions) != 2) {
    stop("conditions must be a list of **two** condition names.")
  }

  condA <- conditions[[1]]
  condB <- conditions[[2]]
  shrink_coef <- paste0(condition_col, "_", condB, "_vs_", condA)

  message("Shrinkage coefficient will be: ", shrink_coef)

  #### Libraries ####
  suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(apeglm)
    library(stringr)
  })

  ##### Helper for renaming columns using library IDs #####
  clean_colnames_using_metadata <- function(raw_colnames, library_ids) {
    cleaned <- sapply(raw_colnames, function(col) {
      hits <- library_ids[sapply(library_ids, function(id) grepl(id, col, fixed = TRUE))]
      if (length(hits) == 0) return(col)
      if (length(hits) > 1) stop(paste("ERROR: Multiple library IDs matched column:", col))
      hits[[1]]
    })
    return(cleaned)
  }

  ######################## Step 1: Load metadata & counts ########################
  sample_meta <- read.delim(sample_metadata)
  library_ids <- sample_meta$library.ID
  counts_df <- read.delim(counts_file)

  count_cols <- grep("Tag", colnames(counts_df), value = TRUE)
  new_names <- clean_colnames_using_metadata(count_cols, library_ids)
  colnames(counts_df)[match(count_cols, colnames(counts_df))] <- new_names
  cts_counts_df <- counts_df[, new_names, drop = FALSE]

  # Convert to integer
  cts_counts_df[] <- lapply(cts_counts_df, function(x) {
    if (is.numeric(x) && any(x != floor(x), na.rm = TRUE)) as.integer(round(x)) else x
  })

  message("==== COUNTS MATRIX ====")
  print(dim(cts_counts_df))
  print(head(cts_counts_df))

  ######################## Step 2: Construct colData ########################
  samples <- colnames(cts_counts_df)
  meta_sub <- sample_meta[sample_meta$library.ID %in% samples, ]
  meta_sub <- meta_sub[match(samples, meta_sub$library.ID), ]
  coldata <- data.frame(
    row.names = samples,
    biorep = factor(meta_sub$Biorep),
    antibody = factor(meta_sub$IP),
    condition = factor(meta_sub$Condition),
    techrep = factor(meta_sub$TechRep)
  )

  message("==== COLDATA MATRIX ====")
  print(dim(coldata))
  print(colnames(coldata))

  ######################## Step 3: DESeq ########################
  dds <- DESeqDataSetFromMatrix(countData = cts_counts_df, colData = coldata, design = ~ condition)
  
  # Default DESeq2
  dds <- DESeq(dds)
  message("Default DESeq2 size factors:")
  print(sizeFactors(dds))

  ######################## Step 4: Custom size factors ########################
  default_size_factors <- sizeFactors(dds)
  if (length(spike_genomes) == 2) {
    norm_col <- "dual.normfactor.ipeff.adj"
  } else if (length(spike_genomes) == 1) {
    norm_col <- paste0(spike_genomes[[1]], ".normfactor.ipeff.adj")
  } else stop("spike_genomes must be length 1 or 2")

  norm_meta_file <- file.path(dirname(sample_metadata), "sample_metadata.norm.tsv")
  if (!file.exists(norm_meta_file)) stop(paste("Normalization metadata not found:", norm_meta_file))
  norm_meta <- read.delim(norm_meta_file)
  seqstats_with_sf <- merge(norm_meta, data.frame(default_size_factors, library.ID = names(default_size_factors)), by = "library.ID")
  custom_size_factors <- seqstats_with_sf[[norm_col]] * seqstats_with_sf$default_size_factors
  names(custom_size_factors) <- seqstats_with_sf$library.ID
  custom_size_factors <- custom_size_factors[names(default_size_factors)]
  sizeFactors(dds) <- custom_size_factors / median(custom_size_factors)

  message("Custom size factors to be used for spike-in normalization:")
  print(sizeFactors(dds))

  ######################## Step 5: Spike-in normalized DESeq ########################
  dds_spikenorm <- DESeq(dds)

  ######################## Step 6: Differential expression ########################
  res_spikenorm <- results(dds_spikenorm, contrast = c("condition", condB, condA))
  res_shr <- lfcShrink(dds_spikenorm, coef = shrink_coef, type = "apeglm")
  res_shr_df <- as.data.frame(res_shr)

  ### Count significant peaks
  up_peaks <- sum(res_shr_df$log2FoldChange > 1 & res_shr_df$padj < 0.05, na.rm = TRUE)
  down_peaks <- sum(res_shr_df$log2FoldChange < -1 & res_shr_df$padj < 0.05, na.rm = TRUE)
  message("Significant UP peaks: ", up_peaks)
  message("Significant DOWN peaks: ", down_peaks)

  return(list(dds = dds, dds_spikenorm = dds_spikenorm, res = res_spikenorm, res_shr = res_shr_df,
              up_peaks = up_peaks, down_peaks = down_peaks))
}

###############################################
### Command-line interface (CLI)
###############################################

args <- commandArgs(trailingOnly = TRUE)
log_file <- "DESeq2_pipeline.log"
log_con <- file(log_file, open = "wt")

# --- Logging: stdout to console+file, messages to file ---
sink(log_con, split = TRUE)          # stdout
sink(log_con, type = "message")      # messages

if (sys.nframe() == 0 && length(args) > 0) {
  tryCatch({
    suppressPackageStartupMessages(library(optparse))
    parser <- OptionParser(option_list = list(
      make_option(c("-c","--counts"), type="character"),
      make_option(c("-m","--metadata"), type="character"),
      make_option(c("-d","--conditions"), type="character"),
      make_option(c("-s","--spike_genomes"), type="character"),
      make_option(c("-o","--outprefix"), type="character", default="DESeq2_results")
    ))
    opts <- parse_args(parser)

    conds <- strsplit(opts$conditions, ",")[[1]]
    spike_list <- strsplit(opts$spike_genomes, ",")[[1]]

    results <- deseq_with_chipwrangler(opts$counts, opts$metadata, conds, spike_genomes = spike_list)

    # Save output tables
    write.table(results$res, file = paste0(opts$outprefix,"_DESeq2.tsv"), sep="\t", quote=FALSE, row.names=TRUE)
    write.table(results$res_shr, file = paste0(opts$outprefix,"_DESeq2_shrunken.tsv"), sep="\t", quote=FALSE, row.names=TRUE)
    write.table(data.frame(UP=results$up_peaks, DOWN=results$down_peaks),
                file = paste0(opts$outprefix,"_DESeq2_significant_summary.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)

    message("DESeq2 finished successfully")

    # Close logging
    sink(type="message")
    sink()
    close(log_con)
    quit(status = 0)

  }, error = function(e){
    message("ERROR: ", e$message)
    sink(type="message")
    sink()
    close(log_con)
    quit(status = 1)
  })
}
