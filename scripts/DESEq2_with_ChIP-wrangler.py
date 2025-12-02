deseq_with_chipwrangler <- function(
    counts_file,
    sample_metadata,
    conditions,             # list of 2 strings, e.g. list("TRP_4hr","TRP_0hr")
    condition_col = "condition"
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
      if (length(hits) > 1) {
        stop(paste("ERROR: Multiple library IDs matched column:", col,
                   "\nMatches:", paste(hits, collapse=", ")))
      }

      id <- hits[[1]]

      after_id <- sub(paste0(".*(", id, ")"), "\\1", col)
      after_clean <- sub("(Tag.*$)", "", after_id)
      final <- gsub("[\\.\\s_]+$", "", after_clean)

      return(id)
    })

    return(cleaned)
  }

  ######################## Step 1: Load metadata & counts ########################
  sample_meta <- read.delim(sample_metadata)
  library_ids <- sample_meta$library.ID

  counts_df <- read.delim(counts_file)

  ### Identify count columns (anything with "Tag") ###
  count_cols <- grep("Tag", colnames(counts_df), value = TRUE)

  ### Clean count column names ###
  new_names <- clean_colnames_using_metadata(count_cols, library_ids)
  colnames(counts_df)[match(count_cols, colnames(counts_df))] <- new_names

  ### Extract counts matrix (only cleaned names) ###
  cts_counts_df <- counts_df[, new_names, drop = FALSE]

  ### Convert to integer counts ###
  cts_counts_df[] <- lapply(cts_counts_df, function(x) {
    if (is.numeric(x) && any(x != floor(x), na.rm = TRUE)) {
      as.integer(round(x))
    } else x
  })


  ######################## Step 2: Construct colData ########################
  samples <- colnames(cts_counts_df)
  coldata <- data.frame(row.names = samples)

  #### YOUR REGEX EXTRACTION — preserved exactly ####
  coldata$treatment <- factor(str_match(samples, '([^_]+)(?:_[^_]+){4}$')[,2])
  coldata$timepoint <- factor(str_match(samples, '([^_]+)(?:_[^_]+){3}$')[,2])
  coldata$biorep    <- factor(str_match(samples, '([^_]+)(?:_[^_]+){2}$')[,2])
  coldata$antibody  <- factor(str_match(samples, '([^_]+)(?:_[^_]+){1}$')[,2])

  #### Your repgroup logic ####
  coldata$repgroup <- ifelse(coldata$biorep %in% c("1","2"), "rep12", "rep34")
  coldata$repgroup <- factor(coldata$repgroup)

  #### Combined condition ####
  coldata$condition <- factor(paste0(coldata$treatment, "_", coldata$timepoint))

  ### QC: must match ###
  if (!all(rownames(coldata) == colnames(cts_counts_df))) {
    stop("ERROR: coldata rownames do NOT match count matrix column names.")
  }


  ######################## Step 3: Build DESeq object ########################
  dds <- DESeqDataSetFromMatrix(
    countData = cts_counts_df,
    colData = coldata,
    design = ~ biorep + condition
  )

  ######################## Step 4: Default DESeq ########################
  dds <- DESeq(dds)

  ######################## Step 5: QC ########################
  message("QC: dispersions summary:")
  print(summary(mcols(dds)$dispersion))

  message("Available DESeq coefficients:")
  print(resultsNames(dds))


  ######################## Step 6: Custom size factors ########################
  default_size_factors <- sizeFactors(dds)
  seqstats <- sample_meta

  default_size_factors_df <- data.frame(
    default_size_factors,
    library.ID = names(default_size_factors)
  )

  seqstats_with_sf <- merge(
    seqstats,
    default_size_factors_df,
    by = "library.ID"
  ) %>%
    mutate(dual_size_factors = dual.normfactor.ipeff.adj * default_size_factors)

  custom_size_factors <- seqstats_with_sf$dual_size_factors
  names(custom_size_factors) <- seqstats_with_sf$library.ID

  custom_size_factors_reordered <- custom_size_factors[names(default_size_factors)]

  ### QC checks ###
  if (!all(names(custom_size_factors_reordered) == names(default_size_factors))) {
    stop("ERROR: Custom size factors are not aligned.")
  }

  ### Apply relative size factors (median normalized) ###
  sizeFactors(dds) <- custom_size_factors_reordered /
                      median(custom_size_factors_reordered)

  ######################## Step 7: Spike-in normalized DESeq ########################
  dds_spikenorm <- DESeq(dds)

  ######################## Step 8: Differential expression ########################
  res_spikenorm <- results(dds_spikenorm,
                           contrast = c("condition", condB, condA))

  ### Shrinkage using defined coefficient ###
  res_shr <- lfcShrink(
    dds_spikenorm,
    coef = shrink_coef,
    type = "apeglm"
  )

  return(list(
    dds = dds,
    dds_spikenorm = dds_spikenorm,
    res = res_spikenorm,
    res_shr = as.data.frame(res_shr)
  ))
}

###############################################
### Command-line interface (CLI)
###############################################

if (sys.nframe() == 0) {
    suppressPackageStartupMessages(library(optparse))

    option_list <- list(
        make_option(c("-c", "--counts"),
                    type = "character", help = "Counts file (TSV)"),

        make_option(c("-m", "--metadata"),
                    type = "character", help = "Sample metadata file (TSV)"),

        make_option(c("-d", "--conditions"),
                    type = "character",
                    help = "Two conditions separated by a comma, e.g. 'TRP_0hr,TRP_4hr'"),

        make_option(c("-o", "--outprefix"),
                    type = "character", default = "DESeq2_results",
                    help = "Prefix for output files [default: %default]")
    )

    parser <- OptionParser(
        usage = "Usage: %prog -c raw_counts.tsv -m metadata.tsv -d condA,condB",
        option_list = option_list
    )

    opts <- parse_args(parser)

    ### Required inputs ###
    if (is.null(opts$counts) || is.null(opts$metadata) || is.null(opts$conditions)) {
        print_help(parser)
        stop("\nERROR: --counts, --metadata, and --conditions are required.\n")
    }

    ### Parse conditions ###
    conds <- strsplit(opts$conditions, ",")[[1]]
    if (length(conds) != 2) {
        stop("ERROR: --conditions must contain exactly two values separated by a comma.")
    }

    message("Running DESeq2 with ChIP-wrangler normalization...")
    message("Counts file:    ", opts$counts)
    message("Metadata file:  ", opts$metadata)
    message("Conditions:     ", paste(conds, collapse = " vs "))
    message("Output prefix:  ", opts$outprefix)

    ### Run the function ###
    results <- deseq_with_chipwrangler(
        counts_file = opts$counts,
        sample_metadata = opts$metadata,
        conditions = conds
    )

    ### Save the results ###
    out_base <- opts$outprefix

    write.table(results$res,
                file = paste0(out_base, "_DESeq2.tsv"),
                sep = "\t", quote = FALSE, row.names = TRUE)

    write.table(results$res_shr,
                file = paste0(out_base, "_DESeq2_shrunken.tsv"),
                sep = "\t", quote = FALSE, row.names = TRUE)

    message("\n✔ Output written:")
    message("  - ", out_base, "_DESeq2.tsv")
    message("  - ", out_base, "_DESeq2_shrunken.tsv")
    message("Done.\n")
}

