#' Preprocess GWAS Summary Statistics for MRcare
#'
#' This function preprocesses GWAS summary statistics to prepare them for use with MRcare.
#' It standardizes column names, performs quality control, and saves the processed data
#' to a file for later use.
#'
#' @param gwas_data Path to a GWAS summary statistics file or a data frame containing GWAS summary statistics
#' @param output_file Path where the processed GWAS data will be saved
#' @param data_type Either "exposure" or "outcome" to determine processing logic
#' @param A1 Column name for effect allele (default: NULL, will try to detect)
#' @param A2 Column name for other allele (default: NULL, will try to detect)
#' @param SNP Column name for SNP ID (default: NULL, will try to detect)
#' @param CHR Column name for chromosome (default: NULL, will try to detect)
#' @param POS Column name for position (default: NULL, will try to detect)
#' @param BETA Column name for effect size (default: NULL, will try to detect)
#' @param SE Column name for standard error (default: NULL, will try to detect)
#' @param Pval Column name for p-value (default: NULL, will try to detect)
#' @param MAF Column name for minor allele frequency (default: NULL, will try to detect)
#' @param N Column name for sample size (default: NULL, will try to detect)
#' @param Zscore Column name for Z-score (default: NULL, will try to detect)
#' @param Ninput Sample size if not available in data (default: NULL)
#' @param maf_filter Minimum MAF threshold to include SNPs (default: 0.01)
#' @param remove_ambiguous Whether to remove strand-ambiguous SNPs (default: TRUE)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Path to the processed GWAS data file
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process exposure data
#' exposure_processed <- preprocess_gwas_data(
#'   gwas_data = "path/to/raw_exposure_gwas.txt",
#'   output_file = "path/to/processed_exposure_gwas.rds",
#'   data_type = "exposure",
#'   Ninput = 10000
#' )
#'
#' # Process outcome data
#' outcome_processed <- preprocess_gwas_data(
#'   gwas_data = "path/to/raw_outcome_gwas.txt",
#'   output_file = "path/to/processed_outcome_gwas.rds",
#'   data_type = "outcome",
#'   Ninput = 15000
#' )
#' }
preprocess_gwas_data <- function(gwas_data, 
                               output_file, 
                               data_type = c("exposure", "outcome"),
                               A1 = NULL, 
                               A2 = NULL, 
                               SNP = NULL, 
                               CHR = NULL, 
                               POS = NULL, 
                               BETA = NULL, 
                               SE = NULL, 
                               Pval = NULL, 
                               MAF = NULL, 
                               N = NULL, 
                               Zscore = NULL,
                               Ninput = NULL,
                               maf_filter = 0.01,
                               remove_ambiguous = TRUE,
                               verbose = TRUE) {
  
  data_type <- match.arg(data_type)
  
  # Start timing for performance tracking
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n===================================================================\n")
    cat(" Preprocessing GWAS Summary Statistics for MRcare\n")
    cat("===================================================================\n\n")
    cat("Processing started at:", format(start_time), "\n\n")
  }
  
  # Read input data
  if (is.character(gwas_data) && file.exists(gwas_data)) {
    if (verbose) cat("Reading GWAS data from file:", gwas_data, "\n")
    
    # Determine file type and read accordingly
    if (grepl("\\.rds$", gwas_data, ignore.case = TRUE)) {
      raw_data <- readRDS(gwas_data)
      if (verbose) cat("Read RDS file successfully\n")
    } else {
      # Assume it's a tab-delimited file
      raw_data = data.table::fread(gwas_data)
      raw_data = as.data.frame(raw_data)
    }
  } else if (is.data.frame(gwas_data)) {
    if (verbose) cat("Using provided data frame\n")
    raw_data <- gwas_data
  } else {
    stop("gwas_data must be either a data frame or a valid file path")
  }
  
  if (verbose) cat("GWAS data contains", nrow(raw_data), "variants\n")
  
  # Standardize column names
  if (verbose) cat("Standardizing column names...\n")
  
  header.inner <- colnames(raw_data)
  header.inner <- tolower(header.inner)
  
  # SNP
  if(is.null(SNP)) {
    try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps", "variant", "variant_id", "varid", "id")
  } else {
    try.snp <- tolower(SNP)
  }
  header.inner[header.inner %in% try.snp] <- "SNP"
  
  # A1
  if(is.null(A1)) {
    try.a1 <- c("a1", "allele1", "allele_1", "effect_allele", "reference_allele", "inc_allele", "ea", "alt", "a1lele1", "al1ele1")
  } else {
    try.a1 <- tolower(A1)
  }
  header.inner[header.inner %in% try.a1] <- "A1"
  
  # A2
  if(is.null(A2)) {
    try.a2 <- c("a2", "allele2", "allele_2", "other_allele", "non_effect_allele", "dec_allele", "nea", "ref", "a0")
  } else {
    try.a2 <- tolower(A2)
  }
  header.inner[header.inner %in% try.a2] <- "A2"
  
  # Z-score
  if(is.null(Zscore)) {
    try.z <- c("zscore", "z-score", "gc_zscore", "z", "tstat", "t-statistic", "tstatistic")
  } else {
    try.z <- tolower(Zscore)
  }
  header.inner[header.inner %in% try.z] <- "Z"
  
  # Chromosome
  if(is.null(CHR)) {
    try.chromosome <- c("chrom", "ch", "chr", "chromosome", "#chr", "chr_id", "chromosome_id")
  } else {
    try.chromosome <- tolower(CHR)
  }
  header.inner[header.inner %in% try.chromosome] <- "CHR"
  
  # Position
  if(is.null(POS)) {
    try.pos <- c("pos", "bp", "position", "base_pair", "basepair", "base_position", "position_bp")
  } else {
    try.pos <- tolower(POS)
  }
  header.inner[header.inner %in% try.pos] <- "POS"
  
  # P-value
  if(is.null(Pval)) {
    try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p", "all_inv_var_meta_p", "p.value", "p-value", "pvalue_gc")
  } else {
    try.p <- tolower(Pval)
  }
  header.inner[header.inner %in% try.p] <- "P"
  
  # Beta
  if(is.null(BETA)) {
    try.beta <- c("b", "beta", "effects", "effect", "all_inv_var_meta_beta", "beta_gc", "effect_size", "effectsize")
  } else {
    try.beta <- tolower(BETA)
  }
  header.inner[header.inner %in% try.beta] <- "BETA"
  
  # Standard error
  if(is.null(SE)) {
    try.se <- c("se", "sebeta", "beta_se", "all_inv_var_meta_sebeta", "standard_error", "stderr", "se_beta")
  } else {
    try.se <- tolower(SE)
  }
  header.inner[header.inner %in% try.se] <- "SE"
  
  # MAF/EAF
  if(is.null(MAF)) {
    try.maf <- c("eaf", "frq", "maf", "frq_u", "f_u", "af1", "all_meta_af", "af", "allele_freq", "a1freq", "freq")
  } else {
    try.maf <- tolower(MAF)
  }
  header.inner[header.inner %in% try.maf] <- "EAF"
  
  # Sample size
  if(is.null(N)) {
    try.n <- c("n", "nsample", "nsum", "samplesize", "sample_size", "totaln", "total_n", "n_samples")
  } else {
    try.n <- tolower(N)
  }
  header.inner[header.inner %in% try.n] <- "N"
  
  # Update column names
  colnames(raw_data) <- header.inner
  

  list.coerce <- c("Z", "BETA", "ODDS_RATIO", "LOG_ODDS", "SE", "AF1","N")

  for (i in 1:length(header.inner)) {
      if (header.inner[i] %in% list.coerce) {
          if (class(raw_data[, header.inner[i]]) != "numeric") {
                class(raw_data[, header.inner[i]]) <- "numeric"
                cat(paste0("Column ", header.inner[i], " has wrong class and has been coerced to numeric."), sep = "\n")
                cat("=============================================================================================================", sep = "\n")
          }
      }
  }

  # Check required columns
  required_cols <- c("SNP")
  missing_cols <- required_cols[!required_cols %in% colnames(raw_data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required column: ", paste(missing_cols, collapse = ", "))
  }
  
  # Handle various combinations of available effect columns
  if (!("BETA" %in% colnames(raw_data))) {
    if ("Z" %in% colnames(raw_data) && "SE" %in% colnames(raw_data)) {
      raw_data$BETA <- raw_data$Z * raw_data$SE
      if (verbose) cat("Calculated BETA from Z and SE\n")
    } else if ("LOG_ODDS" %in% colnames(raw_data)) {
      raw_data$BETA <- raw_data$LOG_ODDS
      if (verbose) cat("Used LOG_ODDS as BETA\n")
    } else if ("ODDS_RATIO" %in% colnames(raw_data)) {
      raw_data$BETA <- log(raw_data$ODDS_RATIO)
      if (verbose) cat("Calculated BETA from ODDS_RATIO\n")
    } else {
      stop("Cannot determine effect size (BETA). Please provide column mapping.")
    }
  }
  
  # Calculate Z-score if missing
  if (!("Z" %in% colnames(raw_data))) {
    if ("BETA" %in% colnames(raw_data) && "SE" %in% colnames(raw_data)) {
      raw_data$Z <- raw_data$BETA / raw_data$SE
      if (verbose) cat("Calculated Z score from BETA and SE\n")
    } else if ("BETA" %in% colnames(raw_data) && "P" %in% colnames(raw_data)) {
      raw_data$Z <- sign(raw_data$BETA) * abs(qnorm(raw_data$P / 2))
      if (verbose) cat("Calculated Z score from BETA and P\n")
    } else {
      if (verbose) cat("Warning: Could not calculate Z score. A Z score column will be missing.\n")
    }
  }
  
  # Calculate P-value if missing
  if (!("P" %in% colnames(raw_data)) && "Z" %in% colnames(raw_data)) {
    raw_data$P <- pnorm(abs(raw_data$Z), lower.tail = FALSE) * 2
    if (verbose) cat("Calculated P values from Z scores\n")
  }
  
  # Handle missing SE
  if (!("SE" %in% colnames(raw_data)) && "BETA" %in% colnames(raw_data) && "Z" %in% colnames(raw_data)) {
    raw_data$SE <- abs(raw_data$BETA / raw_data$Z)
    if (verbose) cat("Calculated SE from BETA and Z\n")
  }
  
  # Add sample size if provided
  if (is.null(raw_data$N) && !is.null(Ninput)) {
    raw_data$N <- Ninput
    if (verbose) cat("Added sample size N =", Ninput, "\n")
  }
  
  # Quality control
  if (verbose) cat("Performing quality control on variants...\n")
  
  # Original count
  original_count <- nrow(raw_data)
  if (verbose) cat("Original variants:", original_count, "\n")
  
  # 1. Remove variants with missing values in key columns
  key_cols <- intersect(c("SNP", "A1", "A2", "BETA", "SE", "P", "Z"), colnames(raw_data))
  complete_cases <- complete.cases(raw_data[, key_cols])
  raw_data <- raw_data[complete_cases, ]
  if (verbose) cat("After removing variants with missing values:", nrow(raw_data), 
                  "(removed", original_count - nrow(raw_data), "variants)\n")
  
  # 2. Remove strand-ambiguous SNPs if requested
  if (remove_ambiguous && all(c("A1", "A2") %in% colnames(raw_data))) {
    # Make sure alleles are uppercase
    raw_data$A1 <- toupper(raw_data$A1)
    raw_data$A2 <- toupper(raw_data$A2)
    
    # Identify strand-ambiguous SNPs
    a1 <- raw_data$A1
    a2 <- raw_data$A2
    
    # Strand-ambiguous pairs: A/T, T/A, C/G, G/C
    ambiguous <- (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | 
                 (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
    
    # Also check that the alleles are single characters (not indels)
    valid_alleles <- nchar(a1) == 1 & nchar(a2) == 1
    
    # Keep only non-ambiguous, valid SNPs
    keep <- !ambiguous & valid_alleles
    
    pre_filter <- nrow(raw_data)
    raw_data <- raw_data[keep, ]
    if (verbose) cat("After removing strand-ambiguous SNPs and indels:", nrow(raw_data), 
                    "(removed", pre_filter - nrow(raw_data), "variants)\n")
  }
  
  # 3. MAF filter if available
  if ("EAF" %in% colnames(raw_data) && !is.null(maf_filter)) {
    pre_filter <- nrow(raw_data)
    
    # First remove any rows with missing MAF
    raw_data <- raw_data[!is.na(raw_data$EAF), ]
    
    # Check for values outside [0,1]
    invalid_eaf <- raw_data$EAF < 0 | raw_data$EAF > 1
    if (any(invalid_eaf)) {
      if (verbose) cat("Warning:", sum(invalid_eaf), "variants had EAF outside [0,1]. These will be removed.\n")
      raw_data <- raw_data[!invalid_eaf, ]
    }
    
    # Calculate MAF by taking the minimum of EAF and (1-EAF)
    maf <- pmin(raw_data$EAF, 1 - raw_data$EAF)
    
    # Apply MAF filter
    raw_data <- raw_data[maf >= maf_filter, ]
    
    if (verbose) cat("After MAF filtering (MAF greater", maf_filter, "):", nrow(raw_data), 
                    "(removed", pre_filter - nrow(raw_data), "variants)\n")
  }
  
  # 4. Remove duplicated SNPs
  pre_filter <- nrow(raw_data)
  duplicated_snps <- duplicated(raw_data$SNP)
  raw_data <- raw_data[!duplicated_snps, ]
  if (verbose) cat("After removing duplicated SNPs:", nrow(raw_data), 
                  "(removed", pre_filter - nrow(raw_data), "variants)\n")
  
  # Ensure all standard columns exist (create them if needed)
  std_columns <- c("SNP", "A1", "A2", "BETA", "SE", "P", "Z", "EAF", "N", "CHR", "POS")
  for (col in std_columns) {
    if (!(col %in% colnames(raw_data))) {
      raw_data[[col]] <- NA
      if (verbose) cat("Added empty column for", col, "\n")
    }
  }
  
  # Determine which columns to keep
  keep_cols <- unique(c(std_columns, colnames(raw_data)))
  raw_data <- raw_data[, keep_cols]
  
  # Add data type identifier
  raw_data$data_type <- data_type
  
  # Save processed data
  if (verbose) cat("Saving processed data to", output_file, "\n")
  
  write.table(raw_data,output_file,row.names=FALSE,quote=FALSE)
  
  # End timing
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  if (verbose) {
    cat("\n===================================================================\n")
    cat(" Preprocessing completed successfully\n")
    cat("===================================================================\n")
    cat("Output file:", output_file, "\n")
    cat("Final variant count:", nrow(raw_data), "\n")
    cat("Processing time:", round(as.numeric(duration), 2), "minutes\n\n")
  }

  cat("Note: the effect allele is A1 and the non effect allele is A2.\n")
  cat("Please ensure that this has been coded correctly.")
  
  # Return the path to the processed file
  return(output_file)
}
