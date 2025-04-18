#' Format CARE Analysis Results
#'
#' Formats the output from CARE analysis into a standardized structure
#' for easier interpretation and use.
#'
#' @param care_results Results from CARE2_boot function
#' @param harmonized_data The harmonized exposure and outcome data
#' @param rivw_result Results from RIVW function
#'
#' @return A list containing formatted results for the BIC method
#'
#' @keywords internal
format_care_results <- function(care_results, harmonized_data, rivw_result) {
  # Initialize results list
  formatted_results <- list()
  
  # Process BIC results only
  formatted_results$BIC <- list(
    estimate = care_results$res$BIC["tilde_theta"],
    se = care_results$res$BIC["se"],
    p = care_results$res$BIC["p"],
    raw_estimate = care_results$res$BIC["hat_theta"],
    raw_p = care_results$res$BIC["hat_p"],
    efron_se = care_results$res$BIC["Efron_se"],
    efron_p = care_results$res$BIC["Efron_p"],
    valid_count = care_results$res$IV["BIC"]
  )
  
  # Add RIVW results
  formatted_results$RIVW <- rivw_result
  
  # Add setting information
  formatted_results$setting <- care_results$res$setting
  
  # Add invalid IV indicators (just for BIC now)
  formatted_results$invalid_ivs <- list(
    BIC = care_results$BIC_IV
  )
  
  # Add harmonized data for reference
  formatted_results$data <- harmonized_data
  
  # Add bootstrap results for further analysis if needed
  #formatted_results$bootstrap_results <- care_results$MRcML_boot
  
  return(formatted_results)
}


#' Extract Independent Variants (Updated to use configuration)
#'
#' Extracts independent variants from GWAS summary statistics using clumping with
#' support for re-randomization and configuration.
#'
#' @param exposure_data Exposure GWAS summary statistics
#' @param p_threshold P-value threshold for selecting significant SNPs
#' @param etamean Mean of the random variable for re-randomization
#' @param ref_panel Path to reference panel for LD clumping
#' @param clump_r2 Clumping r-squared threshold
#' @param clump_kb Clumping window size in kb
#' @param temp_dir Directory for temporary files
#' @param verbose Whether to print detailed progress messages
#' @param rerand Whether to use re-randomization for SNP selection
#' @param plink_path Path to plink executable
#'
#' @return A vector of independent SNP IDs
#'
#' @importFrom utils download.file unzip read.table write.table packageVersion setTxtProgressBar txtProgressBar
#' @keywords internal
extract_independent_variants <- function(exposure_data,
                                       p_threshold = 5e-5,
                                       etamean = 0.5,
                                       ref_panel = NULL,
                                       clump_r2 = 0.001,
                                       clump_kb = 10000,
                                       temp_dir = NULL,
                                       verbose = TRUE,
                                       rerand = TRUE,
                                       plink_path = NULL) {
  
  # Check if configuration should be used
  if (is.null(plink_path) || is.null(ref_panel)) {
    mrcare_config <- get_mrcare_config(create_if_missing = TRUE, verbose = FALSE)
    if (is.null(plink_path)) plink_path <- mrcare_config$plink_path
    if (is.null(ref_panel)) ref_panel <- mrcare_config$ref_panel_path
  }
  
  # Make sure we have plink path
  if (is.null(plink_path)) {
    stop("PLINK executable path not provided and not found in configuration. Run configure_mrcare() first.")
  }
  
  # Make sure we have reference panel
  if (is.null(ref_panel)) {
    stop("Reference panel path not provided and not found in configuration. Run configure_mrcare() first.")
  }
  
  # Create temporary directory if not provided
  if (is.null(temp_dir)) {
    temp_dir <- file.path(tempdir(), "MRcare_temp")
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE)
    }
    if (verbose) cat("Created temporary directory at:", temp_dir, "\n")
  }
  
  # Create a unique identifier for this run to avoid file conflicts
  run_id <- paste0("run_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample(1000:9999, 1))
  
  # Process each chromosome separately to manage memory
  pruned_snps <- c()
  total_sig_snps <- 0
  
  if (verbose) {
    cat("\nExtracting independent variants by chromosome:\n")
    cat("-----------------------------------------------\n")
  }
  
  for (chr_id in 1:22) {
    if (verbose) cat(sprintf("Processing chromosome %d...", chr_id))
    
    # Filter data by chromosome
    chr_data <- exposure_data[exposure_data$CHR == chr_id, ]
    
    if (nrow(chr_data) == 0) {
      if (verbose) cat(" No variants found\n")
      next
    }
    

    # Generate random values for re-randomization if needed
    if (rerand) {
      W <- rnorm(nrow(chr_data), 0, etamean)
    }
    
    # Select significant SNPs using the provided approach
    C_sel <- qnorm(p_threshold/2, lower.tail = FALSE)
    
    # Calculate z-scores
    z_scores <- chr_data$BETA / chr_data$SE
    
    # Apply custom selection logic
    if(rerand) {
      # Re-randomized p-value approach
      z_scores_rerand <- z_scores + W
      pvals_rerand <- pnorm(abs(z_scores_rerand/sqrt(1 + etamean^2)), lower.tail = FALSE) * 2
      
      # Select SNPs that meet both original and re-randomized thresholds
      sig_indices <- which(abs(z_scores) >= C_sel - etamean & abs(z_scores_rerand) >= C_sel)
    } else {
      # Standard approach
      sig_indices <- which(abs(z_scores) >= C_sel)
    }
    
    if (length(sig_indices) == 0) {
      if (verbose) cat(" No significant SNPs found\n")
      next
    }
    
    total_sig_snps <- total_sig_snps + length(sig_indices)
    
    # Select the significant SNPs
    sig_data <- chr_data[sig_indices, ]
    sig_data <- as.data.frame(sig_data)
    
    # Get standard errors for selected SNPs
    tmpSE <- sig_data$SE
    
    
    # Use your custom P-value method; this is to avoid bias induced by clumping
    if (rerand) {
      sig_data$P <- 1/(1 + exp(-tmpSE)) - 0.5
    } else {
      sig_data$P <- pnorm(abs(z_scores), lower.tail = FALSE) * 2
    }
    
    if (verbose) cat(sprintf(" Found %d significant SNPs\n", length(sig_indices)))
    
    temp_dir <-  path.expand(temp_dir)
    # Write SNPs and p-values to files for plink
    snp_file <- file.path(temp_dir, paste0(run_id, "_chr", chr_id, "_snps.txt"))
    p_val_file <- file.path(temp_dir, paste0(run_id, "_chr", chr_id, "_pvals.txt"))
    
    #snp_file <- path.expand(snp_file)
    #p_val_file <- path.expand(p_val_file)

    # Make sure we have SNP and P columns for clumping
    p_vals_df <- data.frame(SNP = sig_data$SNP, P = sig_data$P)
    
    # Write to files
    write.table(sig_data$SNP, file = snp_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(p_vals_df, file = p_val_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    # Improved reference panel path detection logic
    # First try to determine the structure of the reference panel directory
    possible_file_patterns <- c(
      # Structure 1: /path/to/ref_panel/chr1.bed
      file.path(ref_panel, paste0("chr", chr_id, ".bed")),
      # Structure 2: /path/to/ref_panel/chr1/chr1.bed
      file.path(ref_panel, paste0("chr", chr_id), paste0("chr", chr_id, ".bed")),
      # Structure 3: /path/to/ref_panel/1000G.EUR.ALLSNP.QC.CHR1.bed
      file.path(ref_panel, paste0("1000G.", "EUR", ".ALLSNP.QC.CHR", chr_id, ".bed")),
      # Structure 3: /path/to/ref_panel/1000G.AFR.ALLSNP.QC.CHR1.bed
      file.path(ref_panel, paste0("1000G.", "AFR", ".ALLSNP.QC.CHR", chr_id, ".bed")),
      # Structure 3: /path/to/ref_panel/1000G.EAS.ALLSNP.QC.CHR1.bed
      file.path(ref_panel, paste0("1000G.", "EAS", ".ALLSNP.QC.CHR", chr_id, ".bed")),
      # Structure 4: /path/to/ref_panel/1.bed
      file.path(ref_panel, paste0(chr_id, ".bed"))
    )
    
    chr_ref_panel <- NULL
    
    for (pattern in possible_file_patterns) {
      if (file.exists(pattern)) {
        # Found a matching file, extract the base path (without .bed extension)
        chr_ref_panel <- substr(pattern, 1, nchar(pattern) - 4)
        break
      }
    }

    if (is.null(chr_ref_panel)) {
      if (verbose) cat(" Warning: Could not find reference panel files for chr", chr_id, "\n")
    }

    chr_ref_panel <- path.expand(chr_ref_panel)
    # Run clumping
    if (verbose) cat(sprintf("   Running clumping or revised clumping on chr %d...", chr_id))
    

    # TODO: revise this and improve; right now, there is a bug. We may just use my original codes.
    # Not sure what happend.

    # Construct plink command #      " --silent",
    clump_cmd <- paste0(
      shQuote(plink_path),
      " --bfile ", shQuote(chr_ref_panel),
      " --clump ", shQuote(p_val_file),
      " --clump-p1 1",
      " --clump-kb ", clump_kb,
      " --clump-r2 ", clump_r2,
      " --clump-snp-field SNP",
      " --clump-field P",
      " --out ", shQuote(file.path(temp_dir, paste0(run_id, "_chr", chr_id)))
    )
    
    # Run the command, capturing output to avoid printing to console
    system_result <- suppressWarnings(
      system(clump_cmd, intern = TRUE, ignore.stderr = TRUE)
    )
    
    # Read clumped SNPs
    clumped_file <- paste0(file.path(temp_dir, paste0(run_id, "_chr", chr_id)), ".clumped")
    chr_clumped_snps <- c()
    
    if (file.exists(clumped_file)) {
      clumped_data <- tryCatch({
        data.table::fread(clumped_file, header = TRUE)
      }, error = function(e) {
        if (verbose) cat(" Error reading clumped file:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(clumped_data)) {
        clumped_data <- as.data.frame(clumped_data)
        chr_clumped_snps <- clumped_data$SNP
        pruned_snps <- c(pruned_snps, chr_clumped_snps)
        if (verbose) cat(sprintf(" Selected %d SNPs\n", length(chr_clumped_snps)))
      } else {
        if (verbose) cat(" No SNPs selected\n")
      }
    } else {
      if (verbose) cat(" No clumped SNPs found\n")
    }
    
    # Clean up temporary files
    cleanup_files <- list.files(temp_dir, pattern = paste0(run_id, "_chr", chr_id), full.names = TRUE)
    if (length(cleanup_files) > 0) {
      suppressWarnings(file.remove(cleanup_files))
    }
  }
  
  # Summarize results
  if (verbose) {
    cat("-----------------------------------------------\n")
    cat(sprintf("Total significant SNPs across all chromosomes: %d\n", total_sig_snps))
    cat(sprintf("Final number of independent variants: %d\n\n", length(pruned_snps)))
  }
  
  # Return unique SNPs to avoid any potential duplicates
  return(unique(pruned_snps))
}


#' Summarize MRcare Analysis Results
#'
#' Provides a comprehensive summary of the results from an MRcare analysis,
#' including causal estimates, confidence intervals, and diagnostic information.
#'
#' @param object Results object from mr_care function
#' @param ci_level Confidence level for intervals. Default is 0.95
#' @param digits Number of digits to round results to. Default is 3
#'
#' @return A list of class "summary.mrcare" containing formatted summary information
#'
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' \dontrun{
#' results <- mr_care(exposure_data, outcome_data)
#' summary(results)
#' }
summary.mrcare <- function(object, ..., ci_level = 0.95, digits = 3) {
  # Calculate z-value for CI
  z <- qnorm(1 - (1 - ci_level) / 2)
  
  # Create summary object
  summary_obj <- list()
  
  # General information
  summary_obj$general <- list(
    timestamp = object$metadata$timestamp,
    duration_mins = round(object$metadata$duration_mins, 1),
    exposure = object$metadata$exposure_name,
    outcome = object$metadata$outcome_name,
    parameters = object$metadata$parameters
  )
  
  # SNP information
  summary_obj$snp_info <- list(
    total_snps = nrow(object$data),
    valid_snps = object$BIC$valid_count,
    invalid_snps = sum(object$invalid_ivs$BIC != 0, na.rm = TRUE),
    f_statistic = round(object$setting["F"], digits)
  )
  
  # BIC results
  if ("BIC" %in% names(object)) {
    bic_estimate <- object$BIC$estimate
    bic_se <- object$BIC$efron_se
    bic_p <- object$BIC$efron_p
    bic_lower <- bic_estimate - z * bic_se
    bic_upper <- bic_estimate + z * bic_se
    
    summary_obj$BIC <- list(
      estimate = round(bic_estimate, digits),
      se = round(bic_se, digits),
      p_value = format(bic_p, digits = digits, scientific = TRUE),
      ci_lower = round(bic_lower, digits),
      ci_upper = round(bic_upper, digits),
      valid_instruments = object$BIC$valid_count
    )
  }
  
  # RIVW results
  if ("RIVW" %in% names(object)) {
    rivw_estimate <- object$RIVW$beta
    rivw_se <- object$RIVW$sd
    rivw_p <- object$RIVW$p
    rivw_lower <- rivw_estimate - z * rivw_se
    rivw_upper <- rivw_estimate + z * rivw_se
    
    summary_obj$RIVW <- list(
      estimate = round(rivw_estimate, digits),
      se = round(rivw_se, digits),
      p_value = format(rivw_p, digits = digits, scientific = TRUE),
      ci_lower = round(rivw_lower, digits),
      ci_upper = round(rivw_upper, digits),
      nIV = object$RIVW$nIV
    )
  }
  
  # Create formatted strings for printing
  summary_obj$formatted <- list(
    bic_result = paste0(
      "CARE (BIC): ", round(bic_estimate, digits),
      " (", round(bic_lower, digits), " to ", round(bic_upper, digits), ")",
      ", p = ", format(bic_p, digits = digits, scientific = TRUE)
    ),
    rivw_result = paste0(
      "RIVW: ", round(rivw_estimate, digits),
      " (", round(rivw_lower, digits), " to ", round(rivw_upper, digits), ")",
      ", p = ", format(rivw_p, digits = digits, scientific = TRUE)
    )
  )
  
  class(summary_obj) <- "summary.mrcare"
  return(summary_obj)
}

#' Print method for summary.mrcare objects
#'
#' @param x A summary.mrcare object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the summary object
#' @export
print.mrcare <- function(x, ...) {
  cat("=======================================================\n")
  cat(" MRcare Analysis Summary\n")
  cat("=======================================================\n\n")
  
  # Print general information
  cat("Analysis timestamp:", x$general$timestamp, "\n")
  cat("Analysis duration:", x$general$duration_mins, "minutes\n")
  cat("Exposure:", x$general$exposure, "\n")
  cat("Outcome:", x$general$outcome, "\n\n")
  
  # Print SNP information
  cat("Genetic instrument information:\n")
  cat("  Total SNPs used:", x$snp_info$total_snps, "\n")
  cat("  Valid instruments (BIC):", x$snp_info$valid_snps, "\n")
  cat("  Invalid instruments (BIC):", x$snp_info$invalid_snps, "\n")
  cat("  Average F-statistic:", x$snp_info$f_statistic, "\n\n")
  
  # Print causal effect estimates
  cat("Causal effect estimates:\n")
  cat("=======================================================\n")
  cat("Method      Estimate   95% CI            P-value\n")
  cat("-------------------------------------------------------\n")
  
  # BIC result
  if (!is.null(x$BIC)) {
    cat(sprintf("CARE (BIC)  %-10s(%s to %s)  %s\n",
               x$BIC$estimate,
               x$BIC$ci_lower,
               x$BIC$ci_upper,
               x$BIC$p_value))
  }
  
  # RIVW result
  if (!is.null(x$RIVW)) {
    cat(sprintf("RIVW        %-10s(%s to %s)  %s\n",
               x$RIVW$estimate,
               x$RIVW$ci_lower,
               x$RIVW$ci_upper,
               x$RIVW$p_value))
  }
  
  cat("=======================================================\n")
  cat("\nAnalysis parameters:\n")
  cat("  P-threshold:", x$general$parameters$p_threshold, "\n")
  cat("  Clumping r2:", x$general$parameters$clump_r2, "\n")
  cat("  Clumping window:", x$general$parameters$clump_kb, "kb\n")
  cat("  Bias correction:", x$general$parameters$bias_correction, "\n")
  cat("  Algorithm:", x$general$parameters$algorithm, "\n")
  cat("  Bootstrap repetitions:", x$general$parameters$nrep, "\n\n")
  
  cat("Use plot_care(), plot_scatter(), or plot_funnel() for visualizations.\n")
  
  invisible(x)
}
