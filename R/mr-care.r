#' CARE Mendelian Randomization Analysis
#'
#' Performs Mendelian randomization analysis using the CARE (Causal Analysis Robust to plEiotropy) method
#' which is designed to be robust to horizontal pleiotropy and winner's curse.
#'
#' @param exposure_data A data frame containing exposure GWAS summary statistics, or a path to a file containing exposure data
#' @param outcome_data A data frame containing outcome GWAS summary statistics, or a path to a file containing outcome data 
#' @param clump_r2 Clumping r-squared threshold. Default is 0.001
#' @param clump_kb Clumping window size in kb. Default is 10000
#' @param p_threshold P-value threshold for selecting significant SNPs. Default is 5e-5
#' @param ref_panel Path to reference panel for LD clumping. If NULL, will use configured reference panel
#' @param bias_correction Method for correcting winner's curse bias: "no", "direct", or "rerand" (re-randomization). Default is "rerand"
#' @param algorithm Algorithm for analysis: "Lasso" or "CD" (coordinate descent). Default is "CD"
#' @param nrep Number of bootstrap repetitions. Default is 5000
#' @param random_start Number of random starting points for optimization. Default is 0
#' @param etamean Mean of the random variable for re-randomization. Default is 0.5
#' @param temp_dir Directory for temporary files. Default is automatically created in a system temp location
#' @param verbose Whether to print detailed progress messages. Default is TRUE
#' @param plink_path Path to plink executable. If NULL, will use configured plink path
#' @param config Use configuration from configure_mrcare()? Default is TRUE
#'
#' @return A list containing the CARE results with various bias correction methods and model selection approaches
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using file paths with configured environment
#' care_results <- mr_care(
#'   exposure_data = "path/to/exposure_gwas.txt",
#'   outcome_data = "path/to/outcome_gwas.txt"
#' )
#'
#' # Using data frames with custom configuration
#' care_results <- mr_care(
#'   exposure_data = exposure_df,
#'   outcome_data = outcome_df,
#'   config = FALSE,
#'   ref_panel = "/path/to/reference/panel",
#'   plink_path = "/path/to/plink"
#' )
#' }
mr_care <- function(exposure_data, 
                   outcome_data,
                   clump_r2 = 0.001,
                   clump_kb = 10000,
                   p_threshold = 5e-5,
                   ref_panel = NULL,
                   bias_correction = c("rerand", "direct", "no"),
                   algorithm = c("CD", "Lasso"),
                   nrep = 5000,
                   random_start = 0,
                   etamean = 0.5,
                   temp_dir = NULL,
                   verbose = TRUE,
                   plink_path = NULL,
                   config = TRUE) {
  
  # Start timing for performance tracking
  start_time <- Sys.time()
  
  # Banner for user feedback
  if (verbose) {
    cat("\n===================================================================\n")
    cat(" MRcare: Causal Analysis Robust to plEiotropy for Mendelian Randomization\n")
    cat("===================================================================\n\n")
    cat("Analysis started at:", format(start_time), "\n\n")
  }
  
  # Parameter validation
  bias_correction <- match.arg(bias_correction)
  algorithm <- match.arg(algorithm)
  
  # Get configuration if requested
  if (config) {
    if (verbose) cat("Loading MRcare configuration...\n")
    tryCatch({
      mrcare_config <- get_mrcare_config(create_if_missing = FALSE, verbose = FALSE)
      
      if (is.null(mrcare_config)) {
        stop("MRcare configuration not found. Please run configure_mrcare() before running mr_care().")
      }
      
      # Rest of your configuration code
    }, error = function(e) {
      stop("Configuration error: ", e$message,
           "\nPlease run configure_mrcare() before using mr_care().")
    })
    
    # Check if we got a valid configuration
    if (is.null(mrcare_config)) {
      if (verbose) cat("No MRcare configuration found. Will use user-provided settings or defaults.\n")
    }

    # Use configuration values if not explicitly provided
    if (is.null(plink_path)) plink_path <- mrcare_config$plink_path
    if (is.null(ref_panel)) ref_panel <- mrcare_config$ref_panel_path
    
    if (verbose) {
      cat("Using configuration:\n")
      cat(" - PLINK: ", plink_path, "\n")
      cat(" - Reference panel: ", ref_panel, "\n")
    }
  } else {
    if (verbose) cat("Using user-provided settings instead of saved configuration.\n")
  }

  # Verify essential parameters
  if (is.null(plink_path)) {
    stop("PLINK path is required. Either provide plink_path parameter or run configure_mrcare() first.")
  }
  
  if (is.null(ref_panel)) {
    stop("Reference panel path is required. Either provide ref_panel parameter or run configure_mrcare() first.")
  }
  
  # Check if plink and reference panel exist
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }
  
  if (!dir.exists(ref_panel)) {
    stop("Reference panel directory not found at: ", ref_panel)
  }
  
  # Create temporary directory if not provided
  if (is.null(temp_dir)) {
    temp_dir <- file.path(tempdir(), "MRcare_temp", format(Sys.time(), "%Y%m%d_%H%M%S"))
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE)
    }
    if (verbose) cat("Created temporary directory at:", temp_dir, "\n")
  } else {
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE)
      if (verbose) cat("Created temporary directory at:", temp_dir, "\n")
    }
  }
  
  # Process input data - determine if it's a path or a data frame
  # TODO: let's improve and revise this part. I think we may just use the existing ones.
  # another logic is ask the users to preprocess and store the processed data.
  # then here, we only read the processed version for the following steps.
  # I think here, we may overly complicated. Let's revise and try if the package works.
  # Load preprocessed data files
  if (verbose) cat("Loading preprocessed exposure data...\n")
  
  if (is.character(exposure_data)) {
    if (!file.exists(exposure_data)) {
      stop("Exposure data file not found: ", exposure_data)
    }
    if (!grepl("\\.txt$", exposure_data, ignore.case = TRUE)) {
      stop("Exposure data must be a txt file created with preprocess_gwas_data()")
    }
    exposure_processed <- as.data.frame(data.table::fread(exposure_data))
  } else if (is.data.frame(exposure_data)) {
    warning("Providing exposure data as a data frame is deprecated. Please use preprocess_gwas_data() first.")
    exposure_processed <- exposure_data
  } else {
    stop("Exposure data must be a path to a txt file created with preprocess_gwas_data()")
  }
  
  if (verbose) cat("Loading preprocessed outcome data...\n")
  
  if (is.character(outcome_data)) {
    if (!file.exists(outcome_data)) {
      stop("Outcome data file not found: ", outcome_data)
    }
    if (!grepl("\\.txt$", outcome_data, ignore.case = TRUE)) {
      stop("Outcome data must be a txt file created with preprocess_gwas_data()")
    }
    outcome_processed <- as.data.frame(data.table::fread(outcome_data))
  } else if (is.data.frame(outcome_data)) {
    warning("Providing outcome data as a data frame is deprecated. Please use x`() first.")
    outcome_processed <- outcome_data
  } else {
    stop("Outcome data must be a path to a txt file created with preprocess_gwas_data()")
  }
  

  
  # Get sample sizes
  if ("N" %in% colnames(exposure_processed)) {
    # Add logic if most N is NA:
    if(sum(is.na(exposure_processed$N)) > nrow(exposure_processed) * 0.5) {
      stop("More than 50% of the sample size for exposure data is NA. Consider providing the sample size directly and reprocessing the data.")   
    }

    nx <- floor(max(exposure_processed$N, na.rm = TRUE))
  } else {
    stop("Sample size not found in exposure data. Using NA.")
    nx <- NA
  }
  
  if ("N" %in% colnames(outcome_processed)) {
    if(sum(is.na(outcome_processed$N)) > nrow(outcome_processed) * 0.5) {
      stop("More than 50% of the sample size for outcome data is NA. Consider providing the sample size directly and reprocessing the data.")   
    }
    ny <- floor(max(outcome_processed$N, na.rm = TRUE))
  } else {
    stop("Sample size not found in outcome data. Using NA.")
    ny <- NA
  }


  if (verbose) {
    cat(paste("Exposure sample size:", format(nx, big.mark=","), "\n"))
    cat(paste("Outcome sample size:", format(ny, big.mark=","), "\n"))
    cat(paste("Total SNPs in exposure data:", format(nrow(exposure_processed), big.mark=","), "\n"))
    cat(paste("Total SNPs in outcome data:", format(nrow(outcome_processed), big.mark=","), "\n\n"))
    
    cat("Extracting independent variants...\n")
  }
  
  # Extract independent variants
  independent_vars <- extract_independent_variants(
    exposure_data = exposure_processed,
    p_threshold = p_threshold,
    etamean = etamean,
    ref_panel = ref_panel,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    temp_dir = temp_dir,
    verbose = verbose,
    rerand = (bias_correction == "rerand"),
    plink_path = plink_path
  )
  
  if (length(independent_vars) < 3) {
    stop("Less than 3 independent SNPs identified. CARE analysis requires at least 3 SNPs.")
  }
  
  if (verbose) {
    cat(paste("Selected", length(independent_vars), "independent variants\n\n"))
    cat("Harmonizing exposure and outcome data...\n")
  }
    # Prepare data for harmonization
  exposure_subset <- exposure_processed[exposure_processed$SNP %in% independent_vars, ]
  
  rm(exposure_processed)
  
  # Create data frames in TwoSampleMR format
  exposure_formatted <- data.frame(
    SNP = exposure_subset$SNP,
    beta.exposure = exposure_subset$BETA,
    se.exposure = exposure_subset$SE,
    effect_allele.exposure = exposure_subset$A1,
    other_allele.exposure = exposure_subset$A2,
    eaf.exposure = exposure_subset$EAF,
    pval.exposure = exposure_subset$P,
    id.exposure = "exposure",
    exposure = "exposure",
    chr.exposure = exposure_subset$CHR,
    pos.exposure = exposure_subset$POS,
    samplesize.exposure = exposure_subset$N,
    stringsAsFactors = FALSE
  )
  
  outcome_formatted <- data.frame(
    SNP = outcome_processed$SNP,
    beta.outcome = outcome_processed$BETA,
    se.outcome = outcome_processed$SE,
    effect_allele.outcome = outcome_processed$A1,
    other_allele.outcome = outcome_processed$A2,
    eaf.outcome = outcome_processed$EAF,
    pval.outcome = outcome_processed$P,
    id.outcome = "outcome",
    outcome = "outcome",
    chr.outcome = outcome_processed$CHR,
    pos.outcome = outcome_processed$POS,
    samplesize.outcome = outcome_processed$N,
    stringsAsFactors = FALSE
  )
  
  rm(outcome_processed)

  # Harmonize the exposure and outcome data using TwoSampleMR
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop("TwoSampleMR package is required for data harmonization. Please install it.")
  }
  
  harmonized_data <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_formatted,
    outcome_dat = outcome_formatted
  )
  
  # Filter for harmonized SNPs that passed the checks
  if ("mr_keep" %in% colnames(harmonized_data)) {
    harmonized_data <- harmonized_data[harmonized_data$mr_keep, ]
  }
  
  if (nrow(harmonized_data) < 3) {
    stop("Less than 3 SNPs remain after harmonization. CARE analysis requires at least 3 SNPs.")
  }
  
  if (verbose) {
    cat(paste("Harmonized dataset contains", nrow(harmonized_data), "SNPs\n\n"))
    cat("Performing CARE analysis...\n")
  }
  
  # Extract the necessary vectors for analysis
  beta_exp <- harmonized_data$beta.exposure
  beta_out <- harmonized_data$beta.outcome
  se_exp <- harmonized_data$se.exposure
  se_out <- harmonized_data$se.outcome
   

  # Run CARE2_boot function with appropriate parameters
  if (verbose) {
    cat("\nStarting bootstrap analysis with", nrep, "replications...\n")
    cat("Using", algorithm, "algorithm with", 
        ifelse(bias_correction == "rerand", "re-randomization", 
               ifelse(bias_correction == "direct", "direct correction", "no correction")), 
        "bias correction...\n")
  }
  
  # Run CARE analysis
  care_results <- CARE2_boot(
    gamma.exp_sel = beta_exp,
    gamma.out_sel = beta_out,
    se.exp_sel = se_exp,
    se.out_sel = se_out,
    nx = nx,
    ny = ny,
    nrep = nrep,
    algorithm = algorithm,
    random_start = random_start,
    biascorrect = bias_correction,
    etamean = etamean,
    pthr = p_threshold
  )
  
  if (verbose) {
    cat("\nCARE analysis complete!\n")
  }
  
  # Run RIVW as an additional method
  if (verbose) cat("\nCalculating RIVW estimate...\n")
  rivw_result <- RIVW(
    gamma.exp_sel = beta_exp, 
    gamma.out_sel = beta_out, 
    se.exp_sel = se_exp, 
    se.out_sel = se_out,
    etamean = etamean,
    sel.pthr = p_threshold
  )
  
  # Format and return results
  if (verbose) cat("\nFormatting results...\n")
  formatted_results <- format_care_results(care_results, harmonized_data, rivw_result)
  
  # Complete timing
  end_time <- Sys.time()
  analysis_duration <- difftime(end_time, start_time, units = "mins")
  
  if (verbose) {
    cat("\n===================================================================\n")
    cat(" Analysis completed successfully\n")
    cat("===================================================================\n")
    cat("Total time:", round(as.numeric(analysis_duration), 2), "minutes\n")
    
    # Print a brief summary of results
    cat("\nCAUSAL EFFECT ESTIMATES:\n")
    cat("Method\t\tEstimate\tS.E.\t\tP-value\n")
    cat("----------------------------------------------------------------\n")
    cat(sprintf("CARE (BIC):\t%.4f\t\t%.4f\t\t%.4e\n", 
               formatted_results$BIC$estimate, 
               formatted_results$BIC$efron_se, 
               formatted_results$BIC$efron_p))
    cat(sprintf("RIVW:\t\t%.4f\t\t%.4f\t\t%.4e\n", 
               formatted_results$RIVW$beta, 
               formatted_results$RIVW$sd, 
               formatted_results$RIVW$p))
    cat("----------------------------------------------------------------\n")
    cat("For detailed results, examine the returned object\n\n")
  }
  
  # Add analysis metadata
  formatted_results$metadata <- list(
    exposure_name = ifelse(is.character(exposure_data) && file.exists(exposure_data), 
                          basename(exposure_data), "exposure"),
    outcome_name = ifelse(is.character(outcome_data) && file.exists(outcome_data), 
                         basename(outcome_data), "outcome"),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    duration_mins = as.numeric(analysis_duration),
    parameters = list(
      p_threshold = p_threshold,
      clump_r2 = clump_r2,
      clump_kb = clump_kb,
      bias_correction = bias_correction,
      algorithm = algorithm,
      etamean = etamean,
      nrep = nrep
    )
  )
  
  # Class for print method
  class(formatted_results) <- c("mrcare", class(formatted_results))
  if (!inherits(formatted_results, "mrcare")) {
    warning("Class 'mrcare' not properly set. Summary function may not work correctly.")
  }
  
  # Clean up temporary files if created automatically
  if (is.null(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
    if (verbose) cat("Cleaned up temporary files\n")
  }
  
  return(formatted_results)
}
