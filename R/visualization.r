#' Plot MR-CARE Results
#'
#' Creates a forest plot of MR-CARE causal estimates (BIC and RIVW) with confidence intervals
#'
#' @param care_results Results object from mr_care function
#' @param exponentiate Whether to exponentiate estimates (e.g., for odds ratios). Default is FALSE
#' @param ci_level Confidence level for intervals. Default is 0.95
#' @param xlab X-axis label. Default is "Causal effect"
#' @param title Plot title. Default is "MRcare Results"
#' @param colors Vector of colors for different methods. Default is NULL (uses preset colors)
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' results <- mr_care(exposure_data, outcome_data)
#' plot_care(results)
#' 
#' # Odds ratios
#' plot_care(results, exponentiate = TRUE, xlab = "Odds ratio")
#' }
plot_care <- function(care_results, 
                     exponentiate = FALSE,
                     ci_level = 0.95,
                     xlab = "Causal effect",
                     title = "MRcare Results",
                     colors = NULL) {
  
  # Calculate z-value for CI
  z <- qnorm(1 - (1 - ci_level) / 2)
  
  # Prepare data frame for plotting
  plot_data <- data.frame(
    Method = character(),
    Estimate = numeric(),
    Lower = numeric(),
    Upper = numeric(),
    P = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Extract BIC results
  if ("BIC" %in% names(care_results)) {
    estimate <- care_results$BIC$estimate
    se <- care_results$BIC$efron_se  # Use Efron bootstrap SE
    p_value <- care_results$BIC$efron_p
    
    lower <- estimate - z * se
    upper <- estimate + z * se
    
    # Exponentiate if requested (for odds ratios)
    if (exponentiate) {
      estimate <- exp(estimate)
      lower <- exp(lower)
      upper <- exp(upper)
    }
    
    plot_data <- rbind(plot_data, data.frame(
      Method = "CARE (BIC)",
      Estimate = estimate,
      Lower = lower,
      Upper = upper,
      P = p_value
    ))
  }
  
  # Extract RIVW results
  if ("RIVW" %in% names(care_results)) {
    estimate <- care_results$RIVW$beta
    se <- care_results$RIVW$sd
    p_value <- care_results$RIVW$p
    
    lower <- estimate - z * se
    upper <- estimate + z * se
    
    # Exponentiate if requested (for odds ratios)
    if (exponentiate) {
      estimate <- exp(estimate)
      lower <- exp(lower)
      upper <- exp(upper)
    }
    
    plot_data <- rbind(plot_data, data.frame(
      Method = "RIVW",
      Estimate = estimate,
      Lower = lower,
      Upper = upper,
      P = p_value
    ))
  }
  
  # Set colors
  if (is.null(colors)) {
    colors <- c("CARE (BIC)" = "#D95F02", "RIVW" = "#1B9E77")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Method, y = Estimate, 
                               color = Method)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = ifelse(exponentiate, 1, 0), 
                       linetype = "dashed", color = "gray50", linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = title, x = "", y = xlab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # Add p-values as text
  p <- p + ggplot2::geom_text(ggplot2::aes(y = Upper + 0.1 * (max(Upper) - min(Lower)),
                               label = sprintf("p = %.3g", P)), vjust = -0.5)
  
  return(p)
}

#' Plot MR Scatter Plot
#'
#' Creates a scatter plot of genetic associations with exposure and outcome,
#' with causal effect estimate lines for BIC and RIVW methods
#'
#' @param care_results Results object from mr_care function
#' @param highlight_invalid Whether to highlight invalid instruments. Default is TRUE
#' @param xlab X-axis label. Default is "Genetic association with exposure"
#' @param ylab Y-axis label. Default is "Genetic association with outcome"
#' @param title Plot title. Default is "MR Scatter Plot"
#' @param colors Vector of colors for valid and invalid instruments. Default is NULL (uses preset colors)
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' results <- mr_care(exposure_data, outcome_data)
#' plot_scatter(results)
#' }
plot_scatter <- function(care_results, 
                        highlight_invalid = TRUE,
                        xlab = "Genetic association with exposure",
                        ylab = "Genetic association with outcome",
                        title = "MR Scatter Plot",
                        colors = NULL) {
  
  # Extract the harmonized data
  data <- care_results$data
  
  # Create the plot data
  plot_data <- data.frame(
    SNP = data$SNP,
    beta_exp = data$beta.exposure,
    beta_out = data$beta.outcome,
    se_exp = data$se.exposure,
    se_out = data$se.outcome
  )
  
  # Add invalid status if available and requested
  if (highlight_invalid && !is.null(care_results$invalid_ivs$BIC)) {
    invalid_vec <- care_results$invalid_ivs$BIC
    plot_data$invalid <- invalid_vec != 0
  } else {
    plot_data$invalid <- FALSE
  }
  
  # Set colors
  if (is.null(colors)) {
    colors <- c("Valid" = "#0072B2", "Invalid" = "#D55E00")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = beta_exp, y = beta_out, 
                                 color = if(highlight_invalid) invalid else "Valid")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta_out - 1.96 * se_out, 
                               ymax = beta_out + 1.96 * se_out), 
                         alpha = 0.5, width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta_exp - 1.96 * se_exp, 
                                xmax = beta_exp + 1.96 * se_exp), 
                          alpha = 0.5, height = 0) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = colors, 
                               labels = c("Valid", "Invalid"),
                               name = "Instrument status") +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # Add BIC estimate line
  if ("BIC" %in% names(care_results)) {
    bic_estimate <- care_results$BIC$estimate
    p <- p + ggplot2::geom_abline(slope = bic_estimate, intercept = 0, 
                                linetype = "solid", color = "#D95F02", linewidth = 1)
    
    # Add the causal estimate as text
    p <- p + ggplot2::annotate("text", 
                            x = max(plot_data$beta_exp) * 0.8, 
                            y = max(plot_data$beta_out) * 0.9,
                            label = sprintf("CARE (BIC) estimate: %.3g", bic_estimate),
                            color = "#D95F02")
  }
  
  # Add RIVW estimate line
  if ("RIVW" %in% names(care_results)) {
    rivw_estimate <- care_results$RIVW$beta
    p <- p + ggplot2::geom_abline(slope = rivw_estimate, intercept = 0, 
                                linetype = "dashed", color = "#1B9E77", linewidth = 1)
    
    # Add the causal estimate as text
    p <- p + ggplot2::annotate("text", 
                            x = max(plot_data$beta_exp) * 0.8, 
                            y = max(plot_data$beta_out) * 0.8,
                            label = sprintf("RIVW estimate: %.3g", rivw_estimate),
                            color = "#1B9E77")
  }
  
  # Add null effect reference line
  p <- p + ggplot2::geom_abline(slope = 0, intercept = 0, 
                              linetype = "dotted", color = "gray50")
  
  return(p)
}

#' Funnel Plot for MR Analysis
#'
#' Creates a funnel plot to visualize potential pleiotropy in MR analysis
#'
#' @param care_results Results object from mr_care function
#' @param highlight_invalid Whether to highlight invalid instruments. Default is TRUE
#' @param xlab X-axis label. Default is "MR estimate"
#' @param ylab Y-axis label. Default is "Precision (1/SE)"
#' @param title Plot title. Default is "MR Funnel Plot"
#' @param colors Vector of colors for valid and invalid instruments. Default is NULL (uses preset colors)
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' results <- mr_care(exposure_data, outcome_data)
#' plot_funnel(results)
#' }
plot_funnel <- function(care_results, 
                       highlight_invalid = TRUE,
                       xlab = "MR estimate",
                       ylab = "Precision (1/SE)",
                       title = "MR Funnel Plot",
                       colors = NULL) {
  
  # Extract the harmonized data
  data <- care_results$data
  
  # Calculate individual SNP causal estimates and their precision
  plot_data <- data.frame(
    SNP = data$SNP,
    beta_exp = data$beta.exposure,
    beta_out = data$beta.outcome,
    se_exp = data$se.exposure,
    se_out = data$se.outcome
  )
  
  # Calculate Wald ratios for each SNP
  plot_data$mr_estimate <- plot_data$beta_out / plot_data$beta_exp
  
  # Calculate standard errors for Wald ratios
  plot_data$mr_se <- sqrt((plot_data$se_out^2 / plot_data$beta_exp^2) + 
                         (plot_data$beta_out^2 * plot_data$se_exp^2 / plot_data$beta_exp^4))
  
  # Calculate precision (1/SE)
  plot_data$precision <- 1 / plot_data$mr_se
  
  # Add invalid status if available and requested
  if (highlight_invalid && !is.null(care_results$invalid_ivs$BIC)) {
    invalid_vec <- care_results$invalid_ivs$BIC
    plot_data$invalid <- invalid_vec != 0
  } else {
    plot_data$invalid <- FALSE
  }
  
  # Set colors
  if (is.null(colors)) {
    colors <- c("Valid" = "#0072B2", "Invalid" = "#D55E00")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mr_estimate, y = precision, 
                                 color = if(highlight_invalid) invalid else "Valid")) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = colors, 
                               labels = c("Valid", "Invalid"),
                               name = "Instrument status") +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # Add BIC estimate line
  if ("BIC" %in% names(care_results)) {
    bic_estimate <- care_results$BIC$estimate
    p <- p + ggplot2::geom_vline(xintercept = bic_estimate, 
                              linetype = "solid", color = "#D95F02", linewidth = 1)
    
    # Add the causal estimate as text
    p <- p + ggplot2::annotate("text", 
                            x = quantile(plot_data$mr_estimate, 0.9), 
                            y = max(plot_data$precision) * 0.9,
                            label = sprintf("CARE (BIC) estimate: %.3g", bic_estimate),
                            color = "#D95F02")
  }
  
  # Add RIVW estimate line
  if ("RIVW" %in% names(care_results)) {
    rivw_estimate <- care_results$RIVW$beta
    p <- p + ggplot2::geom_vline(xintercept = rivw_estimate, 
                              linetype = "dashed", color = "#1B9E77", linewidth = 1)
    
    # Add the causal estimate as text
    p <- p + ggplot2::annotate("text", 
                            x = quantile(plot_data$mr_estimate, 0.9), 
                            y = max(plot_data$precision) * 0.8,
                            label = sprintf("RIVW estimate: %.3g", rivw_estimate),
                            color = "#1B9E77")
  }
  
  return(p)
}

#' Plot Combined MRcare Analysis Results
#'
#' Creates a combined visualization with forest plot, scatter plot, and funnel plot
#'
#' @param care_results Results object from mr_care function
#' @param highlight_invalid Whether to highlight invalid instruments. Default is TRUE
#' @param exponentiate Whether to exponentiate estimates in forest plot. Default is FALSE
#' @param title Main title for the combined plot. Default is "MRcare Analysis Results"
#'
#' @return A grid-based combined plot
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' \dontrun{
#' results <- mr_care(exposure_data, outcome_data)
#' plot_combined(results)
#' }
plot_combined <- function(care_results,
                         highlight_invalid = TRUE,
                         exponentiate = FALSE,
                         title = "MRcare Analysis Results") {
  
  # Create individual plots
  forest_plot <- plot_care(care_results, exponentiate = exponentiate)
  scatter_plot <- plot_scatter(care_results, highlight_invalid = highlight_invalid)
  funnel_plot <- plot_funnel(care_results, highlight_invalid = highlight_invalid)
  
  # Combine plots using gridExtra
  combined_plot <- gridExtra::grid.arrange(
    forest_plot, scatter_plot, funnel_plot,
    ncol = 1,
    heights = c(1, 2, 2),
    top = title
  )
  
  return(combined_plot)
}
