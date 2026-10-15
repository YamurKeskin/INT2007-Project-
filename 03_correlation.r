# Pseudomonas aeruginosa MIC-BPC Project
# Correlation and uncertainty analysis


# 1. Run the data preparation script.
# This creates the clean_data object containing the 30 isolates.

source("01_prepare_data.R")

# 2. Store the MIC and BPC measurements in shorter variables.

mic_values <- clean_data$MIC_tobramycin_ug_ml
bpc_values <- clean_data$BPC_tobramycin_ug_ml

# 3. Check the data before performing the analysis.

stopifnot(length(mic_values) == length(bpc_values))
stopifnot(nrow(clean_data) == 30)
stopifnot(!any(is.na(mic_values)))
stopifnot(!any(is.na(bpc_values)))

# 4. Calculate the strength and direction of the association.
# Spearman correlation is suitable because the concentration values
# are ordered, discrete and contain tied observations.

observed_rho <- suppressWarnings(
  cor(
    mic_values,
    bpc_values,
    method = "spearman"
  )
)

# Stop the analysis if the correlation cannot be calculated.
stopifnot(is.finite(observed_rho))

# 5. Check that the boot package is available.
# This package is used to calculate bootstrap uncertainty intervals.

if (!requireNamespace("boot", quietly = TRUE)) {
  stop(
    "The boot package is required. Run install.packages('boot') once."
  )
}


# 6. Define the statistic that will be calculated in every
# bootstrap sample.

# Complete isolate rows are resampled. This keeps the MIC and BPC
# measurements from each isolate paired together.

spearman_statistic <- function(dataset, sampled_rows) {

  sampled_data <- dataset[sampled_rows, ]

  sampled_rho <- suppressWarnings(
    cor(
      sampled_data$MIC_tobramycin_ug_ml,
      sampled_data$BPC_tobramycin_ug_ml,
      method = "spearman"
    )
  )

  return(sampled_rho)
}

# 7. Set a seed so that the bootstrap analysis is reproducible.
set.seed(123)

# 8. Define the number of bootstrap samples.
bootstrap_repetitions <- 10000

# 9. Generate the paired bootstrap samples.
bootstrap_result <- boot::boot(
  data = clean_data,
  statistic = spearman_statistic,
  R = bootstrap_repetitions
)

# 10. Extract the bootstrap correlation values.
bootstrap_values <- as.numeric(bootstrap_result$t)

# Remove undefined bootstrap values, if any.
bootstrap_values <- bootstrap_values[
  is.finite(bootstrap_values)
]

# Stop the analysis if no valid bootstrap values were produced.
stopifnot(length(bootstrap_values) > 0)

# 11. Calculate summary statistics for the bootstrap distribution.
bootstrap_mean <- mean(bootstrap_values)
bootstrap_median <- median(bootstrap_values)
bootstrap_bias <- bootstrap_mean - observed_rho
bootstrap_standard_error <- sd(bootstrap_values)

# 12. Calculate the 95% percentile and BCa bootstrap intervals.
bootstrap_intervals <- boot::boot.ci(
  bootstrap_result,
  conf = 0.95,
  type = c("perc", "bca")
)

# Extract the lower and upper limits of the percentile interval.
percentile_lower <- bootstrap_intervals$percent[4]
percentile_upper <- bootstrap_intervals$percent[5]

# Extract the lower and upper limits of the BCa interval.
bca_lower <- bootstrap_intervals$bca[4]
bca_upper <- bootstrap_intervals$bca[5]

# 13. Store the main results in a one-row table.
# No p-value is included because the project focuses on estimation,
# association strength and uncertainty.

correlation_results <- data.frame(
  n = nrow(clean_data),
  method = "Spearman correlation",
  spearman_rho = observed_rho,
  bootstrap_method = "BCa paired bootstrap",
  bootstrap_repetitions = bootstrap_repetitions,
  bootstrap_standard_error = bootstrap_standard_error,
  bootstrap_bias = bootstrap_bias,
  percentile_95_lower = percentile_lower,
  percentile_95_upper = percentile_upper,
  bca_95_lower = bca_lower,
  bca_95_upper = bca_upper
)


# 14. Create a second table for validating the bootstrap results.

bootstrap_validation_results <- data.frame(
  statistic = c(
    "Observed Spearman rho",
    "Bootstrap mean rho",
    "Bootstrap median rho",
    "Bootstrap bias",
    "Bootstrap standard error",
    "Percentile 95% lower",
    "Percentile 95% upper",
    "BCa 95% lower",
    "BCa 95% upper"
  ),

  value = c(
    observed_rho,
    bootstrap_mean,
    bootstrap_median,
    bootstrap_bias,
    bootstrap_standard_error,
    percentile_lower,
    percentile_upper,
    bca_lower,
    bca_upper
  )
)


# 15. Save the main correlation and uncertainty results.
write.csv(
  correlation_results,
  "correlation_results.csv",
  row.names = FALSE
)

# 16. Save all bootstrap correlation values.
# This file will be used to create the bootstrap histogram in Python.
write.csv(
  data.frame(spearman_rho = bootstrap_values),
  "bootstrap_correlation_values.csv",
  row.names = FALSE
)

# 17. Save the bootstrap validation results.
write.csv(
  bootstrap_validation_results,
  "bootstrap_validation_results.csv",
  row.names = FALSE
)

# 18. Print the final results in the R terminal.
print(correlation_results)
print(bootstrap_validation_results)

