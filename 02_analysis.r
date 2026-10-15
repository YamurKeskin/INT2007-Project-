# Pseudomonas aeruginosa MIC-BPC Project
# Statistical analysis
# Research question:
# Is there an association between MIC and BPC values for tobramycin
# among Pseudomonas aeruginosa isolates?


# 1. Run the data preparation script.
# This creates the cleaned dataset called clean_data.
source("01_prepare_data.R")

# 2. Store the MIC and BPC columns in shorter variables.
# This makes the following calculations easier to read.
mic_values <- clean_data$MIC_tobramycin_ug_ml
bpc_values <- clean_data$BPC_tobramycin_ug_ml

# 3. Check that both variables are numeric.
stopifnot(is.numeric(mic_values))
stopifnot(is.numeric(bpc_values))

# 4. Check for missing values.
stopifnot(!any(is.na(mic_values)))
stopifnot(!any(is.na(bpc_values)))

# 5. Calculate descriptive statistics for MIC.
mic_mean <- mean(mic_values)
mic_median <- median(mic_values)
mic_sd <- sd(mic_values)
mic_se <- mic_sd / sqrt(length(mic_values))
mic_iqr <- IQR(mic_values)
mic_min <- min(mic_values)
mic_max <- max(mic_values)

# 6. Calculate descriptive statistics for BPC
bpc_mean <- mean(bpc_values)
bpc_median <- median(bpc_values)
bpc_sd <- sd(bpc_values)
bpc_se <- bpc_sd / sqrt(length(bpc_values))
bpc_iqr <- IQR(bpc_values)
bpc_min <- min(bpc_values)
bpc_max <- max(bpc_values)


# 7. Combine the results into one summary table.
summary_table <- data.frame(
  Measurement = c("MIC", "BPC"),
  Sample_size = c(length(mic_values), length(bpc_values)),
  Mean = c(mic_mean, bpc_mean),
  Median = c(mic_median, bpc_median),
  Standard_deviation = c(mic_sd, bpc_sd),
  Standard_error = c(mic_se, bpc_se),
  IQR = c(mic_iqr, bpc_iqr),
  Minimum = c(mic_min, bpc_min),
  Maximum = c(mic_max, bpc_max)
)


# 8. Print the summary table in the R terminal.
print(summary_table)

# 9. Calculate the BPC-to-MIC ratio for every isolate.
# A ratio above 1 means that the isolate required a higher
# tobramycin concentration to prevent biofilm formation than to
# inhibit visible planktonic growth.

clean_data$BPC_MIC_ratio <-
  clean_data$BPC_tobramycin_ug_ml /
  clean_data$MIC_tobramycin_ug_ml


# 10. Display the isolate-level ratios.
print(clean_data)

# 11. Summarise the BPC-to-MIC ratio
ratio_summary <- data.frame(
  Mean_ratio = mean(clean_data$BPC_MIC_ratio),
  Median_ratio = median(clean_data$BPC_MIC_ratio),
  SD_ratio = sd(clean_data$BPC_MIC_ratio),
  Minimum_ratio = min(clean_data$BPC_MIC_ratio),
  Maximum_ratio = max(clean_data$BPC_MIC_ratio)
)

print(ratio_summary)


# 12. Save the statistical summaries as CSV files.
write.csv(
  summary_table,
  "descriptive_statistics.csv",
  row.names = FALSE
)

write.csv(
  ratio_summary,
  "ratio_summary.csv",
  row.names = FALSE
)

write.csv(
  clean_data,
  "pseudomonas_MIC_BPC_with_ratio.csv",
  row.names = FALSE
)


