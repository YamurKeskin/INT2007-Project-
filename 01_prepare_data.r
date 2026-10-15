# Pseudomonas aeruginosa MIC-BPC project
# This script creates one clean row per clinical isolate.

# 1. Read the two original public-data files.
# check.names = FALSE keeps the original column names unchanged.
mic_raw <- read.csv("MIC_IMC.csv", check.names = FALSE)
bpc_raw <- read.csv("BPC_IMC.csv", check.names = FALSE)

# 2. Confirm that the expected target columns are present.
stopifnot("MIC tobramycin" %in% names(mic_raw))

# The public BPC file may retain an older internal column label.
# Identify its tobramycin column without carrying that legacy label into the
# cleaned dataset.
bpc_tobramycin_columns <- grep(
  "tobramycin",
  names(bpc_raw),
  ignore.case = TRUE,
  value = TRUE
)

stopifnot(length(bpc_tobramycin_columns) == 1)
bpc_tobramycin_column <- bpc_tobramycin_columns[[1]]

# 3. Keep only isolate IDs and the two measurements needed for this project.
mic <- data.frame(
  isolate_id = mic_raw[[1]],
  MIC_tobramycin_ug_ml = mic_raw[["MIC tobramycin"]]
)

bpc <- data.frame(
  isolate_id = bpc_raw[[1]],
  BPC_tobramycin_ug_ml = bpc_raw[[bpc_tobramycin_column]]
)

# 4. Check the raw structure: 60 rows should represent 30 isolate IDs.
stopifnot(nrow(mic) == 60, nrow(bpc) == 60)
stopifnot(length(unique(mic$isolate_id)) == 30)
stopifnot(length(unique(bpc$isolate_id)) == 30)

# 5. Check that technical replicates have the same target value.
# If any isolate has conflicting values, the script stops instead of silently
# deleting a potentially important disagreement.
mic_values_per_isolate <- tapply(
  mic$MIC_tobramycin_ug_ml,
  mic$isolate_id,
  function(x) length(unique(x))
)

bpc_values_per_isolate <- tapply(
  bpc$BPC_tobramycin_ug_ml,
  bpc$isolate_id,
  function(x) length(unique(x))
)

stopifnot(all(mic_values_per_isolate == 1))
stopifnot(all(bpc_values_per_isolate == 1))

# 6. Keep one row per isolate only after the replicate checks pass.
mic_unique <- mic[!duplicated(mic$isolate_id), ]
bpc_unique <- bpc[!duplicated(bpc$isolate_id), ]

# 7. Join by isolate ID rather than relying on row order.
clean_data <- merge(mic_unique, bpc_unique, by = "isolate_id")

# 8. Restore the natural isolate order: Cl-1, Cl-2, ..., Cl-30.
isolate_number <- as.integer(sub("Cl-", "", clean_data$isolate_id))
clean_data <- clean_data[order(isolate_number), ]
rownames(clean_data) <- NULL

# 9. Final quality checks.
stopifnot(nrow(clean_data) == 30)
stopifnot(!anyDuplicated(clean_data$isolate_id))
stopifnot(!anyNA(clean_data))

# 10. Save the analysis-ready dataset.
write.csv(
  clean_data,
  "pseudomonas_MIC_BPC_clean.csv",
  row.names = FALSE
)

# Display the result and a short summary in the console.
print(clean_data)
summary(clean_data)
