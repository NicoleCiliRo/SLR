# ==============================================================================
# DATA VALIDATION FOR META-ANALYSIS
# Checks for issues that cause "missing value where TRUE/FALSE needed" error
# ==============================================================================

cat("==============================================================================\n")
cat("DATA VALIDATION FOR META-ANALYSIS\n")
cat("==============================================================================\n\n")

# Load data
cat("Loading niri_meta.csv...\n")
meta_data <- tryCatch({
  read.csv("niri_meta.csv", stringsAsFactors = FALSE)
}, error = function(e) {
  cat("ERROR: Could not read niri_meta.csv\n")
  cat("Make sure the file is in the current directory.\n")
  cat("Error message:", e$message, "\n")
  stop()
})

cat("✓ File loaded successfully\n")
cat(sprintf("  Number of rows: %d\n", nrow(meta_data)))
cat(sprintf("  Number of columns: %d\n\n", ncol(meta_data)))

# ==============================================================================
# CHECK 1: Column names
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("CHECK 1: Column Names\n")
cat("----------------------------------------------------------------------\n")

required_cols <- c("Study_ID", "Setting", "Sensitivity", "Specificity", 
                   "AUC", "n_diseased", "n_non_diseased", "Device", "Reference_Test")

missing_cols <- setdiff(required_cols, names(meta_data))

if (length(missing_cols) > 0) {
  cat("✗ MISSING COLUMNS:\n")
  for (col in missing_cols) {
    cat(sprintf("  - %s\n", col))
  }
  cat("\nCurrent columns:\n")
  print(names(meta_data))
  stop("Fix missing columns before continuing")
} else {
  cat("✓ All required columns present\n\n")
}

# ==============================================================================
# CHECK 2: Missing values
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("CHECK 2: Missing Values\n")
cat("----------------------------------------------------------------------\n")

has_missing <- FALSE
for (col in required_cols) {
  n_missing <- sum(is.na(meta_data[[col]]))
  if (n_missing > 0) {
    cat(sprintf("✗ %s: %d missing values\n", col, n_missing))
    has_missing <- TRUE
  }
}

if (!has_missing) {
  cat("✓ No missing values\n\n")
} else {
  cat("\nRows with missing values:\n")
  missing_rows <- which(apply(is.na(meta_data[required_cols]), 1, any))
  print(meta_data[missing_rows, c("Study_ID", required_cols)])
  stop("\nFix missing values before continuing")
}

# ==============================================================================
# CHECK 3: Data types and ranges
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("CHECK 3: Data Types and Ranges\n")
cat("----------------------------------------------------------------------\n")

issues <- c()

# Sensitivity (should be 0-1)
if (any(meta_data$Sensitivity < 0 | meta_data$Sensitivity > 1, na.rm = TRUE)) {
  if (max(meta_data$Sensitivity, na.rm = TRUE) > 1 && 
      max(meta_data$Sensitivity, na.rm = TRUE) <= 100) {
    cat("✗ Sensitivity appears to be in percentage format (0-100)\n")
    cat("  Convert to proportions by dividing by 100\n")
    issues <- c(issues, "Sensitivity format")
  } else {
    cat("✗ Sensitivity out of valid range (0-1)\n")
    cat(sprintf("  Range: %.3f to %.3f\n", 
                min(meta_data$Sensitivity, na.rm = TRUE),
                max(meta_data$Sensitivity, na.rm = TRUE)))
    issues <- c(issues, "Sensitivity range")
  }
} else {
  cat(sprintf("✓ Sensitivity valid (%.3f to %.3f)\n",
              min(meta_data$Sensitivity, na.rm = TRUE),
              max(meta_data$Sensitivity, na.rm = TRUE)))
}

# Specificity (should be 0-1)
if (any(meta_data$Specificity < 0 | meta_data$Specificity > 1, na.rm = TRUE)) {
  if (max(meta_data$Specificity, na.rm = TRUE) > 1 && 
      max(meta_data$Specificity, na.rm = TRUE) <= 100) {
    cat("✗ Specificity appears to be in percentage format (0-100)\n")
    cat("  Convert to proportions by dividing by 100\n")
    issues <- c(issues, "Specificity format")
  } else {
    cat("✗ Specificity out of valid range (0-1)\n")
    cat(sprintf("  Range: %.3f to %.3f\n",
                min(meta_data$Specificity, na.rm = TRUE),
                max(meta_data$Specificity, na.rm = TRUE)))
    issues <- c(issues, "Specificity range")
  }
} else {
  cat(sprintf("✓ Specificity valid (%.3f to %.3f)\n",
              min(meta_data$Specificity, na.rm = TRUE),
              max(meta_data$Specificity, na.rm = TRUE)))
}

# Sample sizes
if (any(meta_data$n_diseased < 1, na.rm = TRUE)) {
  cat("✗ n_diseased contains values < 1\n")
  issues <- c(issues, "n_diseased invalid")
} else {
  cat(sprintf("✓ n_diseased valid (%d to %d)\n",
              min(meta_data$n_diseased, na.rm = TRUE),
              max(meta_data$n_diseased, na.rm = TRUE)))
}

if (any(meta_data$n_non_diseased < 1, na.rm = TRUE)) {
  cat("✗ n_non_diseased contains values < 1\n")
  issues <- c(issues, "n_non_diseased invalid")
} else {
  cat(sprintf("✓ n_non_diseased valid (%d to %d)\n",
              min(meta_data$n_non_diseased, na.rm = TRUE),
              max(meta_data$n_non_diseased, na.rm = TRUE)))
}

# ==============================================================================
# CHECK 4: Calculate 2x2 table values and check for issues
# ==============================================================================

cat("\n----------------------------------------------------------------------\n")
cat("CHECK 4: 2x2 Contingency Tables (THIS IS THE CRITICAL CHECK!)\n")
cat("----------------------------------------------------------------------\n\n")

meta_data$TP <- round(meta_data$Sensitivity * meta_data$n_diseased)
meta_data$FN <- meta_data$n_diseased - meta_data$TP
meta_data$TN <- round(meta_data$Specificity * meta_data$n_non_diseased)
meta_data$FP <- meta_data$n_non_diseased - meta_data$TN

cat("Study-by-study validation:\n\n")

problem_studies <- c()

for (i in 1:nrow(meta_data)) {
  study <- meta_data[i, ]
  study_id <- study$Study_ID
  
  cat(sprintf("Study %d: %s\n", i, study_id))
  cat(sprintf("  Disease+: %d (TP=%d, FN=%d)\n", 
              study$n_diseased, study$TP, study$FN))
  cat(sprintf("  Disease-: %d (TN=%d, FP=%d)\n",
              study$n_non_diseased, study$TN, study$FP))
  
  # Check for negative values
  if (study$TP < 0 || study$FN < 0 || study$TN < 0 || study$FP < 0) {
    cat("  ✗ PROBLEM: Negative cell count!\n")
    problem_studies <- c(problem_studies, study_id)
  }
  
  # Check for all zeros in diseased
  if (study$TP == 0 && study$FN == 0) {
    cat("  ✗ PROBLEM: Zero diseased patients (TP=0, FN=0)\n")
    problem_studies <- c(problem_studies, study_id)
  }
  
  # Check for all zeros in non-diseased
  if (study$TN == 0 && study$FP == 0) {
    cat("  ✗ PROBLEM: Zero non-diseased patients (TN=0, FP=0)\n")
    problem_studies <- c(problem_studies, study_id)
  }
  
  # Check for zero cells (will need correction)
  zero_cells <- sum(c(study$TP, study$FN, study$TN, study$FP) == 0)
  if (zero_cells > 0) {
    cat(sprintf("  ⚠ WARNING: %d zero cell(s) - will apply 0.5 correction\n", zero_cells))
  }
  
  # Check for perfect sensitivity or specificity
  if (study$Sensitivity == 1.0) {
    cat("  ⚠ NOTE: Perfect sensitivity (FN=0)\n")
  }
  if (study$Specificity == 1.0) {
    cat("  ⚠ NOTE: Perfect specificity (FP=0)\n")
  }
  if (study$Sensitivity == 0.0) {
    cat("  ✗ PROBLEM: Zero sensitivity (TP=0)\n")
    problem_studies <- c(problem_studies, study_id)
  }
  if (study$Specificity == 0.0) {
    cat("  ✗ PROBLEM: Zero specificity (TN=0)\n")
    problem_studies <- c(problem_studies, study_id)
  }
  
  cat("\n")
}

# ==============================================================================
# CHECK 5: Subgroup sample sizes
# ==============================================================================

cat("----------------------------------------------------------------------\n")
cat("CHECK 5: Subgroup Sample Sizes\n")
cat("----------------------------------------------------------------------\n\n")

settings <- unique(meta_data$Setting)
cat("Settings found:\n")
for (setting in settings) {
  n_studies <- sum(meta_data$Setting == setting)
  cat(sprintf("  %s: %d studies", setting, n_studies))
  if (n_studies < 2) {
    cat(" ⚠ WARNING: Need at least 2 studies per subgroup")
  }
  cat("\n")
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("VALIDATION SUMMARY\n")
cat("==============================================================================\n\n")

if (length(problem_studies) > 0) {
  cat("✗ CRITICAL PROBLEMS FOUND!\n\n")
  cat("Problem studies:\n")
  for (study in unique(problem_studies)) {
    cat(sprintf("  - %s\n", study))
  }
  cat("\nThese studies have issues that will cause the meta-analysis to fail.\n")
  cat("\nCommon fixes:\n")
  cat("  1. Check for data entry errors\n")
  cat("  2. Remove studies with impossible values (e.g., Sen=0 or Spec=0)\n")
  cat("  3. Verify sensitivity/specificity are proportions (0-1, not 0-100)\n")
  cat("  4. Ensure n_diseased and n_non_diseased are correct\n\n")
  
  cat("Detailed data for problem studies:\n")
  problem_idx <- which(meta_data$Study_ID %in% unique(problem_studies))
  print(meta_data[problem_idx, c("Study_ID", "Sensitivity", "Specificity", 
                                  "n_diseased", "n_non_diseased", 
                                  "TP", "FN", "TN", "FP")])
  
  stop("\n\nFIX THE PROBLEMS ABOVE BEFORE RUNNING META-ANALYSIS")
  
} else if (length(issues) > 0) {
  cat("⚠ WARNINGS FOUND:\n")
  for (issue in issues) {
    cat(sprintf("  - %s\n", issue))
  }
  cat("\nFix these warnings before proceeding.\n")
  
} else {
  cat("✓✓✓ ALL CHECKS PASSED! ✓✓✓\n\n")
  cat("Your data is valid and ready for meta-analysis.\n\n")
  cat("Next step: Run meta_analysis.R\n")
  cat("  In RStudio: Open meta_analysis.R and click 'Source'\n")
  cat("  In terminal: Rscript meta_analysis.R\n\n")
}

cat("==============================================================================\n")
