# ==============================================================================
# COMPLETE DIAGNOSTIC TEST ACCURACY META-ANALYSIS
# Includes: Sensitivity, Specificity, PLR, NLR, DOR with forest plots
# ==============================================================================

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos='http://cran.rstudio.com/')
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

cat("Installing/loading required packages...\n")
install_if_missing("meta")
install_if_missing("mada")
install_if_missing("metafor")

cat("\n==============================================================================\n")
cat("DIAGNOSTIC TEST ACCURACY META-ANALYSIS - COMPLETE\n")
cat("==============================================================================\n\n")

# Load data
cat("Loading data from niri_meta.csv...\n")
meta_data <- read.csv("niri_meta.csv", stringsAsFactors = FALSE)
cat("Studies loaded:", nrow(meta_data), "\n")
cat("\nStudies by setting:\n")
print(table(meta_data$Setting))

# Calculate TP, FN, FP, TN
meta_data$TP <- round(meta_data$Sensitivity * meta_data$n_diseased)
meta_data$FN <- meta_data$n_diseased - meta_data$TP
meta_data$TN <- round(meta_data$Specificity * meta_data$n_non_diseased)
meta_data$FP <- meta_data$n_non_diseased - meta_data$TN

# Calculate Likelihood Ratios and DOR
meta_data$PLR <- meta_data$Sensitivity / (1 - meta_data$Specificity)
meta_data$NLR <- (1 - meta_data$Sensitivity) / meta_data$Specificity
meta_data$DOR <- meta_data$PLR / meta_data$NLR

cat("\nData preparation complete!\n")

# ==============================================================================
# BIVARIATE META-ANALYSIS
# ==============================================================================

cat("\n==============================================================================\n")
cat("PERFORMING BIVARIATE META-ANALYSIS (REITSMA MODEL)\n")
cat("==============================================================================\n\n")

fit_reitsma <- reitsma(meta_data, correction = 0.5, correction.control = "all")
summary_reitsma <- summary(fit_reitsma)

# Extract pooled estimates
pooled_sens_logit <- summary_reitsma$coefficients[1, 1]
pooled_spec_logit <- summary_reitsma$coefficients[2, 1]
se_sens <- summary_reitsma$coefficients[1, 2]
se_spec <- summary_reitsma$coefficients[2, 2]

pooled_sens <- plogis(pooled_sens_logit)
pooled_spec <- plogis(pooled_spec_logit)
sens_ci_lower <- plogis(pooled_sens_logit - 1.96 * se_sens)
sens_ci_upper <- plogis(pooled_sens_logit + 1.96 * se_sens)
spec_ci_lower <- plogis(pooled_spec_logit - 1.96 * se_spec)
spec_ci_upper <- plogis(pooled_spec_logit + 1.96 * se_spec)

# Calculate pooled PLR, NLR, DOR
pooled_PLR <- pooled_sens / (1 - pooled_spec)
pooled_NLR <- (1 - pooled_sens) / pooled_spec
pooled_DOR <- pooled_PLR / pooled_NLR

cat("OVERALL POOLED ESTIMATES:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("Pooled Sensitivity: %.3f (95%% CI: %.3f - %.3f) = %.1f%%\n",
            pooled_sens, sens_ci_lower, sens_ci_upper, pooled_sens*100))
cat(sprintf("Pooled Specificity: %.3f (95%% CI: %.3f - %.3f) = %.1f%%\n",
            pooled_spec, spec_ci_lower, spec_ci_upper, pooled_spec*100))
cat(sprintf("Pooled PLR: %.2f\n", pooled_PLR))
cat(sprintf("Pooled NLR: %.2f\n", pooled_NLR))
cat(sprintf("Pooled DOR: %.2f\n", pooled_DOR))
cat("----------------------------------------------------------------------\n")

# ==============================================================================
# FOREST PLOTS - SENSITIVITY & SPECIFICITY
# ==============================================================================

cat("\n==============================================================================\n")
cat("CREATING FOREST PLOTS\n")
cat("==============================================================================\n\n")

# Sensitivity
cat("Creating sensitivity forest plot...\n")
sens_meta <- metaprop(
  event = TP, n = n_diseased, studlab = Study_ID,
  data = meta_data, sm = "PLOGIT", method.tau = "ML",
  subgroup = Setting, random = TRUE, common = FALSE
)

png("forest_plot_sensitivity.png", width=14, height=10, units="in", res=300)
meta::forest(sens_meta,
            rightcols = c("effect", "ci"),
            rightlabs = c("Sensitivity", "95% CI"),
            leftcols = c("studlab", "event", "n"),
            leftlabs = c("Study", "Events", "Total"),
            xlab = "Sensitivity", xlim = c(0, 1),
            main = "Forest Plot - Sensitivity by Setting",
            common = FALSE, random = TRUE,
            col.diamond = "maroon", fontsize = 10, spacing = 1.2)

dev.off()
cat("  ✓ Saved: forest_plot_sensitivity.png\n")

# Specificity
cat("Creating specificity forest plot...\n")
spec_meta <- metaprop(
  event = TN, n = n_non_diseased, studlab = Study_ID,
  data = meta_data, sm = "PLOGIT", method.tau = "ML",
  subgroup = Setting, random = TRUE, common = FALSE
)

png("forest_plot_specificity.png", width=14, height=10, units="in", res=300)
meta::forest(spec_meta,
            rightcols = c("effect", "ci"),
            rightlabs = c("Specificity", "95% CI"),
            leftcols = c("studlab", "event", "n"),
            leftlabs = c("Study", "Events", "Total"),
            xlab = "Specificity", xlim = c(0, 1),
            main = "Forest Plot - Specificity by Setting",
            common = FALSE, random = TRUE,
            col.diamond = "maroon", fontsize = 10, spacing = 1.2)

dev.off()
cat("  ✓ Saved: forest_plot_specificity.png\n")

# ==============================================================================
# FOREST PLOTS - POSITIVE LIKELIHOOD RATIO (PLR)
# ==============================================================================

cat("Creating Positive Likelihood Ratio (PLR) forest plot...\n")

# Meta-analysis of PLR using log transformation
plr_meta <- metagen(
  TE = log(PLR),
  seTE = sqrt(1/TP + 1/(n_diseased - TP) + 1/FP + 1/(n_non_diseased - FP)),
  studlab = Study_ID,
  data = meta_data,
  sm = "PLR",
  method.tau = "ML",
  subgroup = Setting,
  random = TRUE,
  common = FALSE
)

png("forest_plot_PLR.png", width=14, height=10, units="in", res=300)
meta::forest(plr_meta,
            rightcols = c("effect", "ci"),
            rightlabs = c("PLR", "95% CI"),
            leftcols = c("studlab"),
            leftlabs = c("Study"),
            xlab = "Positive Likelihood Ratio",
            main = "Positive Likelihood Ratio (PLR)",
            common = FALSE, random = TRUE,
            col.diamond = "maroon", fontsize = 10, spacing = 1.2,
            xlim = c(0.1, 100), at = c(0.1, 1, 10, 100))

dev.off()
cat("  ✓ Saved: forest_plot_PLR.png\n")

# ==============================================================================
# FOREST PLOTS - NEGATIVE LIKELIHOOD RATIO (NLR)
# ==============================================================================

cat("Creating Negative Likelihood Ratio (NLR) forest plot...\n")

# Meta-analysis of NLR using log transformation
nlr_meta <- metagen(
  TE = log(NLR),
  seTE = sqrt(1/(n_diseased - TP) + 1/TP + 1/TN + 1/(n_non_diseased - TN)),
  studlab = Study_ID,
  data = meta_data,
  sm = "NLR",
  method.tau = "ML",
  subgroup = Setting,
  random = TRUE,
  common = FALSE
)

png("forest_plot_NLR.png", width=14, height=10, units="in", res=300)
meta::forest(nlr_meta,
            rightcols = c("effect", "ci"),
            rightlabs = c("NLR", "95% CI"),
            leftcols = c("studlab"),
            leftlabs = c("Study"),
            xlab = "Negative Likelihood Ratio",
            main = "Negative Likelihood Ratio (NLR)",
            common = FALSE, random = TRUE,
            col.diamond = "maroon", fontsize = 10, spacing = 1.2,
            xlim = c(0.01, 10), at = c(0.01, 0.1, 1, 10))

dev.off()
cat("  ✓ Saved: forest_plot_NLR.png\n")

# ==============================================================================
# FOREST PLOTS - DIAGNOSTIC ODDS RATIO (DOR)
# ==============================================================================

cat("Creating Diagnostic Odds Ratio (DOR) forest plot...\n")

# Calculate DOR directly from 2x2 table
# DOR = (TP * TN) / (FP * FN)
# Apply 0.5 continuity correction for zero cells
meta_data$TP_corr <- ifelse(meta_data$TP == 0, 0.5, meta_data$TP)
meta_data$TN_corr <- ifelse(meta_data$TN == 0, 0.5, meta_data$TN)
meta_data$FP_corr <- ifelse(meta_data$FP == 0, 0.5, meta_data$FP)
meta_data$FN_corr <- ifelse(meta_data$FN == 0, 0.5, meta_data$FN)

meta_data$DOR_calc <- (meta_data$TP_corr * meta_data$TN_corr) /
                      (meta_data$FP_corr * meta_data$FN_corr)

dor_meta <- metagen(
  TE = log(DOR_calc),
  seTE = sqrt(1/TP_corr + 1/TN_corr + 1/FP_corr + 1/FN_corr),
  studlab = Study_ID,
  data = meta_data,
  sm = "OR",
  method.tau = "DL",
  subgroup = Setting,
  random = TRUE,
  common = FALSE
)

png("forest_plot_DOR.png", width=14, height=10, units="in", res=300)
meta::forest(dor_meta,
            rightcols = c("effect", "ci"),
            rightlabs = c("DOR", "95% CI"),
            leftcols = c("studlab"),
            leftlabs = c("Study"),
            xlab = "Diagnostic Odds Ratio",
            main = "Diagnostic Odds Ratio (DOR)",
            common = FALSE, random = TRUE,
            col.diamond = "maroon", fontsize = 10, spacing = 1.2,
            xlim = c(0.1, 1000), at = c(0.1, 1, 10, 100, 1000))

dev.off()
cat("  ✓ Saved: forest_plot_DOR.png\n")

# ==============================================================================
# SROC CURVE
# ==============================================================================

cat("\nCreating SROC curve...\n")

png("sroc_curve.png", width=12, height=10, units="in", res=300)
par(mar=c(5,5,4,2))

plot(fit_reitsma,
  level    = 0.95,
  plotsumm = TRUE,
  predict  = TRUE,       # <-- 95% prediction region
  sroclwd  = 2,
  sroclty  = 1,
  predlty  = 3,
  predlwd  = 1.2,
  pch      = 21,
  cex      = 1.2,
  main     = "Summary ROC with 95% confidence and prediction regions",
)
settings <- unique(meta_data$Setting)
colors <- c("blue", "red")
for (i in 1:length(settings)) {
  setting_data <- meta_data[meta_data$Setting == settings[i], ]
  fpr <- 1 - setting_data$Specificity
  tpr <- setting_data$Sensitivity
  points(fpr, tpr, pch = 19, col = colors[i], cex = 1.5)
}

legend("bottomright", legend = settings, col = colors, pch = 19, cex = 1.2, title = "Setting")
abline(a = 0, b = 1, lty = 2, col = "gray")
grid(col = "lightgray", lty = "dotted")

dev.off()
cat("  ✓ Saved: sroc_curve.png\n")

# ==============================================================================
# SUMMARY TABLE WITH ALL METRICS
# ==============================================================================

cat("\nCreating comprehensive summary table...\n")

meta_data$PPV <- meta_data$TP / (meta_data$TP + meta_data$FP)
meta_data$NPV <- meta_data$TN / (meta_data$TN + meta_data$FN)
meta_data$Accuracy <- (meta_data$TP + meta_data$TN) /
                      (meta_data$TP + meta_data$TN + meta_data$FP + meta_data$FN)

summary_table <- meta_data[, c("Study_ID", "Setting", "Device", "Reference_Test",
                               "Sensitivity", "Specificity", "AUC",
                               "n_diseased", "n_non_diseased",
                               "TP", "FN", "FP", "TN",
                               "PLR", "NLR", "DOR",
                               "PPV", "NPV", "Accuracy")]

write.csv(summary_table, "summary_table_complete.csv", row.names = FALSE)
cat("  ✓ Saved: summary_table_complete.csv\n")

# ==============================================================================
# HETEROGENEITY ANALYSIS
# ==============================================================================

cat("\n==============================================================================\n")
cat("HETEROGENEITY ANALYSIS\n")
cat("==============================================================================\n\n")

cat("SENSITIVITY:\n")
cat(sprintf("  I² = %.1f%%\n", sens_meta$I2 * 100))
cat(sprintf("  τ² = %.4f\n", sens_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", sens_meta$Q, sens_meta$pval.Q))

cat("SPECIFICITY:\n")
cat(sprintf("  I² = %.1f%%\n", spec_meta$I2 * 100))
cat(sprintf("  τ² = %.4f\n", spec_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", spec_meta$Q, spec_meta$pval.Q))

cat("POSITIVE LIKELIHOOD RATIO (PLR):\n")
cat(sprintf("  I² = %.1f%%\n", plr_meta$I2 * 100))
cat(sprintf("  τ² = %.4f\n", plr_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", plr_meta$Q, plr_meta$pval.Q))

cat("NEGATIVE LIKELIHOOD RATIO (NLR):\n")
cat(sprintf("  I² = %.1f%%\n", nlr_meta$I2 * 100))
cat(sprintf("  τ² = %.4f\n", nlr_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", nlr_meta$Q, nlr_meta$pval.Q))

cat("DIAGNOSTIC ODDS RATIO (DOR):\n")
cat(sprintf("  I² = %.1f%%\n", dor_meta$I2 * 100))
cat(sprintf("  τ² = %.4f\n", dor_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", dor_meta$Q, dor_meta$pval.Q))

# ==============================================================================
# FINAL COMPREHENSIVE SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("==============================================================================\n\n")

cat("Generated files:\n")
cat("  • forest_plot_sensitivity.png    - Sensitivity by subgroup\n")
cat("  • forest_plot_specificity.png    - Specificity by subgroup\n")
cat("  • forest_plot_PLR.png            - Positive Likelihood Ratio\n")
cat("  • forest_plot_NLR.png            - Negative Likelihood Ratio\n")
cat("  • forest_plot_DOR.png            - Diagnostic Odds Ratio\n")
cat("  • sroc_curve.png                 - Summary ROC curve\n")
cat("  • summary_table_complete.csv     - All metrics\n\n")

cat("FINAL POOLED ESTIMATES:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("Sensitivity:  %.3f (%.1f%%)\n", pooled_sens, pooled_sens*100))
cat(sprintf("Specificity:  %.3f (%.1f%%)\n", pooled_spec, pooled_spec*100))
cat(sprintf("PLR:          %.2f\n", pooled_PLR))
cat(sprintf("NLR:          %.3f\n", pooled_NLR))
cat(sprintf("DOR:          %.2f\n", pooled_DOR))
cat("----------------------------------------------------------------------\n\n")

cat("INTERPRETATION:\n")
cat("----------------------------------------------------------------------\n")
if (pooled_PLR > 10) {
  cat("PLR > 10: Large increase in probability of disease when test positive\n")
} else if (pooled_PLR > 5) {
  cat("PLR 5-10: Moderate increase in probability of disease when test positive\n")
} else {
  cat("PLR < 5: Small increase in probability of disease when test positive\n")
}

if (pooled_NLR < 0.1) {
  cat("NLR < 0.1: Large decrease in probability of disease when test negative\n")
} else if (pooled_NLR < 0.2) {
  cat("NLR 0.1-0.2: Moderate decrease in probability of disease when test negative\n")
} else {
  cat("NLR > 0.2: Small decrease in probability of disease when test negative\n")
}

if (pooled_DOR > 100) {
  cat("DOR > 100: Excellent discriminatory test performance\n")
} else if (pooled_DOR > 20) {
  cat("DOR 20-100: Good discriminatory test performance\n")
} else if (pooled_DOR > 10) {
  cat("DOR 10-20: Moderate discriminatory test performance\n")
} else {
  cat("DOR < 10: Poor discriminatory test performance\n")
}
cat("----------------------------------------------------------------------\n\n")

cat("All outputs are ready for your living systematic review!\n")
cat("==============================================================================\n")
