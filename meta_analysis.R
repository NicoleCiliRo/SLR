# ==============================================================================
# DIAGNOSTIC TEST ACCURACY META-ANALYSIS - PURE R VERSION
# ==============================================================================
# This is a standalone R script that performs the complete meta-analysis
# No Python required! Just run in R or RStudio
# ==============================================================================

# Install required packages (first time only)
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
install_if_missing("dplyr")

cat("\n")
cat("==============================================================================\n")
cat("DIAGNOSTIC TEST ACCURACY META-ANALYSIS\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

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

cat("\nData preparation complete!\n")

# ==============================================================================
# 2. OVERALL BIVARIATE META-ANALYSIS
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("PERFORMING BIVARIATE META-ANALYSIS (REITSMA MODEL)\n")
cat("==============================================================================\n\n")

# Fit Reitsma bivariate model
fit_reitsma <- reitsma(meta_data, 
                       correction = 0.5,
                       correction.control = "all")

# Get summary
summary_reitsma <- summary(fit_reitsma)

# Extract pooled estimates
pooled_sens_logit <- summary_reitsma$coefficients[1, 1]
pooled_spec_logit <- summary_reitsma$coefficients[2, 1]
se_sens <- summary_reitsma$coefficients[1, 2]
se_spec <- summary_reitsma$coefficients[2, 2]

# Transform from logit scale
pooled_sens <- plogis(pooled_sens_logit)
pooled_spec <- plogis(pooled_spec_logit)

# Calculate 95% CI
sens_ci_lower <- plogis(pooled_sens_logit - 1.96 * se_sens)
sens_ci_upper <- plogis(pooled_sens_logit + 1.96 * se_sens)
spec_ci_lower <- plogis(pooled_spec_logit - 1.96 * se_spec)
spec_ci_upper <- plogis(pooled_spec_logit + 1.96 * se_spec)

cat("OVERALL POOLED ESTIMATES:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("Pooled Sensitivity: %.3f (95%% CI: %.3f - %.3f)\n", 
            pooled_sens, sens_ci_lower, sens_ci_upper))
cat(sprintf("                    %.1f%% (95%% CI: %.1f%% - %.1f%%)\n\n",
            pooled_sens*100, sens_ci_lower*100, sens_ci_upper*100))
cat(sprintf("Pooled Specificity: %.3f (95%% CI: %.3f - %.3f)\n",
            pooled_spec, spec_ci_lower, spec_ci_upper))
cat(sprintf("                    %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
            pooled_spec*100, spec_ci_lower*100, spec_ci_upper*100))
cat("----------------------------------------------------------------------\n")

# ==============================================================================
# 3. SUBGROUP ANALYSIS
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("SUBGROUP ANALYSIS BY SETTING\n")
cat("==============================================================================\n\n")

settings <- unique(meta_data$Setting)

for (setting in settings) {
  cat(sprintf("\n%s STUDIES:\n", toupper(setting)))
  cat("----------------------------------------------------------------------\n")
  
  subset_data <- meta_data[meta_data$Setting == setting, ]
  cat(sprintf("Number of studies: %d\n", nrow(subset_data)))
  
  tryCatch({
    fit_subset <- reitsma(subset_data,
                         correction = 0.5,
                         correction.control = "all")
    summary_subset <- summary(fit_subset)
    
    sens_logit <- summary_subset$coefficients[1, 1]
    spec_logit <- summary_subset$coefficients[2, 1]
    se_s <- summary_subset$coefficients[1, 2]
    se_p <- summary_subset$coefficients[2, 2]
    
    sens <- plogis(sens_logit)
    spec <- plogis(spec_logit)
    sens_ci_l <- plogis(sens_logit - 1.96 * se_s)
    sens_ci_u <- plogis(sens_logit + 1.96 * se_s)
    spec_ci_l <- plogis(spec_logit - 1.96 * se_p)
    spec_ci_u <- plogis(spec_logit + 1.96 * se_p)
    
    cat(sprintf("Pooled Sensitivity: %.3f (95%% CI: %.3f - %.3f)\n",
                sens, sens_ci_l, sens_ci_u))
    cat(sprintf("Pooled Specificity: %.3f (95%% CI: %.3f - %.3f)\n",
                spec, spec_ci_l, spec_ci_u))
    
  }, error = function(e) {
    cat(sprintf("Warning: Could not perform meta-analysis for %s\n", setting))
  })
}

# ==============================================================================
# 4. FOREST PLOTS
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("CREATING FOREST PLOTS\n")
cat("==============================================================================\n\n")

# Sensitivity forest plot
cat("Creating sensitivity forest plot...\n")
sens_meta <- metaprop(
  event = TP,
  n = n_diseased,
  studlab = Study_ID,
  data = meta_data,
  sm = "PLOGIT",
  method.tau = "ML",
  subgroup = Setting,
  random = TRUE,
  common = FALSE
)

png("forest_plot_sensitivity.png", width=14, height=10, units="in", res=300)
meta::forest(sens_meta,
       rightcols = c("effect", "ci"),
       rightlabs = c("Sensitivity", "95% CI"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Events", "Total"),
       xlab = "Sensitivity",
       main = "Forest Plot - Sensitivity by Setting",
       common = FALSE,
       random = TRUE,
       col.diamond.lines = "maroon",
       fontsize = 10,
       spacing = 1.2)

dev.off()
cat("  ✓ Saved: forest_plot_sensitivity.png\n")

# Specificity forest plot
cat("Creating specificity forest plot...\n")
spec_meta <- metaprop(
  event = TN,
  n = n_non_diseased,
  studlab = Study_ID,
  data = meta_data,
  sm = "PLOGIT",
  method.tau = "ML",
  subgroup = Setting,
  random = TRUE,
  common = FALSE
)

png("forest_plot_specificity.png", width=14, height=10, units="in", res=300)
meta::forest(spec_meta,
       rightcols = c("effect", "ci"),
       rightlabs = c("Specificity", "95% CI"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Events", "Total"),
       xlab = "Specificity",
       main = "Forest Plot - Specificity by Setting",
       common = FALSE,
       random = TRUE,
       col.diamond.lines = "maroon",
       fontsize = 10,
       spacing = 1.2)

dev.off()
cat("  ✓ Saved: forest_plot_specificity.png\n")

# ==============================================================================
# 5. SROC CURVE (HSROC + 95% CONFIDENCE REGION + 95% PREDICTION REGION)
# ==============================================================================

cat("\nCreating SROC curve with confidence and prediction regions...\n")

# --- 5.1. Obtener valores resumen para la leyenda -----------------------------

# Resumen del modelo Reitsma (sensibilidad / especificidad con IC95%)
fit_sum <- summary(fit_reitsma, level = 0.95)

# La salida de summary.reitsma devuelve matrices 'sens' y 'spec'
# con columnas: est, 2.5%, 97.5% (los nombres exactos pueden variar ligeramente).
sens_est  <- fit_sum$sens[1, 1]
sens_low  <- fit_sum$sens[1, 2]
sens_high <- fit_sum$sens[1, 3]

spec_est  <- fit_sum$spec[1, 1]
spec_low  <- fit_sum$spec[1, 2]
spec_high <- fit_sum$spec[1, 3]

# AUC global del modelo HSROC
auc_est <- AUC(fit_reitsma)$AUC

# Textos bonitos para la leyenda
lab_summary <- sprintf("SENS = %.2f [%.2f – %.2f]\nSPEC = %.2f [%.2f – %.2f]",
                       sens_est, sens_low, sens_high,
                       spec_est, spec_low, spec_high)

lab_auc <- sprintf("AUC = %.2f", auc_est)

# --- 5.2. Punto resumen y región de confianza (ROCellipse) --------------------

# ROCellipse para obtener:
#  - cr$ROCellipse  -> puntos de la región de confianza 95 %
#  - cr$fprsens     -> punto resumen (FPR, Sens)
cr <- ROCellipse(fit_reitsma, level = 0.95)

sp_fpr  <- cr$fprsens["fpr"]
sp_sens <- cr$fprsens["sens"]

# --- 5.3. Crear la figura -----------------------------------------------------

png("sroc_curve.png", width = 8, height = 8, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))

# Plot base de mada:
#   - HSROC
#   - punto resumen
#   - región de confianza alrededor del punto resumen
#   - región de predicción 95 % (predict = TRUE)
plot(
  fit_reitsma,
  level    = 0.95,
  plotsumm = TRUE,
  predict  = TRUE,       # <-- 95% prediction region
  sroclwd  = 2,
  sroclty  = 1,
  predlty  = 3,
  predlwd  = 1.2,
  pch      = 21,
  cex      = 1.2,
  main     = "Summary ROC with 95% confidence and prediction regions"
  # Importante: NO pasar xlab/ylab para evitar el error de 'xlab matched by multiple'
)
# Add study points colored by setting
cols_fill <- c("clinical" = "red", "in_vitro" = "orange")

points(
  1 - meta_data$Specificity,
  meta_data$Sensitivity,
  pch = 21,
  bg  = cols_fill[ meta_data$Setting ],   # relleno según setting
  col = "black",                          # borde
  cex = 1.4
)
points(
  sp_fpr, sp_sens,
  pch = 23,       # rombo rellenable
  bg  = "red",    # relleno
  col = "black",  # borde
  cex = 2
)

# Datos observados: círculos blancos (sensibilidad vs FPR)
points(
  1 - meta_data$Specificity,   # x = False Positive Rate
  meta_data$Sensitivity,       # y = Sensitivity
  pch = 21,
  bg  = "white",
  cex = 1.3
)

# Línea de no discriminación
abline(a = 0, b = 1, lty = 2, col = "grey60")

# Grid suave
grid(col = "lightgray", lty = "dotted")

# --- 5.4. Leyenda estilo figura de artículo -----------------------------------

legend(
  "bottomright",
  bty    = "n",
  cex    = 0.9,
  pt.cex = c(1.3, NA, NA, NA),
  lwd    = c(NA, 2, 1.2, 1.2),
  lty    = c(NA, 1, 2, 3),
  pch    = c(21, NA, NA, NA),
  col    = c("black", "black", "black", "black"),
  legend = c(
    "Observed data",
    paste0("HSROC curve\n", lab_auc),
    "95% confidence region",
    "95% prediction region"
  )
)

dev.off()
cat("  ✓ Saved: sroc_curve.png\n")


# ==============================================================================
# 7. SUMMARY TABLE
# ==============================================================================

cat("\nCreating summary table...\n")

# Calculate additional metrics
meta_data$PPV <- meta_data$TP / (meta_data$TP + meta_data$FP)
meta_data$NPV <- meta_data$TN / (meta_data$TN + meta_data$FN)
meta_data$Accuracy <- (meta_data$TP + meta_data$TN) / 
                      (meta_data$TP + meta_data$TN + meta_data$FP + meta_data$FN)

# Select columns for output
summary_table <- meta_data[, c("Study_ID", "Setting", "Device", "Reference_Test",
                               "Sensitivity", "Specificity", "AUC",
                               "n_diseased", "n_non_diseased",
                               "TP", "FN", "FP", "TN",
                               "PPV", "NPV", "Accuracy")]

write.csv(summary_table, "summary_table.csv", row.names = FALSE)
cat("  ✓ Saved: summary_table.csv\n")

# ==============================================================================
# HETEROGENEITY ANALYSIS
# ==============================================================================

cat("\n==============================================================================\n")
cat("HETEROGENEITY ANALYSIS\n")
cat("==============================================================================\n\n")

cat("SENSITIVITY:\n")
cat(sprintf("  I² = %.1f%%\n", sens_meta$I2))
cat(sprintf("  τ² = %.4f\n", sens_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n", sens_meta$Q, sens_meta$pval.Q))

cat("\nSPECIFICITY:\n")
cat(sprintf("  I² = %.1f%%\n", spec_meta$I2))
cat(sprintf("  τ² = %.4f\n", spec_meta$tau2))
cat(sprintf("  Q = %.2f (p = %.4f)\n\n", spec_meta$Q, spec_meta$pval.Q))


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("==============================================================================\n\n")

cat("Generated files:\n")
cat("  • forest_plot_sensitivity.png\n")
cat("  • forest_plot_specificity.png\n")
cat("  • sroc_curve.png\n")
cat("  • coupled_forest_plot.png\n")
cat("  • summary_table.csv\n\n")

cat("FINAL RESULTS:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("Overall Pooled Sensitivity: %.3f (%.1f%%)\n", pooled_sens, pooled_sens*100))
cat(sprintf("Overall Pooled Specificity: %.3f (%.1f%%)\n", pooled_spec, pooled_spec*100))
cat("----------------------------------------------------------------------\n\n")

cat("All outputs are ready for your living systematic review!\n")
cat("==============================================================================\n")
