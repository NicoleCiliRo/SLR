# ==============================================================================
# 5. SROC CURVE
# ==============================================================================

cat("\nCreating SROC curve...\n")

png("sroc_curve.png", width=12, height=10, units="in", res=300)
par(mar=c(5,5,4,2))

# Plot SROC
plot(fit_reitsma, 
     sroclwd = 2,
     main = "Summary ROC Curve with Confidence Region",
     cex.main = 1.5,
     cex.lab = 1.2)

# Add study points colored by setting
colors <- c("blue", "red")
settings <- unique(meta_data$Setting)

for (i in 1:length(settings)) {
  setting_data <- meta_data[meta_data$Setting == settings[i], ]
  fpr <- 1 - setting_data$Specificity
  tpr <- setting_data$Sensitivity
  points(fpr, tpr, pch = 19, col = colors[i], cex = 1.5)
}

# Legend
legend("bottomright", 
       legend = settings,
       col = colors[1:length(settings)],
       pch = 19,
       cex = 1.2,
       title = "Setting",
       bty = "n")

# Reference line
abline(a = 0, b = 1, lty = 2, col = "gray", lwd = 1)
grid(col = "lightgray", lty = "dotted")

dev.off()
cat("  ✓ Saved: sroc_curve.png\n")
