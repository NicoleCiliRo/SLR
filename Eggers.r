library(metafor)


dat <- read.csv("niri_meta.csv") # o "summary_table_complete.csv"


dat$TP <- round(dat$Sensitivity * dat$n_diseased)
dat$FN <- dat$n_diseased - dat$TP
dat$TN <- round(dat$Specificity * dat$n_non_diseased)
dat$FP <- dat$n_non_diseased - dat$TN

# 3. Corrección de continuidad para evitar ceros
cc <- 0.5
TPc <- dat$TP + cc
FPc <- dat$FP + cc
FNc <- dat$FN + cc
TNc <- dat$TN + cc

# 4. Calcular DOR, log(DOR) y varianza
dat$DOR <- (TPc * TNc) / (FPc * FNc)
dat$logDOR <- log(dat$DOR)
dat$var <- 1 / TPc + 1 / FPc + 1 / FNc + 1 / TNc

# 5. Meta-análisis de log(DOR)
res <- rma(yi = logDOR, vi = var, data = dat, method = "REML")

# 6. Funnel plot
funnel(res, xlab = "log(DOR)")

# 7. Egger's test de asimetría
regtest(res, model = "rma")

png("funnel_plot_NIR.png", width = 2000, height = 2000, res = 300)
funnel(res, xlab = "log(DOR)")
dev.off()

png("funnel_plot_color.png", width = 2000, height = 2000, res = 300) # o pdf(...)
cols <- ifelse(dat$Setting == "in vivo", "#1F78B4", "#E31A1C") # paleta científica

funnel(res,
    xlab = "log(DOR)",
    pch = 19,
    col = cols,
    cex = 1.4,
    main = "Funnel plot with Egger’s Test"
)

legend("topright",
    legend = c("In vivo", "In vitro"),
    col = c("#1F78B4", "#E31A1C"),
    pch = 19, pt.cex = 1.2, bty = "n"
)

dev.off()
