# ==============================================================================
# MODALITY-SPECIFIC DIAGNOSTIC META-ANALYSIS
# Pooled Sensitivity & Specificity by Modality and Setting
# Living Systematic Review – NIR Imaging
#
# Outputs:
#   - pooled_sens_spec_modality_setting.csv
#   - forest_modality_setting/ (optional forest plots)
# ==============================================================================

# ----------------------------
# 0) Packages
# ----------------------------
if (!requireNamespace("mada", quietly = TRUE)) {
    install.packages("mada", repos = "https://cloud.r-project.org")
}
library(mada)
if (!requireNamespace("meta", quietly = TRUE)) {
    install.packages("meta", repos = "https://cloud.r-project.org")
}
library(meta)

# ----------------------------
# 1) Load data
# ----------------------------
cat("\nLoading data from niri_meta.csv...\n")
meta_data <- read.csv("niri_meta.csv", stringsAsFactors = FALSE)

required_cols <- c(
    "Study_ID", "Setting", "Modality",
    "Sensitivity", "Specificity",
    "n_diseased", "n_non_diseased"
)

missing_cols <- setdiff(required_cols, names(meta_data))
if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

cat("Studies loaded:", nrow(meta_data), "\n")

# ----------------------------
# 2) Build 2x2 tables
# ----------------------------
meta_data$TP <- round(meta_data$Sensitivity * meta_data$n_diseased)
meta_data$FN <- meta_data$n_diseased - meta_data$TP
meta_data$TN <- round(meta_data$Specificity * meta_data$n_non_diseased)
meta_data$FP <- meta_data$n_non_diseased - meta_data$TN

if (any(meta_data$TP < 0 | meta_data$FN < 0 |
    meta_data$TN < 0 | meta_data$FP < 0)) {
    stop("Negative values in TP/FN/FP/TN – check input data.")
}

# ----------------------------
# 3) Standardize labels
# ----------------------------
norm_setting <- function(x) {
    x <- tolower(trimws(x))
    ifelse(x %in% c("in vivo", "invivo", "clinical", "in_vivo"), "in vivo",
        ifelse(x %in% c("in vitro", "invitro", "laboratory", "in_vitro"),
            "in vitro", x
        )
    )
}

meta_data$Setting <- norm_setting(meta_data$Setting)
meta_data$Modality <- toupper(trimws(meta_data$Modality))
meta_data$Modality <- ifelse(meta_data$Modality %in%
    c("NIRT", "NILT", "TRANSILLUMINATION"),
"NIRT",
ifelse(meta_data$Modality %in%
    c("NIRR", "REFLECTANCE"),
"NIRR", NA
)
)

if (any(is.na(meta_data$Modality))) {
    stop("Some studies have missing or invalid Modality labels.")
}

cat("\nStudies by setting:\n")
print(table(meta_data$Setting))

cat("\nStudies by modality:\n")
print(table(meta_data$Modality))

# ----------------------------
# 4) Continuity correction
# ----------------------------
zero_cells <- with(meta_data, TP == 0 | FN == 0 | FP == 0 | TN == 0)

d_cc <- meta_data
if (any(zero_cells)) {
    cat("\nZero cells detected – applying 0.5 continuity correction.\n")
    d_cc[, c("TP", "FN", "FP", "TN")] <- d_cc[, c("TP", "FN", "FP", "TN")] + 0.5
}

# ----------------------------
# 5) Helper functions
# ----------------------------
fit_reitsma_safe <- function(df) {
    if (nrow(df) < 2) {
        return(NULL)
    }

    mada::reitsma(
        cbind(TP, FN, FP, TN) ~ 1,
        data = df
    )
}


invlogit <- function(x) {
    x <- as.numeric(x)
    exp(x) / (1 + exp(x))
}


extract_sens_spec <- function(fit) {
    # Pooled means on logit scale
    mu <- as.numeric(fit$mu)

    # Confidence intervals on logit scale
    ci <- as.matrix(confint(fit))

    data.frame(
        sens = invlogit(mu[1]),
        sens_lcl = invlogit(ci[1, 1]),
        sens_ucl = invlogit(ci[1, 2]),
        spec = invlogit(mu[2]),
        spec_lcl = invlogit(ci[2, 1]),
        spec_ucl = invlogit(ci[2, 2])
    )
}


# ----------------------------
# 6) Define comparison cells
# ----------------------------
cells <- rbind(
    data.frame(Modality = "ALL", Setting = "ALL"),
    expand.grid(Modality = c("NIRT", "NIRR"), Setting = "ALL"),
    expand.grid(Modality = "ALL", Setting = c("in vivo", "in vitro")),
    expand.grid(
        Modality = c("NIRT", "NIRR"),
        Setting = c("in vivo", "in vitro")
    )
)

# ----------------------------
# 7) Run pooled analyses
# ----------------------------
results <- list()

for (i in seq_len(nrow(cells))) {
    m <- cells$Modality[i]
    s <- cells$Setting[i]

    df <- d_cc
    if (m != "ALL") df <- df[df$Modality == m, ]
    if (s != "ALL") df <- df[df$Setting == s, ]

    k <- nrow(df)
    fit <- fit_reitsma_safe(df)

    if (is.null(fit)) {
        results[[i]] <- data.frame(
            Modality = m, Setting = s, k = k,
            sens = NA, sens_lcl = NA, sens_ucl = NA,
            spec = NA, spec_lcl = NA, spec_ucl = NA,
            note = "Not estimable (k<2)"
        )
    } else {
        est <- extract_sens_spec(fit)
        results[[i]] <- data.frame(
            Modality = m, Setting = s, k = k,
            sens = est$sens, sens_lcl = est$sens_lcl, sens_ucl = est$sens_ucl,
            spec = est$spec, spec_lcl = est$spec_lcl, spec_ucl = est$spec_ucl,
            note = ifelse(k < 4, "Exploratory (small k)", "")
        )
    }
}

results_table <- do.call(rbind, results)

# Percent format for readability
results_table$sens_pct <- round(100 * results_table$sens, 1)
results_table$spec_pct <- round(100 * results_table$spec, 1)

# ----------------------------
# 8) Save outputs
# ----------------------------
write.csv(results_table,
    "pooled_sens_spec_modality_setting.csv",
    row.names = FALSE
)

cat("\nSaved: pooled_sens_spec_modality_setting.csv\n\n")
print(results_table)

# ----------------------------
# 9) Forest plots (Sensitivity / Specificity) using meta::metaprop
#    (reitsma objects do NOT have a forest() method)
# ----------------------------

dir.create("forest_modality_setting", showWarnings = FALSE)

make_forest_prop <- function(df, outcome = c("sens", "spec"), filename, main_title) {
    outcome <- match.arg(outcome)

    if (nrow(df) < 2) {
        return(invisible(NULL))
    }

    if (outcome == "sens") {
        event <- df$TP
        total <- df$TP + df$FN
        xlab <- "Sensitivity"
    } else {
        event <- df$TN
        total <- df$TN + df$FP
        xlab <- "Specificity"
    }

    m <- meta::metaprop(
        event = event,
        n = total,
        studlab = df$Study_ID,
        data = df,
        sm = "PLOGIT", # logit scale (good for proportions)
        method.tau = "ML",
        comb.fixed = FALSE,
        comb.random = TRUE,
        hakn = TRUE,
        incr = 0.5, # continuity correction
        allincr = TRUE,
        level = 0.95
    )

    png(filename, width = 9, height = 6, units = "in", res = 300)
    meta::forest(
        m,
        comb.fixed = FALSE,
        comb.random = TRUE,
        prediction = TRUE,
        backtransf = TRUE,
        xlab = xlab,
        leftcols = c("studlab", "event", "n"),
        leftlabs = c("Study", "Events", "Total"),
        rightcols = c("effect", "ci", "w.random"),
        rightlabs = c(xlab, "95% CI", "% Weight"),
        print.I2 = TRUE,
        print.tau2 = TRUE,
        print.Q = TRUE,
        main = main_title
    )
    dev.off()
}

# Generate forest plots for each Modality × Setting cell
for (i in seq_len(nrow(cells))) {
    m <- cells$Modality[i]
    s <- cells$Setting[i]

    df <- d_cc
    if (m != "ALL") df <- df[df$Modality == m, ]
    if (s != "ALL") df <- df[df$Setting == s, ]

    if (nrow(df) < 2) next

    tag <- paste0(gsub(" ", "_", m), "__", gsub(" ", "_", s))

    make_forest_prop(
        df, "sens",
        filename = file.path("forest_modality_setting", paste0("forest_sens_", tag, ".png")),
        main_title = paste("Forest plot – Sensitivity |", m, "|", s)
    )

    make_forest_prop(
        df, "spec",
        filename = file.path("forest_modality_setting", paste0("forest_spec_", tag, ".png")),
        main_title = paste("Forest plot – Specificity |", m, "|", s)
    )
}

cat("\nForest plots saved in folder: forest_modality_setting/\n")
