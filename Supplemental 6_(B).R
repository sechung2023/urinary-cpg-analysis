# ============================================================
# COSP: HR Forest Plot — CPG Median Split + CR+NANA, Cases Only
# ============================================================

library(survival)
library(dplyr)
library(ggplot2)

# ── Prep ─────────────────────────────────────────────────────
cases <- COSPFINAL %>% filter(STATUS == 1)

cpg_median <- median(cases$CPG, na.rm = TRUE)
crnana_median <- median(cases$CR + cases$NANA, na.rm = TRUE)

cases <- cases %>%
  mutate(
    CPG_group = factor(ifelse(CPG >= cpg_median, "High", "Low"),
                       levels = c("Low", "High")),
    CRNANA_group = factor(ifelse(CR + NANA >= crnana_median, "High", "Low"),
                          levels = c("Low", "High")),
    SEX      = factor(SEX, levels = c(0, 1), labels = c("Male", "Female")),
    SKCAT    = factor(SKCAT, levels = c(1, 2), labels = c("NSK", "SK")),
    HT       = factor(HT, levels = c(0, 1), labels = c("Adeno", "Squamous")),
    STAGECAT = factor(STAGECAT, levels = c(0, 1), labels = c("Early", "Late")),
    RACE     = factor(RACE, levels = c(0, 1), labels = c("EA", "AA"))
  )

# ── Model ─────────────────────────────────────────────────────
cox_med <- coxph(Surv(MONTHS, SURVIVAL) ~ CPG_group + CRNANA_group + AGE + SEX + SKCAT + HT + STAGECAT + RACE,
                 data = cases)

# ── Extract results ───────────────────────────────────────────
tidy_cox <- function(model) {
  s     <- summary(model)
  coefs <- as.data.frame(s$coefficients)
  ci    <- as.data.frame(s$conf.int)
  data.frame(
    term  = rownames(coefs),
    HR    = ci[, "exp(coef)"],
    lower = ci[, "lower .95"],
    upper = ci[, "upper .95"],
    pval  = coefs[, "Pr(>|z|)"]
  )
}

res <- tidy_cox(cox_med)

# ── Labels ────────────────────────────────────────────────────
label_map <- c(
  "CPG_groupHigh"    = "CPG: High vs Low",
  "CRNANA_groupHigh" = "CR+NANA: High vs Low",
  "AGE"              = "Age",
  "SEXFemale"        = "Sex: Female vs Male",
  "SKCATSK"          = "Smoking: SK vs NSK",
  "HTSquamous"       = "Histology: Squamous vs Adeno",
  "STAGECATLate"     = "Stage: Late vs Early",
  "RACEAA"           = "Race: AA vs EA"
)

res$label <- label_map[res$term]
if (any(is.na(res$label))) {
  warning("Unmapped terms: ", paste(res$term[is.na(res$label)], collapse = ", "))
  res$label[is.na(res$label)] <- res$term[is.na(res$label)]
}

pval_label <- function(p) {
  ifelse(p >= 0.05,   " ns",
         ifelse(p < 0.0001,  "p < 0.0001",
                ifelse(p < 0.001,   "p < 0.001",
                       ifelse(p < 0.01,    "p < 0.01",
                              "p < 0.05"))))
}

res$significant <- ifelse(res$pval < 0.05, "p < 0.05", "p ≥ 0.05")
res$pval_text   <- pval_label(res$pval)
res$label       <- factor(res$label, levels = rev(res$label))

# ── Plot ──────────────────────────────────────────────────────
ggplot(res, aes(x = HR, y = label, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(aes(x = upper,
                label = sprintf("HR: %.2f (%.2f–%.2f)", HR, lower, upper)),
            hjust = -0.04, vjust = -0.1, size = 3, show.legend = FALSE) +
  geom_text(aes(x = upper, label = pval_text),
            hjust = -0.13, vjust = 1.2, size = 2.7, show.legend = FALSE) +
  scale_x_log10(expand = expansion(mult = c(0.15, 0.55))) +
  scale_color_manual(values = c("p < 0.05" = "red3", "p ≥ 0.05" = "navy")) +
  labs(
    title = "Hazard Ratios — CPG, CR+NANA Median Split, Cases Only (Validation Cohort)",
    x     = "Hazard Ratio (log scale)",
    y     = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    plot.margin     = margin(t = 10, r = 80, b = 10, l = 10, unit = "pt")
  )