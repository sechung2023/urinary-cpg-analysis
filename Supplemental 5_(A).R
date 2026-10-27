# ============================================================
# NIHMD: OR Forest Plot — CPG Median Split + CR+NANA, Case vs Control
# ============================================================
library(dplyr)
library(ggplot2)

# ── Prep ─────────────────────────────────────────────────────
nihmd_data <- NIHMD %>%
  mutate(
    SEX  = factor(SEX, levels = c(0, 1), labels = c("Male", "Female")),
    SK   = factor(SK,  levels = c(0, 1), labels = c("NSK", "SK")),
    RACE = factor(RACE, levels = c(0, 1), labels = c("EA", "AA"))
  )

cpg_median <- median(NIHMD$NIHMDCPG[NIHMD$STATUS == 0], na.rm = TRUE)
crnana_median <- median(NIHMD$CR[NIHMD$STATUS == 0] + NIHMD$NANA[NIHMD$STATUS == 0], na.rm = TRUE)

nihmd_data <- nihmd_data %>%
  mutate(
    CPG_group = factor(ifelse(NIHMDCPG >= cpg_median, "High", "Low"),
                       levels = c("Low", "High")),
    CRNANA_group = factor(ifelse(CR + NANA >= crnana_median, "High", "Low"),
                          levels = c("Low", "High"))
  )

# ── Model ─────────────────────────────────────────────────────
glm_med <- glm(STATUS ~ CPG_group + CRNANA_group + AGE + SEX + SK + RACE,
               data = nihmd_data, family = binomial)

# ── Extract results ───────────────────────────────────────────
tidy_glm <- function(model) {
  coefs <- summary(model)$coefficients
  ci    <- confint(model)
  data.frame(
    term  = rownames(coefs)[-1],
    OR    = exp(coefs[-1, "Estimate"]),
    lower = exp(ci[-1, "2.5 %"]),
    upper = exp(ci[-1, "97.5 %"]),
    pval  = coefs[-1, "Pr(>|z|)"]
  )
}

res <- tidy_glm(glm_med)

# ── Labels ────────────────────────────────────────────────────
label_map <- c(
  "CPG_groupHigh"    = "CPG: High vs Low",
  "CRNANA_groupHigh" = "CR+NANA: High vs Low",
  "AGE"              = "Age",
  "SEXFemale"        = "Sex: Female vs Male",
  "SKSK"             = "Smoking: SK vs NSK",
  "RACEAA"           = "Race: AA vs EA"
)

res$label <- label_map[res$term]
if (any(is.na(res$label))) {
  warning("Unmapped terms: ", paste(res$term[is.na(res$label)], collapse = ", "))
  res$label[is.na(res$label)] <- res$term[is.na(res$label)]
}

pval_label <- function(p) {
  ifelse(p >= 0.05,   "  ns",
         ifelse(p < 0.0001,  "p < 0.0001",
                ifelse(p < 0.001,   "p < 0.001",
                       ifelse(p < 0.01,    "p < 0.01",
                              "p < 0.05"))))
}

res$significant <- ifelse(res$pval < 0.05, "p < 0.05", "p ≥ 0.05")
res$pval_text   <- pval_label(res$pval)
res$label       <- factor(res$label, levels = rev(res$label))

# ── Plot ──────────────────────────────────────────────────────
ggplot(res, aes(x = OR, y = label, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(aes(x = upper,
                label = sprintf("OR: %.2f (%.2f–%.2f)", OR, lower, upper)),
            hjust = -0.045, vjust = -0.2, size = 3, show.legend = FALSE) +
  geom_text(aes(x = upper, label = pval_text),
            hjust = -0.13, vjust = 1.2, size = 2.7, show.legend = FALSE) +
  scale_x_log10(expand = expansion(mult = c(0.35, 0.35))) +
  scale_color_manual(values = c("p < 0.05" = "red3", "p ≥ 0.05" = "navy")) +
  labs(
    title = "Odds Ratios — CPG, CR+NANA Median Split, Case vs Control (Exploratory Cohort)",
    x     = "Odds Ratio (log scale)",
    y     = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )