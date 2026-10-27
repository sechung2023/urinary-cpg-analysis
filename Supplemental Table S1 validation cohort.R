# ============================================================
# COSP: Logistic Regression Table (Full Dataset)
# Outcome: STATUS (0 = control, 1 = case)
# CPG and CR+NANA: Low/High, control-group median split
# Covariates: AGE, SEX, RACE, SK
# NOTE: STAGECAT, HT excluded — coded 0/blank for controls
# ============================================================

library(dplyr)
library(broom)

# ── Prep ─────────────────────────────────────────────────────
cpg_median    <- median(COSPFINAL$CPG[COSPFINAL$STATUS == 0], na.rm = TRUE)
crnana_median <- median(COSPFINAL$CR[COSPFINAL$STATUS == 0] +
                          COSPFINAL$NANA[COSPFINAL$STATUS == 0], na.rm = TRUE)

dat <- COSPFINAL %>%
  mutate(
    CPG_cat = factor(ifelse(CPG >= cpg_median, "High", "Low"),
                     levels = c("Low", "High")),
    CR_NANA_combined = factor(ifelse(CR + NANA >= crnana_median, "High", "Low"),
                              levels = c("Low", "High")),
    STATUS  = as.numeric(STATUS),
    SEX     = factor(SEX,  levels = c(0, 1), labels = c("Male", "Female")),
    RACE    = factor(
      case_when(
        RACE == 0         ~ "EA",
        RACE %in% c(1, 2) ~ "AA",
        TRUE              ~ NA_character_
      ), levels = c("EA", "AA")),
    SK      = factor(SK,   levels = c(0, 1), labels = c("NSK", "SK"))
  )

# ── p-value formatter ────────────────────────────────────────
fmt_p <- function(p) {
  case_when(
    p < 0.001 ~ "<0.001",
    p < 0.01  ~ "<0.01",
    p < 0.05  ~ "<0.05",
    TRUE      ~ sprintf("%.3f", p)
  )
}

# ── Row builder ──────────────────────────────────────────────
make_row <- function(factor_label, levels_label, uni_formula, term, multi_model, data = dat) {
  uni <- glm(as.formula(uni_formula), data = data, family = binomial)
  n   <- nrow(model.frame(uni))
  u   <- tidy(uni,         exponentiate = TRUE, conf.int = TRUE) %>% filter(term == !!term)
  m   <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == !!term)
  
  tibble(
    Factor   = factor_label,
    Levels   = levels_label,
    N        = n,
    Uni_OR   = sprintf("%.2f", u$estimate),
    Uni_CI   = sprintf("%.2f, %.2f", u$conf.low, u$conf.high),
    Uni_p    = fmt_p(u$p.value),
    Multi_OR = sprintf("%.2f", m$estimate),
    Multi_CI = sprintf("%.2f, %.2f", m$conf.low, m$conf.high),
    Multi_p  = fmt_p(m$p.value)
  )
}

# ── Multivariable model (both markers mutually adjusted) ─────
multi_full <- glm(
  STATUS ~ CPG_cat + CR_NANA_combined + AGE + SEX + RACE + SK,
  data = dat, family = binomial
)

# ── Diagnostic: analytic sample ──────────────────────────────
cat("\n--- Sample size check ---\n")
cat("Total rows:               ", nrow(dat), "\n")
cat("Rows used by multi model: ", nrow(model.frame(multi_full)), "\n")
cat("\nRaw RACE:\n");       print(table(COSPFINAL$RACE, useNA = "ifany"))
cat("\nCollapsed RACE:\n");  print(table(dat$RACE, useNA = "ifany"))
cat("-------------------------\n\n")

# ── Assemble table ────────────────────────────────────────────
table_logistic <- bind_rows(
  make_row("CPG concentration",                  "Low; High",
           "STATUS ~ CPG_cat",          "CPG_catHigh",          multi_full),
  make_row("CR + NANA combined",                 "Low; High",
           "STATUS ~ CR_NANA_combined", "CR_NANA_combinedHigh", multi_full),
  make_row("Age",                                paste0("n=", nrow(dat)),
           "STATUS ~ AGE",              "AGE",                  multi_full),
  make_row("Sex",                                "Male; Female",
           "STATUS ~ SEX",              "SEXFemale",            multi_full),
  make_row("Race (ref: EA)",                     "AA",
           "STATUS ~ RACE",             "RACEAA",               multi_full),
  make_row("Smoking",                            "NSK; SK",
           "STATUS ~ SK",               "SKSK",                 multi_full)
)

# ── Print ─────────────────────────────────────────────────────
col_names <- c("Factors", "Levels", "N",
               "OR¹", "95% CI¹", "p-value",
               "OR¹", "95% CI¹", "p-value")

cat("\nCOSP — Logistic Regression Table\n")
cat("¹OR = Odds Ratio, CI = Confidence Interval\n\n")
print(knitr::kable(table_logistic, align = "llrllllll", col.names = col_names))

write.csv(table_logistic, "COSP_OR_table.csv", row.names = FALSE)
cat("\nSaved: COSP_OR_table.csv\n")