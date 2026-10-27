# ============================================================
# COSP: Cox Regression Table (Cases Only)
# Outcome: Surv(MONTHS, EVENT)
# SURVIVAL: 0 = alive, 1 = expired
# CPG and CR+NANA: Low/High, case-only median split
# Covariates: AGE, SEX, RACE, SK, STAGECAT, HT
# ============================================================

library(dplyr)
library(broom)
library(survival)

# ── Config: cohort-specific column names ─────────────────────
CPG_COL  <- "CPG"
CR_COL   <- "CR"
NANA_COL <- "NANA"

# ── Prep ─────────────────────────────────────────────────────
cases_tmp <- COSPFINAL %>% filter(STATUS == 1)

cpg_median    <- median(cases_tmp[[CPG_COL]], na.rm = TRUE)
crnana_median <- median(cases_tmp[[CR_COL]] + cases_tmp[[NANA_COL]], na.rm = TRUE)

cases <- COSPFINAL %>%
  filter(STATUS == 1) %>%
  mutate(
    CPG_cat = factor(ifelse(.data[[CPG_COL]] >= cpg_median, "High", "Low"),
                     levels = c("Low", "High")),
    CR_NANA_combined = factor(
      ifelse(.data[[CR_COL]] + .data[[NANA_COL]] >= crnana_median, "High", "Low"),
      levels = c("Low", "High")),
    EVENT   = as.numeric(SURVIVAL == 1),
    SEX     = factor(SEX,  levels = c(0, 1), labels = c("Male", "Female")),
    RACE    = factor(RACE, levels = c(0, 1), labels = c("EA", "AA")),
    SKCAT   = factor(SKCAT, levels = c(1, 2), labels = c("NSK", "SK")),
    STAGECAT = factor(STAGECAT, levels = c(0, 1),
                      labels = c("Early Stage (I & II)", "Late Stage (III & IV)")),
    HT      = factor(HT, levels = c(0, 1),
                     labels = c("Adenocarcinoma", "Squamous cell carcinoma"))
  ) %>%
  filter(!is.na(STAGECAT), !is.na(HT), !is.na(MONTHS), !is.na(EVENT),
         !is.na(CPG_cat), !is.na(CR_NANA_combined))

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
make_row <- function(factor_label, levels_label, uni_formula, term, multi_model, data = cases) {
  uni <- coxph(as.formula(uni_formula), data = data)
  n   <- uni$n
  u   <- tidy(uni,         exponentiate = TRUE, conf.int = TRUE) %>% filter(term == !!term)
  m   <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == !!term)
  
  tibble(
    Factor   = factor_label,
    Levels   = levels_label,
    N        = n,
    Uni_HR   = sprintf("%.2f", u$estimate),
    Uni_CI   = sprintf("%.2f, %.2f", u$conf.low, u$conf.high),
    Uni_p    = fmt_p(u$p.value),
    Multi_HR = sprintf("%.2f", m$estimate),
    Multi_CI = sprintf("%.2f, %.2f", m$conf.low, m$conf.high),
    Multi_p  = fmt_p(m$p.value)
  )
}

# ── Multivariable model (both markers mutually adjusted) ─────
multi_full <- coxph(
  Surv(MONTHS, EVENT) ~ CPG_cat + CR_NANA_combined + AGE + SEX + RACE + SKCAT + STAGECAT + HT,
  data = cases
)

# ── Diagnostic: confirm analytic sample matches the forest plot ──
cat("\n--- Sample size check ---\n")
cat("Cases (STATUS == 1):        ", nrow(cases_tmp), "\n")
cat("After NA filtering:         ", nrow(cases), "\n")
cat("Rows used by multi model:   ", multi_full$n, "\n")
cat("\nRaw SKCAT distribution (cases):\n")
print(table(cases_tmp$SKCAT, useNA = "ifany"))
cat("\nIf SKCAT has values outside 1/2, those patients are being dropped.\n")
cat("-------------------------\n\n")

# ── Assemble table ────────────────────────────────────────────
table_cox <- bind_rows(
  make_row("CPG concentration",                  "Low; High",
           "Surv(MONTHS, EVENT) ~ CPG_cat",          "CPG_catHigh",                   multi_full),
  make_row("CR + NANA combined",                 "Low; High",
           "Surv(MONTHS, EVENT) ~ CR_NANA_combined", "CR_NANA_combinedHigh",          multi_full),
  make_row("Stage (ref: Early)",                 "Early; Late",
           "Surv(MONTHS, EVENT) ~ STAGECAT",         "STAGECATLate Stage (III & IV)", multi_full),
  make_row("Histology",                          "Adenocarcinoma; Squamous cell carcinoma",
           "Surv(MONTHS, EVENT) ~ HT",               "HTSquamous cell carcinoma",     multi_full),
  make_row("Age",                                paste0("n=", nrow(cases)),
           "Surv(MONTHS, EVENT) ~ AGE",              "AGE",                           multi_full),
  make_row("Sex",                                "Male; Female",
           "Surv(MONTHS, EVENT) ~ SEX",              "SEXFemale",                     multi_full),
  make_row("Race",                               "African American; European American",
           "Surv(MONTHS, EVENT) ~ RACE",             "RACEAA",                        multi_full),
  make_row("Smoking",                            "NSK; SK",
           "Surv(MONTHS, EVENT) ~ SKCAT",            "SKCATSK",                       multi_full)
)

# ── Print ─────────────────────────────────────────────────────
col_names <- c("Factors", "Levels", "N",
               "HR¹", "95% CI¹", "p-value",
               "HR¹", "95% CI¹", "p-value")

cat("\nCOSP — Cox Regression Table (Cases Only)\n")
cat("¹HR = Hazard Ratio, CI = Confidence Interval\n\n")
print(knitr::kable(table_cox, align = "llrllllll", col.names = col_names))

write.csv(table_cox, "COSP_HR_table.csv", row.names = FALSE)
cat("\nSaved: COSP_HR_table.csv\n")