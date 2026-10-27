# ============================================================
# NIHMD: Cox Regression Table (Cases Only)
# ============================================================

library(dplyr)
library(broom)
library(survival)

# ── Prep ─────────────────────────────────────────────────────
cases_tmp <- NIHMD %>% filter(STATUS == 1)

cpg_median    <- median(cases_tmp$NIHMDCPG, na.rm = TRUE)
crnana_median <- median(cases_tmp$CR + cases_tmp$NANA, na.rm = TRUE)

cases <- NIHMD %>%
  filter(STATUS == 1) %>%
  mutate(
    logCPG   = log(NIHMDCPG + 1),
    CPG_cat  = factor(ifelse(NIHMDCPG >= cpg_median, "High", "Low"),
                      levels = c("Low", "High")),
    CR_NANA_combined = factor(ifelse(CR + NANA >= crnana_median, "High", "Low"),
                              levels = c("Low", "High")),
    EVENT    = ifelse(SURVIVAL == 1, 1, 0),
    SEX      = factor(SEX,  levels = c(0, 1), labels = c("Male", "Female")),
    RACE     = factor(RACE, levels = c(0, 1), labels = c("EA", "AA")),
    SKCAT_bin = factor(
      case_when(
        SKCAT == 1 ~ "Never",
        SKCAT %in% c(2, 3) ~ "Ever",
        TRUE ~ NA_character_
      ), levels = c("Never", "Ever")),
    STAGECAT_bin = factor(
      case_when(
        STAGECAT == 1 ~ "Early",
        STAGECAT == 2 ~ "Late",
        TRUE ~ NA_character_
      ), levels = c("Early", "Late")),
    HT       = factor(HT, levels = c(1, 2),
                      labels = c("Adenocarcinoma", "Squamous cell carcinoma"))
  ) %>%
  filter(!is.na(SKCAT_bin), !is.na(STAGECAT_bin), !is.na(HT), !is.na(CR_NANA_combined))

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

# ── Multivariable models ──────────────────────────────────────
multi_full <- coxph(
  Surv(DAYS, EVENT) ~ CPG_cat + CR_NANA_combined + AGE + SEX + RACE + SKCAT_bin + STAGECAT_bin + HT,
  data = cases
)

# ── Assemble table ────────────────────────────────────────────
table_cox <- bind_rows(
  
  # CPG
  make_row("CPG concentration",                  "Low; High",
           "Surv(DAYS, EVENT) ~ CPG_cat",        "CPG_catHigh",               multi_full),
  
  # CR + NANA combined
  make_row("CR + NANA combined",                 "Low; High",
           "Surv(DAYS, EVENT) ~ CR_NANA_combined", "CR_NANA_combinedHigh", multi_full),
  
  # Stage
  make_row("Stage (ref: Early)",                 "Late",
           "Surv(DAYS, EVENT) ~ STAGECAT_bin",   "STAGECAT_binLate",          multi_full),
  
  # Histology
  make_row("Histology",                          "Adenocarcinoma; Squamous cell carcinoma",
           "Surv(DAYS, EVENT) ~ HT",             "HTSquamous cell carcinoma", multi_full),
  
  # Age
  make_row("Age",                                paste0("n=", nrow(cases)),
           "Surv(DAYS, EVENT) ~ AGE",            "AGE",                       multi_full),
  
  # Sex
  make_row("Sex",                                "Male; Female",
           "Surv(DAYS, EVENT) ~ SEX",            "SEXFemale",                 multi_full),
  
  # Race
  make_row("Race",                               "African American; European American",
           "Surv(DAYS, EVENT) ~ RACE",           "RACEAA",                    multi_full),
  
  # Smoking
  make_row("Smoking (ref: Never)",               "Ever",
           "Surv(DAYS, EVENT) ~ SKCAT_bin",      "SKCAT_binEver",             multi_full)
)

# ── Print ─────────────────────────────────────────────────────
col_names <- c("Factors", "Levels", "N",
               "HR¹", "95% CI¹", "p-value",
               "HR¹", "95% CI¹", "p-value")

cat("\nNIHMD — Cox Regression Table (Cases Only)\n")
cat("¹HR = Hazard Ratio, CI = Confidence Interval\n\n")
print(knitr::kable(table_cox, align = "llrllllll", col.names = col_names))

write.csv(table_cox, "NIHMD_HR_table.csv", row.names = FALSE)
cat("\nSaved: NIHMD_HR_table.csv\n")