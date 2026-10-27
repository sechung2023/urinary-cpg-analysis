################################################################################
## 2-plot script: Kaplan-Meier  |  NIHMD + COSPFINAL
##
## Subset : NEVER SMOKERS only (cases)
## Split  : CPG dichotomized at the CASE-ONLY median
##          (median of CPG among all cases in the cohort, computed BEFORE the
##           never-smoker filter -- swap noted below for never-smoker-cases-only)
##
## NOTE: xmax trims the DISPLAYED axis only. survfit/survdiff/p are computed on
##       the full data before plotting, so no statistic is affected.
################################################################################

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

## ---- assumes NIHMD and COSPFINAL are already loaded as data frames ----------

outdir <- "plots"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ============================================================================
## COHORT CONFIG  --  edit column names here only
## ============================================================================
cfg_NIHMD <- list(
  name          = "NIHMD",
  data          = NIHMD,
  cpg           = "NIHMDCPG",   # marker (displayed as CPG)
  status        = "STATUS",
  case_level    = 1,            # value of STATUS marking a CASE
  sk            = "SK",         # binary smoking
  never_level   = 0,            # value of SK marking a NEVER smoker
  time_raw      = "DAYS",       # raw survival time column
  time_to_month = 1 / 30.4375,  # DAYS -> months
  event         = "SURVIVAL",   # 1 = deceased
  xmax          = 160,          # <-- display cut at 160 months (was 170)
  xby           = 20
)

cfg_COSP <- list(
  name          = "COSPFINAL",
  data          = COSPFINAL,
  cpg           = "CPG",      # <-- EDIT if column differs
  status        = "STATUS",
  case_level    = 1,
  sk            = "SK",
  never_level   = 0,
  time_raw      = "REMONTHS",     # already in months
  time_to_month = 1,
  event         = "RESURVIVAL",   # 1 = deceased
  xmax          = 140,
  xby           = 20
)

## ---- palette ----------------------------------------------------------------
km_cols <- c("Low CPG" = "#0057E7", "High CPG" = "#D50032")

fmt_p <- function(p) dplyr::case_when(
  is.na(p)   ~ "p = NA",
  p < 0.001  ~ "p < 0.001",
  TRUE       ~ paste0("p = ", formatC(p, format = "f", digits = 3))
)

## ============================================================================
## KM helper: never-smoker cases, split at the CASE-ONLY median
## ============================================================================
make_km_never <- function(cfg) {
  d <- cfg$data
  
  ## ---- CASE-ONLY median (all cases in the cohort) --------------------------
  med <- median(d[[cfg$cpg]][d[[cfg$status]] == cfg$case_level], na.rm = TRUE)
  ## --- alternative: median within NEVER-SMOKER CASES only -------------------
  ## med <- median(d[[cfg$cpg]][d[[cfg$status]] == cfg$case_level &
  ##                            !is.na(d[[cfg$sk]]) &
  ##                            d[[cfg$sk]] == cfg$never_level], na.rm = TRUE)
  
  ## ---- cases, never smokers only -------------------------------------------
  dc <- d[d[[cfg$status]] == cfg$case_level, ]
  dc$.SK <- dc[[cfg$sk]]
  dc <- dc[!is.na(dc$.SK) & dc$.SK == cfg$never_level, ]
  
  dc$.time  <- dc[[cfg$time_raw]] * cfg$time_to_month
  dc$.event <- as.integer(dc[[cfg$event]])
  dc$.grp   <- factor(ifelse(dc[[cfg$cpg]] >= med, "High CPG", "Low CPG"),
                      levels = c("Low CPG", "High CPG"))
  dc <- dc[complete.cases(dc$.time, dc$.event, dc$.grp), ]
  dc$.grp <- droplevels(dc$.grp)
  
  ## ---- GUARD: survdiff needs >= 2 groups -----------------------------------
  if (nlevels(dc$.grp) < 2L) {
    side <- if (nrow(dc) > 0 && all(dc[[cfg$cpg]] >= med)) "above/at" else "below"
    message(sprintf(
      "SKIP KM  %s | never smokers : all %d cases fall %s the case median (@%.3g) -> only 1 CPG group.",
      cfg$name, nrow(dc), side, med))
    return(NULL)
  }
  
  n_low  <- sum(dc$.grp == "Low CPG")
  n_high <- sum(dc$.grp == "High CPG")
  
  ## ---- stats: computed on ALL follow-up time, independent of the axis cut ---
  fit <- survfit(Surv(.time, .event) ~ .grp, data = dc)
  sd  <- survdiff(Surv(.time, .event) ~ .grp, data = dc)
  p   <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  
  km <- ggsurvplot(
    fit, data = dc,
    conf.int      = FALSE,          # no CI shading
    censor        = FALSE,
    linetype      = "solid",        # solid lines
    size          = 1.1,
    palette       = unname(km_cols[levels(dc$.grp)]),
    legend.title  = "",
    legend.labs   = levels(dc$.grp),
    pval          = fmt_p(p),       # custom case_when() p label
    pval.size     = 5,
    risk.table    = TRUE,
    risk.table.height = 0.24,
    xlim          = c(0, cfg$xmax),
    break.time.by = cfg$xby,
    xlab          = "Months",
    ylab          = "Survival probability",
    title         = sprintf("%s  |  Never-smoker cases (case median split @ %.3g)",
                            cfg$name, med),
    subtitle      = sprintf("Low CPG n = %d, High CPG n = %d", n_low, n_high),
    ggtheme       = theme_bw(base_size = 12)  # NO gridlines here (bleeds into table)
  )
  
  ## gridlines on the plot panel ONLY
  ## coord_cartesian() = visual zoom; data are NOT excluded (no scale limits).
  ## clip = "on" keeps curves inside the panel border.
  km$plot <- km$plot +
    coord_cartesian(xlim = c(0, cfg$xmax), ylim = c(0, 1), clip = "on") +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))
  km
}

## ============================================================================
## SAVE helper
## ============================================================================
save_km <- function(km, file, w = 8, h = 7) {
  png(file.path(outdir, file), width = w, height = h, units = "in", res = 300)
  print(km); dev.off()
}

## ============================================================================
## RUN: 2 plots
## ============================================================================
km_nihmd <- make_km_never(cfg_NIHMD)
if (!is.null(km_nihmd)) {
  print(km_nihmd)
  save_km(km_nihmd, "KM_NIHMD_never_caseMedian.png")
}

km_cosp <- make_km_never(cfg_COSP)
if (!is.null(km_cosp)) {
  print(km_cosp)
  save_km(km_cosp, "KM_COSPFINAL_never_caseMedian.png")
}

message("Done. KM PNGs written to '", normalizePath(outdir), "'.")