################################################################################
## 6-plot script: Kaplan-Meier  |  NIHMD + COSPFINAL
##
## Score  : CR + NANA + CPG combined score
##          = sum of control-referenced z-scores (Glass's delta)
##            cr_z(x) = (x - mean_CONTROL) / sd_CONTROL
##          computed per marker on the FULL cohort, then summed.
## Split  : combined score at the CASES median
##          (median among all cases in the cohort, computed BEFORE subsetting)
## Subsets: 1) all cases   2) never-smoker cases   3) smoker cases
##          -> 3 subsets x {NIHMD, COSP} = 6 plots
##
## OUTPUT : 6 individual PNGs (a placeholder tile is written for any subset that
##          collapses to a single group) + one combined 6-panel PNG.
##          ALL 7 are PRINTED to the screen device.
##
## NOTE: xmax trims the DISPLAYED axis only. survfit/survdiff/p are computed on
##       the full follow-up before plotting, so no statistic is affected.
################################################################################

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(patchwork)

## ---- assumes NIHMD and COSPFINAL are already loaded as data frames ----------

outdir <- "plots"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ---- close any stale file devices left open by an interrupted run -----------
## (a dangling pdf()/png() device silently swallows every print() that follows)
while (!is.null(dev.list())) dev.off()

## Optional: ALSO collect the 7 plots into one multi-page PDF.
## This runs as a SECOND pass, after the screen prints -- it can never
## intercept them. Leave FALSE unless you want the PDF.
PRINT_PDF <- FALSE
pdf_file  <- file.path(outdir, "KM_all_7_plots.pdf")

## ============================================================================
## COHORT CONFIG  --  edit column names here only
## ============================================================================
cfg_NIHMD <- list(
  name          = "NIHMD",
  data          = NIHMD,
  cpg           = "NIHMDCPG",   # marker (displayed as CPG)
  cr            = "CR",         # <-- EDIT if column differs
  nana          = "NANA",       # <-- EDIT if column differs
  status        = "STATUS",
  case_level    = 1,            # value of STATUS marking a CASE
  sk            = "SK",         # binary smoking
  never_level   = 0,            # value of SK marking a NEVER smoker
  time_raw      = "DAYS",       # raw survival time column
  time_to_month = 1 / 30.4375,  # DAYS -> months
  event         = "SURVIVAL",   # 1 = deceased
  xmax          = 160,          # display cut at 160 months (set 170 to revert)
  xby           = 20
)

cfg_COSP <- list(
  name          = "COSPFINAL",
  data          = COSPFINAL,
  cpg           = "CPG",        # <-- EDIT if column differs
  cr            = "CR",         # <-- EDIT if column differs
  nana          = "NANA",       # <-- EDIT if column differs
  status        = "STATUS",
  case_level    = 1,
  sk            = "SK",
  never_level   = 0,
  time_raw      = "REMONTHS",   # already in months
  time_to_month = 1,
  event         = "RESURVIVAL", # 1 = deceased
  xmax          = 140,
  xby           = 20
)

## ---- palette ----------------------------------------------------------------
km_lo <- "#0057E7"   # vivid blue  (Low)
km_hi <- "#D50032"   # vivid red   (High)

fmt_p <- function(p) dplyr::case_when(
  is.na(p)   ~ "p = NA",
  p < 0.001  ~ "p < 0.001",
  TRUE       ~ paste0("p = ", formatC(p, format = "f", digits = 3))
)

## control-referenced z-score (Glass's delta): standardize to control mean/SD
cr_z <- function(x, is_control) {
  mu  <- mean(x[is_control], na.rm = TRUE)
  sdv <- sd(x[is_control],   na.rm = TRUE)
  (x - mu) / sdv
}

## ============================================================================
## GUARD: verify configured columns exist before touching the data
## ============================================================================
require_cols <- function(cfg) {
  need <- c(cpg = cfg$cpg, cr = cfg$cr, nana = cfg$nana, status = cfg$status,
            sk = cfg$sk, time_raw = cfg$time_raw, event = cfg$event)
  have <- names(cfg$data)
  miss <- need[!need %in% have]
  if (length(miss)) {
    stop(sprintf(
      "[%s] configured column(s) not found:\n%s\n\nAvailable columns:\n  %s",
      cfg$name,
      paste(sprintf("  cfg$%s = \"%s\"", names(miss), miss), collapse = "\n"),
      paste(have, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

## ============================================================================
## KM helper  (cases only; CR+NANA+CPG combined score, split at CASES median)
## ============================================================================
make_km_plot <- function(cfg, subset_type = c("all", "never", "smoker")) {
  subset_type <- match.arg(subset_type)
  require_cols(cfg)
  d <- as.data.frame(cfg$data)
  
  is_case <- !is.na(d[[cfg$status]]) & d[[cfg$status]] == cfg$case_level
  ctrl    <- !is_case
  
  ## ---- combined score: sum of control-referenced z-scores ------------------
  d$.score  <- cr_z(as.numeric(d[[cfg$cr]]),   ctrl) +
    cr_z(as.numeric(d[[cfg$nana]]), ctrl) +
    cr_z(as.numeric(d[[cfg$cpg]]),  ctrl)
  score_lab <- "CR+NANA+CPG"
  lo_lab    <- paste("Low",  score_lab)
  hi_lab    <- paste("High", score_lab)
  
  ## ---- CASES median (all cases in the cohort, before the smoking subset) ---
  med     <- median(d$.score[is_case], na.rm = TRUE)
  med_lab <- "cases median"
  
  if (!is.finite(med))
    stop(sprintf("[%s] cases median of the combined score is not finite.",
                 cfg$name), call. = FALSE)
  
  ## ---- cases, then smoking subset -----------------------------------------
  dc <- d[is_case, , drop = FALSE]
  dc$.SK <- dc[[cfg$sk]]
  if (subset_type == "never")
    dc <- dc[!is.na(dc$.SK) & dc$.SK == cfg$never_level, , drop = FALSE]
  if (subset_type == "smoker")
    dc <- dc[!is.na(dc$.SK) & dc$.SK != cfg$never_level, , drop = FALSE]
  
  if (nrow(dc) == 0L) {
    message(sprintf("SKIP KM  %s | %s : no cases in this subset.",
                    cfg$name, subset_type))
    return(NULL)
  }
  
  dc$.time  <- dc[[cfg$time_raw]] * cfg$time_to_month
  dc$.event <- as.integer(dc[[cfg$event]])
  dc$.grp   <- factor(ifelse(dc$.score >= med, hi_lab, lo_lab),
                      levels = c(lo_lab, hi_lab))
  dc <- dc[complete.cases(dc$.time, dc$.event, dc$.grp), , drop = FALSE]
  dc$.grp <- droplevels(dc$.grp)
  
  ## ---- GUARD: survdiff needs >= 2 groups -----------------------------------
  if (nlevels(dc$.grp) < 2L) {
    side <- if (all(dc$.score >= med, na.rm = TRUE)) "above/at" else "below"
    message(sprintf(
      "SKIP KM  %s | %s : all %d cases fall %s the cases median (@%.3g) -> only 1 group.",
      cfg$name, subset_type, nrow(dc), side, med))
    return(NULL)
  }
  
  n_low  <- sum(dc$.grp == lo_lab)
  n_high <- sum(dc$.grp == hi_lab)
  
  sub_lab <- c(all = "All cases", never = "Never-smoker cases",
               smoker = "Smoker cases")[subset_type]
  
  ## ---- stats: computed on ALL follow-up, independent of the axis cut -------
  fit <- survfit(Surv(.time, .event) ~ .grp, data = dc)
  sd  <- survdiff(Surv(.time, .event) ~ .grp, data = dc)
  p   <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  
  km <- ggsurvplot(
    fit, data = dc,
    conf.int      = FALSE,          # no CI shading
    censor        = FALSE,
    linetype      = "solid",        # solid lines
    size          = 1.1,
    palette       = c(km_lo, km_hi),
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
    title         = sprintf("%s | %s\n%s split, %s @ %.3g",
                            cfg$name, sub_lab, score_lab, med_lab, med),
    subtitle      = sprintf("Low n = %d, High n = %d", n_low, n_high),
    ggtheme       = theme_bw(base_size = 12)  # NO gridlines here (bleeds into table)
  )
  
  ## gridlines on the plot panel ONLY; clip = "on" keeps curves inside the box
  km$plot <- km$plot +
    coord_cartesian(xlim = c(0, cfg$xmax), ylim = c(0, 1), clip = "on") +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))
  km
}

## placeholder tile for a skipped KM (keeps the combined panel grid intact)
km_placeholder <- function(cfg, subset_type) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("%s | %s\nsingle group after split\n(KM skipped)",
                             cfg$name, subset_type),
             size = 4, lineheight = 1.1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey97", colour = "grey80"))
}

## ============================================================================
## SAVE helpers  -- each opens a device and ALWAYS closes it via on.exit()
## ============================================================================
save_km <- function(km, file, w = 8, h = 7) {          # individual KM (w/ table)
  png(file.path(outdir, file), width = w, height = h, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)   # closes even if print() errors
  print(km)
}

save_tile <- function(g, file, w = 8, h = 7) {         # individual placeholder
  ggsave(file.path(outdir, file), g, width = w, height = h, dpi = 300)
}

## ============================================================================
## BUILD: 6 individual plots + 1 combined
## ============================================================================
cohorts <- list(cfg_NIHMD, cfg_COSP)
subs    <- c("all", "never", "smoker")

km_panels <- list()   # curve panels only, for the combined grid
km_print  <- list()   # everything to be printed, in order

for (cfg in cohorts) {
  for (s in subs) {
    tag  <- sprintf("%s_%s", cfg$name, s)
    file <- sprintf("KM_%s_combined_casesMedian.png", tag)
    km   <- make_km_plot(cfg, s)
    
    if (!is.null(km)) {
      save_km(km, file)                              # -> PNG
      km_panels[[tag]] <- km$plot +                  # -> 6-grid tile
        theme(legend.position = "top",
              plot.title = element_text(size = 9, face = "bold"))
      km_print[[tag]]  <- km                         # -> screen (curve + table)
    } else {
      ph <- km_placeholder(cfg, s)
      save_tile(ph, file)
      km_panels[[tag]] <- ph
      km_print[[tag]]  <- ph
    }
  }
}

## ---- one combined 6-panel: row 1 = NIHMD, row 2 = COSPFINAL ----------------
ordered <- c(paste0("NIHMD_", subs), paste0("COSPFINAL_", subs))

combined <- wrap_plots(km_panels[ordered], ncol = 3, byrow = TRUE) +
  plot_annotation(
    title = "Kaplan-Meier (cases, CR+NANA+CPG combined-score split, cases median)",
    subtitle = "Row 1: NIHMD   Row 2: COSPFINAL   |   Columns: all cases / never smokers / smokers",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(outdir, "KM_6_panel_combined_casesMedian.png"),
       combined, width = 18, height = 11, dpi = 300, limitsize = FALSE)

km_print[["COMBINED_6_panel"]] <- combined
print_order <- c(ordered, "COMBINED_6_panel")

## ============================================================================
## PRINT ALL 7 TO THE SCREEN
## No file device is open here, so these go to the plot pane.
## In RStudio use the plot-history arrows to page back through all 7.
## ============================================================================
stopifnot(length(print_order) == 7L)

for (i in seq_along(print_order)) {
  nm <- print_order[i]
  message("PRINT  [", i, "/7]  ", nm)
  print(km_print[[nm]])
}

## ============================================================================
## OPTIONAL second pass: same 7 plots into one multi-page PDF
## ============================================================================
if (isTRUE(PRINT_PDF)) {
  pdf(pdf_file, width = 8, height = 7, onefile = TRUE)
  on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
  for (nm in print_order) print(km_print[[nm]])
  dev.off()
  message("All 7 also written to: ", normalizePath(pdf_file))
}

message("Done. 6 individual PNGs + 1 combined 6-panel written to '",
        normalizePath(outdir), "'.  All 7 printed to screen.")