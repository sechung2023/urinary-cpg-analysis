################################################################################
## ROC (4 plots, 2 curves each)  |  NIHMD + COSPFINAL
##
## Curves per plot:  1) CR + NANA + CPG   2) CR + NANA
## Subsets:          Never smokers, Smokers
## Cohorts:          NIHMD, COSPFINAL
##  -> 2 subsets x 2 cohorts = 4 ROC plots
##
## UNADJUSTED: no demographic covariates. Same subjects across the two nested
## models for a fair AUC comparison.
##
## Output: 5 print()s -> each of the 4 plots individually, then a 2x2 panel.
################################################################################

library(pROC)
library(ggplot2)
library(dplyr)
library(patchwork)

## ---- assumes NIHMD and COSPFINAL are already loaded as data frames ----------

outdir <- "plots"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ============================================================================
## COHORT CONFIG  --  edit column names here only
## ============================================================================
cfg_NIHMD <- list(
  name        = "NIHMD",
  data        = NIHMD,
  cpg         = "NIHMDCPG",   # marker (displayed as CPG)
  cr          = "CR",         # <-- EDIT if column differs
  nana        = "NANA",       # <-- EDIT if column differs
  status      = "STATUS",
  case_level  = 1,            # value of STATUS marking a CASE
  sk          = "SK",         # binary smoking (logistic convention)
  never_level = 0             # value of SK marking a NEVER smoker
)

cfg_COSP <- list(
  name        = "COSPFINAL",
  data        = COSPFINAL,
  cpg         = "CPG",        # <-- EDIT if column differs
  cr          = "CR",         # <-- EDIT if column differs
  nana        = "NANA",       # <-- EDIT if column differs
  status      = "STATUS",
  case_level  = 1,
  sk          = "SK",
  never_level = 0
)

## ---- palette (2 curves) -----------------------------------------------------
roc_cols <- c("CR + NANA + CPG" = "#E4002B",   # vivid red
              "CR + NANA"       = "#0033A0")   # vivid blue

## ============================================================================
## ROC helper
## ============================================================================
roc_to_df <- function(r, label) {
  df <- data.frame(fpr   = 1 - r$specificities,
                   tpr   = r$sensitivities,
                   model = sprintf("%s (AUC = %.3f)", label, as.numeric(auc(r))))
  df[order(df$fpr, df$tpr), ]
}

make_roc_plot <- function(cfg, subset_type = c("never", "smoker")) {
  subset_type <- match.arg(subset_type)
  d <- cfg$data
  d$.y   <- as.integer(d[[cfg$status]] == cfg$case_level)
  d$.CPG <- d[[cfg$cpg]]
  d$.CR  <- d[[cfg$cr]]
  d$.NA  <- d[[cfg$nana]]
  d$.SK  <- d[[cfg$sk]]
  
  if (subset_type == "never")  d <- d[!is.na(d$.SK) & d$.SK == cfg$never_level, ]
  if (subset_type == "smoker") d <- d[!is.na(d$.SK) & d$.SK != cfg$never_level, ]
  
  ## same subjects across both models (fair AUC comparison)
  d <- d[complete.cases(d$.y, d$.CPG, d$.CR, d$.NA), ]
  
  ## unadjusted models: marker terms only
  m3 <- glm(.y ~ .CR + .NA + .CPG, data = d, family = binomial)
  m2 <- glm(.y ~ .CR + .NA,        data = d, family = binomial)
  
  r3 <- roc(d$.y, predict(m3, type = "response"), quiet = TRUE)
  r2 <- roc(d$.y, predict(m2, type = "response"), quiet = TRUE)
  
  df <- rbind(roc_to_df(r3, "CR + NANA + CPG"),
              roc_to_df(r2, "CR + NANA"))
  ## legend order = rbind order (full, then mid); colour by base label
  df$model <- factor(df$model, levels = unique(df$model))
  base_lab <- sub(" \\(AUC.*$", "", levels(df$model))
  pal <- setNames(roc_cols[base_lab], levels(df$model))
  
  sub_lab <- c(never = "Never smokers", smoker = "Smokers")[subset_type]
  n_case <- sum(d$.y == 1); n_ctrl <- sum(d$.y == 0)
  
  ggplot(df, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 1.05) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous("1 - Specificity", limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous("Sensitivity",     limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed() +
    labs(title = sprintf("%s  |  %s", cfg$name, sub_lab),
         subtitle = sprintf("cases = %d, controls = %d  |  unadjusted",
                            n_case, n_ctrl)) +
    theme_bw(base_size = 12) +
    theme(legend.position   = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
          panel.grid.minor  = element_blank(),
          plot.title        = element_text(face = "bold"))
}

## ============================================================================
## SAVE helper
## ============================================================================
save_gg <- function(p, file, w = 6, h = 6)
  ggsave(file.path(outdir, file), p, width = w, height = h, dpi = 300)

## ============================================================================
## RUN: 4 plots
## ============================================================================
cohorts <- list(cfg_NIHMD, cfg_COSP)
subs    <- c("never", "smoker")

roc_panels <- list()
for (cfg in cohorts) {
  for (s in subs) {
    tag <- sprintf("%s_%s", cfg$name, s)
    rp  <- make_roc_plot(cfg, s)
    roc_panels[[tag]] <- rp
    save_gg(rp, sprintf("ROC_%s.png", tag))
  }
}

## ---- 5 print()s: each of the 4 plots individually, then the 2x2 panel ------
ordered <- c("NIHMD_never",     "NIHMD_smoker",
             "COSPFINAL_never", "COSPFINAL_smoker")

for (tag in ordered) print(roc_panels[[tag]])   # plots 1-4, one each

combined <- wrap_plots(roc_panels[ordered], ncol = 2, byrow = TRUE) +
  plot_annotation(
    title = "Lung-cancer CPG biomarker  |  ROC: CR+NANA+CPG vs CR+NANA (unadjusted)",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))

print(combined)                                 # plot 5: combined panel
ggsave(file.path(outdir, "ALL_4_ROC_panel.png"),
       combined, width = 13, height = 13, dpi = 300, limitsize = FALSE)

message("Done. 4 ROC plots + combined panel written to '",
        normalizePath(outdir), "'.")