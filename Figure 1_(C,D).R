# =============================================================================
# CPG boxplots: control-referenced z-score (Glass's delta)
# Groups: Control-Never, Control-Ever, Case-Never, Case-Ever
# NIHMD and COSPFINAL as SEPARATE plots
# Uses NIHMD and COSPFINAL already loaded in the environment
#
# Comparisons tested (4 only, not all 6):
#   1) Control-Never vs Control-Ever
#   2) Control-Never vs Case-Never
#   3) Control-Ever  vs Case-Ever
#   4) Case-Never    vs Case-Ever
# Non-significant comparisons are dropped from the plot entirely.
#
# Smoking coding (verified against STATUS cross-tabs):
#   NIHMD SK = 0 never / 1 ever                 -- all subjects
#   COSP  SK = 0 never / 1 former / 2 current   -- all subjects; 1 & 2 -> "Ever"
# =============================================================================

library(tidyverse)
library(ggpubr)     # geom_bracket
set.seed(1)         # reproducible jitter

# ------------------------------- CONFIG --------------------------------------
cfg_NIHMD <- list(
  cohort_name = "NIHMD",
  data_obj    = "NIHMD",
  cpg_col     = "NIHMDCPG",
  smoke_col   = "SK",
  smoke_map   = c("0" = "Never", "1" = "Ever"),
  title       = "NIHMD \u2014 CPG z-score",
  out_png     = "CPG_zscore_status_smoking_NIHMD.png"
)

cfg_COSP <- list(
  cohort_name = "COSPFINAL",
  data_obj    = "COSPFINAL",
  cpg_col     = "CPG",
  smoke_col   = "SK",
  smoke_map   = c("0" = "Never", "1" = "Ever", "2" = "Ever"),
  title       = "COSPFINAL \u2014 CPG z-score",
  out_png     = "CPG_zscore_status_smoking_COSPFINAL.png"
)

TEST         <- "kruskal"      # "kruskal" or "wilcox" for the pairwise tests
P_ADJUST     <- "BH"           # "none" for unadjusted p-values
SIG_ALPHA    <- 0.05           # comparisons at or above this are not drawn
DROP_NONPOS  <- FALSE           # drop CPG <= 0 (matches the earlier pipeline)

CAP          <- 2.5            # display cap for whisker ends and jittered points
YLIM         <- c(-1.1, 2.2)   # visible panel range; all data stay in the statistics
Y_BREAKS     <- c(-1, 0, 1, 2)

BRACKET_START<- 2.65           # y of the lowest bracket level (above the panel)
BRACKET_STEP <- 0.40           # vertical gap between bracket levels
TITLE_GAP    <- 158            # pt of space under the title for the brackets

FIG_W        <- 6.5
FIG_H        <- 7.0
SAVE_PNG     <- TRUE

# Default ggplot2 hue colours, matching the original two-group script
PAL <- c("Control" = "#F8766D", "Case" = "#00BFC4")

# Group keys in x-axis order
GRP_KEYS <- c("Control_Never", "Control_Ever", "Case_Never", "Case_Ever")

COMPARISONS <- list(
  c("Control_Never", "Control_Ever"),
  c("Control_Never", "Case_Never"),
  c("Control_Ever",  "Case_Ever"),
  c("Case_Never",    "Case_Ever")
)

# ------------------------------- HELPERS -------------------------------------
format_p <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "p < 0.0001",
    p < 0.001  ~ "p < 0.001",
    p < 0.01   ~ "p < 0.01",
    p < 0.05   ~ "p < 0.05",
    TRUE       ~ "ns"
  )
}

# Control-referenced z-score: mean/SD taken from ALL controls pooled (Glass's delta)
prep_cohort <- function(cfg) {
  df <- as.data.frame(get(cfg$data_obj))
  
  d <- df %>%
    filter(!is.na(.data[[cfg$cpg_col]]), !is.na(STATUS), !is.na(.data[[cfg$smoke_col]])) %>%
    mutate(
      cpg    = as.numeric(.data[[cfg$cpg_col]]),
      Status = factor(ifelse(STATUS == 1, "Case", "Control"),
                      levels = c("Control", "Case")),
      Smoke  = unname(cfg$smoke_map[as.character(.data[[cfg$smoke_col]])])
    ) %>%
    filter(!is.na(Smoke))
  
  if (DROP_NONPOS) {
    n_drop <- sum(d$cpg <= 0, na.rm = TRUE)
    if (n_drop > 0) message(cfg$cohort_name, ": dropping ", n_drop,
                            " subject(s) with CPG <= 0.")
    d <- d %>% filter(cpg > 0)
  }
  
  ctrl <- filter(d, Status == "Control")
  m <- mean(ctrl$cpg); s <- sd(ctrl$cpg)
  
  d %>%
    mutate(
      zscore = (cpg - m) / s,
      key    = factor(paste(Status, Smoke, sep = "_"), levels = GRP_KEYS)
    ) %>%
    filter(!is.na(key))
}

# x-axis labels carrying the group n
make_labels <- function(d) {
  n_tab <- d %>% count(key, .drop = FALSE)
  lab   <- setNames(
    paste0(sub("_", "\n", as.character(n_tab$key)), "\n(n = ", n_tab$n, ")"),
    as.character(n_tab$key)
  )
  lab
}

# Five-number summary from ALL points; only whisker ENDS clamped to +/- cap
box_stats <- function(df, cap = CAP) {
  df %>%
    group_by(key, Status) %>%
    summarise(
      lower  = quantile(zscore, 0.25, type = 7),
      middle = median(zscore),
      upper  = quantile(zscore, 0.75, type = 7),
      iqr    = upper - lower,
      ymin   = min(zscore[zscore >= lower - 1.5 * iqr]),
      ymax   = max(zscore[zscore <= upper + 1.5 * iqr]),
      .groups = "drop"
    ) %>%
    mutate(ymin = pmax(ymin, -cap),
           ymax = pmin(ymax,  cap))
}

# One pairwise test on ALL points in the two groups
pair_test <- function(d, g1, g2) {
  sub <- d %>% filter(key %in% c(g1, g2))
  g   <- droplevels(sub$key)
  if (nlevels(g) < 2) return(NA_real_)
  if (TEST == "kruskal") kruskal.test(x = sub$zscore, g = g)$p.value
  else                   wilcox.test(sub$zscore ~ g)$p.value
}

# Greedy packing: brackets that don't overlap horizontally share a level
assign_levels <- function(cmp) {
  cmp <- cmp[order(cmp$i2 - cmp$i1, cmp$i1), ]
  cmp$level <- NA_integer_
  for (k in seq_len(nrow(cmp))) {
    lv <- 1L
    repeat {
      occ <- cmp[!is.na(cmp$level) & cmp$level == lv, , drop = FALSE]
      clash <- FALSE
      if (nrow(occ) > 0) {
        clash <- any(!(cmp$i2[k] < occ$i1 | cmp$i1[k] > occ$i2))
      }
      if (!clash) break
      lv <- lv + 1L
    }
    cmp$level[k] <- lv
  }
  cmp
}

# ------------------------------- BUILD ---------------------------------------
make_box <- function(cfg) {
  
  d   <- prep_cohort(cfg)
  lab <- make_labels(d)
  
  cat("\n\n==================== ", cfg$cohort_name, " ====================\n", sep = "")
  print(d %>% group_by(key) %>%
          summarise(n = n(), median_z = round(median(zscore), 3),
                    IQR_z = round(IQR(zscore), 3),
                    median_raw = round(median(cpg), 2), .groups = "drop") %>%
          as.data.frame())
  
  kw_all <- kruskal.test(x = d$zscore, g = droplevels(d$key))
  cat("\nOmnibus Kruskal-Wallis across all groups: chi-squared = ",
      signif(unname(kw_all$statistic), 5), ", df = ", unname(kw_all$parameter),
      ", p = ", signif(kw_all$p.value, 4), "\n", sep = "")
  
  # ---- the four requested comparisons ----
  cmp <- do.call(rbind, lapply(COMPARISONS, function(pr) {
    data.frame(g1 = pr[1], g2 = pr[2],
               p_raw = pair_test(d, pr[1], pr[2]),
               stringsAsFactors = FALSE)
  }))
  cmp$p_adj <- p.adjust(cmp$p_raw, method = P_ADJUST)
  cmp$i1    <- match(cmp$g1, GRP_KEYS)
  cmp$i2    <- match(cmp$g2, GRP_KEYS)
  cmp$label <- format_p(cmp$p_adj)
  
  cat("\nPairwise ", TEST, " tests (", P_ADJUST, "-adjusted):\n", sep = "")
  print(cmp %>% transmute(comparison = paste(g1, "vs", g2),
                          p_raw = signif(p_raw, 4),
                          p_adj = signif(p_adj, 4),
                          shown = ifelse(p_adj < SIG_ALPHA, "yes", "no (ns)")) %>%
          as.data.frame())
  
  br <- cmp %>% filter(!is.na(p_adj), p_adj < SIG_ALPHA)
  if (nrow(br) > 0) {
    br <- assign_levels(br)
    br$y      <- BRACKET_START + (br$level - 1) * BRACKET_STEP
    # NOTE: brackets must reference the factor LEVELS (keys), not the display
    # labels -- passing labels trains them onto the x scale as new categories.
    br$xmin_l <- br$g1
    br$xmax_l <- br$g2
  }
  
  n_out <- sum(d$zscore < -CAP | d$zscore > CAP)
  if (n_out > 0) cat("Note: ", n_out, " point(s) beyond +/-", CAP,
                     " are cropped from the display but remain in all tests.\n", sep = "")
  
  stats <- box_stats(d)
  jit   <- filter(d, zscore >= -CAP, zscore <= CAP)   # display crop only
  
  ggplot() +
    geom_jitter(data = jit,
                aes(x = key, y = zscore, color = Status),
                position = position_jitter(width = 0.15, seed = 1),
                size = 1.6, alpha = 0.55) +
    geom_boxplot(data = stats,
                 aes(x = key, ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax, fill = Status),
                 stat = "identity", alpha = 0.4, color = "black", width = 0.6) +
    {if (nrow(br) > 0)
      geom_bracket(data = br,
                   aes(xmin = xmin_l, xmax = xmax_l, label = label, y.position = y),
                   inherit.aes = FALSE, tip.length = 0.02, label.size = 3.4)} +
    scale_x_discrete(labels = lab, drop = FALSE) +
    scale_y_continuous(breaks = Y_BREAKS) +
    scale_fill_manual(values = PAL) +
    scale_color_manual(values = PAL) +
    # clip = "off": brackets render in the title gap above the panel
    coord_cartesian(ylim = YLIM, clip = "off") +
    labs(x = NULL, y = "CPG (control-referenced z-score)", title = cfg$title) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 12,
                                 margin = margin(b = TITLE_GAP)),
      axis.text   = element_text(colour = "black"),
      axis.text.x = element_text(size = 9.5, lineheight = 1.05),
      axis.line   = element_line(linewidth = 0.5),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
    )
}

# ------------------------------- RUN -----------------------------------------
p_nihmd <- make_box(cfg_NIHMD)
p_cosp  <- make_box(cfg_COSP)

if (SAVE_PNG) {
  ggsave(cfg_NIHMD$out_png, p_nihmd, width = FIG_W, height = FIG_H, dpi = 300)
  ggsave(cfg_COSP$out_png,  p_cosp,  width = FIG_W, height = FIG_H, dpi = 300)
  cat("\nSaved: ", cfg_NIHMD$out_png, " and ", cfg_COSP$out_png, "\n", sep = "")
}

# ---- print each (run one at a time to view separately) ----
print(p_nihmd)
print(p_cosp)