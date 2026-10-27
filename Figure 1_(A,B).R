# ============================================================
# CPG boxplots: control-referenced z-score, cases vs controls
# NIHMD and COSPFINAL as SEPARATE plots
# Uses NIHMD and COSPFINAL already loaded in the environment
# ============================================================
library(tidyverse)
library(ggpubr)     # geom_bracket
set.seed(1)         # reproducible jitter

format_p <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "p < 0.0001",
    p < 0.001  ~ "p < 0.001",
    p < 0.01   ~ "p < 0.01",
    p < 0.05   ~ "p < 0.05",
    TRUE       ~ "ns"          # not shown (filtered out below)
  )
}

# add control-referenced z-score (Glass's Δ, scale)
prep_cohort <- function(df, cpg_col) {
  df <- df %>%
    filter(!is.na(.data[[cpg_col]]), !is.na(STATUS), .data[[cpg_col]] > 0) %>%
    mutate(
      Group   = factor(ifelse(STATUS == 1, "Case", "Control"),
                       levels = c("Control", "Case")),
      cpg     = .data[[cpg_col]],
    )
  ctrl <- filter(df, Group == "Control")
  m    <- mean(ctrl$cpg)
  s    <- sd(ctrl$cpg)
  df %>% mutate(zscore = (cpg - m) / s)   # Glass's Δ vs control
}

# five-number summary from ALL points; only whisker ENDS clamped to ±cap
box_stats <- function(df, cap = 2.3) {
  df %>%
    group_by(Group) %>%
    summarise(
      lower  = quantile(zscore, 0.25, type = 7),
      middle = median(zscore),
      upper  = quantile(zscore, 0.75, type = 7),
      iqr    = upper - lower,
      ymin   = min(zscore[zscore >= lower - 1.5 * iqr]),  # Tukey whisker, all data
      ymax   = max(zscore[zscore <= upper + 1.5 * iqr]),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(ymin, -cap),   # cap whisker display only
      ymax = pmin(ymax,  cap)
    )
}

make_box <- function(df, title, bracket_y = 2.7, cap = 2.5) {
  p     <- wilcox.test(df$zscore ~ df$Group)$p.value    # ALL points
  stats <- box_stats(df, cap)
  jit   <- filter(df, zscore >= -cap, zscore <= cap)    # display crop only
  
  ggplot() +
    geom_jitter(data = jit,
                aes(x = Group, y = zscore, color = Group),
                position = position_jitter(width = 0.15, seed = 1),
                size = 2, alpha = 0.7) +
    geom_boxplot(data = stats,
                 aes(x = Group, ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax, fill = Group),
                 stat = "identity", alpha = 0.4, color = "black") +
    geom_bracket(xmin = "Control", xmax = "Case",
                 y.position = bracket_y,
                 label = format_p(p), tip.length = 0.02, label.size = 3.4) +
    scale_y_continuous(breaks = c(-1, 0, 1, 2)) +
    # clip = "off": bracket renders in the title gap above the panel; axis fixed at ±2.5
    coord_cartesian(ylim = c(-1.1, 2.2), clip = "off") +
    labs(x = NULL, y = "CPG (control-referenced z-score)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      # title gap (margin b) is IDENTICAL for both plots -> identical panel length,
      # and the bracket rises into this gap, so the title always sits above it
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 12,
                                 margin = margin(b = 85)),
      axis.text   = element_text(colour = "black"),
      axis.line   = element_line(linewidth = 0.5),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
    )
}

# ---- prep from in-environment data frames ----
nihmd_p <- prep_cohort(NIHMD,     "NIHMDCPG")
cosp_p  <- prep_cohort(COSPFINAL, "CPG")

# ---- two separate plots (NIHMD bracket set higher; same layout otherwise) ----
p_nihmd <- make_box(nihmd_p, "NIHMD — CPG z-score",     bracket_y = 2.7)
p_cosp  <- make_box(cosp_p,  "COSPFINAL — CPG z-score", bracket_y = 2.5)

# ---- save each to its own file (IDENTICAL dimensions) ----
ggsave("CPG_zscore_NIHMD.png",     p_nihmd, width = 6, height = 6.5, dpi = 300)
ggsave("CPG_zscore_COSPFINAL.png", p_cosp,  width = 6, height = 6.5, dpi = 300)

# ---- print each (run one at a time to view separately) ----
p_nihmd
p_cosp