# ============================================================
# COSPFINAL CPG boxplots: control-referenced z-score
# Control vs Non-smoker LCC vs Former Smoker LCC vs Current Smoker LCC
# Uses COSPFINAL already loaded in the environment
# Panel geometry + bracket spacing locked to match the NIHMD script
# ============================================================
library(tidyverse)
library(ggpubr)     # stat_pvalue_manual
set.seed(1)         # reproducible jitter

# ---- config ----
CPG_COL   <- "CPG"
GRP_COL   <- "SKCAT2"                # COSP: 0=control, 1=never, 2=former, 3=current

USE_ASSAY_RANGE <- FALSE   # FALSE = all values; TRUE = restrict to ASSAY_MIN..ASSAY_MAX
ASSAY_MIN <- 70
ASSAY_MAX <- 5000

CAP       <- 2.2                     # display crop for dots + whisker ends
YLIM      <- c(-1.1, 2.5)            # visible panel window
BREAKS    <- c(-1, 0, 1, 2)

# ---- geometry lock (must be IDENTICAL to the NIHMD script) ----
BRK_START  <- 2.9      # lowest bracket, just above the panel ceiling
BRK_STEP   <- 0.45     # FIXED gap between brackets — never compressed
TOP_BLOCK  <- 156      # total height above panel, in pt (NIHMD: ~14 title + 142 margin)
MARGIN_T   <- 15       # plot.margin top when a title is drawn
BASE_TITLE <- 12       # title size when it fits
FIG_W      <- 8.5
FIG_H      <- 7.5

# always report a threshold RANGE (never an exact value)
format_p <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "p < 0.0001",
    p < 0.001  ~ "p < 0.001",
    p < 0.01   ~ "p < 0.01",
    p < 0.05   ~ "p < 0.05",
    TRUE       ~ "ns"
  )
}

# fixed-step bracket stack (identical spacing regardless of how many survive)
bracket_positions <- function(k) {
  if (k <= 0) return(numeric(0))
  BRK_START + BRK_STEP * (seq_len(k) - 1)
}

# shrink/drop the title only as far as needed to clear the bracket stack
title_size_for <- function(k) {
  if (k <= 3)      BASE_TITLE   # same as NIHMD
  else if (k == 4) 8
  else             0            # 0 = no title at all
}

# ---- prep: 4-level smoking group + control-referenced z-score (Glass's Δ) ----
prep_cosp <- function(df, cpg_col, grp_col, use_range = USE_ASSAY_RANGE) {
  df <- df %>%
    mutate(cpg = as.numeric(.data[[cpg_col]])) %>%
    filter(!is.na(cpg), !is.na(.data[[grp_col]]))
  
  if (use_range) {
    df <- df %>% filter(cpg >= ASSAY_MIN, cpg <= ASSAY_MAX)
  }
  
  df <- df %>%
    mutate(
      Group = dplyr::case_when(
        .data[[grp_col]] == 0 ~ "Control",
        .data[[grp_col]] == 1 ~ "Non-smoker LCC",
        .data[[grp_col]] == 2 ~ "Former Smoker LCC",
        .data[[grp_col]] == 3 ~ "Current Smoker LCC"
      ),
      Group = factor(
        Group,
        levels = c("Control", "Non-smoker LCC",
                   "Former Smoker LCC", "Current Smoker LCC")
      )
    ) %>%
    filter(!is.na(Group))
  
  ctrl <- filter(df, Group == "Control")
  m    <- mean(ctrl$cpg)
  s    <- sd(ctrl$cpg)
  
  df %>% mutate(zscore = (cpg - m) / s)
}

# ---- axis labels carry the per-group n actually analyzed ----
group_n_labels <- function(df) {
  n_tab <- df %>% count(Group, .drop = FALSE)
  setNames(
    paste0(as.character(n_tab$Group), "\n(n = ", n_tab$n, ")"),
    as.character(n_tab$Group)
  )
}

# five-number summary from ALL retained points; only whisker ENDS clamped to ±cap
box_stats <- function(df, cap = CAP) {
  df %>%
    group_by(Group) %>%
    summarise(
      lower  = quantile(zscore, 0.25, type = 7),
      middle = median(zscore),
      upper  = quantile(zscore, 0.75, type = 7),
      iqr    = upper - lower,
      ymin   = min(zscore[zscore >= lower - 1.5 * iqr]),
      ymax   = max(zscore[zscore <= upper + 1.5 * iqr]),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(ymin, -cap),
      ymax = pmin(ymax,  cap)
    )
}

# pairwise Wilcoxon on ALL retained points (no display cap in the test)
pw <- function(df, g1, g2) {
  d <- df %>% filter(Group %in% c(g1, g2)) %>% mutate(Group = droplevels(Group))
  if (dplyr::n_distinct(d$Group) < 2) return(NA_real_)
  wilcox.test(zscore ~ Group, data = d)$p.value
}

all_pairs <- tibble::tribble(
  ~group1,              ~group2,
  "Control",            "Non-smoker LCC",
  "Control",            "Former Smoker LCC",
  "Control",            "Current Smoker LCC",
  "Non-smoker LCC",     "Former Smoker LCC",
  "Non-smoker LCC",     "Current Smoker LCC",
  "Former Smoker LCC",  "Current Smoker LCC"
)

make_box <- function(df, title, cap = CAP) {
  stats <- box_stats(df, cap)
  jit   <- filter(df, zscore >= -cap, zscore <= cap)   # display crop only
  labs4 <- group_n_labels(df)
  
  comps <- all_pairs %>%
    mutate(p = purrr::map2_dbl(group1, group2, ~ pw(df, .x, .y))) %>%
    filter(!is.na(p), p < 0.05) %>%           # drop non-significant brackets
    mutate(
      label      = format_p(p),
      y.position = bracket_positions(n())
    )
  
  print(comps)
  
  k     <- nrow(comps)
  t_pt  <- title_size_for(k)
  
  if (k >= 5) {
    message("NOTE: ", k, " significant brackets at fixed spacing — the stack ",
            "will run past the top of the figure. Raise FIG_H and TOP_BLOCK ",
            "by the SAME amount in BOTH scripts to make room without ",
            "changing panel height.")
  }
  
  p <- ggplot() +
    geom_jitter(data = jit,
                aes(x = Group, y = zscore, color = Group),
                position = position_jitter(width = 0.15, seed = 1),
                size = 2, alpha = 0.7) +
    geom_boxplot(data = stats,
                 aes(x = Group, ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax, fill = Group),
                 stat = "identity", alpha = 0.4, color = "black")
  
  if (k > 0) {
    p <- p + stat_pvalue_manual(
      comps, label = "label",
      xmin = "group1", xmax = "group2", y.position = "y.position",
      tip.length = 0.02, label.size = 3.4
    )
  }
  
  p <- p +
    scale_x_discrete(labels = labs4) +
    scale_y_continuous(breaks = BREAKS) +
    coord_cartesian(ylim = YLIM, clip = "off") +
    labs(x = NULL, y = "CPG (control-referenced z-score)",
         title = if (t_pt > 0) title else NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text   = element_text(colour = "black"),
      axis.line   = element_line(linewidth = 0.5)
    )
  
  # total space above the panel is held at MARGIN_T + TOP_BLOCK either way,
  # so shrinking or dropping the title never changes panel height
  if (t_pt > 0) {
    p <- p + theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, size = t_pt,
                                 margin = margin(b = TOP_BLOCK - 1.2 * t_pt)),
      plot.margin = margin(t = MARGIN_T, r = 10, b = 10, l = 10)
    )
  } else {
    p <- p + theme(
      plot.title  = element_blank(),
      plot.margin = margin(t = MARGIN_T + TOP_BLOCK, r = 10, b = 10, l = 10)
    )
  }
  
  p
}

# ---- audit exactly what (if anything) is being dropped ----
cat("rows in COSPFINAL: ", nrow(COSPFINAL), "\n")
print(table(COSPFINAL[[GRP_COL]], useNA = "ifany"))
print(
  COSPFINAL %>%
    mutate(cpg_parsed = as.numeric(.data[[CPG_COL]])) %>%
    filter(is.na(cpg_parsed) | is.na(.data[[GRP_COL]]) |
             !.data[[GRP_COL]] %in% 0:3) %>%
    select(dplyr::all_of(c(CPG_COL, "STATUS", GRP_COL)))
)

# ---- side-by-side: all values vs assay-range subset ----
df_all   <- prep_cosp(COSPFINAL, CPG_COL, GRP_COL, use_range = FALSE)
df_range <- prep_cosp(COSPFINAL, CPG_COL, GRP_COL, use_range = TRUE)

cat("\n--- group n under each basis ---\n")
print(
  full_join(
    count(df_all,   Group, name = "n_all_values"),
    count(df_range, Group, name = "n_assay_range"),
    by = "Group"
  ) %>%
    mutate(dropped_by_range = n_all_values - n_assay_range)
)

cat("\n--- pairwise Wilcoxon p under each basis ---\n")
print(
  all_pairs %>%
    mutate(
      p_all_values  = purrr::map2_dbl(group1, group2, ~ pw(df_all,   .x, .y)),
      p_assay_range = purrr::map2_dbl(group1, group2, ~ pw(df_range, .x, .y)),
      lab_all       = format_p(p_all_values),
      lab_range     = format_p(p_assay_range)
    )
)

# ---- plot on whichever basis USE_ASSAY_RANGE selects ----
cosp_p <- prep_cosp(COSPFINAL, CPG_COL, GRP_COL)

print(count(cosp_p, Group))
print(
  cosp_p %>%
    group_by(Group) %>%
    summarise(
      n            = n(),
      median_z     = median(zscore),
      hidden_above = sum(zscore >  CAP),
      hidden_below = sum(zscore < -CAP),
    )
)
print(kruskal.test(zscore ~ Group, data = cosp_p))

p_cosp <- make_box(cosp_p, "COSPFINAL — CPG z-score by smoking")

ggsave("CPG_zscore_smoking4_COSPFINAL.png", p_cosp,
       width = FIG_W, height = FIG_H, dpi = 300)

print(p_cosp)