# ============================================================
# NIHMD CPG boxplots: control-referenced z-score
# Control vs Non-smoker LCC vs Former Smoker LCC vs Current Smoker LCC
# Uses NIHMD already loaded in the environment
# ============================================================
library(tidyverse)
library(ggpubr)     # stat_pvalue_manual
set.seed(1)         # reproducible jitter

# ---- config ----
CPG_COL   <- "NIHMDCPG"
CAP       <- 2.2                     # display crop for dots + whisker ends
YLIM      <- c(-1.1, 2.5)            # visible panel window
BREAKS    <- c(-1, 0, 1, 2)
BRK_START <- 2.9                     # first bracket height (above panel)
BRK_STEP  <- 0.45                    # vertical spacing between brackets

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

# ---- prep: 4-level smoking group + control-referenced z-score (Glass's Δ) ----
prep_nihmd <- function(df, cpg_col) {
  df <- df %>%
    mutate(cpg = readr::parse_number(as.character(.data[[cpg_col]]))) %>%
    filter(!is.na(cpg), !is.na(SKCAT)) %>%
    mutate(
      Group = dplyr::case_when(
        SKCAT == 0 ~ "Control",
        SKCAT == 1 ~ "Non-smoker LCC",
        SKCAT == 2 ~ "Former Smoker LCC",
        SKCAT == 3 ~ "Current Smoker LCC"
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

# ---- axis labels carry the per-group n (from ALL retained values) ----
group_n_labels <- function(df) {
  n_tab <- df %>% count(Group, .drop = FALSE)
  setNames(
    paste0(as.character(n_tab$Group), "\n(n = ", n_tab$n, ")"),
    as.character(n_tab$Group)
  )
}

# five-number summary from ALL points; only whisker ENDS clamped to ±cap
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

# pairwise Wilcoxon on ALL points (no cap in the test)
pw <- function(df, g1, g2) {
  d <- df %>% filter(Group %in% c(g1, g2)) %>% mutate(Group = droplevels(Group))
  wilcox.test(zscore ~ Group, data = d)$p.value
}

make_box <- function(df, title, cap = CAP) {
  stats <- box_stats(df, cap)
  jit   <- filter(df, zscore >= -cap, zscore <= cap)   # display crop only
  labs4 <- group_n_labels(df)
  
  # all six comparisons, shortest span first so brackets nest cleanly
  comps <- tibble::tribble(
    ~group1,              ~group2,
    "Control",            "Non-smoker LCC",
    "Non-smoker LCC",     "Former Smoker LCC",
    "Former Smoker LCC",  "Current Smoker LCC",
    "Control",            "Former Smoker LCC",
    "Non-smoker LCC",     "Current Smoker LCC",
    "Control",            "Current Smoker LCC"
  ) %>%
    mutate(p = purrr::map2_dbl(group1, group2, ~ pw(df, .x, .y))) %>%
    filter(p < 0.05) %>%                     # drop non-significant brackets
    mutate(
      label      = format_p(p),
      y.position = BRK_START + BRK_STEP * (seq_len(n()) - 1)
    )
  
  print(comps)
  
  # reserve headroom for however many brackets survived
  top_gap <- 40 + 34 * nrow(comps)
  
  p <- ggplot() +
    geom_jitter(data = jit,
                aes(x = Group, y = zscore, color = Group),
                position = position_jitter(width = 0.15, seed = 1),
                size = 2, alpha = 0.7) +
    geom_boxplot(data = stats,
                 aes(x = Group, ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax, fill = Group),
                 stat = "identity", alpha = 0.4, color = "black")
  
  if (nrow(comps) > 0) {
    p <- p + stat_pvalue_manual(
      comps, label = "label",
      xmin = "group1", xmax = "group2", y.position = "y.position",
      tip.length = 0.02, label.size = 3.4
    )
  }
  
  p +
    scale_x_discrete(labels = labs4) +
    scale_y_continuous(breaks = BREAKS) +
    coord_cartesian(ylim = YLIM, clip = "off") +
    labs(x = NULL, y = "CPG (control-referenced z-score)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 12,
                                 margin = margin(b = top_gap)),
      axis.text   = element_text(colour = "black"),
      axis.line   = element_line(linewidth = 0.5),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
    )
}

# ---- prep from in-environment data frame ----
nihmd_p <- prep_nihmd(NIHMD, CPG_COL)

# ---- audit exactly what (if anything) is still being dropped ----
cat("rows in NIHMD:      ", nrow(NIHMD), "\n")
cat("rows in nihmd_p:    ", nrow(nihmd_p), "\n")
print(
  NIHMD %>%
    mutate(cpg_parsed = readr::parse_number(as.character(.data[[CPG_COL]]))) %>%
    filter(is.na(cpg_parsed) | is.na(SKCAT) | !SKCAT %in% 0:3) %>%
    select(dplyr::all_of(CPG_COL), SKCAT)
)

# how many rows made it in, and how many dots the crop hides
print(count(nihmd_p, Group))
print(
  nihmd_p %>%
    group_by(Group) %>%
    summarise(
      n            = n(),
      median_z     = median(zscore),
      hidden_above = sum(zscore >  CAP),
      hidden_below = sum(zscore < -CAP),
      .groups = "drop"
    )
)
print(kruskal.test(zscore ~ Group, data = nihmd_p))

# ---- plot ----
p_nihmd <- make_box(nihmd_p, "NIHMD — CPG z-score by smoking")

ggsave("CPG_zscore_smoking4_NIHMD.png", p_nihmd,
       width = 8.5, height = 7.5, dpi = 300)

print(p_nihmd)