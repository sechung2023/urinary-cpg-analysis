# ============================================================
# NIHMD: SHAP Summary (Beeswarm) Plot
# Outcome: STATUS (0 = control, 1 = case)
# Features: NIHMDCPG, NIHMDCR, NIHMDNANA, SEX, RACE, SK, AGE
# ============================================================
## ── Packages ─────────────────────────────────────────────────
for (p in c("xgboost", "shapviz", "ggplot2", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(xgboost)
library(shapviz)
library(ggplot2)
library(dplyr)

## ── Settings ─────────────────────────────────────────────────
seed     <- 123
nrounds  <- 200
features <- c("NIHMDCPG", "CRNANA", "SEX", "RACE", "SK", "AGE")
xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "logloss",
  max_depth        = 3,
  eta              = 0.1,
  subsample        = 0.8,
  colsample_bytree = 0.8
)

## ── Prep data ─────────────────────────────────────────────────
nihmd_feat <- NIHMD %>%
  mutate(
    NIHMDCPG = as.numeric(NIHMDCPG),
    CRNANA   = as.numeric(CR) + as.numeric(NANA),
    SEX      = as.numeric(SEX),
    RACE     = as.numeric(RACE),
    SK       = as.numeric(SK),
    AGE      = as.numeric(AGE),
    STATUS   = as.numeric(STATUS)
  ) %>%
  filter(STATUS %in% c(0, 1)) %>%
  select(all_of(c(features, "STATUS"))) %>%
  filter(complete.cases(.))

cat("Rows used after cleaning:\n")
cat("  NIHMD:", nrow(nihmd_feat), "(", sum(nihmd_feat$STATUS == 1),
    "cancer /", sum(nihmd_feat$STATUS == 0), "control)\n")

## ── Fit XGBoost ──────────────────────────────────────────────
X <- as.matrix(nihmd_feat[, features])
y <- nihmd_feat$STATUS

set.seed(seed)
dtrain <- xgb.DMatrix(X, label = y)
fit <- xgb.train(params = xgb_params, data = dtrain,
                 nrounds = nrounds, verbose = 0)

## ── SHAP values ──────────────────────────────────────────────
shp <- shapviz(fit, X_pred = X, X = as.data.frame(X))

# Rename for display
colnames(shp$S) <- gsub("NIHMDCPG", "CPG", gsub("CRNANA", "CR+NANA", colnames(shp$S)))
colnames(shp$X) <- gsub("NIHMDCPG", "CPG", gsub("CRNANA", "CR+NANA", colnames(shp$X)))

## ── Beeswarm plot ─────────────────────────────────────────────
# Compute mean |SHAP| per feature for annotations
mean_shap <- sort(colMeans(abs(shp$S)), decreasing = TRUE)
shap_labels <- data.frame(
  feature  = names(mean_shap),
  mean_abs = round(as.numeric(mean_shap), 4),
  y_pos    = seq(length(mean_shap), 1)  # sv_importance orders top-to-bottom
)

p <- sv_importance(shp, kind = "beeswarm") +
  ggtitle("SHAP Summary — Exploratory Cohort") +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(t = 10, r = 60, b = 10, l = 10, unit = "pt")
  ) +
  annotate(
    "text",
    x     = max(abs(shp$S)) * 1.05,
    y     = shap_labels$y_pos,
    label = sprintf("%.4f", shap_labels$mean_abs),
    hjust = 0,
    size  = 3.2,
    color = "gray30"
  )

print(p)
ggsave("SHAP_NIHMD.png", p, width = 8, height = 5, dpi = 200)
ggsave("SHAP_NIHMD.pdf", p, width = 8, height = 5)

## ── Mean |SHAP| ranking ───────────────────────────────────────
cat("\nMean |SHAP| importance ranking:\n")
print(sort(colMeans(abs(shp$S)), decreasing = TRUE))
cat("\nSaved: SHAP_NIHMD.png / SHAP_NIHMD.pdf\n")