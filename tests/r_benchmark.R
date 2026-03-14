# HTA JavaScript Meta-Analysis Benchmark vs R
# Validates HTA engine against metafor and meta packages

library(metafor)
library(meta)

cat("=== HTA JavaScript vs R Benchmark ===\n\n")

# Test Dataset 1: Basic pairwise meta-analysis
cat("--- Test 1: Basic Pairwise Meta-Analysis ---\n")
yi <- c(0.5, 0.3, 0.7, 0.4, 0.6)
sei <- c(0.15, 0.20, 0.12, 0.18, 0.14)

# Random-effects with REML (DerSimonian-Laird)
res_dl <- rma(yi = yi, sei = sei, method = "DL")
cat("R (metafor DL):\n")
cat(sprintf("  Pooled effect: %.4f (95%% CI: %.4f to %.4f)\n",
            coef(res_dl), res_dl$ci.lb, res_dl$ci.ub))
cat(sprintf("  tau^2: %.4f, I^2: %.1f%%\n", res_dl$tau2, res_dl$I2))

# Fixed-effect
res_fe <- rma(yi = yi, sei = sei, method = "FE")
cat("\nR (metafor FE):\n")
cat(sprintf("  Pooled effect: %.4f (95%% CI: %.4f to %.4f)\n",
            coef(res_fe), res_fe$ci.lb, res_fe$ci.ub))

# HKSJ adjustment
res_hksj <- rma(yi = yi, sei = sei, method = "DL", test = "knha")
cat("\nR (metafor DL + HKSJ):\n")
cat(sprintf("  Pooled effect: %.4f (95%% CI: %.4f to %.4f)\n",
            coef(res_hksj), res_hksj$ci.lb, res_hksj$ci.ub))

cat("\n--- Test 2: Publication Bias Tests ---\n")

# Larger dataset for publication bias tests
yi2 <- c(0.8, 0.6, 0.9, 0.7, 0.5, 1.0, 0.85, 0.4, 0.3, 0.75)
sei2 <- c(0.15, 0.18, 0.12, 0.20, 0.22, 0.10, 0.14, 0.25, 0.28, 0.16)

# Egger's test
res2 <- rma(yi = yi2, sei = sei2, method = "DL")
egger <- regtest(res2, model = "lm")
cat("Egger's Test (R):\n")
cat(sprintf("  Intercept: %.4f, t = %.4f, p = %.4f\n",
            egger$est[1], egger$zval[1], egger$pval[1]))

# Begg's test
begg <- ranktest(res2)
cat("\nBegg's Rank Correlation (R):\n")
cat(sprintf("  tau = %.4f, p = %.4f\n", begg$tau, begg$pval))

# Trim-and-fill
tf <- trimfill(res2)
cat("\nTrim-and-Fill (R):\n")
cat(sprintf("  k0 = %d studies imputed\n", tf$k0))
cat(sprintf("  Adjusted effect: %.4f (95%% CI: %.4f to %.4f)\n",
            coef(tf), tf$ci.lb, tf$ci.ub))

cat("\n--- Test 3: Meta-Regression ---\n")

# Meta-regression
yi3 <- c(0.5, 0.3, 0.7, 0.4, 0.6)
sei3 <- c(0.15, 0.20, 0.12, 0.18, 0.14)
mod <- c(1, 2, 3, 2, 3)

res_reg <- rma(yi = yi3, sei = sei3, mods = ~ mod, method = "DL")
cat("Meta-regression (R):\n")
cat(sprintf("  Intercept: %.4f (SE: %.4f)\n",
            coef(res_reg)[1], sqrt(vcov(res_reg)[1,1])))
cat(sprintf("  Slope: %.4f (SE: %.4f)\n",
            coef(res_reg)[2], sqrt(vcov(res_reg)[2,2])))
cat(sprintf("  Residual tau^2: %.4f, R^2: %.1f%%\n",
            res_reg$tau2, max(0, res_reg$R2)))

cat("\n--- Test 4: Subgroup Analysis ---\n")

yi4 <- c(0.5, 0.4, 0.6, 0.8, 0.9, 0.7)
sei4 <- c(0.15, 0.18, 0.12, 0.14, 0.16, 0.20)
group <- c("A", "A", "A", "B", "B", "B")

# Using metafor
res_sub <- rma(yi = yi4, sei = sei4, method = "DL")
res_A <- rma(yi = yi4[group == "A"], sei = sei4[group == "A"], method = "DL")
res_B <- rma(yi = yi4[group == "B"], sei = sei4[group == "B"], method = "DL")

cat("Subgroup Analysis (R):\n")
cat(sprintf("  Group A: effect = %.4f (tau^2 = %.4f)\n", coef(res_A), res_A$tau2))
cat(sprintf("  Group B: effect = %.4f (tau^2 = %.4f)\n", coef(res_B), res_B$tau2))

# Q between
Q_A <- res_A$QE
Q_B <- res_B$QE
Q_within <- Q_A + Q_B
Q_between <- res_sub$QE - Q_within
df_between <- 1  # 2 groups - 1
p_between <- 1 - pchisq(Q_between, df_between)
cat(sprintf("  Q_between: %.4f, df: %d, p: %.4f\n", Q_between, df_between, p_between))

cat("\n--- Test 5: Influence Diagnostics ---\n")

# Influence diagnostics
res5 <- rma(yi = yi, sei = sei, method = "DL")
inf <- influence(res5)
cat("Influence Diagnostics (R):\n")
cat("  Study | Cook's D | DFFITS | Hat\n")
for (i in 1:length(yi)) {
  cat(sprintf("  %d     | %.4f   | %.4f | %.4f\n",
              i, inf$inf$cook.d[i], inf$inf$dffits[i], inf$inf$hat[i]))
}

cat("\n--- Test 6: Heterogeneity Statistics ---\n")

res6 <- rma(yi = yi, sei = sei, method = "DL")
cat("Heterogeneity (R metafor):\n")
cat(sprintf("  Q = %.4f (df=%d, p=%.4f)\n", res6$QE, res6$k-1, res6$QEp))
cat(sprintf("  tau^2 = %.4f (tau = %.4f)\n", res6$tau2, sqrt(res6$tau2)))
cat(sprintf("  I^2 = %.1f%%, H^2 = %.4f\n", res6$I2, res6$H2))

# Prediction interval
pi <- predict(res6)
cat(sprintf("  Prediction interval: %.4f to %.4f\n", pi$cr.lb, pi$cr.ub))

cat("\n--- Test 7: Cumulative Meta-Analysis ---\n")

res7 <- rma(yi = yi, sei = sei, method = "DL")
cum <- cumul(res7)
cat("Cumulative Meta-Analysis (R):\n")
for (i in 1:length(yi)) {
  cat(sprintf("  After %d studies: effect = %.4f (95%% CI: %.4f to %.4f)\n",
              i, cum$estimate[i], cum$ci.lb[i], cum$ci.ub[i]))
}

cat("\n--- Test 8: Leave-One-Out ---\n")

res8 <- rma(yi = yi, sei = sei, method = "DL")
loo <- leave1out(res8)
cat("Leave-One-Out Analysis (R):\n")
for (i in 1:length(yi)) {
  cat(sprintf("  Without study %d: effect = %.4f, tau^2 = %.4f\n",
              i, loo$estimate[i], loo$tau2[i]))
}

cat("\n=== EXPECTED VALUES FOR HTA VALIDATION ===\n")
cat("Copy these expected values into the HTA test suite:\n\n")

cat("const R_BENCHMARK = {\n")
cat("  pairwise: {\n")
cat(sprintf("    dl: { effect: %.4f, ci_lower: %.4f, ci_upper: %.4f, tau2: %.6f, I2: %.1f },\n",
            coef(res_dl), res_dl$ci.lb, res_dl$ci.ub, res_dl$tau2, res_dl$I2))
cat(sprintf("    fe: { effect: %.4f, ci_lower: %.4f, ci_upper: %.4f },\n",
            coef(res_fe), res_fe$ci.lb, res_fe$ci.ub))
cat(sprintf("    hksj: { effect: %.4f, ci_lower: %.4f, ci_upper: %.4f }\n",
            coef(res_hksj), res_hksj$ci.lb, res_hksj$ci.ub))
cat("  },\n")
cat("  publicationBias: {\n")
cat(sprintf("    egger: { intercept: %.4f, pvalue: %.4f },\n",
            egger$est[1], egger$pval[1]))
cat(sprintf("    begg: { tau: %.4f, pvalue: %.4f },\n",
            begg$tau, begg$pval))
cat(sprintf("    trimfill: { k0: %d, adjusted: %.4f }\n",
            tf$k0, coef(tf)))
cat("  },\n")
cat("  regression: {\n")
cat(sprintf("    intercept: %.4f, slope: %.4f, tau2_resid: %.6f\n",
            coef(res_reg)[1], coef(res_reg)[2], res_reg$tau2))
cat("  },\n")
cat("  heterogeneity: {\n")
cat(sprintf("    Q: %.4f, tau2: %.6f, I2: %.1f, H2: %.4f,\n",
            res6$QE, res6$tau2, res6$I2, res6$H2))
cat(sprintf("    prediction_lower: %.4f, prediction_upper: %.4f\n",
            pi$cr.lb, pi$cr.ub))
cat("  }\n")
cat("};\n")

cat("\n=== Benchmark Complete ===\n")
