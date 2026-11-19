# Classical synthetic control (Abadie) implemented manually with CVXR.
# augsynth::multisynth hits a tidyr version incompatibility in this environment,
# so we implement SCM directly:
#   - treated-aggregate = mean log_obs across the 7 treated zones per month
#   - find w \in simplex (w_i >= 0, sum w = 1) over 54 control zones
#     minimizing sum_t (treated_mean_t - sum_i w_i * y_{i,t})^2 across pre-period
#   - apply w to post-period controls to get synthetic counterfactual
#   - ATT_t = treated_mean_t - synthetic_t for post-period months
#
# Also runs a per-treated-zone version (7 separate SCM fits) and averages.
# Reports pre-period RMSPE as a fit diagnostic.

suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(CVXR); library(ggplot2)
}))

BUILD <- "analysis"
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids61     <- readRDS(file.path(BUILD, "zones61_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)
treat_date <- as.Date("2017-01-01")

md <- full_data %>% filter(locationID %in% ids61) %>%
  mutate(log_obs = log(observations)) %>%
  arrange(locationID, month)
bal <- md %>% count(locationID) %>% filter(n==max(n))
md  <- md %>% filter(locationID %in% bal$locationID)

# wide: rows = month, cols = locationID
Y_long <- md %>% select(month, locationID, log_obs)
Y_wide <- Y_long %>% pivot_wider(names_from=locationID, values_from=log_obs) %>%
  arrange(month)
months <- Y_wide$month
Y <- as.matrix(Y_wide %>% select(-month))
colnames(Y) <- colnames(Y_wide)[-1]

treat_col_ids <- as.character(treat_ids)
ctrl_col_ids  <- setdiff(colnames(Y), treat_col_ids)
cat("treated cols:", length(treat_col_ids), "| control cols:", length(ctrl_col_ids),
    "| months:", length(months), "\n")

pre_idx  <- which(months <  treat_date)
post_idx <- which(months >= treat_date)
cat("pre months:", length(pre_idx), " post months:", length(post_idx), "\n")

# --- SCM for aggregate treated mean ---------------------------------------
y_treated <- rowMeans(Y[, treat_col_ids])
Yctrl     <- Y[, ctrl_col_ids]
n_ctrl    <- ncol(Yctrl)

w  <- Variable(n_ctrl)
obj<- Minimize(sum((y_treated[pre_idx] - Yctrl[pre_idx, ] %*% w)^2))
cons <- list(w >= 0, sum(w) == 1)
prob <- Problem(obj, cons)
res  <- solve(prob)
w_hat<- as.numeric(res$getValue(w))
synth<- as.numeric(Yctrl %*% w_hat)

att      <- y_treated - synth
pre_rmspe<- sqrt(mean(att[pre_idx]^2))
post_rmspe<- sqrt(mean(att[post_idx]^2))
avg_att  <- mean(att[post_idx])

cat("\n========== SCM (aggregate, 7 treated -> 54-donor simplex) ==========\n")
cat(sprintf("Pre-period RMSPE:  %0.4f log-pts\n", pre_rmspe))
cat(sprintf("Post-period RMSPE: %0.4f log-pts  (ratio post/pre = %0.2f)\n",
            post_rmspe, post_rmspe/pre_rmspe))
cat(sprintf("Avg post-period ATT: %+0.4f log-pts  (=  %+0.1f%%)\n",
            avg_att, (exp(avg_att)-1)*100))

# Donors with non-trivial weight
w_df <- data.frame(zone=ctrl_col_ids, weight=w_hat) %>% arrange(desc(weight))
cat("\nTop donor weights (weight > 0.01):\n")
print(w_df %>% filter(weight > 0.01) %>% mutate(weight=round(weight,4)))

# Monthly ATT path
att_df <- data.frame(month=months, treated=y_treated, synthetic=synth, att=att,
                     period=ifelse(months < treat_date, "pre", "post"))
write.csv(att_df, file.path(BUILD, "scm_aggregate_att.csv"), row.names=FALSE)

p1 <- ggplot(att_df, aes(month)) +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_line(aes(y=treated, color="Treated mean"), linewidth=1) +
  geom_line(aes(y=synthetic, color="Synthetic control"), linewidth=1, linetype="dashed") +
  scale_color_manual(values=c("Treated mean"="firebrick", "Synthetic control"="steelblue")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Synthetic control (aggregate): treated-mean log trips vs synthetic counterfactual",
       subtitle=sprintf("Pre-period RMSPE %0.4f, post-period mean ATT %+0.4f (%+0.1f%%)",
                        pre_rmspe, avg_att, (exp(avg_att)-1)*100),
       x=NULL, y="log trips (zone-mean)", color=NULL) +
  theme_minimal() + theme(legend.position="top")
ggsave(file.path(BUILD, "scm_aggregate_61.png"), p1, width=10, height=5.5, dpi=130)

p2 <- ggplot(att_df, aes(month, att, color=period)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray40") +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_color_manual(values=c("pre"="steelblue","post"="firebrick")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Synthetic control gap: treated mean minus synthetic",
       subtitle=sprintf("Pre-period should sit near zero if the match is good. Avg post ATT = %+0.3f log-pts.",
                        avg_att),
       x=NULL, y="log-points (treated - synthetic)", color=NULL) +
  theme_minimal() + theme(legend.position="top")
ggsave(file.path(BUILD, "scm_gap_61.png"), p2, width=10, height=5.5, dpi=130)

# --- SCM per treated zone, then average -----------------------------------
cat("\n========== SCM per treated zone ==========\n")
per_zone_atts <- list()
for (tz in treat_col_ids) {
  y_tz <- Y[, tz]
  w  <- Variable(n_ctrl)
  obj<- Minimize(sum((y_tz[pre_idx] - Yctrl[pre_idx, ] %*% w)^2))
  prob <- Problem(obj, list(w>=0, sum(w)==1))
  res  <- solve(prob)
  w_hat<- as.numeric(res$getValue(w))
  synth_tz <- as.numeric(Yctrl %*% w_hat)
  att_tz   <- y_tz - synth_tz
  per_zone_atts[[tz]] <- att_tz
  cat(sprintf("  zone %-5s pre-RMSPE = %0.4f | avg post ATT = %+0.4f (%+0.1f%%)\n",
              tz, sqrt(mean(att_tz[pre_idx]^2)),
              mean(att_tz[post_idx]), (exp(mean(att_tz[post_idx]))-1)*100))
}
per_zone_mat <- do.call(cbind, per_zone_atts)
avg_att_perzone <- mean(rowMeans(per_zone_mat)[post_idx])
cat(sprintf("\nAverage across 7 per-zone SCMs: avg post ATT = %+0.4f log-pts (%+0.1f%%)\n",
            avg_att_perzone, (exp(avg_att_perzone)-1)*100))

# --- Placebo / permutation inference (Abadie) -----------------------------
# Reassign "treatment" to each control zone in turn, fit a synthetic control for
# it from the OTHER controls, and compare the real treated corridor's post/pre
# RMSPE ratio to this placebo distribution. This is the SCM analogue of a p-value:
# a real effect should make the treated unit an outlier among the placebos.
cat("\n========== Placebo (permutation) inference ==========\n")
ratio_treated <- post_rmspe / pre_rmspe
placebo <- lapply(ctrl_col_ids, function(c0) {
  donors <- setdiff(ctrl_col_ids, c0)
  yv  <- Y[, c0]; Yd <- Y[, donors]
  wv  <- Variable(ncol(Yd))
  pr  <- Problem(Minimize(sum((yv[pre_idx] - Yd[pre_idx, ] %*% wv)^2)),
                 list(wv >= 0, sum(wv) == 1))
  wh  <- as.numeric(solve(pr)$getValue(wv))
  g   <- yv - as.numeric(Yd %*% wh)
  data.frame(zone = c0, month = months, gap = g,
             pre_rmspe  = sqrt(mean(g[pre_idx]^2)),
             post_rmspe = sqrt(mean(g[post_idx]^2)),
             avg_post   = mean(g[post_idx]))
})
placebo_df <- bind_rows(placebo)
pl_sum <- placebo_df %>% group_by(zone) %>%
  summarise(pre_rmspe = first(pre_rmspe), post_rmspe = first(post_rmspe),
            ratio = first(post_rmspe)/first(pre_rmspe),
            avg_post = first(avg_post), .groups = "drop")

treated_pre_rmspe <- pre_rmspe   # treated aggregate pre-period fit, for the filter below
p_ratio <- (sum(pl_sum$ratio >= ratio_treated) + 1) / (nrow(pl_sum) + 1)
p_gap   <- (sum(pl_sum$avg_post <= avg_att) + 1) / (nrow(pl_sum) + 1)  # effect is negative
# standard Abadie practice: exclude placebos that cannot even fit their own
# pre-period (pre-RMSPE > 2x the treated fit), so bad fits can't fake a big effect.
well    <- pl_sum %>% filter(pre_rmspe <= 2 * treated_pre_rmspe)
p_ratio_wf <- (sum(well$ratio >= ratio_treated) + 1) / (nrow(well) + 1)
cat(sprintf("Treated post/pre RMSPE ratio = %.2f\n", ratio_treated))
cat(sprintf("Permutation p (RMSPE ratio), all %d controls:          %.3f\n",
            nrow(pl_sum), p_ratio))
cat(sprintf("Permutation p (RMSPE ratio), %d well-fitting placebos:  %.3f\n",
            nrow(well), p_ratio_wf))
cat(sprintf("Permutation p (avg gap as negative as treated):        %.3f\n", p_gap))
write.csv(pl_sum, file.path(BUILD, "scm_placebo_summary.csv"), row.names = FALSE)

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_vline(xintercept = treat_date, linetype = "dashed", color = "red") +
  geom_line(data = placebo_df, aes(month, gap, group = zone),
            color = "gray75", alpha = 0.6, size = 0.3) +
  geom_line(data = data.frame(month = months, gap = att), aes(month, gap),
            color = "firebrick", size = 1.2) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(title = "Placebo (permutation) test: treated gap vs. control-zone placebos",
       subtitle = sprintf(paste0("Red = real treated corridor. Each gray line = a control zone ",
                                 "treated as a fake intervention.\nPermutation p = %.3f (RMSPE ratio, all controls)."),
                          p_ratio),
       x = NULL, y = "gap: observed - synthetic (log-points)") +
  theme_minimal()
ggsave(file.path(BUILD, "scm_placebo_61.png"), p3, width = 10, height = 5.5, dpi = 130)

# --- Leave-one-donor-out (aggregate SCM robustness) -----------------------
cat("\n========== Leave-one-donor-out (aggregate SCM) ==========\n")
pos_donors <- (w_df %>% filter(weight > 0.01))$zone
loo_don <- sapply(pos_donors, function(d) {
  donors <- setdiff(ctrl_col_ids, d); Yd <- Y[, donors]
  wv <- Variable(ncol(Yd))
  pr <- Problem(Minimize(sum((y_treated[pre_idx] - Yd[pre_idx, ] %*% wv)^2)),
                list(wv >= 0, sum(wv) == 1))
  wh <- as.numeric(solve(pr)$getValue(wv))
  mean((y_treated - as.numeric(Yd %*% wh))[post_idx])
})
cat(sprintf("Full-donor avg ATT:   %+0.4f (%+0.1f%%)\n", avg_att, (exp(avg_att)-1)*100))
cat(sprintf("Drop-one-donor range: %+0.4f to %+0.4f  (%+0.1f%% to %+0.1f%%)\n",
            min(loo_don), max(loo_don), (exp(min(loo_don))-1)*100, (exp(max(loo_don))-1)*100))

# --- Leave-one-treated-zone-out (per-zone estimate robustness) -------------
cat("\n========== Leave-one-treated-zone-out (per-zone avg) ==========\n")
tj <- sapply(seq_along(treat_col_ids), function(i)
  mean(rowMeans(per_zone_mat[, -i, drop = FALSE])[post_idx]))
cat(sprintf("Per-zone avg ATT:       %+0.4f (%+0.1f%%)\n",
            avg_att_perzone, (exp(avg_att_perzone)-1)*100))
cat(sprintf("Drop-one-treated range: %+0.4f to %+0.4f  (%+0.1f%% to %+0.1f%%)\n",
            min(tj), max(tj), (exp(min(tj))-1)*100, (exp(max(tj))-1)*100))

# per-zone avg is the primary SCM estimate (targets the same ATT estimand as the
# TWFE DiD); the pooled aggregate is the conservative floor.
summary_tbl <- data.frame(
  method=c("SCM per-zone avg (primary)","SCM aggregate (pooled)","TWFE DiD (reference)","UWS alt-control","Residential core"),
  logpts=c(avg_att_perzone, avg_att, -0.1420, -0.1438, -0.1440),
  pct=round((exp(c(avg_att_perzone, avg_att, -0.1420, -0.1438, -0.1440))-1)*100, 2))
cat("\n--- Comparison with TWFE / other robustness rows ---\n")
print(summary_tbl)
write.csv(summary_tbl, file.path(BUILD, "scm_summary.csv"), row.names=FALSE)

cat("\nSaved: scm_aggregate_att.csv, scm_aggregate_61.png, scm_gap_61.png, scm_summary.csv\n")
cat("================ SYNTHETIC-CONTROL DIAGNOSTICS COMPLETE ================\n")
