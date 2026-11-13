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
