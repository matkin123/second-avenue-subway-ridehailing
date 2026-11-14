# Residual diagnostics for the three headline specifications:
#   (1) Synthetic control, per-zone  (primary SCM, 09_synthetic_control.R)
#   (2) All Manhattan, 61-zone two-way FE  (headline, 05_did_main_61zones.R)
#   (3) Residential core, 58-zone two-way FE  (drop 153/88/186, robustness row)
#
# Produces two figures:
#   resid_over_time_61_58_scm.png  -- residuals by month (parallel-trend /
#                                     serial-correlation check)
#   resid_vs_fitted_61_58_scm.png  -- residual vs fitted (heteroskedasticity /
#                                     outlier check)
#
# Residual definitions differ by method and this is intentional:
#   TWFE  -> residual = log_obs - LSDV fitted (zone FE + month FE + ATT removed).
#            A clean diagnostic should show NO structure and NO pre-trend.
#   SCM   -> residual = observed - synthetic. Pre-period = genuine fit error
#            (should hug 0); post-period = the ATT itself (the -10.6%), so SCM
#            treated residuals DROP after Jan-2017 by construction.
# In resid-vs-fitted the SCM panel is restricted to the pre-period, the only
# window where "residual = fit error" holds for all three.

suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(CVXR); library(ggplot2)
}))

BUILD      <- "analysis"
treat_ids  <- c(75,262,263,140,141,236,237)
treat_date <- as.Date("2017-01-01")
rc_drop    <- c(153,88,186)   # Marble Hill, FiDi-S, Penn Sq

full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids61     <- readRDS(file.path(BUILD, "zones61_ids.rds"))

md <- full_data %>%
  filter(locationID %in% ids61) %>%
  mutate(month     = as.Date(month),
         Treatment = as.integer(locationID %in% treat_ids),
         pre_post  = as.integer(month > as.Date("2016-12-31") & month < as.Date("2017-12-31")),
         log_obs   = log(observations))
bal <- md %>% count(locationID) %>% filter(n == max(n))
md  <- md %>% filter(locationID %in% bal$locationID) %>% arrange(locationID, month)
cat("balanced panel:", length(unique(md$locationID)), "zones x",
    length(unique(md$month)), "months =", nrow(md), "rows\n")

# ---- TWFE residuals via LSDV (identical coef to plm within, clean row order) ----
twfe_resid <- function(data, label) {
  fit <- lm(log_obs ~ factor(locationID) + factor(month) + Treatment:pre_post, data = data)
  b   <- coef(fit)["Treatment:pre_post"]
  cat(sprintf("%-28s DiD coef = %+0.4f (%+0.1f%%), n = %d\n",
              label, b, (exp(b)-1)*100, nrow(data)))
  data %>% transmute(model = label, month, locationID,
                     group  = ifelse(Treatment == 1, "Treated", "Control"),
                     period = ifelse(pre_post == 1, "Post", "Pre"),
                     fitted = as.numeric(fitted(fit)),
                     resid  = as.numeric(residuals(fit)))
}

r61 <- twfe_resid(md,                                  "All Manhattan - 61 zones TWFE")
r58 <- twfe_resid(md %>% filter(!(locationID %in% rc_drop)),
                                                       "Residential core - 58 zones")

# ---- SCM per-zone residuals (replicates 09_synthetic_control.R) ----
Y_wide <- md %>% select(month, locationID, log_obs) %>%
  pivot_wider(names_from = locationID, values_from = log_obs) %>% arrange(month)
months <- Y_wide$month
Y <- as.matrix(Y_wide %>% select(-month)); colnames(Y) <- colnames(Y_wide)[-1]

treat_col <- as.character(treat_ids)
ctrl_col  <- setdiff(colnames(Y), treat_col)
pre_idx   <- which(months <  treat_date)
Yctrl     <- Y[, ctrl_col]; n_ctrl <- ncol(Yctrl)

scm_rows <- list()
for (tz in treat_col) {
  y_tz <- Y[, tz]
  w    <- Variable(n_ctrl)
  prob <- Problem(Minimize(sum((y_tz[pre_idx] - Yctrl[pre_idx, ] %*% w)^2)),
                  list(w >= 0, sum(w) == 1))
  w_hat <- as.numeric(solve(prob)$getValue(w))
  synth <- as.numeric(Yctrl %*% w_hat)
  scm_rows[[tz]] <- data.frame(model = "Synthetic control - per-zone",
                               month = months, locationID = as.integer(tz),
                               group = "Treated",
                               period = ifelse(months < treat_date, "Pre", "Post"),
                               fitted = synth, resid = y_tz - synth)
}
rscm <- bind_rows(scm_rows)
cat(sprintf("%-28s pre-period mean |resid| = %0.4f | post-period mean resid = %+0.4f (%+0.1f%%)\n",
            "SCM per-zone", mean(abs(rscm$resid[rscm$period=="Pre"])),
            mean(rscm$resid[rscm$period=="Post"]),
            (exp(mean(rscm$resid[rscm$period=="Post"]))-1)*100))

model_levels <- c("Synthetic control - per-zone",
                  "All Manhattan - 61 zones TWFE",
                  "Residential core - 58 zones")
all_r <- bind_rows(r61, r58, rscm) %>%
  mutate(model = factor(model, levels = model_levels))

# ================= VIEW 1: residuals over time =================
monthly_mean <- all_r %>% group_by(model, group, month) %>%
  summarise(resid = mean(resid), .groups = "drop")

p1 <- ggplot(all_r, aes(month, resid, color = group)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.3) +
  geom_vline(xintercept = treat_date, color = "red", linetype = "dashed", size = 0.5) +
  geom_point(alpha = 0.18, size = 0.9) +
  geom_line(data = monthly_mean, aes(group = group), size = 1.2) +
  geom_point(data = monthly_mean, aes(group = group), size = 1.6) +
  facet_wrap(~model, nrow = 1) +
  scale_color_manual(values = c("Treated" = "firebrick", "Control" = "steelblue")) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "4 months") +
  labs(title = "Residuals over time, by specification",
       subtitle = paste0("Points = zone-months; bold lines = group monthly mean. ",
                         "TWFE residuals have the treatment effect removed (should be flat, no pre-trend).\n",
                         "SCM residual = treated - synthetic: pre-period is fit error (~0), post-period IS the ATT, so it drops after Jan-2017."),
       x = NULL, y = "residual (log-points)", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
ggsave(file.path(BUILD, "resid_over_time_61_58_scm.png"), p1, width = 12, height = 5, dpi = 130)

# ================= VIEW 2: residuals vs fitted =================
# SCM restricted to pre-period (the only window where residual = fit error).
rvf <- all_r %>%
  filter(!(model == "Synthetic control - per-zone" & period == "Post"))

p2 <- ggplot(rvf, aes(fitted, resid, color = group)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.3) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "black", size = 0.8, formula = y ~ x) +
  facet_wrap(~model, nrow = 1, scales = "free_x") +
  scale_color_manual(values = c("Treated" = "firebrick", "Control" = "steelblue")) +
  labs(title = "Residuals vs fitted values, by specification",
       subtitle = paste0("Flat even band around 0 = homoskedastic, no outliers. ",
                         "Black line = loess. SCM panel is pre-period only (fit residuals)."),
       x = "fitted value (log trips)", y = "residual (log-points)", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
ggsave(file.path(BUILD, "resid_vs_fitted_61_58_scm.png"), p2, width = 12, height = 5, dpi = 130)

cat("\nSaved: resid_over_time_61_58_scm.png, resid_vs_fitted_61_58_scm.png\n")
