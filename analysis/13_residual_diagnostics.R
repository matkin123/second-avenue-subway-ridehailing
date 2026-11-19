# Residual diagnostics for the two headline specifications:
#   (1) Synthetic control, per-zone  (primary SCM, 09_synthetic_control.R)
#   (2) All Manhattan, 61-zone two-way FE  (headline, 05_did_main_61zones.R)
#
# Produces two figures:
#   resid_over_time_61_scm.png  -- residuals by month (parallel-trend /
#                                  serial-correlation check)
#   resid_vs_fitted_61_scm.png  -- residual vs fitted (heteroskedasticity /
#                                  outlier check)
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
  library(plm); library(lmtest); library(sandwich)
}))

BUILD      <- "analysis"
treat_ids  <- c(75,262,263,140,141,236,237)
treat_date <- as.Date("2017-01-01")

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

# ================= FORMAL INFERENCE DIAGNOSTICS (61-zone TWFE) =================
# These replace "eyeball the residual cloud" with actual tests. The cluster-robust
# SE used throughout is the *remedy* for serial correlation and heteroskedasticity;
# the tests below are what justify reaching for it.
cat("\n========== FORMAL INFERENCE DIAGNOSTICS (61-zone TWFE) ==========\n")

## (a) Serial correlation within a zone across the 24 months (Bertrand-Duflo-
##     Mullainathan 2004 is the canonical DiD concern).
pmd   <- pdata.frame(md, index = c("locationID","month"))
sc_fe <- plm(log_obs ~ Treatment:pre_post + factor(month),
             data = pmd, model = "within", effect = "individual")
ww <- tryCatch(pwartest(sc_fe), error = function(e) NULL)  # Wooldridge FE serial-corr
bg <- tryCatch(pbgtest(sc_fe),  error = function(e) NULL)  # Breusch-Godfrey LM
cat("-- Serial correlation (within zone, over months) --\n")
if (!is.null(ww)) cat(sprintf("  Wooldridge (pwartest):     F = %6.2f  p = %s\n",
                              ww$statistic, format.pval(ww$p.value)))
if (!is.null(bg)) cat(sprintf("  Breusch-Godfrey (pbgtest): chisq = %5.2f  p = %s\n",
                              bg$statistic, format.pval(bg$p.value)))

## (b) Cross-sectional / spatial dependence of the residuals (Pesaran CD). Month
##     FE remove common shocks but NOT local shocks shared by neighbouring zones.
w61  <- readRDS(file.path(BUILD, "wts61.rds"))
pr   <- pdata.frame(md %>% mutate(resid_fe = r61$resid) %>%
                      select(locationID, month, resid_fe),
                    index = c("locationID","month"))
cd_gen <- tryCatch(pcdtest(pr$resid_fe, test = "cd"),
                   error = function(e) NULL)
cd_loc <- tryCatch(pcdtest(pr$resid_fe, test = "cd", w = w61$mat),
                   error = function(e) NULL)
cat("-- Cross-sectional dependence of residuals (Pesaran CD) --\n")
if (!is.null(cd_gen)) cat(sprintf("  Global CD (all zone pairs):  z = %6.2f  p = %s\n",
                                  cd_gen$statistic, format.pval(cd_gen$p.value)))
if (!is.null(cd_loc)) cat(sprintf("  Local CD (queen-contiguity): z = %6.2f  p = %s\n",
                                  cd_loc$statistic, format.pval(cd_loc$p.value)))

## (c) Wild cluster bootstrap-t for the DiD coefficient (Cameron-Gelbach-Miller),
##     Rademacher weights, null imposed -- robust to the modest cluster count (61).
wild_cluster_p <- function(data, B = 999, seed = 1) {
  set.seed(seed)
  X  <- model.matrix(~ factor(locationID) + factor(month) + Treatment:pre_post, data)
  y  <- data$log_obs
  cl <- as.character(data$locationID)
  xk <- which(colnames(X) == "Treatment:pre_post")
  XtXinv <- solve(crossprod(X))
  h  <- as.numeric(X %*% XtXinv[, xk])            # beta_xk = h'y; cluster SE uses h
  N  <- nrow(X); k <- ncol(X); G <- length(unique(cl))
  cadj  <- (G/(G-1)) * ((N-1)/(N-k))
  cr_se <- function(res) sqrt(cadj * sum(rowsum(h * res, cl)^2))
  b_hat    <- sum(h * y)
  res_full <- as.numeric(y - X %*% (XtXinv %*% crossprod(X, y)))
  t_hat    <- b_hat / cr_se(res_full)
  Xr   <- X[, -xk, drop = FALSE]                  # restricted fit imposes DiD = 0
  yhr  <- as.numeric(Xr %*% solve(crossprod(Xr), crossprod(Xr, y)))
  resr <- y - yhr
  cl_u <- unique(cl)
  t_star <- numeric(B)
  for (b in seq_len(B)) {
    e  <- setNames(sample(c(-1, 1), G, replace = TRUE), cl_u)
    v  <- resr * e[cl]
    rb <- as.numeric(v - X %*% (XtXinv %*% crossprod(X, v)))    # M_X v
    t_star[b] <- sum(h * v) / cr_se(rb)
  }
  list(t = t_hat, p = mean(abs(t_star) >= abs(t_hat)), B = B)
}
wc <- wild_cluster_p(md)
cat("-- Wild cluster bootstrap (B =", wc$B, ", Rademacher, null imposed) --\n")
cat(sprintf("  DiD t-stat = %.2f  wild-cluster p = %s\n", wc$t, format.pval(wc$p)))

## (d) Conley spatial-HAC SEs: robust to BOTH spatial correlation across zones
##     AND serial correlation within zone (Conley 1999 / Hsiang 2010). This is the
##     direct fix for the cross-sectional dependence the CD test flags -- zone
##     clustering does not address correlation between *different* zones.
suppressWarnings(suppressMessages(library(sf)))
shp  <- st_read(file.path(BUILD, "geo/taxi_zones/taxi_zones/taxi_zones.shp"), quiet = TRUE)
cent <- shp[shp$LocationID %in% ids61, ]
xy   <- st_coordinates(st_centroid(st_geometry(cent))) * 0.0003048   # US-ft -> km
rownames(xy) <- as.character(cent$LocationID)

Xc <- model.matrix(~ factor(locationID) + factor(month) + Treatment:pre_post, md)
yc <- md$log_obs
xkc<- which(colnames(Xc) == "Treatment:pre_post")
XtXic <- solve(crossprod(Xc))
hc <- as.numeric(Xc %*% XtXic[, xkc])            # score weights; b = hc'y
uc <- as.numeric(yc - Xc %*% (XtXic %*% crossprod(Xc, yc)))
ec <- hc * uc                                    # scalar score per obs
Nc <- nrow(Xc); kc <- ncol(Xc); b_hat_c <- sum(hc * yc)

zord <- sort(unique(md$locationID)); mord <- sort(unique(md$month))
E <- matrix(0, length(zord), length(mord), dimnames = list(zord, as.character(mord)))
E[cbind(match(md$locationID, zord), match(md$month, mord))] <- ec
Dkm <- as.matrix(dist(xy[as.character(zord), ]))   # zone-pair distances (km)

conley_se <- function(dist_cut, lag, dofc = Nc/(Nc-kc)) {
  S  <- 1 - Dkm/dist_cut; S[S < 0] <- 0            # Bartlett kernel in space (keeps matrix)
  G  <- t(E) %*% S %*% E                           # month x month
  Kt <- outer(seq_along(mord), seq_along(mord),
              function(a, b) pmax(0, 1 - abs(a - b)/(lag + 1)))  # Bartlett in time
  sqrt(sum(Kt * G) * dofc)
}
cat("-- Conley spatial-HAC SEs for DiD coef (b =", sprintf("%+0.4f", b_hat_c), ") --\n")
cgrid <- expand.grid(dist_km = c(1, 2, 5), lag = c(1, 3, 6))
cres  <- do.call(rbind, lapply(seq_len(nrow(cgrid)), function(i) {
  se <- conley_se(cgrid$dist_km[i], cgrid$lag[i])
  data.frame(dist_km = cgrid$dist_km[i], lag = cgrid$lag[i], se = se,
             t = b_hat_c/se, p = 2*pnorm(-abs(b_hat_c/se)))
}))
for (i in seq_len(nrow(cres)))
  cat(sprintf("  cutoff %2.0f km, lag %d:  SE = %0.4f  t = %6.2f  p = %s\n",
              cres$dist_km[i], cres$lag[i], cres$se[i], cres$t[i], format.pval(cres$p[i])))
write.csv(cres, file.path(BUILD, "resid_conley_se_61.csv"), row.names = FALSE)

## (e) Influence: leave-one-treated-zone-out jackknife on the DiD coefficient.
##     With only 7 treated zones the estimate could hinge on one corridor zone.
did_coef <- function(d) coef(lm(log_obs ~ factor(locationID) + factor(month) +
                                  Treatment:pre_post, d))["Treatment:pre_post"]
b_all <- did_coef(md)
jk    <- sapply(treat_ids, function(z) did_coef(md %>% filter(locationID != z)))
cat("-- Leave-one-treated-zone-out jackknife (influence) --\n")
cat(sprintf("  Full-sample DiD: %+0.4f (%+0.1f%%)\n", b_all, (exp(b_all)-1)*100))
cat(sprintf("  Drop-one range:  %+0.4f to %+0.4f  (%+0.1f%% to %+0.1f%%)\n",
            min(jk), max(jk), (exp(min(jk))-1)*100, (exp(max(jk))-1)*100))
write.csv(data.frame(dropped_zone = treat_ids, did = as.numeric(jk),
                     pct = (exp(as.numeric(jk))-1)*100),
          file.path(BUILD, "resid_influence_jackknife_61.csv"), row.names = FALSE)

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
                  "All Manhattan - 61 zones TWFE")
all_r <- bind_rows(r61, rscm) %>%
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
       subtitle = paste0("Points = zone-months; bold lines = group monthly mean.\n",
                         "TWFE: effect removed, residuals stay flat. ",
                         "SCM (treated - synthetic): ~0 pre-2017, post-2017 drop IS the ATT."),
       x = NULL, y = "residual (log-points)", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
ggsave(file.path(BUILD, "resid_over_time_61_scm.png"), p1, width = 8.5, height = 5, dpi = 130)

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
ggsave(file.path(BUILD, "resid_vs_fitted_61_scm.png"), p2, width = 8.5, height = 5, dpi = 130)

cat("\nSaved: resid_over_time_61_scm.png, resid_vs_fitted_61_scm.png\n")
