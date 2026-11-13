# Pre-trend diagnostics for the 61-zone DiD (parameterization-corrected).
#
# IMPORTANT: plm::plm(..., model="within", effect="twoways") with
# `Treatment:month_factor` silently drops the LAST factor level (here Nov 2017)
# as the implicit reference, NOT the level named in `relevel(..., ref=...)`.
# Nov 2017 is a holiday spike, so coefficients expressed "vs Nov 2017" look
# artificially inflated. The correct reading is to re-reference each coefficient
# to a chosen pre-period base (Dec 2016), propagating cluster SEs via the
# variance-covariance of the two coefficients.
#
# Outputs:
#   event_study_61.png        : canonical event-study plot, ref = Dec 2016
#   trajectories_61.png       : group-mean log trips, de-meaned
#   trajectories_61_raw.png   : group-mean raw trips/zone/month
#   detrended_effect.png      : observed vs pre-trend extrapolation
#   event_study_coefs_61.csv  : full coef table w/ both parameterizations
#   pretrend_summary_61.csv   : summary tests

suppressWarnings(suppressMessages({
  library(dplyr); library(plm); library(lmtest); library(car)
  library(sandwich); library(broom); library(ggplot2); library(tidyr)
}))

BUILD <- "analysis"
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids61     <- readRDS(file.path(BUILD, "zones61_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)

md <- full_data %>%
  filter(locationID %in% ids61) %>%
  mutate(Treatment = as.integer(locationID %in% treat_ids),
         pre_post  = as.integer(month > as.Date("2016-12-31") & month < as.Date("2017-12-31")),
         log_obs   = log(observations))
bal <- md %>% count(locationID) %>% filter(n == max(n))
md  <- md %>% filter(locationID %in% bal$locationID) %>% arrange(locationID, month)

all_months <- sort(unique(md$month))
treat_date <- as.Date("2017-01-01")
ref_date   <- as.Date("2016-12-01")

# ---- 1) event-study estimation (two-way FE; plm internally drops last mo) -----
md$month_factor <- factor(md$month)
es <- plm("log_obs ~ Treatment:month_factor",
          data=md, index=c("locationID","month"),
          model="within", effect="twoways")
cf  <- coef(es)
V   <- vcovHC(es, type="HC1", cluster="group")
estimated_mos <- as.Date(gsub("Treatment:month_factor","", names(cf)))
dropped_mo    <- setdiff(as.character(all_months), as.character(estimated_mos))
stopifnot(length(dropped_mo)==1)
cat("plm implicitly dropped month:", dropped_mo, "(this is the raw reference)\n")

# ---- 2) re-reference each coef to ref_date (= Dec 2016) -----------------------
ref_term <- paste0("Treatment:month_factor", ref_date)
stopifnot(ref_term %in% names(cf))
nmo <- length(all_months)
# Build contrast matrix L: one row per month, each row = (coef[m] - coef[ref])
L <- matrix(0, nrow=nmo, ncol=length(cf))
colnames(L) <- names(cf)
for (i in seq_len(nmo)) {
  mo <- all_months[i]; tr <- paste0("Treatment:month_factor", mo)
  if (as.character(mo) != dropped_mo) L[i, tr] <- 1
  L[i, ref_term] <- L[i, ref_term] - 1  # subtract ref
}
rownames(L) <- as.character(all_months)
# contrast estimates and SEs
cf_rel <- as.numeric(L %*% cf)
Vc_rel <- L %*% V %*% t(L)
se_rel <- sqrt(diag(Vc_rel))

ec <- data.frame(
  month = all_months,
  estimate  = cf_rel,
  std.error = se_rel,
  conf.low  = cf_rel - 1.96*se_rel,
  conf.high = cf_rel + 1.96*se_rel,
  period    = ifelse(all_months < ref_date, "pre",
             ifelse(all_months == ref_date, "ref", "post")),
  stringsAsFactors = FALSE)

# raw (vs Nov 2017) side-by-side for transparency
ec$raw_vs_dropped <- NA_real_
for (i in seq_len(nmo)) {
  mo <- all_months[i]; tr <- paste0("Treatment:month_factor", mo)
  ec$raw_vs_dropped[i] <- if (as.character(mo)==dropped_mo) 0 else unname(cf[tr])
}

write.csv(ec, file.path(BUILD, "event_study_coefs_61.csv"), row.names=FALSE)

# ---- 3) joint F-test: pre-period = 0 (re-referenced, ref=Dec 2016) ------------
pre_rows  <- which(ec$period == "pre")
post_rows <- which(ec$period == "post")
Lpre <- L[pre_rows, , drop=FALSE]
cf_pre <- Lpre %*% cf
V_pre  <- Lpre %*% V %*% t(Lpre)
Wpre   <- as.numeric(t(cf_pre) %*% solve(V_pre) %*% cf_pre)
df1_pre<- nrow(Lpre); df2 <- es$df.residual
Fpre   <- Wpre/df1_pre; Ppre <- pf(Fpre, df1_pre, df2, lower.tail=FALSE)

Lpost <- L[post_rows, , drop=FALSE]
cf_post<- Lpost %*% cf
V_post <- Lpost %*% V %*% t(Lpost)
Wpost  <- as.numeric(t(cf_post) %*% solve(V_post) %*% cf_post)
df1_post<- nrow(Lpost); Fpost <- Wpost/df1_post
Ppost  <- pf(Fpost, df1_post, df2, lower.tail=FALSE)

# ---- 4) linear pre-trend slope (weighted regression on pre coefs) -------------
pre_ec <- ec %>% filter(period == "pre") %>%
  mutate(month_idx = as.numeric(month - ref_date) / 30.44)
trend_fit <- lm(estimate ~ month_idx, data=pre_ec, weights = 1/(std.error^2))
slope    <- unname(coef(trend_fit)["month_idx"])
slope_se <- unname(summary(trend_fit)$coefficients["month_idx","Std. Error"])
slope_p  <- 2*pt(-abs(slope/slope_se), df=nrow(pre_ec)-2)

# ---- 5) FE DiD with group-specific linear trend -------------------------------
md$t_idx  <- as.numeric(md$month - min(md$month)) / 30.44
f_base    <- "log_obs ~ Treatment:pre_post + factor(month)"
f_trend   <- "log_obs ~ Treatment:pre_post + Treatment:t_idx + factor(month)"
m_base  <- plm(f_base,  data=md, index=c("locationID","month"),
               model="within", effect="twoways")
m_trend <- plm(f_trend, data=md, index=c("locationID","month"),
               model="within", effect="twoways")
Vb <- vcovHC(m_base,  type="HC1", cluster="group")
Vt <- vcovHC(m_trend, type="HC1", cluster="group")
b_base  <- coef(m_base)["Treatment:pre_post"]
b_trend <- coef(m_trend)["Treatment:pre_post"]
se_base_c  <- sqrt(Vb["Treatment:pre_post","Treatment:pre_post"])
se_trend_c <- sqrt(Vt["Treatment:pre_post","Treatment:pre_post"])
tr_slope    <- coef(m_trend)["Treatment:t_idx"]
tr_slope_se <- sqrt(Vt["Treatment:t_idx","Treatment:t_idx"])

# ---- 6) sensitivity: drop month-containing-opening (Dec 2016 OR Jan 2017) -----
md_noDec <- md %>% filter(month != as.Date("2016-12-01"))
m_noDec <- plm(f_base, data=md_noDec, index=c("locationID","month"),
               model="within", effect="twoways")
b_noDec <- coef(m_noDec)["Treatment:pre_post"]
se_noDec<- sqrt(vcovHC(m_noDec, type="HC1", cluster="group")["Treatment:pre_post","Treatment:pre_post"])

md_noJan <- md %>% filter(month != as.Date("2017-01-01"))
m_noJan <- plm(f_base, data=md_noJan, index=c("locationID","month"),
               model="within", effect="twoways")
b_noJan <- coef(m_noJan)["Treatment:pre_post"]
se_noJan<- sqrt(vcovHC(m_noJan, type="HC1", cluster="group")["Treatment:pre_post","Treatment:pre_post"])

# ---- 7) summary CSV -----------------------------------------------------------
summary_tbl <- data.frame(
  stat = c("pre_mean_coef","pre_max_abs_coef","pre_n_CIs_excluding_0","pre_n",
           "post_mean_coef","post_max_abs_coef","post_n_negative","post_n",
           "jointF_pre_stat","jointF_pre_pvalue",
           "jointF_post_stat","jointF_post_pvalue",
           "linear_pre_slope_per_month","linear_pre_slope_SE","linear_pre_slope_p",
           "FE_DiD_baseline","FE_DiD_baseline_clusterSE",
           "FE_DiD_with_group_trend","FE_DiD_with_group_trend_clusterSE",
           "treated_linear_trend_per_month","treated_linear_trend_SE",
           "FE_DiD_drop_Dec2016","FE_DiD_drop_Dec2016_clusterSE",
           "FE_DiD_drop_Jan2017","FE_DiD_drop_Jan2017_clusterSE"),
  value = c(mean(ec$estimate[ec$period=="pre"]),
            max(abs(ec$estimate[ec$period=="pre"])),
            sum(ec$conf.low[ec$period=="pre"] > 0 | ec$conf.high[ec$period=="pre"] < 0),
            sum(ec$period=="pre"),
            mean(ec$estimate[ec$period=="post"]),
            max(abs(ec$estimate[ec$period=="post"])),
            sum(ec$estimate[ec$period=="post"] < 0),
            sum(ec$period=="post"),
            Fpre, Ppre, Fpost, Ppost,
            slope, slope_se, slope_p,
            b_base, se_base_c,
            b_trend, se_trend_c,
            unname(tr_slope), unname(tr_slope_se),
            b_noDec, se_noDec, b_noJan, se_noJan))
write.csv(summary_tbl, file.path(BUILD, "pretrend_summary_61.csv"), row.names=FALSE)

cat("\n========== PRE-TREND DIAGNOSTICS (61 zones, ref = Dec 2016) ==========\n\n")
cat(sprintf("Pre-period (n=%d):  mean coef = %+0.4f  max|coef| = %0.4f  CIs excluding 0 = %d\n",
    sum(ec$period=="pre"), mean(ec$estimate[ec$period=="pre"]),
    max(abs(ec$estimate[ec$period=="pre"])),
    sum(ec$conf.low[ec$period=="pre"] > 0 | ec$conf.high[ec$period=="pre"] < 0)))
cat(sprintf("Post-period (n=%d): mean coef = %+0.4f  max|coef| = %0.4f  negative in %d/%d\n\n",
    sum(ec$period=="post"), mean(ec$estimate[ec$period=="post"]),
    max(abs(ec$estimate[ec$period=="post"])),
    sum(ec$estimate[ec$period=="post"] < 0), sum(ec$period=="post")))

cat(sprintf("Joint F pre-period = 0:   F(%d,%d) = %.3f  p = %s\n",
            df1_pre, df2, Fpre, format.pval(Ppre)))
cat(sprintf("Joint F post-period = 0:  F(%d,%d) = %.3f  p = %s\n\n",
            df1_post, df2, Fpost, format.pval(Ppost)))

cat(sprintf("Linear pre-trend slope: %+0.5f log-pts/month (SE %0.5f, p = %s)\n\n",
            slope, slope_se, format.pval(slope_p)))

cat("DiD baseline vs with group-specific linear trend vs month-drop sensitivities:\n")
cat(sprintf("  baseline FE DiD:                     %+0.4f (cluster SE %0.4f)  => %+0.1f%%\n",
            b_base,  se_base_c,  (exp(b_base)-1)*100))
cat(sprintf("  FE DiD + Treatment x linear trend:   %+0.4f (cluster SE %0.4f)  => %+0.1f%%\n",
            b_trend, se_trend_c, (exp(b_trend)-1)*100))
cat(sprintf("     [Treatment trend coef:            %+0.4f (SE %0.4f)]\n",
            tr_slope, tr_slope_se))
cat(sprintf("  FE DiD dropping Dec 2016:            %+0.4f (cluster SE %0.4f)  => %+0.1f%%\n",
            b_noDec, se_noDec, (exp(b_noDec)-1)*100))
cat(sprintf("  FE DiD dropping Jan 2017:            %+0.4f (cluster SE %0.4f)  => %+0.1f%%\n",
            b_noJan, se_noJan, (exp(b_noJan)-1)*100))

# ---- 8) PLOTS -----------------------------------------------------------------
pre_fit_df <- data.frame(month = all_months,
                         month_idx = as.numeric(all_months - ref_date)/30.44)
pre_fit_df$estimate <- predict(trend_fit, newdata=pre_fit_df)

p1 <- ggplot(ec, aes(month, estimate, color=period)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray40") +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=10) +
  geom_point(size=2.8) +
  geom_line(data=pre_fit_df, aes(month, estimate), inherit.aes=FALSE,
            color="steelblue", linetype="dashed", linewidth=0.7) +
  annotate("text", x=as.Date("2016-01-15"), y=0.08,
           label=sprintf("pre-period slope: %+0.4f/mo (p=%s)",
                         slope, format.pval(slope_p)),
           hjust=0, size=3.3, color="steelblue") +
  annotate("text", x=as.Date("2017-07-01"), y=0.08,
           label=sprintf("post mean: %+0.3f  (joint F=%.1f, p<0.001)",
                         mean(ec$estimate[ec$period=="post"]), Fpost),
           hjust=0.5, size=3.3, color="firebrick") +
  scale_color_manual(values=c("pre"="steelblue","ref"="black","post"="firebrick"),
                     breaks=c("pre","ref","post")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Event Study, 61-zone sample — parameterization-corrected",
       subtitle=sprintf("Two-way FE, log outcome, ref = Dec 2016. Pre-trend joint F = %.2f (p = %s).",
                        Fpre, format.pval(Ppre)),
       x=NULL, y="Treated-vs-control gap, log-points (vs Dec 2016)",
       color="period") + theme_minimal()
ggsave(file.path(BUILD, "event_study_61.png"), p1, width=10, height=5.5, dpi=130)

# group-mean trajectories (de-meaned within group)
agg <- md %>%
  group_by(Treatment, month) %>%
  summarise(mean_log = mean(log_obs), mean_obs = mean(observations), .groups="drop") %>%
  group_by(Treatment) %>%
  mutate(norm_log = mean_log - mean(mean_log[month < treat_date])) %>% ungroup() %>%
  mutate(grp = ifelse(Treatment==1, "Treated (7 zones)", "Control (54 zones)"))

p2 <- ggplot(agg, aes(month, norm_log, color=grp)) +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  scale_color_manual(values=c("Treated (7 zones)"="firebrick",
                              "Control (54 zones)"="steelblue")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Group-mean log trips, normalized to pre-period mean (61-zone sample)",
       subtitle="Parallel-trends check: do the two lines move together before Jan 2017?",
       x=NULL, y="log trips - group's pre-period mean", color=NULL) +
  theme_minimal() + theme(legend.position="top")
ggsave(file.path(BUILD, "trajectories_61.png"), p2, width=10, height=5.5, dpi=130)

agg_raw <- md %>%
  group_by(Treatment, month) %>%
  summarise(mean_obs = mean(observations), .groups="drop") %>%
  mutate(grp = ifelse(Treatment==1, "Treated (7 zones)", "Control (54 zones)"))

p3 <- ggplot(agg_raw, aes(month, mean_obs/1000, color=grp)) +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  scale_color_manual(values=c("Treated (7 zones)"="firebrick",
                              "Control (54 zones)"="steelblue")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Mean ridehailing trips per zone per month, by group (61-zone sample)",
       subtitle="Raw levels — shows level gap and 2016 trajectory of each group",
       x=NULL, y="trips / zone / month (thousands)", color=NULL) +
  theme_minimal() + theme(legend.position="top")
ggsave(file.path(BUILD, "trajectories_61_raw.png"), p3, width=10, height=5.5, dpi=130)

# detrended: observed vs pre-trend extrapolation
proj <- data.frame(month = all_months,
                   month_idx = as.numeric(all_months - ref_date)/30.44)
proj$projection <- predict(trend_fit, newdata=proj)
ec2 <- ec %>% left_join(proj %>% select(month, projection), by="month") %>%
  mutate(detrended = estimate - projection)

p4 <- ggplot(ec2, aes(month)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray40") +
  geom_vline(xintercept=treat_date, color="red", linetype="dashed") +
  geom_point(aes(y=estimate, color="Observed event-study coef"), size=2.3) +
  geom_line(aes(y=projection, color="Pre-trend extrapolation"), linewidth=0.9) +
  geom_point(aes(y=detrended, color="Detrended (obs - extrap)"), size=2.3, shape=17) +
  scale_color_manual(values=c("Observed event-study coef"="firebrick",
                              "Pre-trend extrapolation"="steelblue",
                              "Detrended (obs - extrap)"="darkgreen")) +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Pre-trend extrapolation and detrended effect (crude Honest-DiD style bound)",
       subtitle="Blue line: what the treated-vs-control gap WOULD be if the 2016 trajectory simply continued.",
       x=NULL, y="log-points vs Dec 2016", color=NULL) +
  theme_minimal() + theme(legend.position="top")
ggsave(file.path(BUILD, "detrended_effect.png"), p4, width=10, height=5.5, dpi=130)

cat("\nSaved plots: event_study_61.png, trajectories_61.png, trajectories_61_raw.png, detrended_effect.png\n")
cat("Saved tables: event_study_coefs_61.csv, pretrend_summary_61.csv\n")
cat("\n================ PRE-TREND DIAGNOSTICS COMPLETE ================\n")
