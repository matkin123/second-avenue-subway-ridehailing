# HonestDiD (Rambachan-Roth 2023) sensitivity analysis on the 61-zone
# event study. Computes robust confidence sets allowing the post-period
# differential trend to be as large as Mbar * the maximum pre-period
# deviation (relative-magnitudes restriction) and as large as M (log-pts/mo)
# of second-difference curvature (smoothness restriction).
#
# Output: honest_did_rm_61.csv / honest_did_sm_61.csv  (bounds per Mbar / M)
#         honest_did_61_rm.png / honest_did_61_sm.png   (sensitivity plots)

suppressWarnings(suppressMessages({
  library(dplyr); library(plm); library(sandwich); library(ggplot2)
  library(HonestDiD)
}))

BUILD <- "analysis"
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids61     <- readRDS(file.path(BUILD, "zones61_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)

md <- full_data %>%
  filter(locationID %in% ids61) %>%
  mutate(Treatment = as.integer(locationID %in% treat_ids),
         log_obs   = log(observations))
bal <- md %>% count(locationID) %>% filter(n == max(n))
md  <- md %>% filter(locationID %in% bal$locationID) %>% arrange(locationID, month)
md$month_factor <- factor(md$month)

# Fit event study (plm silently drops the LAST month as implicit ref)
es <- plm("log_obs ~ Treatment:month_factor",
          data=md, index=c("locationID","month"),
          model="within", effect="twoways")
cf  <- coef(es)
V   <- vcovHC(es, type="HC1", cluster="group")
all_mos <- sort(unique(md$month))
dropped_mo <- setdiff(as.character(all_mos),
                     gsub("Treatment:month_factor","", names(cf)))
cat("plm-dropped (raw implicit ref):", dropped_mo, "\n")

ref_date <- as.Date("2016-12-01")
ref_term <- paste0("Treatment:month_factor", ref_date)

# Build re-referenced betahat vector (nmo-1 long: all months except ref)
# Ordering required by HonestDiD: pre-period first (increasing time), then post-period
nonref_mos <- all_mos[all_mos != ref_date]
nonref_mos <- sort(nonref_mos)
pre_mos    <- nonref_mos[nonref_mos < ref_date]
post_mos   <- nonref_mos[nonref_mos > ref_date]
ord_mos    <- c(pre_mos, post_mos)

L <- matrix(0, nrow=length(ord_mos), ncol=length(cf)); colnames(L)=names(cf)
for (i in seq_along(ord_mos)) {
  mo <- ord_mos[i]; tr <- paste0("Treatment:month_factor", mo)
  if (as.character(mo) != dropped_mo) L[i, tr] <- 1
  L[i, ref_term] <- L[i, ref_term] - 1
}
betahat <- as.numeric(L %*% cf)
sigma   <- L %*% V %*% t(L)

numPre  <- length(pre_mos)
numPost <- length(post_mos)
cat("betahat length:", length(betahat), " | numPre =", numPre, " | numPost =", numPost, "\n\n")

# Original event-study 95% CIs (treated-vs-control gap in each post month vs Dec 2016)
se_vec  <- sqrt(diag(sigma))
print(data.frame(month=ord_mos, est=round(betahat,4),
                 lo=round(betahat-1.96*se_vec,4),
                 hi=round(betahat+1.96*se_vec,4)))

# ---- 1) sensitivity using ATT on l_vec = average of all post-period months ----
l_vec_avg <- rep(1/numPost, numPost)
mbar_grid <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)

cat("\n========== HonestDiD: relative-magnitudes (Mbar) bounds, avg post-ATT ==========\n")
out_rm <- createSensitivityResults_relativeMagnitudes(
  betahat=betahat, sigma=sigma,
  numPrePeriods=numPre, numPostPeriods=numPost,
  l_vec=l_vec_avg, Mbarvec=mbar_grid)
print(out_rm)

# ---- 2) smoothness restriction (M = max second-diff of counterfactual trend) ----
m_grid <- seq(0, 0.05, by=0.01)
cat("\n========== HonestDiD: smoothness (M) bounds, avg post-ATT ==========\n")
out_sm <- createSensitivityResults(
  betahat=betahat, sigma=sigma,
  numPrePeriods=numPre, numPostPeriods=numPost,
  l_vec=l_vec_avg, Mvec=m_grid)
print(out_sm)

# Original (non-robust) CI
orig_ci <- constructOriginalCS(betahat=betahat, sigma=sigma,
                               numPrePeriods=numPre, numPostPeriods=numPost,
                               l_vec=l_vec_avg)
cat("\nOriginal 95% CI for avg post-ATT: [", round(orig_ci$lb,4), ",", round(orig_ci$ub,4), "]\n")

write.csv(cbind(method="relative_magnitudes", out_rm),
          file.path(BUILD, "honest_did_rm_61.csv"), row.names=FALSE)
write.csv(cbind(method="smoothness", out_sm),
          file.path(BUILD, "honest_did_sm_61.csv"), row.names=FALSE)

# ---- Plots ----
pl_rm <- createSensitivityPlot_relativeMagnitudes(
            robustResults=out_rm, originalResults=orig_ci) +
         labs(title="HonestDiD: Avg post-period ATT under relative-magnitudes restriction",
              subtitle="Robust 95% CI stays negative while Mbar * max|pre-dev| < |post-ATT|. x = Mbar.",
              x="Mbar (max post-period bias as multiple of max pre-period deviation)",
              y="Avg post-period ATT (log-points)")
ggsave(file.path(BUILD, "honest_did_61_rm.png"), pl_rm, width=9, height=5, dpi=130)

pl_sm <- createSensitivityPlot(robustResults=out_sm, originalResults=orig_ci) +
         labs(title="HonestDiD: Avg post-period ATT under smoothness restriction",
              subtitle="M = allowed per-period curvature of the counterfactual trend (log-pts).",
              x="M (smoothness bound)", y="Avg post-period ATT (log-points)")
ggsave(file.path(BUILD, "honest_did_61_sm.png"), pl_sm, width=9, height=5, dpi=130)

# Breakdown Mbar: largest Mbar at which the robust CI still excludes 0
rm_tbl <- as.data.frame(out_rm)
rm_tbl$sig_neg <- rm_tbl$ub < 0
if (any(rm_tbl$sig_neg)) {
  breakdown <- max(rm_tbl$Mbar[rm_tbl$sig_neg])
  cat("\n>>> Breakdown Mbar (largest Mbar at which effect stays negative/significant):",
      breakdown, "<<<\n")
} else cat("\n>>> Even Mbar=0 would not give a significant effect — check data <<<\n")

cat("\nSaved: honest_did_rm_61.csv, honest_did_sm_61.csv, honest_did_61_rm.png, honest_did_61_sm.png\n")
cat("================ HonestDiD DIAGNOSTICS COMPLETE ================\n")
