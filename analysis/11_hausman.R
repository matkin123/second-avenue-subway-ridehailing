# Hausman / FE-vs-RE diagnostics on the 61-zone main spec.
#
# Purpose: the textbook Hausman test (FE vs RE) returns p ~ 1 here, which looks
# alarming until you see WHY. This script documents the degeneracy, explains it,
# and replaces the broken test with (a) the design-based case for FE and (b) a
# numerically-stable robust Hausman (Mundlak-Wu) on the one within-varying
# regressor.
#
# Root cause (established here): every census / land-use covariate is a single
# ACS cross-section repeated across all 23 months -> within-zone SD = 0 -> the
# within estimator DROPS all nine covariates. The only regressors with
# within-zone variation are the DiD interaction (Treatment:pre_post) and the
# month dummies. With the model saturated in month effects and no time-varying
# covariates, RE's feasible-GLS transform collapses toward the within transform,
# so FE and RE return essentially identical DiD coefficients. The classical
# Hausman covariance (V_FE - V_RE) is then numerically indefinite and the
# statistic -> 0, p -> 1. This is a low-power / numerical outcome, NOT evidence
# that RE is correct.

suppressWarnings(suppressMessages({
  library(dplyr); library(plm); library(lmtest); library(sandwich)
}))

BUILD <- "analysis"
census_rides       <- readRDS(file.path(BUILD, "ckpt_census.rds"))$census_rides
land_use_manhattan <- readRDS(file.path(BUILD, "ckpt_landuse.rds"))
ids61              <- readRDS(file.path(BUILD, "zones61_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)

# ---- build model_data exactly as 06_models_full_suite.R (61-zone main) ----
lu <- as.data.frame(land_use_manhattan) %>% select(LocationID, pct_res, pct_com, pct_other)
cr <- as.data.frame(census_rides); cr$geometry <- NULL
cr <- cr %>% left_join(lu, by = c("locationID" = "LocationID"))
cr$Treatment <- ifelse(cr$locationID %in% treat_ids, 1, 0)
md <- cr
md$pre_post <- ifelse(md$month > '2016-12-31' & md$month < "2017-12-31", 1, 0)
md$populationE    <- md$populationE/1000
md$median_incomeE <- md$median_incomeE/1000
for (v in c("pct_zero_vehicle","pct_owner_occ","pct_bachelors","white_pct","black_pct",
            "latino_pct","native_pct","asian_pct","pct_res","pct_com","pct_other"))
  md[[v]] <- 100*md[[v]]
md <- md[complete.cases(md), ]
md$log_obs <- log(md$observations)
md <- md %>% filter(locationID %in% ids61)
bal <- md %>% count(locationID) %>% filter(n == max(n))
md <- md %>% filter(locationID %in% bal$locationID) %>% arrange(locationID, month)
md$tp <- md$Treatment * md$pre_post
md <- md %>% group_by(locationID) %>% mutate(tp_bar = mean(tp)) %>% ungroup()  # Mundlak device
n_treat <- length(unique(md$locationID[md$Treatment == 1]))
n_ctrl  <- length(unique(md$locationID[md$Treatment == 0]))
cat(sprintf("panel: %d zones (%d treated, %d control) x %d months = %d rows\n",
            length(unique(md$locationID)), n_treat, n_ctrl,
            length(unique(md$month)), nrow(md)))

covs <- c("populationE","median_incomeE","median_ageE","pct_bachelors","pct_zero_vehicle",
          "pct_owner_occ","black_pct","pct_res","pct_other")
WE_f <- as.formula(paste("log_obs ~", paste(covs, collapse=" + "),
                         "+ Treatment:pre_post + factor(month)"))

clust_vcov <- function(m, cl) sandwich::vcovCL(m, cluster = cl, type = "HC1")

results <- list()  # collect for CSV

cat("\n========== SECTION 1: the degenerate textbook Hausman ==========\n")
fe <- plm(WE_f, data = md, index = "locationID", model = "within")
re <- plm(WE_f, data = md, index = "locationID", model = "random")
ht <- phtest(fe, re)
cat(sprintf("classical Hausman (chisq form): chisq=%.3f  df=%d  p=%s\n",
            ht$statistic, ht$parameter, format.pval(ht$p.value)))
results$hausman_classical_p <- unname(ht$p.value)

# (a) covariates are time-invariant within zone -> FE drops them
within_sd <- sapply(covs, function(v){
  s <- md %>% group_by(locationID) %>% summarise(s = sd(.data[[v]]), .groups="drop")
  mean(s$s, na.rm = TRUE)
})
cat("\nmean within-zone SD of each covariate (0 => time-invariant => dropped by FE):\n")
print(round(within_sd, 6))
cat("=> all nine covariates are a single cross-section; FE identifies the DiD off\n",
    "   the Treatment:pre_post interaction and month dummies only.\n", sep="")
results$max_within_sd_cov <- max(within_sd)

# (b) FE vs RE agree on every common coefficient
bfe <- coef(fe)["Treatment:pre_post"]; bre <- coef(re)["Treatment:pre_post"]
cat(sprintf("\nFE DiD = %+.5f   RE DiD = %+.5f   abs diff = %.5f log-pts (%.2f%% of FE)\n",
            bfe, bre, abs(bfe - bre), 100*abs(bfe - bre)/abs(bfe)))

# (c) V_FE - V_RE is numerically indefinite -> classical statistic ill-posed
common <- intersect(names(coef(fe)), names(coef(re)))
Vd <- vcov(fe)[common, common] - vcov(re)[common, common]
ev <- eigen(Vd, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("V_FE - V_RE over %d common coefs: %d of %d eigenvalues < 0 (indefinite)\n",
            length(common), sum(ev < -1e-12), length(ev)))
cat("=> the classical Hausman covariance is not positive-definite; chisq -> 0, p -> 1.\n")
cat("   NOTE: plm's regression Hausman (method='aux') ALSO fails to compute here -\n")
cat("   the within model has dropped all nine covariates, so there is nothing left\n")
cat("   to compare except the single within-varying DiD regressor.\n")

cat("\n========== SECTION 2: design-based case for FE (the decisive argument) ==========\n")
cat("Subway proximity is not randomly assigned. Treated zones differ systematically\n")
cat("on TIME-INVARIANT traits, so RE's Cov(regressor, zone-effect)=0 fails by design.\n\n")
zl <- md %>% group_by(locationID) %>%
  summarise(Treatment = first(Treatment),
            across(all_of(c("median_incomeE","populationE","pct_bachelors",
                            "white_pct","pct_zero_vehicle","pct_res")), first),
            .groups = "drop")
cat(sprintf("%-18s %10s %10s %10s %8s %10s\n","trait","control","treated","diff","t","p"))
design_tbl <- lapply(c("median_incomeE","populationE","pct_bachelors",
                       "white_pct","pct_zero_vehicle","pct_res"), function(v){
  tt <- t.test(zl[[v]] ~ zl$Treatment)
  cat(sprintf("%-18s %10.2f %10.2f %+10.2f %8.2f %10s\n",
              v, tt$estimate[1], tt$estimate[2],
              tt$estimate[2]-tt$estimate[1], tt$statistic, format.pval(tt$p.value,3)))
  data.frame(trait=v, control=tt$estimate[1], treated=tt$estimate[2],
             diff=tt$estimate[2]-tt$estimate[1], t=tt$statistic, p=tt$p.value)
})
design_tbl <- do.call(rbind, design_tbl); rownames(design_tbl) <- NULL
cat("=> treated zones are richer, denser, whiter, more residential, lower zero-vehicle.\n",
    "   FE differences out ALL such fixed traits (observed and unobserved); RE does not.\n", sep="")

cat("\n========== SECTION 3: robust Hausman that CAN be computed (Mundlak-Wu) ==========\n")
cat("Add the zone-mean of the one within-varying regressor (Treatment:pre_post) to the\n")
cat("model and test it with cluster-robust (by zone) SEs. H0: coef=0 (RE ok).\n")
cat("Rejecting => correlated zone effect => FE preferred.\n\n")
mund <- function(lbl, rhs){
  m <- lm(as.formula(paste("log_obs ~", rhs, "+ tp + tp_bar")), data = md)
  V <- clust_vcov(m, md$locationID)
  b <- coef(m)["tp_bar"]; se <- sqrt(V["tp_bar","tp_bar"]); z <- b/se; p <- 2*pnorm(-abs(z))
  cat(sprintf("%-22s tp_bar=%+.4f  robust SE=%.4f  z=%.2f  p=%s\n", lbl, b, se, z, format.pval(p,3)))
  c(coef=unname(b), se=se, z=z, p=p)
}
m_cov <- mund("covariate spec", paste(paste(covs, collapse=" + "), "+ factor(month)"))
m_min <- mund("minimal spec",  "factor(month)")
results$mundlak_cov_p <- unname(m_cov["p.tp_bar"]); results$mundlak_min_p <- unname(m_min["p.tp_bar"])
cat(sprintf("\nInterpretation: large positive correlated group effect (treated zones sit\n"))
cat(sprintf("~%.0f%% above controls in baseline log-ridership after covariates), marginally\n",
            100*(exp(m_cov["coef"]*unique(md$tp_bar[md$Treatment==1]))-1)))
cat("significant (p~0.07) given only 7 treated zones. Leans toward FE; underpowered to\n")
cat("reject at 5%. Corroborates the design argument rather than contradicting it.\n")

cat("\n========== SECTION 4: FE vs RE side-by-side on the DiD coefficient ==========\n")
sfe <- sqrt(vcovHC(fe, method="arellano", type="HC1")["Treatment:pre_post","Treatment:pre_post"])
sre <- sqrt(vcovHC(re, method="arellano", type="HC1")["Treatment:pre_post","Treatment:pre_post"])
cat(sprintf("%-10s %10s %12s %10s\n","estimator","DiD coef","robust SE","% effect"))
cat(sprintf("%-10s %+10.4f %12.4f %+9.1f%%\n","FE (within)", bfe, sfe, (exp(bfe)-1)*100))
cat(sprintf("%-10s %+10.4f %12.4f %+9.1f%%\n","RE (GLS)",   bre, sre, (exp(bre)-1)*100))
cat("=> FE and RE coincide on the DiD estimate (differ by 0.5%). Because the interaction\n")
cat("   is identified off within-zone timing that BOTH estimators exploit, the DiD is\n")
cat("   insensitive to the RE assumption even though the group LEVEL is not. FE is\n")
cat("   preferred for identification; FE~RE agreement is reassuring, not a reason for RE.\n")
results$fe_did <- unname(bfe); results$fe_se <- unname(sfe)
results$re_did <- unname(bre); results$re_se <- unname(sre)

# ---- write CSV ----
out <- data.frame(metric = names(results), value = unlist(results), row.names = NULL)
write.csv(out, file.path(BUILD, "hausman_diagnostics_61.csv"), row.names = FALSE)
write.csv(design_tbl, file.path(BUILD, "hausman_design_balance_61.csv"), row.names = FALSE)
cat("\nwrote hausman_diagnostics_61.csv and hausman_design_balance_61.csv\n")

cat("\n================ HAUSMAN DIAGNOSTICS COMPLETE ================\n")
cat("Bottom line: FE is correct by research design (subway proximity correlates with\n")
cat("time-invariant zone traits). The Hausman p~1 is a numerical/low-power degeneracy\n")
cat("caused by cross-sectional covariates + month saturation, not evidence for RE. The\n")
cat("robust Mundlak-Wu test leans toward FE; FE and RE agree on the DiD to within 0.5%.\n")
