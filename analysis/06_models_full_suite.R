# Full model suite on the covariate panel — MAIN SPEC = 61 Manhattan zones.
# Builds model_data from census_rides + land_use_manhattan, aligns the sample to
# the 61-zone residential/commercial set (66 Manhattan zones with ride data minus
# the 5 parks/islands), then runs OLS / two-way FE / pFtest / Hausman / event
# study / clustered SEs / comparison table / three SEM spatial models + robustness.
#
# The 52-zone contiguous subset is kept as a single robustness row.
#
# NOTE: census/land-use covariate VALUES come from interpolation (tidycensus
# population-weighted interp + TIGER-2010 areal interp), so OLS & SEM coefficients
# are approximate. The FE DiD is immune (within estimator absorbs time-invariant
# covariates).

suppressWarnings(suppressMessages({
  library(dplyr); library(spdep); library(plm); library(lmtest)
  library(sandwich); library(splm); library(broom); library(ggplot2)
}))

BUILD <- "analysis"
census_rides        <- readRDS(file.path(BUILD, "ckpt_census.rds"))$census_rides
land_use_manhattan  <- readRDS(file.path(BUILD, "ckpt_landuse.rds"))
ids61               <- readRDS(file.path(BUILD, "zones61_ids.rds"))
ids52               <- readRDS(file.path(BUILD, "zones52_ids.rds"))
w61                 <- readRDS(file.path(BUILD, "wts61.rds"))

`%||%` <- function(a,b) if (is.null(a)) b else a
modelsummary <- function(models, vcov=NULL, coef_map=NULL, ...) {
  terms <- if (is.null(coef_map)) "Treatment:pre_post" else names(coef_map)
  cols <- lapply(seq_along(models), function(i){
    m <- models[[i]]; V <- if(!is.null(vcov)) vcov[[i]] else NULL
    cf <- coef(m)
    se <- tryCatch(if(is.null(V)) sqrt(diag(as.matrix(stats::vcov(m)))) else sqrt(diag(V)),
                   error=function(e) setNames(rep(NA_real_, length(cf)), names(cf)))
    sapply(terms, function(tn) if(tn %in% names(cf)) sprintf("%.4f (%.4f)", cf[tn], se[tn]) else "")
  })
  tab <- as.data.frame(cols, stringsAsFactors=FALSE, optional=TRUE)
  names(tab) <- names(models)
  rownames(tab) <- if(is.null(coef_map)) terms else unname(coef_map[terms])
  nobs_row <- sapply(models, function(m) tryCatch(stats::nobs(m), error=function(e) NA))
  r2_row <- sapply(models, function(m) tryCatch({s<-summary(m)
      if(!is.null(s$r.squared)) round(unname(s$r.squared[1]),3) else NA}, error=function(e) NA))
  tab <- rbind(tab, "Num.Obs."=as.character(nobs_row), "R2"=as.character(r2_row))
  cat("\n---- model comparison (coef (SE)) ----\n"); print(tab); cat("\n"); invisible(tab)
}

# ---- Step 5: build model_data ----
land_use_df <- as.data.frame(land_use_manhattan) %>% select(LocationID, pct_res, pct_com, pct_other)
census_rides <- census_rides %>% left_join(land_use_df, by=c("locationID"="LocationID"))
census_rides$Treatment <- ifelse(census_rides$locationID %in% c(75,262,263,140,141,236,237), 1, 0)

model_data <- census_rides %>% select(-geometry)
model_data$pre_post <- ifelse(model_data$month > '2016-12-31' & model_data$month < "2017-12-31", 1, 0)
model_data$populationE    <- model_data$populationE/1000
model_data$median_incomeE <- model_data$median_incomeE/1000
for (v in c("pct_zero_vehicle","pct_owner_occ","pct_bachelors","white_pct","black_pct",
            "latino_pct","native_pct","asian_pct","pct_res","pct_com","pct_other"))
  model_data[[v]] <- 100*model_data[[v]]
model_data <- model_data[complete.cases(model_data), ]
cc_zones <- length(unique(model_data$locationID))
model_data$log_obs <- log(model_data$observations)

cat("complete-cases zones:", cc_zones,
    "| overlap with 61-main:", length(intersect(unique(model_data$locationID), ids61)),
    "| overlap with 52-robust:", length(intersect(unique(model_data$locationID), ids52)), "\n")

# MAIN SPEC: 61-zone balanced panel
model_data <- model_data %>% filter(locationID %in% ids61)
bal <- model_data %>% count(locationID) %>% filter(n==max(n))
model_data <- model_data %>% filter(locationID %in% bal$locationID) %>%
  arrange(locationID, month)
zones <- sort(unique(model_data$locationID))
cat("analysis panel (MAIN, 61):", length(zones), "zones x",
    length(unique(model_data$month)), "months =", nrow(model_data), "rows\n")
mean_baseline_treatment <- mean(model_data$observations[model_data$Treatment==1 & model_data$pre_post==0])

# spatial weights aligned to the 61 zones (queen contiguity, built from tz_shp)
stopifnot(identical(zones, sort(as.integer(w61$ids))))
tz_listw <- w61$listw
tz_mat   <- w61$mat

cat("\n========== SECTION 1: OLS (corrected, + Treatment main effect) ==========\n")
OLS_f3 <- "log_obs ~ populationE + median_incomeE + median_ageE + pct_zero_vehicle + black_pct + latino_pct + asian_pct + pct_com + pct_other + Treatment + Treatment:pre_post"
OLS_model3_log <- lm(OLS_f3, data=model_data, na.action=na.exclude)
cat("OLS DiD (Treatment:pre_post):", round(coef(OLS_model3_log)["Treatment:pre_post"],4),
    "| Treatment main:", round(coef(OLS_model3_log)["Treatment"],4),
    "| R2:", round(summary(OLS_model3_log)$r.squared,3), "\n")

cat("\n========== SECTION 2: FIXED EFFECTS (two-way within) ==========\n")
WE_f <- "log_obs ~ populationE + median_incomeE + median_ageE + pct_bachelors + pct_zero_vehicle + pct_owner_occ + black_pct + pct_res + pct_other + Treatment:pre_post + factor(month)"
WE_model_log <- plm(WE_f, data=model_data, index="locationID", model="within")
b_did <- coef(WE_model_log)["Treatment:pre_post"]
cat("FE DiD (log-pts):", round(b_did,4), "| % change:", round((exp(b_did)-1)*100,1),
    "% | trips/zone/mo:", round((exp(b_did)-1)*mean_baseline_treatment), "\n")

OLS_for_ftest <- plm(WE_f, data=model_data, index="locationID", model="pooling")
cat("pFtest (FE vs pooled) p =", format.pval(pFtest(WE_model_log, OLS_for_ftest)$p.value), "\n")
RE_model_log <- plm(WE_f, data=model_data, index="locationID", model="random")
# p~1 is expected, not a bug: covariates are time-invariant so FE drops them and
# FE~RE coincide => V_FE-V_RE is indefinite. See 11_hausman.R for the full explanation.
cat("Hausman (FE vs RE) p =", format.pval(phtest(WE_model_log, RE_model_log)$p.value), "\n")

cat("\n========== SECTION 3: SPATIAL DEPENDENCE TESTS ==========\n")
model_data$residuals_log <- residuals(WE_model_log)
cat("CD test on log_obs    z =", round(pcdtest(pdata.frame(model_data)$log_obs, index="month", w=tz_mat)$statistic,3), "\n")
cat("CD test on residuals  z =", round(pcdtest(pdata.frame(model_data)$residuals_log, index="month", w=tz_mat)$statistic,3), "\n")

cat("\n========== SECTION 4: EVENT STUDY ==========\n")
model_data$month_factor <- relevel(factor(model_data$month), ref="2016-12-01")
es <- plm("log_obs ~ Treatment:month_factor", data=model_data,
          index=c("locationID","month"), model="within", effect="twoways")
ec <- tidy(es, conf.int=TRUE); ec <- ec[grepl("Treatment:month_factor", ec$term),]
ec$month <- as.Date(gsub("Treatment:month_factor","",ec$term))
pre <- ec$estimate[ec$month < as.Date("2017-01-01")]
post<- ec$estimate[ec$month >= as.Date("2017-01-01")]
cat("pre-period mean coef:", round(mean(pre),4), "| post-period mean coef:", round(mean(post),4),
    "| post negative:", sum(post<0), "of", length(post), "\n")

cat("\n========== SECTION 5: CLUSTER-ROBUST SEs ==========\n")
Vc <- vcovHC(WE_model_log, type="HC1", cluster="group")
se_o <- sqrt(vcov(WE_model_log)["Treatment:pre_post","Treatment:pre_post"])
se_c <- sqrt(Vc["Treatment:pre_post","Treatment:pre_post"])
cat("DiD:", round(b_did,4), "| original SE:", round(se_o,4),
    "(p", format.pval(2*pt(-abs(b_did/se_o),Inf)), ") | clustered SE:", round(se_c,4),
    "(p", format.pval(2*pt(-abs(b_did/se_c),Inf)), ")\n")

cat("\n========== SECTION 6: MODEL COMPARISON TABLE ==========\n")
WE_loc_only <- plm(WE_f, data=model_data, index="locationID", model="within", effect="individual")
modelsummary(
  list("OLS"=OLS_model3_log, "FE (loc only)"=WE_loc_only,
       "FE (two-way)"=WE_model_log, "FE (cluster SE)"=WE_model_log),
  vcov=list(NULL,NULL,NULL,Vc),
  coef_map=c("Treatment:pre_post"="DiD: Treatment x Post","Treatment"="Treatment (main)"))

cat("\n========== SECTION 7: SPATIAL ERROR MODELS (splm::spreml) ==========\n")
SEM_data <- pdata.frame(model_data, index=c("locationID"))
run_sem <- function(lbl, f) {
  r <- tryCatch({
    m <- spreml(f, data=SEM_data, w=tz_listw, lag=FALSE, errors="sem")
    co <- m$coefficients
    did <- co[grep("Treatment", names(co))]
    rho <- m$errcomp["rho"] %||% NA
    cat(sprintf("%-26s DiD=%s | rho=%s\n", lbl,
                paste(round(did,4), collapse=","),
                if(is.na(rho)) "NA" else round(rho,4)))
  }, error=function(e) cat(sprintf("%-26s FAILED: %s\n", lbl, conditionMessage(e))))
  invisible(r)
}
run_sem("SEM1 loc+time FE + DiD", "log_obs ~ locationID + factor(month) + Treatment:pre_post")
run_sem("SEM2 loc FE + DiD",      "log_obs ~ locationID + Treatment:pre_post")
run_sem("SEM3 time FE + covars",  "log_obs ~ factor(month) + populationE + median_incomeE + median_ageE + pct_bachelors + pct_zero_vehicle + pct_owner_occ + black_pct + pct_com + pct_other + Treatment:pre_post")

cat("\n========== SECTION 8: ROBUSTNESS ==========\n")
nine <- model_data %>% filter(pre_post==0 | (pre_post==1 & month <= as.Date("2017-09-01")))
b9 <- coef(plm(WE_f, data=nine, index="locationID", model="within"))["Treatment:pre_post"]
alt_ids <- c(238,239,142,143,151,24,166)
alt <- model_data %>% filter(Treatment==1 | locationID %in% alt_ids)
balt <- coef(plm(WE_f, data=alt, index="locationID", model="within"))["Treatment:pre_post"]
rc <- model_data %>% filter(!(locationID %in% c(153,88,186)))
brc <- coef(plm(WE_f, data=rc, index="locationID", model="within"))["Treatment:pre_post"]
md52 <- model_data %>% filter(locationID %in% ids52)
b52  <- coef(plm(WE_f, data=md52, index="locationID", model="within"))["Treatment:pre_post"]
cat(sprintf("%-32s %+8.4f %+7.1f%%\n", "Main (61-zone)", b_did, (exp(b_did)-1)*100))
cat(sprintf("%-32s %+8.4f %+7.1f%%\n", "Residential core (58-zone)", brc, (exp(brc)-1)*100))
cat(sprintf("%-32s %+8.4f %+7.1f%%\n", "9-month window", b9,   (exp(b9)-1)*100))
cat(sprintf("%-32s %+8.4f %+7.1f%% (%d ctrl)\n", "UWS alt-control", balt, (exp(balt)-1)*100,
            length(unique(alt$locationID[alt$Treatment==0]))))
cat(sprintf("%-32s %+8.4f %+7.1f%%\n", "52-zone subset", b52, (exp(b52)-1)*100))

cat("\n================ FULL CENSUS/LAND-USE PIPELINE COMPLETE ================\n")
