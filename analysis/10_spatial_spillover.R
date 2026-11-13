# ============================================================================
# Spatial spillover estimation for the Second Avenue Subway DiD (61-zone main).
#
# Motivation:
#   A spatial error model (SEM) fixes spatially-correlated ERRORS, not treatment
#   INTERFERENCE (SUTVA). The baseline TWFE (-13.2%) treats every non-treated
#   zone as a clean control, but control zones bordering the new line may be affected
#   (riders walking in from just outside; drivers repositioning one zone over),
#   contaminating the counterfactual.
#
# Methodology (best-in-class for DiD + spatial spillover):
#   (A) Ring / "donut" spillover-robust DiD. Partition controls into a first-order
#       adjacency RING (queen-contiguous to >=1 treated zone) vs FAR controls.
#       - Donut: drop the ring, identify the direct effect off clean far controls.
#       - Butts (2023, "DiD Estimation with Spatial Spillovers"): keep all, add a
#         Ring x Post term; the Treatment x Post coef is the spillover-robust
#         direct effect (vs far controls) and the Ring x Post coef estimates the
#         spillover onto neighbors. Cf. Clarke (2017).
#   (B) SLX / spatial-lag-of-treatment (Halleck Vega & Elhorst 2015; LeSage &
#       Pace 2009): add W(Treatment x Post), W row-standardised queen contiguity.
#       Direct = own Treatment x Post; indirect/spillover = the WT x Post coef.
#   Design-based ring + reduced-form SLX are preferred over SAR/SDM per Gibbons &
#   Overman (2012, "Mostly Pointless Spatial Econometrics?"); the full SEM
#   specification (see 06_models_full_suite.R) also hit a singularity here,
#   reinforcing the choice.
# ============================================================================

suppressWarnings(suppressMessages({
  library(dplyr); library(spdep); library(plm); library(lmtest); library(sandwich)
}))

BUILD <- "analysis"
census_rides       <- readRDS(file.path(BUILD, "ckpt_census.rds"))$census_rides
land_use_manhattan <- readRDS(file.path(BUILD, "ckpt_landuse.rds"))
ids61              <- readRDS(file.path(BUILD, "zones61_ids.rds"))
w61                <- readRDS(file.path(BUILD, "wts61.rds"))
treat_ids <- c(75, 262, 263, 140, 141, 236, 237)

# ---- rebuild model_data exactly as 06_models_full_suite.R ------------------
land_use_df <- as.data.frame(land_use_manhattan) %>% select(LocationID, pct_res, pct_com, pct_other)
census_rides <- census_rides %>% left_join(land_use_df, by = c("locationID" = "LocationID"))
census_rides$Treatment <- ifelse(census_rides$locationID %in% treat_ids, 1, 0)
model_data <- census_rides %>% select(-geometry)
model_data$pre_post <- ifelse(model_data$month > "2016-12-31" & model_data$month < "2017-12-31", 1, 0)
model_data <- model_data[complete.cases(model_data[, c(
  "populationE","median_incomeE","median_ageE","pct_zero_vehicle","pct_owner_occ",
  "pct_bachelors","observations")]), ]
model_data$log_obs <- log(model_data$observations)
model_data <- model_data %>% filter(locationID %in% ids61)
bal <- model_data %>% count(locationID) %>% filter(n == max(n))
model_data <- model_data %>% filter(locationID %in% bal$locationID) %>% arrange(locationID, month)
zones <- sort(unique(model_data$locationID))
stopifnot(identical(zones, sort(as.integer(w61$ids))))
mean_base <- mean(model_data$observations[model_data$Treatment == 1 & model_data$pre_post == 0])
pct <- function(b) (exp(b) - 1) * 100

clustered <- function(m) {
  V <- vcovHC(m, type = "HC1", cluster = "group")
  ct <- coeftest(m, vcov = V)
  ct
}
row_of <- function(ct, term) ct[rownames(ct) == term, , drop = FALSE]

cat("panel:", length(zones), "zones x", length(unique(model_data$month)),
    "months =", nrow(model_data), "rows | treated baseline mean =",
    round(mean_base), "trips/zone/mo\n\n")

# ---- 0. Baseline TWFE (reproduce -0.142) -----------------------------------
m0 <- plm(log_obs ~ Treatment:pre_post + factor(month),
          data = model_data, index = "locationID", model = "within")
c0 <- clustered(m0); b0 <- coef(m0)["Treatment:pre_post"]
cat("========== 0. BASELINE TWFE (all controls) ==========\n")
cat(sprintf("Treatment x Post: %.4f  (%+.1f%%)  clustered SE %.4f  p=%.4f\n\n",
            b0, pct(b0), row_of(c0,"Treatment:pre_post")[2], row_of(c0,"Treatment:pre_post")[4]))

# ---- Identify the first-order adjacency ring -------------------------------
ids <- w61$ids; nb <- w61$nb
is_treat <- ids %in% treat_ids
treat_pos <- which(is_treat)
is_adj <- sapply(seq_along(ids), function(i) any(nb[[i]] %in% treat_pos)) & !is_treat
adj_ids <- ids[is_adj]; far_ids <- ids[!is_treat & !is_adj]
model_data$Adjacent <- as.integer(model_data$locationID %in% adj_ids)
# explicit interaction columns (robust to formula-expansion quirks in plm)
model_data$treat_post <- model_data$Treatment * model_data$pre_post
model_data$ring_post  <- model_data$Adjacent  * model_data$pre_post
cat("ring partition:", sum(is_treat), "treated |", length(adj_ids),
    "adjacent (ring) |", length(far_ids), "far controls\n")
cat("  ring zone IDs:", paste(adj_ids, collapse = ","), "\n\n")

# ---- A1. Donut: drop the ring, direct effect off clean far controls --------
donut <- model_data %>% filter(Adjacent == 0)   # treated + far controls only
mD <- plm(log_obs ~ Treatment:pre_post + factor(month),
          data = donut, index = "locationID", model = "within")
cD <- clustered(mD); bD <- coef(mD)["Treatment:pre_post"]
cat("========== A1. DONUT DiD (ring dropped) ==========\n")
cat(sprintf("Direct (vs far controls): %.4f  (%+.1f%%)  SE %.4f  p=%.4f\n\n",
            bD, pct(bD), row_of(cD,"Treatment:pre_post")[2], row_of(cD,"Treatment:pre_post")[4]))

# ---- A2. Butts (2023) ring spillover-explicit DiD --------------------------
mB <- plm(log_obs ~ treat_post + ring_post + factor(month),
          data = model_data, index = "locationID", model = "within")
cB <- clustered(mB)
bBd <- coef(mB)["treat_post"]; bBs <- coef(mB)["ring_post"]
cat("========== A2. BUTTS (2023) RING SPILLOVER-EXPLICIT ==========\n")
cat(sprintf("Direct  (Treatment x Post): %.4f  (%+.1f%%)  SE %.4f  p=%.4f\n",
            bBd, pct(bBd), row_of(cB,"treat_post")[2], row_of(cB,"treat_post")[4]))
cat(sprintf("Spill   (Ring x Post):      %.4f  (%+.1f%%)  SE %.4f  p=%.4f\n\n",
            bBs, pct(bBs), row_of(cB,"ring_post")[2], row_of(cB,"ring_post")[4]))

# ---- B. SLX / spatial-lag-of-treatment -------------------------------------
# W row-standardised (nb2mat style="W"); WTreat_i = fraction of i's neighbours treated.
Wt <- as.numeric(w61$mat %*% as.numeric(is_treat)); names(Wt) <- ids
model_data$WTreat <- Wt[as.character(model_data$locationID)]
model_data$WTreatPost <- model_data$WTreat * model_data$pre_post
mS <- plm(log_obs ~ Treatment:pre_post + WTreatPost + factor(month),
          data = model_data, index = "locationID", model = "within")
cS <- clustered(mS)
bSd <- coef(mS)["Treatment:pre_post"]; bSs <- coef(mS)["WTreatPost"]
cat("========== B. SLX SPATIAL-LAG-OF-TREATMENT ==========\n")
cat(sprintf("Direct (Treatment x Post):   %.4f  (%+.1f%%)  SE %.4f  p=%.4f\n",
            bSd, pct(bSd), row_of(cS,"Treatment:pre_post")[2], row_of(cS,"Treatment:pre_post")[4]))
cat(sprintf("Spill  (W.Treat x Post):     %.4f  (%+.1f%% if fully ringed by treatment)  SE %.4f  p=%.4f\n\n",
            bSs, pct(bSs), row_of(cS,"WTreatPost")[2], row_of(cS,"WTreatPost")[4]))

# ---- summary ----------------------------------------------------------------
cat("================= SUMMARY =================\n")
cat(sprintf("%-34s %+7.1f%%\n", "Baseline TWFE (all controls)", pct(b0)))
cat(sprintf("%-34s %+7.1f%%\n", "Donut (ring dropped)", pct(bD)))
cat(sprintf("%-34s %+7.1f%%\n", "Butts direct (vs far controls)", pct(bBd)))
cat(sprintf("%-34s %+7.1f%%\n", "Butts ring spillover", pct(bBs)))
cat(sprintf("%-34s %+7.1f%%\n", "SLX direct", pct(bSd)))
cat(sprintf("%-34s %+7.1f%%\n", "SLX neighbour-treatment spillover", pct(bSs)))
saveRDS(list(b0=b0,bD=bD,bBd=bBd,bBs=bBs,bSd=bSd,bSs=bSs,
             adj_ids=adj_ids, far_ids=far_ids,
             c0=c0,cD=cD,cB=cB,cS=cS),
        file.path(BUILD, "spatial_spillover_out.rds"))
cat("\nsaved spatial_spillover_out.rds\n")
