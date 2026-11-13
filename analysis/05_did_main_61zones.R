# Main core DiD on the 61 Manhattan zones (headline specification).
# Uses the full 61-zone residential/commercial set (66 Manhattan zones with
# ride data minus 5 park/island zones with ~0 pop). The 52-zone contiguous
# subset is run as a robustness row in 12_did_robustness_52zones.R.

suppressWarnings(suppressMessages({
  library(dplyr); library(plm); library(lmtest)
  library(sandwich); library(broom); library(ggplot2)
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

balance <- md %>% count(locationID) %>% filter(n == max(n))
md <- md %>% filter(locationID %in% balance$locationID) %>% arrange(locationID, month)

cat("zeros in observations:", sum(md$observations==0), "\n")
cat("zones in 61-panel:", length(unique(md$locationID)),
    "| treated:", length(unique(md$locationID[md$Treatment==1])),
    "| months:", length(unique(md$month)), "| rows:", nrow(md), "\n")
dropped <- setdiff(ids61, unique(md$locationID))
if (length(dropped)) cat("zones dropped by balance filter:", paste(dropped, collapse=","), "\n")
mean_baseline <- mean(md$observations[md$Treatment==1 & md$pre_post==0])
cat("mean baseline trips (treatment zones, pre-2017):", round(mean_baseline), "\n\n")

f <- "%s ~ Treatment:pre_post + factor(month)"

cat("========== RAW two-way FE (61 zones) ==========\n")
raw_fe <- plm(sprintf(f,"observations"), data=md, index=c("locationID","month"),
              model="within", effect="twoways")
b_raw <- coef(raw_fe)["Treatment:pre_post"]
ci_raw <- confint(raw_fe)["Treatment:pre_post",]
cat("raw FE DiD coef (trips/zone/month):", round(b_raw),
    " 95% CI [", round(ci_raw[1]), ",", round(ci_raw[2]), "]\n\n")

cat("========== MAIN HEADLINE: log-outcome two-way FE (61 zones) ==========\n")
log_fe <- plm(sprintf(f,"log_obs"), data=md, index=c("locationID","month"),
              model="within", effect="twoways")
b_log <- coef(log_fe)["Treatment:pre_post"]
cat("log FE DiD coef (log-points):", round(b_log,4),
    "| % change:", round((exp(b_log)-1)*100,1),
    "% | trips/zone/mo:", round((exp(b_log)-1)*mean_baseline), "\n\n")

cat("========== pFtest & Hausman ==========\n")
pool   <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="pooling")
we_loc <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="within")
cat("pFtest p =", format.pval(pFtest(we_loc, pool)$p.value), "\n")
re <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="random")
# p~1 is expected (degenerate): see 11_hausman.R for the full explanation.
cat("Hausman p =", format.pval(phtest(we_loc, re)$p.value), "\n\n")

cat("========== Cluster-robust SEs (zone-level, 61 zones) ==========\n")
V    <- vcovHC(we_loc, type="HC1", cluster="group")
se_o <- sqrt(vcov(we_loc)["Treatment:pre_post","Treatment:pre_post"])
se_c <- sqrt(V["Treatment:pre_post","Treatment:pre_post"])
bL   <- coef(we_loc)["Treatment:pre_post"]
cat("DiD coef:", round(bL,4),
    "| original SE:", round(se_o,4), "(p", format.pval(2*pt(-abs(bL/se_o), df=Inf)), ")",
    "| clustered SE:", round(se_c,4), "(p", format.pval(2*pt(-abs(bL/se_c), df=Inf)), ")\n\n")

cat("========== Robustness summary ==========\n")
nine <- md %>% filter(pre_post==0 | (pre_post==1 & month <= as.Date("2017-09-01")))
b9 <- coef(plm(sprintf(f,"log_obs"), data=nine, index="locationID", model="within"))["Treatment:pre_post"]

alt_ids <- c(238,239,142,143,151,24,166)
alt <- md %>% filter(Treatment==1 | locationID %in% alt_ids)
balt <- coef(plm(sprintf(f,"log_obs"), data=alt, index="locationID", model="within"))["Treatment:pre_post"]

# residential-core subset: drop atypical zones 153 (Marble Hill, geographically Bronx),
# 88 (FiDi-S, commercial/tourist), 186 (Penn Sq, commercial/tourist)
rc <- md %>% filter(!(locationID %in% c(153,88,186)))
brc <- coef(plm(sprintf(f,"log_obs"), data=rc, index="locationID", model="within"))["Treatment:pre_post"]

ids52 <- readRDS(file.path(BUILD, "zones52_ids.rds"))
md52  <- md %>% filter(locationID %in% ids52)
b52   <- coef(plm(sprintf(f,"log_obs"), data=md52, index="locationID", model="within"))["Treatment:pre_post"]

cat(sprintf("%-32s %10s %10s %14s\n","spec","log-pts","% change","trips/zone/mo"))
for (nm in c("Main (61-zone)","Residential core (58-zone)","9-month post window","UWS alt-control","52-zone subset")) {
  b <- switch(nm,
              "Main (61-zone)"=b_log,
              "Residential core (58-zone)"=brc,
              "9-month post window"=b9,
              "UWS alt-control"=balt,
              "52-zone subset"=b52)
  cat(sprintf("%-32s %+10.4f %+9.1f%% %14.0f\n", nm, b, (exp(b)-1)*100, (exp(b)-1)*mean_baseline))
}

cat("\n================ 61-ZONE MAIN-SPEC VALIDATION COMPLETE ================\n")
