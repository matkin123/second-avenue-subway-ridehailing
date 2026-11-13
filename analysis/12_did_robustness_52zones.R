# Robustness DiD on a 52-zone contiguous Manhattan subset.
# The 52 LocationIDs (zones52_ids.rds) form a self-contained queen-contiguous
# block (no out-of-set neighbors); all 52 are Manhattan zones and all 7 treatment
# zones are present. Re-estimating the FE / event-study / Hausman / robustness
# numbers on this tighter, fully-contiguous subsample is a spatial robustness
# check on the 61-zone headline. No census data is needed here (the within
# estimator absorbs the time-invariant covariates anyway).

suppressWarnings(suppressMessages({
  library(dplyr); library(plm); library(lmtest)
  library(sandwich); library(broom); library(ggplot2)
}))

BUILD <- "analysis"
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data
ids52     <- readRDS(file.path(BUILD, "zones52_ids.rds"))
treat_ids <- c(75,262,263,140,141,236,237)

md <- full_data %>%
  filter(locationID %in% ids52) %>%
  mutate(Treatment = as.integer(locationID %in% treat_ids),
         pre_post  = as.integer(month > as.Date("2016-12-31") & month < as.Date("2017-12-31")),
         log_obs   = log(observations))

# keep the full 23-month balanced panel for clean comparison
balance <- md %>% count(locationID) %>% filter(n == max(n))
md <- md %>% filter(locationID %in% balance$locationID) %>% arrange(locationID, month)

cat("zeros in observations:", sum(md$observations==0), "\n")
cat("zones in 52-panel:", length(unique(md$locationID)),
    "| treated:", length(unique(md$locationID[md$Treatment==1])),
    "| months:", length(unique(md$month)), "| rows:", nrow(md), "\n")
dropped <- setdiff(ids52, unique(md$locationID))
if (length(dropped)) cat("zones dropped by balance filter:", paste(dropped, collapse=","), "\n")
mean_baseline <- mean(md$observations[md$Treatment==1 & md$pre_post==0])
cat("mean baseline trips (treatment zones, pre-2017):", round(mean_baseline), "\n\n")

f <- "%s ~ Treatment:pre_post + factor(month)"

cat("========== raw-outcome two-way FE (52 zones) ==========\n")
raw_fe <- plm(sprintf(f,"observations"), data=md, index=c("locationID","month"),
              model="within", effect="twoways")
b_raw <- coef(raw_fe)["Treatment:pre_post"]
ci_raw <- confint(raw_fe)["Treatment:pre_post",]
cat("raw FE DiD coef (trips/zone/month):", round(b_raw), "\n")
cat("  95% CI:", round(ci_raw[1]), "to", round(ci_raw[2]), "\n\n")

cat("========== log-outcome two-way FE (52 zones) ==========\n")
log_fe <- plm(sprintf(f,"log_obs"), data=md, index=c("locationID","month"),
              model="within", effect="twoways")
b_log <- coef(log_fe)["Treatment:pre_post"]
cat("log FE DiD coef (log-points):", round(b_log,4), "\n")
cat("  % change:", round((exp(b_log)-1)*100,1), "%\n")
cat("  abs trips/zone/month:", round((exp(b_log)-1)*mean_baseline), "\n")
cat("  annual across 7 zones x 11 months:", round((exp(b_log)-1)*mean_baseline*7*11), "\n\n")

cat("========== pFtest (FE vs pooled) & Hausman (FE vs RE) ==========\n")
pool   <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="pooling")
we_loc <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="within")
cat("pFtest p =", format.pval(pFtest(we_loc, pool)$p.value), "\n")
re <- plm(sprintf(f,"log_obs"), data=md, index="locationID", model="random")
cat("Hausman p =", format.pval(phtest(we_loc, re)$p.value), "\n\n")

cat("========== Event study (parallel-trends test, 52 zones) ==========\n")
md$month_factor <- relevel(factor(md$month), ref="2016-12-01")
es <- plm("log_obs ~ Treatment:month_factor", data=md,
          index=c("locationID","month"), model="within", effect="twoways")
ec <- tidy(es, conf.int=TRUE)
ec <- ec[grepl("Treatment:month_factor", ec$term),]
ec$month <- as.Date(gsub("Treatment:month_factor","",ec$term))
pre  <- ec$estimate[ec$month < as.Date("2017-01-01")]
post <- ec$estimate[ec$month >= as.Date("2017-01-01")]
pre_span0 <- all(ec$conf.low[ec$month<as.Date("2017-01-01")] < 0 &
                 ec$conf.high[ec$month<as.Date("2017-01-01")] > 0)
cat("pre-period coefs: mean", round(mean(pre),4), "| max|coef|", round(max(abs(pre)),4),
    "| all pre CIs span 0:", pre_span0, "\n")
cat("post-period coefs: mean", round(mean(post),4), "| n negative:",
    sum(post<0), "of", length(post), "\n")
ref <- data.frame(term="ref", estimate=0, std.error=0, statistic=NA, p.value=NA,
                  conf.low=0, conf.high=0, month=as.Date("2016-12-01"))
ecp <- rbind(ec[,names(ref)], ref); ecp <- ecp[order(ecp$month),]
p <- ggplot(ecp, aes(month, estimate)) +
  geom_point(size=2.3) + geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=10) +
  geom_hline(yintercept=0, linetype="dashed", color="gray40") +
  geom_vline(xintercept=as.Date("2017-01-01"), color="red", linetype="dashed") +
  scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
  labs(title="Event Study (52-zone contiguous subset): Second Avenue Subway effect",
       subtitle="Two-way FE, ref = Dec 2016. Pre-period near zero = parallel trends.",
       x=NULL, y="log-points vs Dec 2016") + theme_minimal()
ggsave(file.path(BUILD,"event_study_real52.png"), p, width=9, height=5, dpi=120)
cat("event study plot saved: event_study_real52.png\n\n")

cat("========== Cluster-robust SEs (zone level, 52 zones) ==========\n")
V    <- vcovHC(we_loc, type="HC1", cluster="group")
se_o <- sqrt(vcov(we_loc)["Treatment:pre_post","Treatment:pre_post"])
se_c <- sqrt(V["Treatment:pre_post","Treatment:pre_post"])
bL   <- coef(we_loc)["Treatment:pre_post"]
cat("DiD coef:", round(bL,4),
    "| original SE:", round(se_o,4), "(p", format.pval(2*pt(-abs(bL/se_o), df=Inf)), ")",
    "| clustered SE:", round(se_c,4), "(p", format.pval(2*pt(-abs(bL/se_c), df=Inf)), ")\n\n")

cat("========== Robustness: 9-month window & UWS alt-control ==========\n")
nine <- md %>% filter(pre_post==0 | (pre_post==1 & month <= as.Date("2017-09-01")))
b9 <- coef(plm(sprintf(f,"log_obs"), data=nine, index="locationID", model="within"))["Treatment:pre_post"]
alt_ids <- c(238,239,142,143,151,24,166)
alt <- md %>% filter(Treatment==1 | locationID %in% alt_ids)
balt <- coef(plm(sprintf(f,"log_obs"), data=alt, index="locationID", model="within"))["Treatment:pre_post"]
cat(sprintf("%-26s %10s %10s %14s\n","spec","log-pts","% change","trips/zone/mo"))
for (nm in c("Main (52-zone MH controls)","9-month post window","UWS alt-control")) {
  b <- switch(nm, "Main (52-zone MH controls)"=b_log, "9-month post window"=b9, "UWS alt-control"=balt)
  cat(sprintf("%-26s %+10.4f %+9.1f%% %14.0f\n", nm, b, (exp(b)-1)*100, (exp(b)-1)*mean_baseline))
}
cat("\nUWS alt-control uses", length(unique(alt$locationID[alt$Treatment==0])), "control zones present in the 52-set\n")
cat("\n================ 52-ZONE SUBSET ROBUSTNESS COMPLETE ================\n")
