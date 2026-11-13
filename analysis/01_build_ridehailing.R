# ============================================================================
# Rebuilds ckpt_step1.rds (full_data) from the raw monthly ride-hailing counts,
# correcting a confirmed data-quality defect in the NYC TLC source:
#
#   The published Nov-2017 FHV trip-record file duplicated whole daily batches,
#   inflating every taxi zone by ~1.45x (26.0M rows vs 18.0M trips actually
#   dispatched per the FHV Base Aggregate Report). Oct and Dec match their
#   aggregates exactly, so the raw Nov-2017 counts are inflated.
#
# Fix:
#   - De-duplicate Nov-2017: multiply each zone's Nov-2017 count by a per-zone
#     factor = distinct_trips / total_rows measured from the source trip records
#     (nov2017_dedup_factors.csv). Recovers the true de-duplicated volume.
#   - KEEP Dec-2017: December is CLEAN (19.97M records ~ 19.92M dispatched) and
#     is the normal holiday peak, so it is retained.
#
# Run from the repo root:  Rscript analysis/01_build_ridehailing.R
# ============================================================================
suppressWarnings(suppressMessages({ library(dplyr); library(lubridate) }))

BUILD <- "analysis"
stopifnot(file.exists("data/2016sum.csv"), file.exists("data/2017sum.csv"))

raw <- rbind(read.csv("data/2016sum.csv"), read.csv("data/2017sum.csv"))
raw$month <- ymd(raw$month)
raw <- raw %>% filter(locationID > 0) %>% select(-MoM_change)

# ---- corrected build: de-dupe Nov-2017, keep Dec-2017 ----------------------
fac <- read.csv(file.path(BUILD, "nov2017_dedup_factors.csv"))  # locationID, factor
full_data <- raw %>%
  left_join(fac, by = "locationID") %>%
  mutate(factor = ifelse(is.na(factor), 1, factor),
         observations = ifelse(month == ymd("2017-11-01"),
                               round(observations * factor), observations)) %>%
  select(-factor) %>%
  arrange(locationID, month) %>%
  group_by(locationID) %>%
  mutate(MoM_change = (observations/lag(observations) - 1) * 100) %>% ungroup()
full_data$MoM_change[is.na(full_data$MoM_change)] <- 0
full_data <- full_data %>% filter(observations > 500)   # Dec-2017 retained

saveRDS(list(full_data = full_data), file.path(BUILD, "ckpt_step1.rds"))

cat(sprintf("\nWROTE corrected ckpt_step1.rds: rows=%d, zones=%d, months=%s..%s\n",
            nrow(full_data), length(unique(full_data$locationID)),
            min(full_data$month), max(full_data$month)))
chkm <- full_data %>% group_by(month) %>% summarise(tot = sum(observations)) %>%
  filter(month >= ymd("2017-09-01"))
cat("citywide totals (corrected), Sep-2017 onward:\n"); print(as.data.frame(chkm))
