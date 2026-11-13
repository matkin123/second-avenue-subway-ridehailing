# Build the 61-zone main analytical sample + spatial weights.
#
# The 61 = 66 Manhattan taxi zones with ride data, minus the 5 park/island zones
# with ~0 population (12,43,120,128,194). This is what complete.cases() on the
# covariate pipeline returns, and it is the main analytical sample.
#
# Spatial models need a contiguity structure, so we build queen-contiguity on
# tz_shp restricted to the 61 IDs.
#
# Outputs:
#   analysis/zones61_ids.rds : integer vector of 61 LocationIDs
#   analysis/wts61.rds       : list(nb=poly2nb result, listw=nb2listw, mat=nb2mat)

suppressWarnings(suppressMessages({
  library(dplyr); library(sf); library(spdep)
}))

BUILD <- "analysis"
ck_c  <- readRDS(file.path(BUILD, "ckpt_census.rds"))
cr    <- ck_c$census_rides
tz    <- ck_c$tz_shp

# per-zone covariate NA flag (first row per zone)
cr_df <- cr; cr_df$geometry <- NULL
pz <- cr_df %>% group_by(locationID) %>% slice(1) %>% ungroup() %>%
  mutate(any_na = is.na(populationE) | is.na(median_incomeE) |
                   is.na(med_gross_rentE) | is.na(median_ageE) |
                   is.na(pct_zero_vehicle) | is.na(pct_owner_occ) |
                   is.na(pct_bachelors) | is.na(pop_density))

ids61 <- sort(pz$locationID[!pz$any_na])
cat("All Manhattan zones with ride data:", length(unique(cr_df$locationID)), "\n")
cat("NA-on-covariate (parks/islands) dropped:",
    paste(sort(pz$locationID[pz$any_na]), collapse=","), "\n")
cat("zones61:", length(ids61), "- IDs:\n"); print(ids61)

# sanity: all 7 treatment zones present
treat_ids <- c(75,262,263,140,141,236,237)
stopifnot(all(treat_ids %in% ids61))

saveRDS(ids61, file.path(BUILD, "zones61_ids.rds"))

tz61 <- tz[tz$LocationID %in% ids61, ]
tz61 <- tz61[match(ids61, tz61$LocationID), ]
stopifnot(identical(as.integer(tz61$LocationID), as.integer(ids61)))

nb   <- poly2nb(tz61, queen=TRUE, row.names = as.character(ids61))
listw<- nb2listw(nb, zero.policy=TRUE, style="W")
mat  <- nb2mat(nb, zero.policy=TRUE)

n_iso <- sum(card(nb)==0)
deg   <- card(nb)
cat("\nqueen-contiguity: n=", length(nb),
    " | isolates:", n_iso, " | max degree:", max(deg),
    " | mean degree:", round(mean(deg),2), "\n")
if (n_iso) cat("  isolates (IDs):", paste(ids61[deg==0], collapse=","), "\n")

saveRDS(list(nb=nb, listw=listw, mat=mat, ids=ids61),
        file.path(BUILD, "wts61.rds"))
cat("saved zones61_ids.rds and wts61.rds\n")
