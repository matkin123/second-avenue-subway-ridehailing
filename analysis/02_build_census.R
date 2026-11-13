# Build census demographics -> taxi zones (population-weighted interpolation).
# Uses the TLC taxi_zones shapefile (with LocationID attributes) in analysis/geo
# and the CENSUS_API_KEY from the environment. tidycensus pulls tract/block-group
# geometry automatically.

suppressWarnings(suppressMessages({
  library(tidycensus); library(sf); library(units); library(dplyr)
}))
options(tigris_use_cache = TRUE)
stopifnot(nzchar(Sys.getenv("CENSUS_API_KEY")))

BUILD <- "analysis"
GEO   <- file.path(BUILD, "geo")
full_data <- readRDS(file.path(BUILD, "ckpt_step1.rds"))$full_data

variables_to_get <- c(
  population="B01003_001", median_income="B19013_001",
  zero_vehicle_count="B08201_002", total_hh_vehicle="B08201_001",
  owner_occ_hh="B08137_002", total_hh_occ="B08137_001",
  med_gross_rent="B25031_001", median_age="B01002_001",
  bachelors="B15003_022", white="B02001_002", black="B02001_003",
  native="B02001_004", asian="B02001_005", pi="B02001_006",
  latino="B03003_003")

cat("fetching tract ACS 2016 (wide)...\n")
nyc_acs <- get_acs(geography="tract", variables=variables_to_get, state="NY",
                   county="New York", survey="acs5",
                   key=Sys.getenv("CENSUS_API_KEY"),
                   geometry=TRUE, output="wide", year=2016, progress_bar=FALSE)
nyc_acs <- nyc_acs %>% select(!ends_with("M"))

cat("reading canonical taxi_zones (Manhattan)...\n")
tz_shp <- st_read(file.path(GEO, "taxi_zones/taxi_zones"), quiet=TRUE)
tz_shp <- tz_shp[tz_shp$borough == "Manhattan", ]
tz_shp <- st_transform(tz_shp, 4269)

cat("fetching block-group population 2016 (interp weights)...\n")
pop_2016_acs <- get_acs(geography="block group", variables=c(population="B01003_001"),
                        state="NY", county="New York", survey="acs5",
                        key=Sys.getenv("CENSUS_API_KEY"),
                        geometry=TRUE, output="wide", year=2016, progress_bar=FALSE)
pop_2016_acs <- pop_2016_acs[!st_is_empty(pop_2016_acs), , drop=FALSE]
pop_2016_acs <- pop_2016_acs %>% select(GEOID, populationE, geometry) %>%
  rename(pop16 = populationE)

cat("population-weighted interpolation (sums + averages)...\n")
sums_df <- nyc_acs %>% select(-c(median_incomeE, med_gross_rentE, median_ageE))
interp_sums <- interpolate_pw(sums_df, tz_shp, to_id="LocationID",
                              extensive=TRUE, weights=pop_2016_acs, weight_column="pop16")
averages_df <- nyc_acs %>% select(median_incomeE, med_gross_rentE, median_ageE)
interp_avgs <- interpolate_pw(averages_df, tz_shp, to_id="LocationID",
                              extensive=FALSE, weights=pop_2016_acs, weight_column="pop16")

interp_avgs_df <- as.data.frame(interp_avgs) %>% select(-geometry)
interp_sums_df <- as.data.frame(interp_sums)
pop_interp <- merge(interp_sums_df, interp_avgs_df, by="LocationID")
pop_interp <- st_as_sf(pop_interp)

pop_interp <- pop_interp %>% mutate(
  pct_zero_vehicle = zero_vehicle_countE/total_hh_vehicleE,
  pct_owner_occ    = owner_occ_hhE/total_hh_occE,
  pct_bachelors    = bachelorsE/populationE,
  pop_density      = as.numeric(set_units(populationE/st_area(.), "1/km2")),
  white_pct=whiteE/populationE, black_pct=blackE/populationE,
  native_pct=nativeE/populationE, asian_pct=asianE/populationE,
  latino_pct=latinoE/populationE)

census_rides <- merge(full_data, pop_interp, by.x="locationID", by.y="LocationID")
census_rides <- census_rides %>% select(-c(zero_vehicle_countE, total_hh_vehicleE,
  owner_occ_hhE, total_hh_occE, bachelorsE, whiteE, blackE, nativeE, asianE, piE, latinoE))

cat("census_rides rows:", nrow(census_rides),
    "| Manhattan zones:", length(unique(census_rides$locationID)), "\n")
saveRDS(list(census_rides=census_rides, tz_shp=tz_shp), file.path(BUILD, "ckpt_census.rds"))
cat("saved ckpt_census.rds\n")
