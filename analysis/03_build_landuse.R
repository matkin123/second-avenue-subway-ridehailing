# Build PLUTO land-use -> census tract -> taxi zone (areal interpolation).
# Uses TIGER 2010 tracts (NAME10 as ctlabel) for New York County, which join 1:1
# to PLUTO CT2010. tz_shp comes from the census checkpoint (same Manhattan
# taxi-zone geometry, EPSG 4269).

suppressWarnings(suppressMessages({
  library(tigris); library(sf); library(dplyr)
}))
options(tigris_use_cache = TRUE)

BUILD <- "analysis"
tz_shp <- readRDS(file.path(BUILD, "ckpt_census.rds"))$tz_shp

pluto <- read.csv("data/MN.csv")
pluto <- pluto %>% select(CT2010, LandUse, LotArea, BldgArea, ComArea, ResArea) %>%
  mutate(largest_area = ifelse(BldgArea > LotArea, BldgArea, LotArea))
pluto_summary <- pluto %>% group_by(CT2010) %>%
  summarise(total_area=sum(largest_area), total_res=sum(ResArea), total_com=sum(ComArea),
            .groups="drop")
pluto_summary$pct_res   <- pluto_summary$total_res / pluto_summary$total_area
pluto_summary$pct_com   <- pluto_summary$total_com / pluto_summary$total_area
pluto_summary$pct_other <- (pluto_summary$total_area -
  (pluto_summary$total_com + pluto_summary$total_res)) / pluto_summary$total_area

cat("fetching TIGER 2010 census tracts (New York County)...\n")
ct <- tracts(state="NY", county="New York", year=2010, progress_bar=FALSE)
ct <- st_transform(ct, st_crs(tz_shp))
ct_as_df <- as.data.frame(ct) %>% transmute(ctlabel = as.numeric(NAME10), geometry)
ct_as_df <- ct_as_df[!duplicated(ct_as_df$ctlabel), ]

land_use_w_geom <- st_as_sf(left_join(pluto_summary, ct_as_df, by=c("CT2010"="ctlabel")))
land_use_w_geom <- land_use_w_geom[!st_is_empty(land_use_w_geom), ]
land_use_w_geom <- land_use_w_geom %>% select(pct_res, pct_com, pct_other, geometry)

cat("areal interpolation CT -> TZ...\n")
suppressWarnings({
  land_use_interpolate <- st_interpolate_aw(land_use_w_geom, tz_shp, extensive=FALSE)
  land_use_manhattan   <- st_join(tz_shp, land_use_interpolate, join=st_equals, left=TRUE)
})

ldf <- as.data.frame(land_use_manhattan)
cat("land_use_manhattan zones:", nrow(ldf),
    "| non-NA pct_res:", sum(!is.na(ldf$pct_res)), "\n")
cat("  pct_res mean:", round(mean(ldf$pct_res, na.rm=TRUE),3),
    "| pct_com mean:", round(mean(ldf$pct_com, na.rm=TRUE),3),
    "| pct_other mean:", round(mean(ldf$pct_other, na.rm=TRUE),3), "\n")
saveRDS(land_use_manhattan, file.path(BUILD, "ckpt_landuse.rds"))
cat("saved ckpt_landuse.rds\n")
