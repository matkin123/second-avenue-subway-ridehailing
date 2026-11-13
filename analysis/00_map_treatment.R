#!/usr/bin/env Rscript
# Treatment-vs-control map of Manhattan taxi zones.
# Run from the repo root:  Rscript analysis/00_map_treatment.R
suppressPackageStartupMessages({ library(sf); library(ggplot2) })

shp <- "analysis/geo/taxi_zones/taxi_zones/taxi_zones.shp"
out <- "analysis/map_treatment.png"

treated <- c(75, 140, 141, 236, 237, 262, 263)

z <- st_read(shp, quiet = TRUE)
mn <- z[z$borough == "Manhattan", ]
mn$grp <- ifelse(mn$LocationID %in% treated, "Treatment", "Control")

col_treat <- "#B4232A"
col_ctrl  <- "#DCDCDE"

p <- ggplot(mn) +
  geom_sf(aes(fill = grp), colour = "white", lwd = 0.18) +
  scale_fill_manual(values = c("Control" = col_ctrl, "Treatment" = col_treat),
                    breaks = c("Treatment", "Control")) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.12, 0.74),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 14, colour = "#333333"),
    legend.key.size = unit(0.42, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(out, p, width = 4.2, height = 5.6, dpi = 220, bg = "white")
cat("map written:", out, "\n")
cat("treated zones found:", sum(mn$grp == "Treatment"), "of 7\n")
