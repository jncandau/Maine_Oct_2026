# =============================================================================
# 02_describe_trap_data.R
# -----------------------------------------------------------------------------
# Describe the Maine light-trap series produced by R/01_read_trap_data.R:
# sampling effort per year and the geography of the trap network.
#
# INPUT   data/processed/trap_daily.rds   (one row per sampled trap-night)
# OUTPUT  output/tables/traps_per_year.csv    effort summary, one row per year
#         output/figures/trap_sites_map.png   map of the trap localities
#
# NOTE    trap_daily.rds holds only the records with a daily series (535 of the
#         643 Maine light-trap records); the 108 season-total-only records are
#         absent by construction. Counts below therefore describe the traps
#         usable for phenology, not every trap the Maine Forest Service ran.
#
# REQUIREMENTS  dplyr, ggplot2, maps, ggrepel
#   renv::install(c("maps", "ggrepel")); renv::snapshot()
#
# USAGE   source("R/02_describe_trap_data.R")   # from the project root
# =============================================================================

library(dplyr)
library(ggplot2)

dir.create(file.path("output", "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("output", "figures"), showWarnings = FALSE, recursive = TRUE)

trap_daily <- readRDS(file.path("data", "processed", "trap_daily.rds"))


## ---- 1. Traps per year ------------------------------------------------------
## A "trap" is one trapping record (flightID). `n_sites` differs from
## `n_records` only where two records share a locality in the same year.

traps_per_year <- trap_daily %>%
  group_by(year) %>%
  summarise(
    n_traps         = n_distinct(flightID),
    n_sites         = n_distinct(locality),
    n_nights        = sum(status == "counted"),
    nights_per_trap = round(n_nights / n_distinct(flightID), 1),
    first_doy       = min(doy),
    last_doy        = max(doy),
    moths_total     = sum(moths, na.rm = TRUE),
    .groups         = "drop"
  )

write.csv(traps_per_year, file.path("output", "tables", "traps_per_year.csv"),
          row.names = FALSE)


## ---- 2. Trap localities -----------------------------------------------------
## Coordinates are constant within a locality in this dataset (checked), so the
## first pair is representative.

sites <- trap_daily %>%
  group_by(locality) %>%
  summarise(
    lat      = first(latitude),
    lon      = first(longitude),
    n_years  = n_distinct(year),
    n_nights = sum(status == "counted"),
    .groups  = "drop"
  )

## Sites worth naming on the map: the longest-running series and the two
## latitudinal extremes. Labelling all 75 would be unreadable.
label_sites <- bind_rows(
  slice_max(sites, n_years, n = 1),
  slice_max(sites, lat,     n = 1),
  slice_min(sites, lat,     n = 1)
) %>%
  distinct(locality, .keep_all = TRUE) %>%
  ## Short map labels: the workbook localities carry ", Maine", township codes
  ## and alternative place names, all of which are too long to sit on a map.
  mutate(label = locality %>%
           sub(",? *Maine$", "", .) %>%
           sub(" Township.*$", "", .) %>%
           sub(",? *T[0-9]+ R[0-9]+.*$", "", .) %>%
           sub("^Dennistown / Jackman.*$", "Jackman", .))


## ---- 3. Map -----------------------------------------------------------------

## Base layers: US states for Maine and its neighbours, plus the Canadian
## border to the north and east (the trap network follows that border closely).
## ggplot2::map_data() reads the boundary databases of the `maps` package.
states <- map_data("state", region = c("maine", "new hampshire", "vermont",
                                       "massachusetts"))
canada <- map_data("world", region = "Canada")

xlim <- range(sites$lon) + c(-0.9, 0.9)
ylim <- range(sites$lat) + c(-0.4, 0.4)

map_plot <- ggplot() +
  geom_polygon(data = canada, aes(long, lat, group = group),
               fill = "grey96", colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = states, aes(long, lat, group = group),
               fill = "white", colour = "grey60", linewidth = 0.3) +
  geom_point(data = sites, aes(lon, lat, fill = n_years),
             shape = 21, size = 2.6, colour = "grey20", stroke = 0.3) +
  ggrepel::geom_text_repel(
    data = label_sites, aes(lon, lat, label = label),
    size = 2.5, colour = "grey15", segment.colour = "grey45",
    segment.size = 0.3, min.segment.length = 0, box.padding = 1.2,
    point.padding = 0.4, force = 12, max.overlaps = Inf,
    ## Push the labels off the trap network, towards the map margins.
    nudge_x = c(-0.9, 0.9, -0.7, -0.6), nudge_y = c(-0.35, 0.35, 0.35, -0.3),
    seed = 1
  ) +
  scale_fill_viridis_c(
    name   = "Years\nsampled",
    option = "viridis", direction = -1,
    breaks = c(1, 10, 20, 29), limits = c(1, max(sites$n_years))
  ) +
  annotate("text", x = xlim[1] + 0.15, y = ylim[2] - 0.12,
           ## Plain ASCII: the default PNG device does not render the accent.
           label = "Quebec", hjust = 0, size = 2.6, colour = "grey55") +
  annotate("text", x = xlim[2] - 0.15, y = ylim[2] - 0.12,
           label = "New Brunswick", hjust = 1, size = 2.6, colour = "grey55") +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    ## Claim checked against `sites` before rendering: 4 of the 75 localities
    ## have n_years >= 15, and the two 29-year series sit on opposite borders.
    title    = sprintf(
      "Traps covered the whole state, but only %d of %d localities ran 15 years or more",
      sum(sites$n_years >= 15), nrow(sites)),
    subtitle = "Light-trap localities with at least one daily capture series, 1961-1989 and 2016",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid       = element_blank(),
    plot.title       = element_text(size = 9, face = "plain", hjust = 0),
    plot.subtitle    = element_text(size = 8, colour = "grey35"),
    legend.key.width = unit(0.35, "cm"),
    legend.position  = "right"
  )

ggsave(file.path("output", "figures", "trap_sites_map.png"),
       plot = map_plot, width = 6.0, height = 5.6, dpi = 300)

message(sprintf("Described %d traps at %d localities over %d years (%d trap-nights).",
                n_distinct(trap_daily$flightID), nrow(sites),
                nrow(traps_per_year), sum(trap_daily$status == "counted")))
