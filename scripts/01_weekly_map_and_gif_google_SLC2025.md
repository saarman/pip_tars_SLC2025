Weekly **Cx. pipiens** and **Cx. tarsalis** abundance: SLC 2025 field
season
================
Norah Saarman
2026-03-13

- [Create Weekly Maps (Square-Root
  Scaling)](#create-weekly-maps-square-root-scaling)
  - [Google Zoom 11](#google-zoom-11)
  - [Google Zoom 12](#google-zoom-12)
- [Create GIFS for Zoom 11 and 12](#create-gifs-for-zoom-11-and-12)

Weekly maps of **Culex pipiens** and **Culex tarsalis** abundance across
SLC for the 2025 field season, overlay onto Google satellite image.

Load libraries:

``` r
library(ggmap)
```

    ## Loading required package: ggplot2

    ## ℹ Google's Terms of Service: <https://mapsplatform.google.com>
    ##   Stadia Maps' Terms of Service: <https://stadiamaps.com/terms-of-service>
    ##   OpenStreetMap's Tile Usage Policy: <https://operations.osmfoundation.org/policies/tiles>
    ## ℹ Please cite ggmap if you use it! Use `citation("ggmap")` for details.

``` r
library(ggplot2)
library(sf)
```

    ## Linking to GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1; sf_use_s2() is TRUE

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(gifski)
#library(magick)
```

# Create Weekly Maps (Square-Root Scaling)

Loop over all weeks with square-root scaling

Notes:

Weekly file contained multiple observations per site and species (from
different trap types or sampling dates), so added code to combine
multiple observations from different trap types and sampling dates per
week.

Also applied a square-root transformation to the variable before it is
mapped to point size.

## Google Zoom 11

``` r
# Weeks to plot
weeks <- 18:40

# 1. Register your Google Maps API key
register_google(key = "AIzaSyB6IptHVZgQeRsvUYtG2kTdHxYaGYOxPag")

# 2. Build a bounding box
lon_min <- -112.8
lon_max <- -111.2
lat_min <- 40.6
lat_max <- 41

# 3. Download the satellite map using the bounding box
bg <- get_map(
  location = c(
    left   = lon_min,
    bottom = lat_min,
    right  = lon_max,
    top    = lat_max
  ),
  maptype = "satellite",
  zoom = 11   #zoom in with zoom = 12
)
```

    ## ! Bounding box given to Google - spatial extent only approximate.

    ## ℹ <https://maps.googleapis.com/maps/api/staticmap?center=40.8,-112&zoom=11&size=640x640&scale=2&maptype=satellite&language=en-EN&key=xxx>

``` r
# 4. Find the global max across all weeks so the size scale is identical
global_max <- max(sapply(weeks, function(w) {
  path <- paste0("../data/d_weeks_2025/spatial_d_week_", w, ".csv")
  spatial_data <- read.csv(path)
  max(spatial_data$females, na.rm = TRUE)
}))
print(paste("global max:",global_max))
```

    ## [1] "global max: 10271"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  path <- paste0("../data/d_weeks_2025/spatial_d_week_", w, ".csv")
  spatial_data_raw <- read.csv(path)
  
  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
  group_by(site_code, name, longitude, latitude, species) %>%
  summarise(females = sum(females, na.rm = TRUE), .groups = "drop")

  # Reorder species so pipiens is plotted last (on top)
  spatial_data$species <- factor(
    spatial_data$species,
    levels = c("Culex tarsalis", "Culex pipiens")
  )

  map <- ggmap(bg) +

    # 1. Draw tarsalis FIRST (underneath)
    geom_point(
      data = subset(spatial_data, species == "Culex tarsalis"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        shape = species,
        alpha = species
      )
    ) +

    # 2. Draw pipiens SECOND (on top)
    geom_point(
      data = subset(spatial_data, species == "Culex pipiens"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        shape = species,
        alpha = species
      )
    ) +

    scale_color_manual(
      values = c(
        "Culex pipiens"  = "#00E5FF",
        "Culex tarsalis" = "#FF00AA"
      ),
      name = "Species",
      drop = FALSE
    ) +

    scale_shape_manual(
      values = c(
        "Culex pipiens"  = 16,  # circle
        "Culex tarsalis" = 16   # circle
      ),
      name = "Species",
      drop = FALSE
    ) +

    scale_alpha_manual(
      values = c(
        "Culex pipiens"  = 0.6,
        "Culex tarsalis" = 0.6
      ),
      guide = "none",
      drop = FALSE
    ) +

    scale_size_continuous(
    trans = "sqrt",
    range = c(3, 20),
    limits = c(1, global_max),
    breaks = c(1, 10, 50, 500, 5000),
    name = "Count"
    ) +
    
    #scale_size( # calls continuous automatically
    #  range = c(3, 20),
    #  limits = c(1, global_max),
    #  breaks = c(1, 10, 50, 100, 500, 5000),
    #  name = "Count"
    #) +

    labs(
      title = paste("Week", w),
      x = "Longitude",
      y = "Latitude"
    ) +

    theme_minimal()

  print(map)

  file_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/week_",
    w,
    "_zoom11.png"
  )
  
  # save one high res image
  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)


  gif_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/week_",
    w,
    "_zoom11_gif.png"
  )
  
  # save one low res image for gif
  ggsave(filename = gif_name, plot = map, width = 5, height = 5, dpi = 140)
}
```

![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-1.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-2.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-3.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-4.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-5.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-6.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-7.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-8.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-9.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-10.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-11.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-12.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-13.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-14.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-15.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-16.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-17.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-18.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-19.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-20.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-21.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-22.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-23.png)<!-- -->

## Google Zoom 12

``` r
# Weeks to plot
weeks <- 18:40

# 1. Register your Google Maps API key
register_google(key = "AIzaSyB6IptHVZgQeRsvUYtG2kTdHxYaGYOxPag")

# 2. Build a bounding box
lon_min <- -112.8
lon_max <- -111.2
lat_min <- 40.6
lat_max <- 41

# 3. Download the satellite map using the bounding box
bg <- get_map(
  location = c(
    left   = lon_min,
    bottom = lat_min,
    right  = lon_max,
    top    = lat_max
  ),
  maptype = "satellite",
  zoom = 12   #zoom out with zoom = 11
)
```

    ## ! Bounding box given to Google - spatial extent only approximate.

    ## ℹ <https://maps.googleapis.com/maps/api/staticmap?center=40.8,-112&zoom=12&size=640x640&scale=2&maptype=satellite&language=en-EN&key=xxx>

``` r
# 4. Find the global max across all weeks so the size scale is identical
global_max <- max(sapply(weeks, function(w) {
  path <- paste0("../data/d_weeks_2025/spatial_d_week_", w, ".csv")
  spatial_data <- read.csv(path)
  max(spatial_data$females, na.rm = TRUE)
}))
print(paste("global max:",global_max))
```

    ## [1] "global max: 10271"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  path <- paste0("../data/d_weeks_2025/spatial_d_week_", w, ".csv")
  spatial_data_raw <- read.csv(path)
  
  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
  group_by(site_code, name, longitude, latitude, species) %>%
  summarise(females = sum(females, na.rm = TRUE), .groups = "drop")

  # Reorder species so pipiens is plotted last (on top)
  spatial_data$species <- factor(
    spatial_data$species,
    levels = c("Culex tarsalis", "Culex pipiens")
  )

  map <- ggmap(bg) +

    # 1. Draw tarsalis FIRST (underneath)
    geom_point(
      data = subset(spatial_data, species == "Culex tarsalis"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        shape = species,
        alpha = species
      )
    ) +

    # 2. Draw pipiens SECOND (on top)
    geom_point(
      data = subset(spatial_data, species == "Culex pipiens"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        shape = species,
        alpha = species
      )
    ) +

    scale_color_manual(
      values = c(
        "Culex pipiens"  = "#00E5FF",
        "Culex tarsalis" = "#FF00AA"
      ),
      name = "Species",
      drop = FALSE
    ) +

    scale_shape_manual(
      values = c(
        "Culex pipiens"  = 16,  # circle
        "Culex tarsalis" = 16   # circle
      ),
      name = "Species",
      drop = FALSE
    ) +

    scale_alpha_manual(
      values = c(
        "Culex pipiens"  = 0.6,
        "Culex tarsalis" = 0.6
      ),
      guide = "none",
      drop = FALSE
    ) +

    scale_size_continuous(
    trans = "sqrt",
    range = c(3, 20),
    limits = c(1, global_max),
    breaks = c(1, 10, 50, 500, 5000),
    name = "Count"
    ) +
    
    #scale_size( # calls continuous automatically
    #  range = c(3, 20),
    #  limits = c(1, global_max),
    #  breaks = c(1, 10, 50, 100, 500, 5000),
    #  name = "Count"
    #) +

    labs(
      title = paste("Week", w),
      x = "Longitude",
      y = "Latitude"
    ) +

    theme_minimal()

  print(map)

  file_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/week_",
    w,
    "_zoom12.png"
  )

  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)
  
  
    gif_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/week_",
    w,
    "_zoom12_gif.png"
  )
  
  # save one low res image for gif
  ggsave(filename = gif_name, plot = map, width = 5, height = 5, dpi = 140)
}
```

![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-1.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-2.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-3.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-4.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-5.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-6.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-7.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-8.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-9.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-10.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-11.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-12.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-13.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-14.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-15.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-16.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-17.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-18.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-19.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-20.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-21.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-22.png)<!-- -->![](01_weekly_map_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-23.png)<!-- -->

# Create GIFS for Zoom 11 and 12

``` r
library(gifski)

# Create GIF for zoom 11
gif_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11"

gif_imgs <- list.files(
  gif_dir,
  pattern = "week_.*_zoom11_gif\\.png$",
  full.names = TRUE
)

gif_weeks <- as.numeric(sub(".*week_([0-9]+).*", "\\1", basename(gif_imgs)))
gif_imgs <- gif_imgs[order(gif_weeks)]

gifski(
  png_files = gif_imgs,
  gif_file  = file.path(gif_dir, "mosquito_2025_zoom11.gif"),
  width     = 700,
  height    = 700,
  delay     = 0.6
)
```

    ## [1] "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/mosquito_2025_zoom11.gif"

``` r
# Create GIF for zoom 12
gif_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12"

gif_imgs <- list.files(
  gif_dir,
  pattern = "week_.*_zoom12_gif\\.png$",
  full.names = TRUE
)

gif_weeks <- as.numeric(sub(".*week_([0-9]+).*", "\\1", basename(gif_imgs)))
gif_imgs <- gif_imgs[order(gif_weeks)]

gifski(
  png_files = gif_imgs,
  gif_file  = file.path(gif_dir, "mosquito_2025_zoom12.gif"),
  width     = 700,
  height    = 700,
  delay     = 0.6
)
```

    ## [1] "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/mosquito_2025_zoom12.gif"
