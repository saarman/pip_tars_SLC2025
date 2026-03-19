Weekly Maps **Cx. pipiens** and **Cx. tarsalis** abundance: SLC 2025
field season
================
Norah Saarman
2026-03-19

- [Setup](#setup)
- [SQRT Weekly Maps (Square-Root
  Scaling)](#sqrt-weekly-maps-square-root-scaling)
  - [Zoom 11](#zoom-11)
  - [Zoom 12](#zoom-12)
  - [SQRT weekly GIFS for Zoom 11 and
    12](#sqrt-weekly-gifs-for-zoom-11-and-12)
- [SQRT Seasonal Seasonal averages (Square-Root
  Scaling)](#sqrt-seasonal-seasonal-averages-square-root-scaling)
  - [Zoom 11](#zoom-11-1)
  - [Zoom 12](#zoom-12-1)
- [Raw Weekly Maps (Raw Scaling)](#raw-weekly-maps-raw-scaling)
  - [Zoom 11](#zoom-11-2)
  - [Zoom 12](#zoom-12-2)
  - [Raw weekly GIFs for raw, zoom 11 and
    12](#raw-weekly-gifs-for-raw-zoom-11-and-12)

# Setup

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

# SQRT Weekly Maps (Square-Root Scaling)

Loop over all weeks with square-root scaling

Notes:

Weekly file contained multiple observations per site and species (from
different trap types or sampling dates), so added code to combine
multiple observations from different trap types and sampling dates per
week.

Also applied a square-root transformation to the variable before it is
mapped to point size.

### Zoom 11

``` r
# Weeks to plot
weeks <- 15:40

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025)

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
global_max <- all_data %>%
  filter(disease_week %in% weeks) %>%
  group_by(disease_week, site_code, site_name, longitude, latitude, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 17797"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  spatial_data_raw <- all_data %>%
    filter(disease_week == w)
  
  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

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
        size = count,
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
        size = count,
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

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-1.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-2.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-3.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-4.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-5.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-6.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-7.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-8.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-9.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-10.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-11.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-12.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-13.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-14.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-15.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-16.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-17.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-18.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-19.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-20.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-21.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-22.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-23.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-24.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-25.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom11-26.png)<!-- -->

### Zoom 12

``` r
# Weeks to plot
weeks <- 15:40

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025)

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
global_max <- all_data %>%
  filter(disease_week %in% weeks) %>%
  group_by(disease_week, site_code, site_name, longitude, latitude, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 17797"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  spatial_data_raw <- all_data %>%
    filter(disease_week == w)
  
  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

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
        size = count,
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
        size = count,
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

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-1.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-2.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-3.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-4.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-5.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-6.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-7.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-8.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-9.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-10.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-11.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-12.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-13.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-14.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-15.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-16.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-17.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-18.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-19.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-20.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-21.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-22.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-23.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-24.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-25.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-weeks-zoom12-26.png)<!-- -->

### SQRT weekly GIFS for Zoom 11 and 12

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

# SQRT Seasonal Seasonal averages (Square-Root Scaling)

``` r
library(dplyr)

# read data
pipiens  <- read.csv("../data/pipiens_2025.csv")
tarsalis <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens, tarsalis)

# define seasons from week numbers
all_data <- all_data %>%
  mutate(
    season_calc = case_when(
      disease_week %in% 15:22 ~ "early",
      disease_week %in% 23:31 ~ "mid",
      disease_week %in% 32:40 ~ "late",
      TRUE ~ NA_character_
    )
  )

# compute mean counts across weeks
seasonal_avg <- all_data %>%
  filter(!is.na(season_calc)) %>%
  group_by(season_calc, site_code, site_name, longitude, latitude, species, disease_week) %>%
  summarise(weekly_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(season_calc, site_code, site_name, longitude, latitude, species) %>%
  summarise(mean_count = mean(weekly_count, na.rm = TRUE), .groups = "drop")

seasonal_avg
```

    ## # A tibble: 333 × 7
    ##    season_calc site_code site_name         longitude latitude species mean_count
    ##    <chr>           <int> <chr>                 <dbl>    <dbl> <chr>        <dbl>
    ##  1 early               1 RAC Soccer            -112.     40.8 Culex …       5.88
    ##  2 early               1 RAC Soccer            -112.     40.8 Culex …      93   
    ##  3 early               3 Hinckley              -112.     40.8 Culex …      10   
    ##  4 early               3 Hinckley              -112.     40.8 Culex …     160.  
    ##  5 early               4 Rudy                  -112.     40.8 Culex …       6.25
    ##  6 early               4 Rudy                  -112.     40.8 Culex …      83.4 
    ##  7 early               5 Amelia Earhart (…     -112.     40.8 Culex …       6   
    ##  8 early               5 Amelia Earhart (…     -112.     40.8 Culex …      49.7 
    ##  9 early               9 Drechsel              -112.     40.8 Culex …       2.5 
    ## 10 early               9 Drechsel              -112.     40.8 Culex …     151.  
    ## # ℹ 323 more rows

### Zoom 11

``` r
# Seasons to plot
seasons <- c("early", "mid", "late")

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025) %>%
  mutate(
    season_plot = case_when(
      disease_week %in% 15:22 ~ "early",
      disease_week %in% 23:31 ~ "mid",
      disease_week %in% 32:40 ~ "late",
      TRUE ~ NA_character_
    )
  )

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
# 4. Find the global max across all seasons so the size scale is identical
global_max <- all_data %>%
  filter(!is.na(season_plot)) %>%
  group_by(season_plot, site_code, site_name, longitude, latitude, species, disease_week) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(season_plot, site_code, site_name, longitude, latitude, species) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(mean_count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 6115.55555555556"

``` r
# 5. Loop over all seasons
for (s in seasons) {

  spatial_data_raw <- all_data %>%
    filter(season_plot == s)
  
  # Combine multiple observations from different trap types and sampling dates per week,
  # then average across weeks within season
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species, disease_week) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

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
        size = mean_count,
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
        size = mean_count,
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
      name = "Mean count"
    ) +

    labs(
      title = paste(tools::toTitleCase(s), "season"),
      x = "Longitude",
      y = "Latitude"
    ) +

    theme_minimal()

  print(map)

  file_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/",
    s,
    "_zoom11.png"
  )
  
  # save one high res image
  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)
}
```

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom11-1.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom11-2.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom11-3.png)<!-- -->

### Zoom 12

``` r
# Seasons to plot
seasons <- c("early", "mid", "late")

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025) %>%
  mutate(
    season_plot = case_when(
      disease_week %in% 15:22 ~ "early",
      disease_week %in% 23:31 ~ "mid",
      disease_week %in% 32:40 ~ "late",
      TRUE ~ NA_character_
    )
  )

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
# 4. Find the global max across all seasons so the size scale is identical
global_max <- all_data %>%
  filter(!is.na(season_plot)) %>%
  group_by(season_plot, site_code, site_name, longitude, latitude, species, disease_week) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(season_plot, site_code, site_name, longitude, latitude, species) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(mean_count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 6115.55555555556"

``` r
# 5. Loop over all seasons
for (s in seasons) {

  spatial_data_raw <- all_data %>%
    filter(season_plot == s)
  
  # Combine multiple observations from different trap types and sampling dates per week,
  # then average across weeks within season
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species, disease_week) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

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
        size = mean_count,
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
        size = mean_count,
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
      name = "Mean count"
    ) +

    labs(
      title = paste(tools::toTitleCase(s), "season"),
      x = "Longitude",
      y = "Latitude"
    ) +

    theme_minimal()

  print(map)

  file_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/",
    s,
    "_zoom12.png"
  )

  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)
}
```

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom12-1.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom12-2.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/all-seasons-zoom12-3.png)<!-- -->

# Raw Weekly Maps (Raw Scaling)

Loop over all weeks with raw scaling.

Notes:

Weekly file contained multiple observations per site and species (from
different trap types or sampling dates), so added code to combine
multiple observations from different trap types and sampling dates per
week.

Point size is mapped directly to the raw count value for visualization
only.

### Zoom 11

``` r
# Weeks to plot
weeks <- 15:40

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025)

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
  zoom = 11
)
```

    ## ! Bounding box given to Google - spatial extent only approximate.

    ## ℹ <https://maps.googleapis.com/maps/api/staticmap?center=40.8,-112&zoom=11&size=640x640&scale=2&maptype=satellite&language=en-EN&key=xxx>

``` r
# 4. Find the global max across all weeks so the size scale is identical
global_max <- all_data %>%
  filter(disease_week %in% weeks) %>%
  group_by(disease_week, site_code, site_name, longitude, latitude, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 17797"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  spatial_data_raw <- all_data %>%
    filter(disease_week == w)

  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

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
        size = count,
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
        size = count,
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
        "Culex pipiens"  = 16,
        "Culex tarsalis" = 16
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
      range = c(3, 20),
      limits = c(1, global_max),
      breaks = c(1, 10, 50, 500, 5000),
      name = "Count"
    ) +

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
    "_zoom11_raw.png"
  )

  # save one high res image
  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)

  gif_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/week_",
    w,
    "_zoom11_gif_raw.png"
  )

  # save one low res image for gif
  ggsave(filename = gif_name, plot = map, width = 5, height = 5, dpi = 140)
}
```

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-1.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-2.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-3.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-4.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-5.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-6.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-7.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-8.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-9.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-10.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-11.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-12.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-13.png)<!-- -->

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-14.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-15.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-16.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-17.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-18.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-19.png)<!-- -->

    ## Warning: Removed 4 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 4 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 4 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-20.png)<!-- -->![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-21.png)<!-- -->

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-22.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-23.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-24.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-25.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z11-raw-26.png)<!-- -->

### Zoom 12

``` r
# Weeks to plot
weeks <- 15:40

# Read the two species files once and combine them
pipiens_2025  <- read.csv("../data/pipiens_2025.csv")
tarsalis_2025 <- read.csv("../data/tarsalis_2025.csv")

all_data <- bind_rows(pipiens_2025, tarsalis_2025)

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
  zoom = 12
)
```

    ## ! Bounding box given to Google - spatial extent only approximate.

    ## ℹ <https://maps.googleapis.com/maps/api/staticmap?center=40.8,-112&zoom=12&size=640x640&scale=2&maptype=satellite&language=en-EN&key=xxx>

``` r
# 4. Find the global max across all weeks so the size scale is identical
global_max <- all_data %>%
  filter(disease_week %in% weeks) %>%
  group_by(disease_week, site_code, site_name, longitude, latitude, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  summarise(global_max = max(count, na.rm = TRUE)) %>%
  pull(global_max)

print(paste("global max:", global_max))
```

    ## [1] "global max: 17797"

``` r
# 5. Loop over all weeks
for (w in weeks) {

  spatial_data_raw <- all_data %>%
    filter(disease_week == w)

  # Combine multiple observations from different trap types and sampling dates per week
  spatial_data <- spatial_data_raw %>%
    group_by(site_code, site_name, longitude, latitude, species) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

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
        size = count,
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
        size = count,
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
        "Culex pipiens"  = 16,
        "Culex tarsalis" = 16
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
      range = c(3, 20),
      limits = c(1, global_max),
      breaks = c(1, 10, 50, 500, 5000),
      name = "Count"
    ) +

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
    "_zoom12_raw.png"
  )

  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)

  gif_name <- paste0(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/week_",
    w,
    "_zoom12_gif_raw.png"
  )

  # save one low res image for gif
  ggsave(filename = gif_name, plot = map, width = 5, height = 5, dpi = 140)
}
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-1.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-2.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-3.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-4.png)<!-- -->

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-5.png)<!-- -->

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-6.png)<!-- -->

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 5 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-7.png)<!-- -->

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-8.png)<!-- -->

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-9.png)<!-- -->

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-10.png)<!-- -->

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-11.png)<!-- -->

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-12.png)<!-- -->

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-13.png)<!-- -->

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-14.png)<!-- -->

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-15.png)<!-- -->

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-16.png)<!-- -->

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-17.png)<!-- -->

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 10 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-18.png)<!-- -->

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-19.png)<!-- -->

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-20.png)<!-- -->

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-21.png)<!-- -->

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 16 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-22.png)<!-- -->

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-23.png)<!-- -->

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 19 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-24.png)<!-- -->

    ## Warning: Removed 13 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 13 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 13 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 17 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-25.png)<!-- -->

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_maps_and_gif_google_SLC2025_files/figure-gfm/z12-raw-26.png)<!-- -->

### Raw weekly GIFs for raw, zoom 11 and 12

``` r
library(gifski)

# Create GIF for zoom 11
gif_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11"

gif_imgs <- list.files(
  gif_dir,
  pattern = "week_.*_zoom11_gif\\_raw.png$",
  full.names = TRUE
)

gif_weeks <- as.numeric(sub(".*week_([0-9]+).*", "\\1", basename(gif_imgs)))
gif_imgs <- gif_imgs[order(gif_weeks)]

gifski(
  png_files = gif_imgs,
  gif_file  = file.path(gif_dir, "mosquito_2025_zoom11_raw.gif"),
  width     = 700,
  height    = 700,
  delay     = 0.6
)
```

    ## [1] "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom11/mosquito_2025_zoom11_raw.gif"

``` r
# Create GIF for zoom 12
gif_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12"

gif_imgs <- list.files(
  gif_dir,
  pattern = "week_.*_zoom12_gif\\_raw.png$",
  full.names = TRUE
)

gif_weeks <- as.numeric(sub(".*week_([0-9]+).*", "\\1", basename(gif_imgs)))
gif_imgs <- gif_imgs[order(gif_weeks)]

gifski(
  png_files = gif_imgs,
  gif_file  = file.path(gif_dir, "mosquito_2025_zoom12_raw.gif"),
  width     = 700,
  height    = 700,
  delay     = 0.6
)
```

    ## [1] "/uufs/chpc.utah.edu/common/home/saarman-group1/pip_tars_SLC2025/figures/d_weeks_2025_maps/zoom12/mosquito_2025_zoom12_raw.gif"
