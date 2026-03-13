#Author:
#Contact
#Email:
#Date
#purpose of this script

###########################################################


library(ggmap)
library(ggplot2)


for (w in 18:20){
  
path = paste0("C:/Users/a02386149/Box/katie_spatial/spatial_d_week_",w,".csv")
spatial_data <- read.csv(path)


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


# Reorder species so pipiens is plotted last (on top)
spatial_data$species <- factor(
  spatial_data$species,
  levels = c("Culex tarsalis", "Culex pipiens")
)

map <- ggmap(bg) +
  
  # 1. Draw tarsalis FIRST (triangles, underneath)
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
  
  # 2. Draw pipiens SECOND (circles, on top)
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
    name = "Species"
  ) +
  scale_shape_manual(
    values = c(
      "Culex pipiens"  = 16,  # circle
      "Culex tarsalis" = 16   # triangle
    ),
    name = "Species"
  ) +
  scale_alpha_manual(
    values = c(
      "Culex pipiens"  = 0.6,
      "Culex tarsalis" = 0.6
    ),
    guide = "none"
  ) +
  scale_size(
    range = c(3, 20),
    name = "Count"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

map

  file_name = paste0("C:/Users/a02386149/Box/katie_spatial/images/week_",w,".png")


  ggsave( filename = file_name, plot = map, width = 8, height = 8, dpi = 300 )
  
}
w=18





####################################################################
library(ggmap)
library(ggplot2)
library(sf)

# 1. Register API key once
register_google(key = "AIzaSyB6IptHVZgQeRsvUYtG2kTdHxYaGYOxPag")

# 2. Define a single bounding box for ALL weeks
lon_min <- -112.8   
lon_max <- -111.2   
lat_min <- 40.6
lat_max <- 41

# 3. Download the satellite map ONCE
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


# 5. Loop through weeks
for (w in 18:20) {
  
  path <- paste0("C:/Users/a02386149/Box/katie_spatial/spatial_d_week_", w, ".csv")
  spatial_data <- read.csv(path)
  
  # reorder species so pipiens is drawn on top
  spatial_data$species <- factor(
    spatial_data$species,
    levels = c("Culex tarsalis", "Culex pipiens")
  )
  
  map <- ggmap(bg) +
   
    # tarsalis first (underneath)
    geom_point(
      data = subset(spatial_data, species == "Culex tarsalis"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        alpha = species
      ),
      shape = 16
    ) +
    
    # pipiens second (on top)
    geom_point(
      data = subset(spatial_data, species == "Culex pipiens"),
      aes(
        x = longitude,
        y = latitude,
        size = females,
        color = species,
        alpha = species
      ),
      shape = 16
    ) +
    
    scale_color_manual(
      values = c(
        "Culex pipiens"  = "#00E5FF",
        "Culex tarsalis" = "#FF00AA"
      ),
      name = "Species"
    ) +
    scale_alpha_manual(
      values = c("Culex pipiens" = 0.6, "Culex tarsalis" = 0.6),
      guide = "none"
    ) +
    scale_size(
      range = c(3, 20),
      name = "Count"
    ) +
    labs(
      title = paste("Week", w),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
 
  
  # save image
  file_name <- paste0("C:/Users/a02386149/Box/katie_spatial/images/week_", w, ".png")
  ggsave(filename = file_name, plot = map, width = 8, height = 8, dpi = 300)
}




