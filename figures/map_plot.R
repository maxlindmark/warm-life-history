# Map plot
library(rnaturalearth)
library(rnaturalearthdata)
library(broom)
library(rgdal)
library(ggmap)
library(tidyverse)
library(sf)
library(png)
library(patchwork)
library(cowplot)

theme_set(theme_light(base_size = 12))

# Read UTM function
LongLatToUTM <- function(x, y, zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}


## BIG map
# Specify ranges for big map
ymin = 53; ymax = 69; xmin = 8; xmax = 32

map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", continent = "europe")

# Crop the polygon for plotting and efficiency:
# st_bbox(map_data) # find the rough coordinates
swe_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax))))

# Transform our map into UTM 33 coordinates, which is the equal-area projection we fit in:
utm_zone33 <- 32633
swe_coast_proj <- sf::st_transform(swe_coast, crs = utm_zone33)

# Add point to Forsmark & Biotest
df <- data.frame(lon =  18.1, lat = 60.4)

# Add UTM coords
utm_coords <- LongLatToUTM(df$lon, df$lat, zone = 33)
df$X <- utm_coords$X
df$Y <- utm_coords$Y

xmin <- 303379.1; xmax <- 958492.4; xrange <- xmax - xmin
ymin <- 5983578; ymax <- 6450163; yrange <- ymax - ymin

p1 <-
ggplot(swe_coast_proj) +
  geom_sf() +
  geom_point(data = df, aes(x = X, y = Y), size = 5, shape = 15, color = "red") +
  annotate("text", label = "Sweden", x = xmin + 0.33*xrange, y = ymin + 1.4*yrange,
           color = "black", size = 6) +
  labs(x = "Longitude", y = "Latitude") +
  geom_segment(aes(x = xmin + 0.58*xrange, y = ymin + 1.535*yrange,
                   xend = xmin + 1*xrange, yend = ymin + 1.73*yrange),
               arrow = arrow(length = unit(0.5, "cm")), color = "red",
               lineend = "butt", linejoin = "bevel", size = 1.5) +
  NULL

p1
  
## SMALL map (inset)
bt <- readPNG("figures/maps/biotest_map.png")

p2 <- ggdraw() +
  draw_image(bt)

p2

p1 + inset_element(p2, left = 0.45, bottom = 0.45, right = 1, top = 1)

ggsave("figures/map.png", width = 6.5, height = 6.5, dpi = 600)
