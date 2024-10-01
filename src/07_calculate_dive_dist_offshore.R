# ADD DISTANCE TO SHORE TO DIVES
# Add Distance to Shore + Bathymetry
library(geosphere)
library(rnaturalearth)
library(sf)

setwd('~/seal_telemetry/')

# Load Data
dives = read_csv(file = './data/L3/dive/Hg_2019-2023_DiveTypeClassification.csv')

# convert to sf object
proj = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'

dfsf <- dives %>% st_as_sf(coords = c('meanX_km','meanY_km')) %>% 
  st_set_crs(value = st_crs(proj)) %>% st_transform(4326)

coastline = rnaturalearth::ne_coastline(scale = 'medium')

# Define bounding box for the US East Coast and Canada
east_coast_bbox <- st_as_sfc(st_bbox(c(xmin = -76, ymin = 37, xmax = -59, ymax = 46), crs = 4326))

# Subset coastline data to include only the US East Coast and Canada
east_coast <- st_intersection(x = coastline, east_coast_bbox) %>%
  st_union() 
plot(east_coast)

# Convert the sf MULTILINESTRING object to a SpatialLines object
east_coast_sp <- as(east_coast, "Spatial")

# use dist2Line from geosphere - only works for WGS84 to get the distance to shore of each location in a trip
dist <- geosphere::dist2Line(p = st_coordinates(dfsf), line = east_coast_sp)

#combine initial data with distance to coastline
dist_ = as.data.frame(dist)
dist_ = dist_ %>% rename(LANDLON = lon, LANDLAT = lat, distanttoshore = distance) %>% dplyr::select(-c(ID))

# Bind distance information to original dataset
dives = cbind(dives, dist_)
write_csv(x = dives, file = './data/L3/dive/Hg_2019-2023_DiveTypeClassification_withDistanceToShore.csv')