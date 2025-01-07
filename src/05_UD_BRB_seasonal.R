#EIH
#2024-06-07

# Kernel Density Utilization Distributions - Biased Brownian Bridge
rm(list = ls())
gc()
library(adehabitatHR)
library(adehabitatLT)
library(marmap)
library(raster)
library(sf)
library(tidyr)
library(dplyr)
library(readr)


setwd("~/seal_telemetry/")

source('./READ-PSB-MoveSeals/src/fxns_helper.R')
rm(list = setdiff(ls(), 'volras'))

# Create land for plotting
# Load the land spatial polygons dataframe and subset based on relevant landmasses
land = rnaturalearth::ne_states()
land = land[land$admin %in% c("Canada", 'United States of America'),]
land = land[which(land$postal %in% c("ME", "NH", "VT", "RI", "CT", "MA", "PA", "NY", "DE", "NJ", "VA", "MD", "NC", 'OH', 'WV', 'SC', 'KY', 'TN', 'IL', 'WI', 'MN', 'IA', 'IN', 'MI') | land$region == "Eastern Canada"),]
plot(land['featurecla'])

# Now minor islands (islands less than 2 square km)
minor_islands = st_read("./data/shapefiles/ne_10m_minor_islands/ne_10m_minor_islands.shp") %>%
  mutate(name = featurecla) %>%
  dplyr::select(name, geometry)

#get the bounding box of the land and subset minor islands based on this bounding box
bounding_box <- st_bbox(land) %>% st_as_sfc()
minor_islands = st_intersection(minor_islands, bounding_box)

# Combine the land and minor islands
hires_land = bind_rows(land, minor_islands)
rm(minor_islands, land)

# DO WE HAVE ENOUGH INDIVIDUALS FOR POPULATION LEVEL INFERENCE
# Shimada et al. 2021 - e.g. rarefaction curves. If you keep adding more individuals and the overlap index is high, you likely have sufficient individuals for population UD
# Methods in SDLfilter 


df = read_csv("./data/L2/SSM/Hg-2019-2023-SSM-Tracks-CompleteTrips-Reroute.csv") %>% 
  arrange(ptt, date) %>%
  dplyr::select(ptt, TripID, id, date, lon, lat, x, y, x.se, y.se)


hist(df$x.se, breaks = 500, xlim = c(0,30))
hist(df$y.se, breaks = 200, xlim = c(0,30))
quantile(df$x.se, 0.99)
quantile(df$y.se, 0.99)

df = df %>% arrange(id, date)
dupes = which(duplicated(x = df[,c('ptt', 'date')]))

df = df[-dupes,]

# Subset to desired time period
df = df[df$date < as.POSIXct(x = '2023-06-01 00:00:00', tz = 'UTC'), ]

#PROJECT TO World Mercator in KM (animotum spits out x, y)
lon0 = mean(df$lon)
lat0 = mean(df$lat)
lat1 = min(df$lat)
lat2 = max(df$lat)
ogproj = "+proj=merc +datum=WGS84 +units=km +no_defs"
# need equal area projection for kernel density analysis
projAEA = sprintf("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                  lat1, lat2, lat0, lon0)

dfsf = sf::st_as_sf(x = df, coords = c('x', 'y')) %>% 
  st_set_crs(CRS(ogproj)) %>% st_transform(projAEA) %>%
  mutate(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2]) %>% 
  arrange(ptt, date) %>% 
  mutate(Month = as.numeric(format(date, '%m')))

dfsf$ptt = as.character(dfsf$ptt)
uid = unique(dfsf$ptt)

#Plot gut check 
library(ggplot2)
ggplot() +
  geom_sf(data = hires_land) +
  geom_sf(data = dfsf)

## Data range with 40km buffer
buff <- 40000
xmin <- min(dfsf$x) - buff; xmax <- max(dfsf$x) + buff
ymin <- min(dfsf$y) - buff; ymax <- max(dfsf$y) + buff

## Make a grid layer
cell.size <- 1000 # (1000 m = 1km)
x <- seq(xmin, xmax, cell.size)
y <- seq(ymin, ymax, cell.size)
xy.df <- expand.grid(x = x, y = y)
xy.coords <- SpatialPixels(SpatialPoints(xy.df), proj4string = projAEA)
xy.sp <- SpatialPoints(xy.coords)
z <- rep(1, nrow(xy.df))
xyz <- cbind(xy.df, z)
grid_spdf <- SpatialPixelsDataFrame(xy.coords, xyz)

rm(list = setdiff(ls(), c('dfsf', 'grid_spdf', 'projAEA', 'uid', 'raster.standardize', 'volras', 'hires_land')))
gc()


# Loop over each season, get individual ids in that season and calculate individual UDs
# Calculate with an 15 day smoothing window
# Get sequence of Julian days to calculate 

## Creates an object of class Itraj
## UD per seal
season_L = list()
seasons = data.frame(Season = c(rep('Winter', 2), rep('Spring', 3), rep('Summer', 3), rep('Fall', 3)), Month = c(1,2,3,4,5,6,7,8,9,10,11))
sid = unique(seasons$Season)

# Get the month that corresponds to the maximum number of positions... assign accordingly to season
dfsf = dfsf %>% group_by(TripID) %>% 
  mutate(TripMonth = names(which.max(table(Month))))

for(i in 1:length(sid)){
  
  ## Tracking data
  months = seasons$Month[seasons$Season == sid[i]]
  # assign to season based on the month that contains the maximum number of positions for a given trip
  seal.data <- dfsf[dfsf$TripMonth %in% months, ]
  
  seal.data = seal.data[order(seal.data$TripID, seal.data$date), ]
  table(seal.data$Month)
  
  xy = data.frame(x = seal.data$x, y = seal.data$y)
  DateTime = seal.data$date
  id = seal.data$TripID
  #burst = seal.data$id
  data.ltraj <- as.ltraj(xy = xy, date = DateTime, 
                         id = id,
                         proj4string = CRS(projAEA))
  
  ## Parameters for BRB
  ### A BIT ABOUT biased random bridge BRB: 
  ## Similar to the Brownian bridge approach (see ?kernelbb), with several noticeable improvements. 
  ## Brownian bridge approach supposes that the animal moves in a purely random fashion from the starting relocation and reaches the next relocation randomly. 
  # The BRB approach goes further by adding an advection component (i.e., a "drift") to the purely diffusive movement: 
  # is is supposed that the animal movement is governed by a drift component (a general tendency to move in the direction of the next relocation) 
  # and a diffusion component (tendency to move in other directions than the direction of the drift).
  
  # An important aspect of the BRB approach is that the drift component is allowed to change in direction and strength from one step to the other, 
  # but should remain constant during each of them. For this reason, it is required to set an upper time threshold Tmax. Steps characterized by a longer duration are not taken into account into the estimation of the pdf. 
  # This upper threshold should be based on biological grounds.
  tmax = 2*3600 # 
  Lmin = 100 # smallest distance below which we would consider the animal not to be moving
  hmin = 20000 # 20 km - this bandwidth smoother reflects the 99th percentile of the x and y standard errors from SSM 
  
  # Estimate the diffusion coefficient with maximum likelihood
  vv <- BRB.likD(data.ltraj, Tmax=tmax, Lmin=Lmin)
  ud = BRB(ltr = data.ltraj, D = vv, Tmax = tmax, Lmin = Lmin, hmin = hmin, grid = grid_spdf, type = 'UD', filtershort = TRUE)
  
  ## Convert the UD to raster
  r = lapply(X = ud, FUN = function(x){
    tmp = raster::raster(x)
    crs(tmp) = CRS(projAEA)
    return(tmp)
  })
  names(r) = unique(seal.data$TripID)
  season_L[[i]] <- r
  names(season_L)[[i]] = sid[i]
}


sex = read_csv(file = './data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv') %>% dplyr::select(ptt, sex)
FID = sex$ptt[sex$sex == "F"]
MID = sex$ptt[sex$sex == 'M']
Fpoly = list() 
Mpoly = list()
FStackL = list()
MStackL = list()

for (i in 1:length(sid)){
  curseason = sid[i]
  # subset the season
  curL = season_L[[curseason]]
  
  # Get list of ptts in season i
  TripID = names(curL)
  ptt = as.numeric(sapply(strsplit(TripID, '-'), '[[', 1))
  
  # Get the FIDs 
  fidx = which(ptt %in% FID)
  FStackR = raster::stack(curL[c(fidx)])
  names(FStackR) = TripID[which(ptt %in% FID)]
  
  # Get the MIDs
  midx = which(ptt %in% MID)
  MStackR = raster::stack(curL[c(midx)])
  names(MStackR) = TripID[which(ptt %in% MID)]
  
  # Save the F and M stacks
  raster::writeRaster(x = FStackR, filename = paste('./data/L3/UD/Hg_PreConstruction_FemaleUDs_', curseason, '.grd', sep = ''), overwrite = T)
  raster::writeRaster(x = MStackR, filename = paste('./data/L3/UD/Hg_PreConstruction_MaleUDs_', curseason, '.grd', sep = ''), overwrite = T)
  
}


for (i in 1:length(sid)){
  
  curseason = sid[i]
  
  FStackR = raster::stack(x = paste('./data/L3/UD/UD_rasters/Hg_PreConstruction_FemaleUDs_', curseason, '.grd', sep = ''))
  MStackR = raster::stack(x = paste('./data/L3/UD/UD_rasters/Hg_PreConstruction_MaleUDs_', curseason, '.grd', sep = ''))
  
  Fptts = gsub(pattern = 'X', replacement = '', x = sapply(strsplit(names(FStackR), '\\.'), '[[', 1))
  Mptts = gsub(pattern = 'X', replacement = '', x = sapply(strsplit(names(MStackR), '\\.'), '[[', 1))
  
  
  # Calculate mean (this is the 'collective UD' for each sex for season i)
  curCollectiveFUD = raster::calc(x = FStackR, fun = mean)
  curCollectiveMUD = raster::calc(x=MStackR, fun = mean)
  
  perc = c(95, 50)
  res = lapply(X = 1:length(perc), FUN = function(x) {
    tmp = volras(x = curCollectiveFUD, percent = perc[x], simplify = TRUE, returnwhat = 'both')
    return(tmp)
  })
  
  resF = bind_rows(res[[1]][[1]], res[[2]][[1]])
  resF$area_sqkm[1:2] = c(res[[1]][[2]], res[[2]][[2]])
  resF$PercentVolume = as.character(perc)
  resF$Sex = 'Female'
  resF$Season = curseason
  resF$N = length(unique(Fptts))
  resF$Ntrips = nlayers(FStackR)
  resF$type = 'collective'
  Fpoly[[i]] = resF
  
  res = lapply(X = 1:length(perc), FUN = function(x) {
    tmp = volras(x = curCollectiveMUD, percent = perc[x], simplify = TRUE, returnwhat = 'both')
    return(tmp)
  })
  
  resM = bind_rows(res[[1]][[1]], res[[2]][[1]])
  resM$area_sqkm[1:2] = c(res[[1]][[2]], res[[2]][[2]])
  resM$PercentVolume = as.character(perc)
  resM$Sex = 'Male'
  resM$Season = curseason
  resM$N = length(unique(Mptts))
  resM$Ntrips = nlayers(MStackR)
  resM$type = 'collective'
  Mpoly[[i]] = resM
}

FUD = dplyr::bind_rows(Fpoly) %>% sf::st_transform(4326)
MUD = dplyr::bind_rows(Mpoly) %>% sf::st_transform(4326)

st_write(obj = FUD, dsn = './data/L3/UD/colleHg_PreConstruction_FemaleCollectiveUDs_AllSeasonPolygons.shp', append = F)

st_write(obj = MUD, dsn = './data/L3/UD/Hg_PreConstruction_MaleCollectiveUDs_AllSeasonPolygons.shp', append = F)
