library(adehabitatHR)
library(adehabitatLT)
library(raster)
library(sf)
library(stars)
library(tidyr)
library(dplyr)
library(readr)
# 3D Kernel 
# https://cran.r-project.org/web/packages/splancs/splancs.pdf 

setwd("~/seal_telemetry/")


# DO WE HAVE ENOUGH INDIVIDUALS FOR POPULATION LEVEL INFERENCE
# Shimada et al. 2021 - e.g. rarefaction curves. If you keep adding more individuals and the overlap index is high, you likely have sufficient individuals for population UD
# Methods in SDLfilter 
library(SDLfilter)


df = read_csv("./data/L2/SSM/Hg-2019-2023-SMM-Tracks-CompleteTrips-Reroute.csv") %>% arrange(ptt, date) %>% dplyr::select(ptt, id, date, lon, lat, x, y)

dupes = which(duplicated(x = df[,c('ptt', 'date')]))

df = df[-dupes,]

df = df[df$date < as.POSIXct('2023-06-01 00:00:00', tz = 'UTC'), ]

#PROJECT TO Mercator 
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

rm(list = setdiff(ls(), c('dfsf', 'grid_spdf', 'projAEA', 'uid', 'volras', 'hires_land')))
gc()

## Creates an object of class Itraj
## UD per seal
ud_raster <- list()

for(i in 1:length(uid)){
  
  ## Tracking data
  ID = uid[i]
  seal.data <- df[df$ptt == ID, ]
  seal.data = seal.data[order(seal.data$date), ]
  xy = data.frame(x = seal.data$x, y = seal.data$y)
  DateTime = seal.data$date
  id = seal.data$ptt
  #burst = seal.data$id
  data.ltraj <- as.ltraj(xy = xy, date = DateTime, id = id, proj4string = CRS(projMERC))
  
  ## Parameters for BRB
  tmax = 2*3600
  Lmin = 100 # smallest distance below which we would consider the animal not to be moving
  hmin = 20000 # km - the 99% of the standard error estimates of the SSM
  vv <- BRB.likD(data.ltraj, Tmax=tmax, Lmin=Lmin)
  ud = BRB(ltr = data.ltraj, D = vv, Tmax = tmax, Lmin = Lmin, hmin = hmin, grid = grid_spdf, type = 'UD')
  
  
  ## Convert the UD to raster
  r = raster::raster(ud)
  ud_raster[[i]] <- r
  names(ud_raster)[[i]] = ID
}

# Convert to stars 
udstars = lapply(ud_raster, st_as_stars)


############## PARTITION THE ANALYSIS BY SEX ###################################
# there is evidence that males and females have different dispersal patterns, with males tending to disperse more widely
# Should probably run these analyses partitioning out males and females 
meta = read_csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv")

males = meta$ptt[which(meta$sex == 'M')]
females = meta$ptt[which(meta$sex == 'F')]

Fstars = udstars[which(names(udstars) %in% females)]
Mstars = udstars[which(names(udstars) %in% males)]

# Calculate the overlap 
R = 3000 # approx # of individuals * 100
overlapF <- boot_overlap(Fstars, R = 3000, method = "PHR", percent = 100)
save(list = 'overlapF', file = './data/L3/Hg-Pre-Construction-OverlapProbability100PVC_BRB1kmres_Female.RData')

# Calculate the overlap 
overlapM <- boot_overlap(Mstars, R = 3000, method = "PHR", percent = 100)
save(list = 'overlapM', file = './data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Male.RData')


####################### TOTAL ###################################
# Calculate the overlap - takes about 11 hours
overlap <- boot_overlap(udstars, R = 6000, method = "PHR", percent = 100)
save('overlap', file = './data/L3/HgPre-Construction_OverlapProbability100PVC_BRB1kmres-ALL.RData')

########################################################################################################################
################################# SCRATCH ##############################################################################
################################# SEASONAL BREAK OUT ###################################################################
# FWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_FemaleUDs_Winter.grd')
# MWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_MaleUDs_Winter.grd')
# 
# WUDs = c(FWUDs, MWUDs)
# 
# WUDstar = lapply(X = WUDs, FUN = function(x){
#   # Convert the layer to a stars object
#   stars_obj <- st_as_stars(x)
#   return(stars_obj)
# })
# 
# Woverlap = boot_overlap(WUDstar, R = 5000, method = "PHR", percent = 100)
# save(list = 'Woverlap', file = './data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Winter.RData')
# 
# ################################# SPRING
# FWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_FemaleUDs_Spring.grd')
# MWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_MaleUDs_Spring.grd')
# 
# WUDs = c(FWUDs, MWUDs)
# 
# WUDstar = lapply(X = WUDs, FUN = function(x){
#   # Convert the layer to a stars object
#   stars_obj <- st_as_stars(x)
#   return(stars_obj)
# })
# 
# Woverlap = boot_overlap(WUDstar, R = 5000, method = "PHR", percent = 100)
# save(list = 'Woverlap', file = './data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Spring.RData')
# 
# ####################### Summer
# FWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_FemaleUDs_Summer.grd')
# MWUDs = terra::rast(x = './data/L3/UD/Hg_PreConstruction_MaleUDs_Summer.grd')
# 
# WUDs = c(FWUDs, MWUDs)
# 
# WUDstar = lapply(X = WUDs, FUN = function(x){
#   # Convert the layer to a stars object
#   stars_obj <- st_as_stars(x)
#   return(stars_obj)
# })
# 
# Woverlap = boot_overlap(WUDstar, R = 5000, method = "PHR", percent = 100)
# save(list = 'Woverlap', file = './data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Summer.RData')




