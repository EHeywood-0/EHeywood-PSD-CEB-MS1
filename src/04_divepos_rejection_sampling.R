### DIVE LOCATION REJECTION SAMPLING 
# EIH
# 2024-06-05

######################################################### SETUP #########################################################################
rm(list = ls())
gc()

library(dplyr)
library(readr)
library(tidyr)
library(sf)
library(stars)
library(aniMotum)

#install.packages("aniMotum", repos = c("https://ianjonsen.r-universe.dev", "https://cloud.r-project.org"))

setwd('~/seal_telemetry/')

# Source necessary functions for rejection sampling
source('./READ-PSB-MoveSeals/src/fxns_rejection_sampling.R') 
rm(list = setdiff(ls(), c('rmvnorm_prec', 'sample_path_reject', 'impute_bathy')))
######################################################### LOAD VARS #####################################################################
data = read_csv("./data/L2/SSM/Hg-2019-2023-SSM-DivePositions.csv")
dives = data %>% dplyr::select(SegID, diveID, date, DiveDur, IDI, Depth)
data_sf = st_as_sf(data, coords = c('lon', 'lat')) %>% st_set_crs(4326)

# Get bathymetry data at 15 arc second resolution (~350m @ 40 degrees latitude)
# resolution of the grid, in minutes (default is 4)
# Minimum resolution of the ETOPO 2022 dataset is 15 arc-second == 0.25 minutes
coords = st_coordinates(data_sf)
bathy = marmap::getNOAA.bathy(lon1 = min(coords[,1])-2,
                              lon2 = max(coords[,1])+2,
                              lat1 = min(coords[,2])-2,
                              lat2 = max(coords[,2])+2,
                              resolution = 0.25) # 15 arc second resolution


# Convert to raster and then to stars object
bathymetry = marmap::as.raster(bathy) %>% st_as_stars()

# Load SSM Objects + Random Effects
load(file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_paramobjects.RData')

# Already done
fdone = list.files('./data/L3/dive/imputeddivepos/', pattern = '*.csv', full.names = F)
finishedIDs = sapply(strsplit(fdone, '_'), '[[', 1)

fits = parres[which(!names(parres) %in% finishedIDs)]

######################################################### REJECTION SAMPLING ################################################################
library(future)
library(furrr)


niter = 100
tol <- 3
stablei <- 100

#Detect the number of logical cores
num_cores <- parallelly::availableCores(logical = TRUE)
num_workers = round(0.85*num_cores, 0)

# Set up the future plan for parallel execution
future::plan(multisession, workers = num_workers)  # or any other plan like multicore, cluster, etc.

result2 = furrr::future_map_dfr(.x = seq_along(fits), .options = furrr_options(seed=TRUE),.f = function(.x){
  
  
  res = sample_path_reject(x = fits[[.x]],
                           bathymetry = bathymetry,
                           dives = dives,
                           n = niter,
                           tolerance = tol,
                           stableiter = stablei,
                           progress = FALSE)
  
  numdives = nrow(fits[[.x]]$predicted)
  res2 = res %>% st_drop_geometry() %>% mutate(x = sample_x, y = sample_y) %>%
    st_as_sf(coords = c('sample_x', 'sample_y')) %>%
    st_set_crs(st_crs(fits[[.x]]$predicted)) %>% st_transform(4326) %>%
    
    mutate(lat = st_coordinates(.)[,2], 
           lon = st_coordinates(.)[,1], 
           imputeiter = sample) %>% 
    st_drop_geometry() %>%
    group_by(sample) %>% 
    mutate(n_dives_imputed = n(), n_dives_lost = numdives - n()) %>%
    ungroup() %>%
    dplyr::select(id, diveID, imputeiter, n_dives_imputed, n_dives_lost,
                  date, lon, lat, x, y, 
                  IDI, DiveDur, Depth, bathydepth) %>%
    arrange(date, diveID)
  write_csv(x = res2, file = sprintf('./data/L3/dive/imputeddivepos/%s_imputeddivepos.csv', fits[[.x]]$predicted$id[1]))
  
  
  return(res)
    
  })


