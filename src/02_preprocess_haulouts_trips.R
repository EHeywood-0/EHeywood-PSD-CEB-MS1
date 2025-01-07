# Deploy Haulout Detection Methods in HauloutDetection.Rmd on entire seal telemetry dataset/
# EIH
# updated: 2024-09-27
# inputs: fxns_helper.R
#         Hg_2019-2023_Prefiltered_Locs.csv
#         Hg_2019-2023_PercentDry+Behavior_HauloutMethods_Merged.csv
#         Hg_2019-2023_PercentDry_HauloutAtSea_Events.csv

# set up
rm(list = ls())
gc()
setwd("~/seal_telemetry/")

# Data manipulation and document knitting
library(readr)
library(dplyr)
library(knitr)
library(tidyr)

# Spatial libs
library(sf)
library(stars)

# Modeling
library(aniMotum)
library(marmap)
# Library for US and Canada Land Polygons and Minor Island Polygons
library(rnaturalearthhires)

# Define haulout duration threshold in hours (Nowak et al. 2023) - 33.6 +/- 29 hours
ho_dur_max = 24*7

distthreshold = 7500 # 7.5 km - if the mean fitted / imputed distance from shore falls outside this threshold it will be classified as not hauled out


################################################################################################################################
#################################################### LOAD DATA #################################################################
# Load pre-filtered location data
locs = read_csv(file = './data/L1/locs/Hg_2019-2023_Prefiltered_Locs.csv') %>% 
  dplyr::select(id, datetime, lc, lat, lon, error.radius, smaj, smin, eor)

target_ptt = unique(locs$id)

# READ IN PERCENT DRY HAULOUT AND HAULOUT BEHAVIOR PROCESSED DATA WITH OVERLAP REMOVED
haulouts = read_csv("./data/L1/haulout/Hg_2019-2023_PercentDry+Behavior_HauloutMethods_Merged.csv") %>%
  arrange(ptt, startho)

atsea = read_csv("./data/L1/haulout/Hg_2019-2023_PercentDry_HauloutAtSea_Events.csv") %>% 
  filter(What == 'AtSea') %>% arrange(id, Start)
################################################################################################################################
### ASSIGN HAULOUT PERIODS #### 

# STEP 1 ##############################################################################
#  Use the haulout periods in the haulout methods merged csv to assign locations as hauled out
# This csv contains haulouts identified from both beh haulout csvs and the histos percent dry timeline
locsL = split(locs, locs$id)

# Iterate over each ptt id to assign haulout behavior column
storageL = list()

for(i in 1:length(locsL)){
  tmp_loc = locsL[[i]]
  pttid = unique(tmp_loc$id)
  
  HO = haulouts[haulouts$ptt == pttid,]
  AS = atsea[atsea$id == pttid, ]
  
  tmp_loc$haulout = NA
  tmp_loc$hotype = NA
  # Eliminate at sea time periods that occur while the BEH csv said haulout
  badASidxStart = unlist(sapply(X = AS$Start, FUN = function(x){
    if (any(x > HO$startho & x < HO$endho)){
      return(which(AS$Start == x))
    }else{return(NULL)}
  }, simplify = TRUE))
  badASidxend = unlist(sapply(X = AS$End, FUN = function(x){
    if (any(x < HO$endho & x > HO$startho)){
      return(which(AS$Start == x))
    }else{return(NULL)}
  }, simplify = TRUE))
  
  badtot = c(badASidxStart, badASidxend)
  if (length(badtot)>0){
    AS = AS[-badtot, ]
    }
  
  if (nrow(HO)>0){
    
    for (j in 1:nrow(HO)){
      s = HO$startho[j]
      e = HO$endho[j]
      hoidx = which(tmp_loc$datetime >= s & tmp_loc$datetime <= e)
      tmp_loc$haulout[hoidx] = TRUE
      tmp_loc$hotype[hoidx] = HO$type[j]
      } # end for loop
    }else{tmp_loc$haulout = NA}
  
  if (nrow(AS)>0){
    
    for (k in 1:nrow(AS)){
      s = AS$Start[k]
      e = AS$End[k]
      atseaidx = which(tmp_loc$datetime >= s & tmp_loc$datetime <= e)
      tmp_loc$haulout[atseaidx] = FALSE
      tmp_loc$hotype[atseaidx] = 'PercentDry AtSea'
    }
  }else{tmp_loc}
  storageL[[i]] = tmp_loc
}

data = bind_rows(storageL) %>% arrange(id, datetime)

rm(list = setdiff(ls(), c('data', 'distthreshold', 'ho_dur_max')))
########################################################################################################
############################################################################################################################################################################
###### STEP 2 and 2A ######
### 2. As a second pass, model the data to account for location error and then intersect track points with a land shapefile containing minor islands (islands < 2 sq km.)
### 2A. "On land" bouts with > 4 consecutive points get assigned as haulout behavior. Use these to supplement the incomplete haulout record CSVs. 

# Load the land spatial polygons dataframe and subset based on relevant landmasses
land = rnaturalearthhires::states10
land = land[land$admin %in% c("Canada", 'United States of America'),]
land = land[which(land$postal %in% c("ME", "NH", "VT", "RI", "CT", "MA", "PA", "NY", "DE", "NJ", "VA", "MD", "NC", 'OH', 'WV', 'SC', 'KY', 'TN', 'IL', 'WI', 'MN', 'IA', 'IN', 'MI') | land$region == "Eastern Canada"),]
#plot(land['featurecla'])

# Now minor islands (islands less than 2 square km)
minor_islands = st_read("./READ-PSB-MoveSeals/ne_10m_minor_islands/ne_10m_minor_islands.shp") %>%
  mutate(name = featurecla) %>%
  dplyr::select(name, geometry)

#get the bounding box of the land and subset minor islands based on this bounding box
bounding_box <- st_bbox(land) %>% st_as_sfc()
minor_islands_subset = st_intersection(minor_islands, bounding_box)

# Combine the land and minor islands
hires_land = bind_rows(land, minor_islands_subset)

#########################################################################################################################################################
##################################################################### FIT SSM ########################################################################### 

# Format anilocs dataframe for animotum expected columns and column names
anilocs = data %>% 
  dplyr::select(id, datetime, lc, lon, lat, smaj, smin, eor) %>% arrange(id, datetime) %>% rename(date = datetime)

anilocs$id = as.character(anilocs$id)

ids = unique(data$id)


# Source EIH modified sample path function from J Hatch
source("./READ-PSB-MoveSeals/src/fxns_helper.R")

# Set seed (for reproducibility) and fit
set.seed(29)
fit = fit_ssm(x = anilocs, 
              vmax = 3.5, 
              time.step = data.frame(id = anilocs$id, date = anilocs$date), 
              model = "crw")


# Look at standard error estimations that were not calculated
idsnose = sapply(X = 1:nrow(fit), FUN = function(x){
  tmppar = as.data.frame(fit$ssm[[x]]$par)
  
  if (any(is.na(tmppar$`Std. Error`))){
    idforrefit = fit$id[x]
  }else{idforrefit = NA}
  
  return(idforrefit)
  
  
  
})

idsnose = idsnose[!is.na(idsnose)]

# GREAT NO CONVERGENCE ISSUES or SE estimation issues! 
nonconvergeids = fit$id[fit$converged == FALSE]

idsforrefit = unique(c(idsnose, nonconvergeids))


# EIH Modified Josh Hatch's 'sample_paths' to sample fitted paths
pathsims = sample_paths(x = fit, n = 100)

# Get the predictions from the fitted model (the predicted locations at the original fitted date time stamps)
predicted = grab(fit, what = "predicted")
# give the best fit data a sample number of 0
predicted$sample = 0

# remove all velocity related columns and their associated standard errors
predicted = predicted %>% dplyr::select(-c(u, v, u.se, v.se, s, s.se))

pathsims = pathsims %>% st_transform(4326) %>% 
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  rename(x = sample_x, y = sample_y) %>% dplyr::select(id, date, lon, lat, x, y, sample)

# bind them together
allsims = bind_rows(predicted, pathsims)


# Convert fitted locs to a sf object
sim_sf = st_as_sf(allsims, coords = c('lon', 'lat')) %>% st_set_crs(value = 4326)

rm(pathsims)


# Extract bathy data at 15 arcsecond resolution
bathy = marmap::getNOAA.bathy(lon1 = min(data$lon)-0.5, lon2 = max(data$lon)+0.5,
                              lat1 = min(data$lat)-0.5, lat2 = max(data$lat)+0.5, resolution = 0.25)
bathy_raster = marmap::as.raster(bathy) %>% st_as_stars()

bathy_depths = st_extract(bathy_raster, sim_sf)

allsims$bathydepth = bathy_depths$layer

# Get mean estimates for bathymetry all sims...
samplebathy = allsims %>% dplyr::select(id, date, sample, bathydepth) %>% pivot_wider(names_from = 'sample', values_from = 'bathydepth', names_prefix = 'sample-') 
bathymeans = apply(samplebathy[,-c(1,2)], MARGIN = 1, FUN = mean)
samplebathy$meanbathy_m = bathymeans

indices_above5m = which(samplebathy$meanbathy_m >= -5)

samplebathy$HauloutBathy = FALSE
samplebathy$HauloutBathy[indices_above5m] = TRUE
samplebathy = samplebathy %>% dplyr::select(id, date, meanbathy_m, HauloutBathy)

# NOW DISTANCE 
hires_land_union = st_union(hires_land)

# Calculate great circle distance for each point (to closest)
dist = st_distance(sim_sf, hires_land_union, which = 'Great Circle')
allsims$dist_m = as.numeric(dist)
sampledist = allsims %>% dplyr::select(id, date, sample, dist_m) %>% pivot_wider(names_from = 'sample', values_from = 'dist_m', names_prefix = 'sample-') 
distmeans = apply(sampledist[,-c(1,2)], MARGIN = 1, FUN = mean)

sampledist$meandist_m = distmeans

distthreshold = 7500
indices_further_10km = which(sampledist$meandist_m >= distthreshold)
sampledist$TooFar = FALSE
sampledist$TooFar[indices_further_10km] = TRUE

# Merge back
predicted_ = left_join(predicted, samplebathy,by = c("id","date")) %>% dplyr::select(c(id, date, meanbathy_m, HauloutBathy))
predicted_ = left_join(predicted_, sampledist, by = c('id', 'date')) %>% dplyr::select(c(id, date, meanbathy_m, HauloutBathy, meandist_m, TooFar))

data$id = as.character(data$id)
# Merge the fitted with the filt_locs
data_ = left_join(data, predicted_, join_by(id, datetime == date))

write.csv(x = data_, file = './data/L1/locs/HG_2019-2023_PreTripHauloutClassification.csv', row.names = F)

######################################### ASSIGN FINAL HAULOUT DETERMINATION ####################################
#################################################################################################################

V1 = data_
# HAULOUT ASSIGNMENT LOGIC
### if any of the three lines of evidence say hauled out, then assign hauled out...
V1$final_haulout = FALSE

V1$final_haulout[V1$haulout == TRUE] = TRUE
V1$final_haulout[V1$HauloutBathy == TRUE] = TRUE

#test177036 = V1[V1$id == '177036',]

# If either percent dry or bathy says it is hauled out but it is too far from land - say it isn't hauled out
V1$final_haulout[which(V1$HauloutBathy == TRUE & V1$TooFar == TRUE)] = FALSE
V1$final_haulout[which(V1$hotype == 'percent dry' & V1$TooFar == TRUE)] = FALSE

# 177036 has issues we solve here - this one clearly hauled out but percent dry faulted
V1$final_haulout[which(V1$id == 177036 & V1$HauloutBathy == TRUE)] = TRUE

# RECLASS HOTYPE
V1$hotype[which(V1$final_haulout & is.na(V1$hotype))] = paste(V1$hotype[which(V1$final_haulout & is.na(V1$hotype))], 'Spatial', sep = '-')


############################################################################################################################################################################
####################################################### TRIP/HAULOUT ID ASSIGNMENT ##############################################################################################
data_V2 = V1

# Assign a haulout ID for each haulout 'bout'
trips = data_V2 %>% group_by(id) %>% arrange(datetime) %>%
  mutate(HauloutCumSum = cumsum(ifelse(!final_haulout, 1, 0))) %>%
  arrange(desc(id), datetime) %>% 
  mutate(HauloutID = ifelse(!final_haulout, NA, paste(id, HauloutCumSum, sep = "-"))) %>%
  ungroup() %>%
  group_by(id, HauloutID) %>%
  # get the mean lat and long of each haulout event ID
  mutate(MeanHauloutLat = mean(lat, na.rm=TRUE), MeanHauloutLon = mean(lon, na.rm = TRUE),
         Haulout_Start = min(datetime), Haulout_End = max(datetime),
         Nhauloutlocs = n()) %>%
  ungroup() %>%
  # get trip id
  group_by(id) %>%
  mutate(TripID = cumsum(ifelse(!final_haulout, 0, 1))) %>%
  arrange(desc(id), datetime) %>% 
  mutate(TripID = ifelse(!final_haulout, paste(id, TripID, sep = "-"), NA)) %>%
  ungroup()


# THIS JUST GETS NEW TRIP AND HAULOUT IDS THAT START AT 1 for each PTT
tripIDS = trips %>% group_by(id) %>% distinct(TripID) %>% filter(!is.na(TripID)) %>%
  mutate(NewTripID = paste(id, row_number(), sep = "-")) %>% ungroup() %>% dplyr::select(-c(id))

hauloutIDs = trips %>% group_by(id) %>% distinct(HauloutID) %>% filter(!is.na(HauloutID)) %>%
  mutate(NewHauloutID = paste(id, row_number(), sep = "-")) %>% ungroup() %>% dplyr::select(-c(id))

tmp = merge(trips, tripIDS, "TripID", all.x = TRUE) %>% arrange(desc(id), datetime)

newtrips = merge(tmp, hauloutIDs, "HauloutID", all.x = TRUE) %>% arrange(desc(id), datetime) %>%
  mutate(TripID = NewTripID, HauloutID = NewHauloutID) %>% dplyr::select(-c(NewTripID, NewHauloutID))


###############################################################################################################################################################
###################################################################### SEGMENT TRIPS ##########################################################################
# Assign Segments
# Use Theo Michelot's split at gap function in helper_functions.R
# Rename TripID to ID as an ID column is expected for the functiona and we want to segment the trip data
tripsseg = newtrips %>% arrange(id, datetime) %>% rename(ID = TripID, time = datetime)

# Get maximum gap time in minutes
maxgap = 24*60

tripsseg = split_at_gap(data = tripsseg, max_gap = maxgap)

# SELECT SOME FINAL COLUMNS TO KEEP 
tripsforwrite = tripsseg %>% rename(SegID = ID, TripID = ID_old, datetime = time, MeanImputedDistFromShore_M = meandist_m, MeanImputedBathy_m = meanbathy_m) %>%
  dplyr::select(id, TripID, SegID, HauloutID, datetime, lc, lon, lat, error.radius, smaj, 
                smin, eor, haulout, hotype, HauloutBathy, final_haulout,
                MeanImputedBathy_m, MeanImputedDistFromShore_M)




write_csv(tripsforwrite, file = "./data/L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv")





