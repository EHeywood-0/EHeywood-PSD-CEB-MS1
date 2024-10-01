# EIH
# Rejection Sampling Sensitivity Analysis:
# Modified: 2024-09-30
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

rm(list = setdiff(ls(), c('bathymetry', 'impute_bathy', 'rmvnorm_prec', 'sample_path_reject', 'dives')))


# Load SSM Objects + Random Effects
load(file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_paramobjects.RData')

##################################################################################################################
################################################ TESTING #########################################################
##################################################################################################################
# look at the stability and precision of the N = 2-100 imputed locations and their extracted bathymetry
#`225840-6-1`

# one indivdiual trip with progress bar (ndives = 426)
timehighN = system.time({

progressr::with_progress({
 result = sample_path_reject(x = parres[[1]],
                         bathymetry = bathymetry,
                          dives = dives,
                          n = 100,
                          tolerance = 3,
                          stableiter = 100,
                          progress = TRUE)
 })
 })

# summarise the results of the imputation
resultsum = result %>% st_drop_geometry() %>% group_by(id, diveID) %>%
   summarise(MeanBathy = mean(bathydepth, na.rm=T), SDBathy = sd(bathydepth, na.rm=T),
             divedepth = unique(Depth), N = n(), meanx = mean(sample_x), sdx = sd(sample_x))

# get the dive ids that have greater than 190 imputations
ids = unique(resultsum$diveID[resultsum$N > 190])

resfilt = result %>% filter(diveID %in% ids) %>% st_drop_geometry()

# # # Initialize an empty list to store the results
boot_list <- list()

library(data.table)
resfilt = data.table(resfilt)
write_csv(x = resfilt, "./data/L3/dive/sensitivity_testing_imputed_dive_pos_200.csv")
# Loop over values of N from 10 to 100 by 2
set.seed(2024-05-31) # anchor the RNG again, as we use `sample` again
Ns = seq(10, 200, 2)

for (i in 1:length(Ns)) {
  N = Ns[i]
  
  # create 200 bootstrapped datasets for each N
  dt_sample <- replicate(200, resfilt[, .SD[sample(x = .N, size = N, replace = TRUE)], by = diveID], simplify = FALSE)

  dt_sample = bind_rows(dt_sample, .id = 'replicate')

  dt_sample$N = N
  
  # calculate the mean, sd of x, y coords, and depth of each replicate at each level of N for each diveID
  summ = dt_sample %>% group_by(N, replicate, diveID) %>%
    summarise(ux = mean(sample_x),
              uy = mean(sample_y),
              sdx = sd(sample_x),
              sdy = sd(sample_y),
              ubathy = mean(bathydepth),
              sdbathy = sd(bathydepth),
              se_ux = sd(sample_x) / sqrt(n()),
              se_uy = sd(sample_y) / sqrt(n()),
              se_ub = sd(bathydepth) / sqrt(n()))


  boot_list[[i]] = summ


}

# row bind the bootstrapped results
bootres = bind_rows(boot_list)
write_csv(x = bootres, file = './data/L3/dive/sensitivity_testing_bootstrapped_data.csv')


