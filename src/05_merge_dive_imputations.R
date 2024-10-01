library(readr)
library(dplyr)
library(tidyr)
setwd("~/seal_telemetry")

path = paste(getwd(), '/data/L3/dive/imputeddivepos/', sep = '')

fnames = list.files(path = path,
                    pattern = '*.csv$', all.files = TRUE,
                    recursive = TRUE, full.names = TRUE)

data = purrr::map_dfr(.x = fnames, .f = function(.x){

  tmp = read_csv(.x) %>% group_by(diveID) %>%
    mutate(Nimputes = n()) %>%
    rename(bathy_m = bathydepth)
  return(tmp)

})

data = data %>%
  mutate(totaldivespossible = n_dives_lost + n_dives_imputed) %>%
  mutate(percloss = n_dives_lost / totaldivespossible)


datmerge = data[data$bathy_m < 0, ]
View(head(datmerge))
write_csv(x = datmerge, file = './data/L3/dive/Hg_2019-2023_DivePosEstimates_Imputations_BestEst.csv')


###################################################################
################### SCRATCH #######################################
###################################################################
# Merge in original 'best' estimates
# og = read_csv("./data/L2/SSM/Hg-2019-2023-SSM-DivePositions.csv")
# og$imputeiter = 0
# og$id = og$SegID
# og$date = og$mid_dive_dt
# og$Depth = -og$Depth
# 
# bathy = marmap::getNOAA.bathy(lon1 = -78, lon2 = -59, lat1 = 37, lat2 = 48, resolution = 0.004)
# bathy = marmap::as.raster(bathy)
# 
# ogsf = og %>% sf::st_as_sf(coords = c('lon', 'lat')) %>% sf::st_set_crs(4326)
# 
# ogsf$bathy_m = raster::extract(x = bathy, y = ogsf)
# og = ogsf %>% mutate(lat = sf::st_coordinates(.)[,2], lon = sf::st_coordinates(.)[,1]) %>% sf::st_drop_geometry()
# 
# og = og[,c('id', 'diveID', 'imputeiter', 'date', 'lon', 'lat', 'x', 'y', 'x.se', 'y.se', 'IDI', 'DiveDur', 'Depth', 'bathy_m')]
# 
# diveIDs = unique(data$diveID)
# og = og[og$diveID %in% diveIDs, ]
# og = og[!is.na(og$bathy_m),]
#datmerge = bind_rows(og, data)