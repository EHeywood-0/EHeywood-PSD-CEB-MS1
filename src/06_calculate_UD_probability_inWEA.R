# GET PROBABILITY OF OVERLAP WITH WEA 

library(raster)
library(sf)
library(readr)

raster01 = function(r){
  
  # get the min max values
  minmax_r = range(values(r), na.rm=TRUE) 
  
  # rescale 
  r01 = (r-minmax_r[1]) / (diff(minmax_r))
  return(r01)
}


seasons = c('Winter', "Spring", 'Summer', 'Fall')

prj = crs(x = raster::stack(x = paste0('./data/L3/UD/UD_rasters/Hg_PreConstruction_FemaleUDs_', seasons[1],'.grd')))

# read in wind energy lease areas - subsetted, unioned, and dissolved for those in the northeast
wea = st_read("./data/shapefiles/Wind_Lease_Boundaries__BOEM__2025-01-01/Offshore_Wind_Lease_Outlines_Northeast_Dissolve.shp") %>%
  st_transform(prj)


collectL = list()
# GET AREAS FOR FEMALES
for (season in 1:length(seasons)){
  
  FStackR = raster::stack(x = paste0('./data/L3/UD/UD_rasters/Hg_PreConstruction_FemaleUDs_', seasons[season],'.grd'))
  MStackR = raster::stack(x = paste0('./data/L3/UD/UD_rasters/Hg_PreConstruction_MaleUDs_', seasons[season],'.grd'))
  
  # # Transform probabilites so they sum to one 
  # FStackR_T = raster::stack(lapply(X = 1:nlayers(FStackR), FUN = function(x){
  #   return(raster01(FStackR[[x]]))
  # }))
  
  TotF = cellStats(x = FStackR, stat = 'sum', na.rm=T)
  TotM = cellStats(x = MStackR, stat = 'sum', na.rm=T)
  
  # Mask the raster stacks by the weas
  Fmask = raster::mask(x = FStackR, mask = wea)
  Mmask = raster::mask(x = MStackR, mask = wea)
  
  
  # Sum the values within the masked area 
  Fwea_sum = cellStats(Fmask, stat = 'sum', na.rm = T)
  
  Mwea_sum = cellStats(Mmask, stat = 'sum', na.rm = T)
  

  
  Fres = data.frame(Sex = 'Female', 
                    TripID = gsub(pattern = 'X', replacement = '', gsub(pattern = '\\.', '-', names(FStackR))),
                    Season = seasons[season], 
                    Probability_in_Wea = Fwea_sum/TotF)
  
  Mres = data.frame(Sex = 'Male', 
                    TripID = gsub(pattern = 'X', replacement = '', gsub(pattern = '\\.', '-', names(MStackR))),
                    Season = seasons[season], 
                    Probability_in_Wea = Mwea_sum/TotM)
  allres = bind_rows(Fres, Mres)
  
  collectL[[season]] = allres
  
}

final = bind_rows(collectL)

write_csv(final, './data/L4/Hg-PreConstruction-TripUDs-ProbWEA.csv')
