# plot trip and haulout locations in leaflet for validation
# Purpose: A visual gutcheck of the trip / haulout event segmentation process
# EIH
# inputs: HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv
# updated: 2024-09-27



# EIH
# 2024-03-29

# Dynamic mapping
library(leaflet)
library(leaflegend)
library(htmlwidgets)
# Spatial libraries
library(sf)

library(readr)
library(dplyr)
library(sftime)

setwd('~/seal_telemetry/data/')
#################################################################################################################################
############################################### PLOTTING HAULOUT LOCS ###########################################################

df = read_csv("./L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv")
df$Year = format(df$datetime, "%Y")


# plot haulout locations to check for reasonableness
haulout_locs = df[df$final_haulout == TRUE, ]

haulout_locs$hotype[haulout_locs$hotype == 'PercentDry AtSea'] = 'Spatial'
haulout_locs$hotype[haulout_locs$hotype == 'NA-Spatial'] = 'Spatial'

m = leaflet(haulout_locs) %>%
  setView(lng = mean(haulout_locs$lon), lat = mean(haulout_locs$lat), zoom = 5) %>%
  addProviderTiles(providers$Esri.WorldImagery)

pal <- colorFactor(palette = c("orange", "green", 'blue', 'purple'), levels = unique(haulout_locs$hotype))

for (i in 1:length(unique(df$Year))){
  year = unique(haulout_locs$Year)[i]
  tmp = haulout_locs[haulout_locs$Year == year,]
  
  m = m %>% 
    addCircleMarkers(data = tmp,
                     lng = ~lon, 
                     lat = ~lat, 
                     color = ~pal(hotype), 
                     stroke = TRUE,
                     fillOpacity = 0.5, 
                     radius = 5,
                     group = as.character(unique(tmp$Year))) %>%
    addLegend("bottomleft", pal = pal, 
              title=as.character(unique(tmp$Year)),
              values = haulout_locs$hotype,
              group = as.character(unique(tmp$Year)))
}

# Add layers control
m = m %>% # Layers control
  
  addLayersControl(overlayGroups = c(unique(haulout_locs$Year)),
                   position = c("bottomright"),
                   options = layersControlOptions(collapsed = FALSE)
  )
# Save out plot
saveWidget(m, file="../plots/ssm_validation/Hg-2019-2023_classified_hauloutlocs_visual_validation_map.html")


############################################### PLOTTING TRIP LOCS ###########################################################
tr = df[!df$final_haulout, ]
ptts = unique(tr$id) # get ptts


for (ptt in ptts){
  
  # Get this ptt animal
  tptt = tr[tr$id == ptt, ]
  
  # We're only going to plot the trips with >= 20 locs - group by segID (ID) and filter
#  t = tptt %>% group_by(SegID) %>% mutate(N = n()) %>% filter(N >= 50)
  
  sft_lines = tptt %>% dplyr::select(SegID, lon, lat, datetime) %>% arrange(SegID, datetime) %>%
    st_as_sftime(time_column_name = 'datetime', coords = c('lon', 'lat')) %>% 
    st_set_crs(4326) %>% group_by(SegID) %>%
    dplyr::summarize(do_union=FALSE) %>%  # do_union=FALSE doesn't work as well
    st_cast("LINESTRING") 
  palL = colorFactor(palette = 'magma', domain = unique(tptt$SegID))
  t = tptt
  m = leaflet(t) %>%
    setView(lng = mean(t$lon), lat = mean(t$lat), zoom = 5) %>%
    addProviderTiles(providers$Esri.WorldImagery)
  
  segids = unique(t$SegID)
  
  for (id in segids){
    tmp = t[t$SegID == id, ]
    tmp = tmp %>% arrange(datetime)
    l = sft_lines[sft_lines$SegID ==id,]
    tmp$datenum = as.numeric(tmp$datetime)
    pal = colorNumeric(palette = 'viridis', domain = sort(tmp$datenum), alpha = 0.75)
    m = m %>% 
      addCircleMarkers(data = tmp[which(tmp$datenum == min(tmp$datenum)),],
                       lng = ~lon, 
                       lat = ~lat, 
                       stroke = TRUE,
                       color = 'green',
                       label = paste('Start Time', tmp$datetime[which(tmp$datetime == min(tmp$datetime))]),
                       radius = 10,
                       group = as.character(unique(tmp$SegID))) %>%
      addCircleMarkers(data = tmp[which(tmp$datenum == max(tmp$datenum)),],
                       lng = ~lon, 
                       lat = ~lat, 
                       stroke = TRUE,
                       color = 'red',
                       label = paste('End Time', tmp$datetime[which(tmp$datetime == max(tmp$datetime))]),
                       radius = 10,
                       group = as.character(unique(tmp$SegID))) %>%
      
      addCircleMarkers(data = tmp,
                       lng = ~lon, 
                       lat = ~lat, 
                       stroke = TRUE,
                       color = ~pal(tmp$datenum),
                       fillOpacity = 0.5, 
                       label = paste('Average Distance to shore (m):', tmp$MeanImputedDistFromShore_M, 'Mean Bathy', tmp$MeanImputedBathy_m, 'Time: ', tmp$datetime),
                       radius = 5,
                       group = as.character(unique(tmp$SegID))) %>%
      addPolylines(data = l, color = ~palL(unique(tmp$SegID)), group = as.character(unique(tmp$SegID)))
  }
  
  # Add layers control
  m = m %>% # Layers control
    
    addLayersControl(overlayGroups = c(unique(t$SegID)),
                     position = c("bottomright"),
                     options = layersControlOptions(collapsed = TRUE)
    )
  
  # Save out plot
  fn = sprintf(fmt = "../plots/ssm_validation/Hg-%s_classified_trip_visual_validation_map.html", ptt)
  saveWidget(m, file= fn)
  
  
}

###############################################################################################################################################################
