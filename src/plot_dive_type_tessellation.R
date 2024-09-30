# EIH
# 2024-06-24
# Gridded Probability of Benthic Diving .. . . and potentially more
# st_make_grid() - hexagonal grid

setwd('~/seal_telemetry/data/')

library(readr)
library(dplyr)
library(sf)

library(ggplot2)
library(scico)
library(metR)
#library(tidyterra)
#library(ggOceanMaps)
#library(mapdata)
#library(viridis)

####################################### LOAD NECESSARY DATA SOURCES ###################################################
# High resolution land
# Load the land spatial polygons dataframe and subset based on relevant landmasses
land = rnaturalearthhires::states10
land = land[land$admin %in% c("Canada", 'United States of America'),]
land = land[which(land$postal %in% c("ME", "NH", "VT", "RI", "CT", "MA", "PA", "NY", "DE", "NJ", "VA", "MD", "NC", 'OH', 'WV', 'SC', 'KY', 'TN', 'IL', 'WI', 'MN', 'IA', 'IN', 'MI') | land$region == "Eastern Canada"),]
plot(land['featurecla'])

# Now minor islands (islands less than 2 square km)
minor_islands = st_read("../READ-PSB-MoveSeals/ne_10m_minor_islands/ne_10m_minor_islands.shp") %>%
  mutate(name = featurecla) %>%
  dplyr::select(name, geometry)

#get the bounding box of the land and subset minor islands based on this bounding box
bounding_box <- st_bbox(land) %>% st_as_sfc()
minor_islands = st_intersection(minor_islands, bounding_box)

# Combine the land and minor islands
hires_land = bind_rows(land, minor_islands)
rm(minor_islands, land)
gc()

# Read in L3 dive type classification data
# convert to sf object
proj = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'

data = read_csv("./L3/dive/Hg_2019-2023_DiveTypeClassification.csv") %>%
  st_as_sf(coords = c('meanX_km','meanY_km')) %>% st_set_crs(st_crs(proj)) %>%
  st_transform(4326) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])

# Filter to pre-construction
data = data[data$date < as.POSIXct(x = '2023-06-01', tz = 'UTC'), ]


# BATHYMETRY
bathy = marmap::getNOAA.bathy(lon1 = min(data$lon) - 1, lon2 = max(data$lon) + 1, 
                              lat1 = min(data$lat) - 1, lat2 = max(data$lat) + 1, resolution = 1)

bf = marmap::fortify.bathy(bathy)


# WEAs
we = st_read(".//shapefiles/BOEM_Wind_Leases_8_30_2023.shp")
we = st_intersection(we, bounding_box)
plot(we['STATE'])

we = st_union(we)

gc()
####################################################### END DATA LOAD #######################################################

############################################# SEASONALLY ################################################
## GRIDDED HEXAGONALLY ################################################################

data$Month = as.numeric(format(data$date, '%m'))
data$Monb = format(data$date, '%b')

data$Season = NA
data$Season[data$Month %in% c(1,2)] = 'Winter'
data$Season[data$Month %in% c(3,4,5)] = 'Spring'

data$Season[data$Month %in% c(6,7,8)] = 'Summer'
data$Season[data$Month %in% c(9,10,11,12)] = 'Fall'
seasons = unique(data$Season)


box = st_bbox(obj = data)

lat0 = mean(data$lat)
lon0 = mean(data$lon)

# Specify 2 standard parallels at 1/6 latitude range of study from the top and bottom lats

# 1/6 south of northern extreme
lat1 = as.numeric(box[4] - ((box[4] - box[2]) * 0.16666))

# 1/6 north of southern extreme
lat2 = as.numeric(box[2] + ((box[4] - box[2]) * 0.16666))

# Custom albers equal area (need equal area for hexagonal tessellation)
aea_custom = sprintf('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f +ellps=WGS84 +datum=WGS84', 
                       lat1, lat2, lat0, lon0)


# Classify dives
data$BenthDem = ifelse(data$PWCProb_Benthic > 0.5, 1, 0)
data = data %>% filter(meanbathy_m < -20) %>% dplyr::select(sex, ptt, id, Depth, DiveDur, IDI, BenthDem, Season)

# Calculate st_grids for each season 
storeL = list()

for (season in seasons){
  
  # Subset data to season, taking those dives where the mean bathy (mean of imputed data) is deeper than 20m
  dat_sf = data %>% filter(Season == season) %>% 
    st_transform(aea_custom)
  
  # Create a honeycomb grid
  area_honeycomb_grid <- st_make_grid(dat_sf, c(10000, 10000), what = "polygons", square = FALSE)
  
  # Convert the grid to an sf object and add grid ID
  honeycomb_grid_sf <- st_sf(geometry = st_sfc(area_honeycomb_grid)) %>%
    mutate(grid_id = 1:length(area_honeycomb_grid))
  
  # Join the data with the grid to calculate mean PWCProb_Benthic for each grid cell
  dat_with_grid <- st_join(dat_sf, honeycomb_grid_sf, join = st_intersects)
  
  # Summarize mean PWCProb_Benthic for each grid cell
  honeycomb_summary <- dat_with_grid %>%
    group_by(grid_id) %>%
    summarise(
      Season = unique(Season),
      n_dives = n(),
      n_ptt = length(unique(ptt)),
      PropBenthic = sum(BenthDem)/n(),
      mean_DiveDepth = mean(Depth, na.rm = TRUE),
      mean_DiveDur = mean(DiveDur, na.rm=TRUE),
      ProbF = length(which(sex == 'F')) / n()
    ) %>% st_drop_geometry()
  
  # Join the summary back to the grid for plotting
  honeycomb_count <- left_join(honeycomb_grid_sf, honeycomb_summary, by = 'grid_id') %>%
    filter(n_dives > 1)
  
  # Plot the number of dives and mean PWCProb_Benthic
  #plot(honeycomb_count["n_dives"], main = "Number of Dives")
  #plot(honeycomb_count["mean_PWCProb_Benthic"], main = "Mean Probability of Benthic Dive")
  #plot(honeycomb_count["n_ptt"], main = "Numbr of Individuals")
  #plot(honeycomb_count["mean_DiveDur"][honeycomb_count$mean_DiveDur < 300,], main = "Mean Dive Duration")
  #plot(honeycomb_count["ProbF"], main = "Proportion of Dives - Female")
  
  storeL[[season]] = honeycomb_count
  
}

# Bind rows
res = bind_rows(storeL) %>% filter(n_dives > 2, Season != "Fall") 

# spatially transform to wgs84 for plotting in ggplot2
resl = res %>% st_transform(4326)

resl$Season = factor(resl$Season, levels = c("Winter", "Spring", "Summer"))
#resl$What = factor(resl$What, levels = c('mean_PWCProb_Benthic', 'n_dives'))
# Color palette
# Extract the 'oleron' color palette
oleron_palette <- scico(256, palette = "oleron")

# Split the palette into two halves
half_length <- length(oleron_palette) / 2
first_half <- oleron_palette[1:half_length]
first_half = first_half[1:115]

second_half <- colorRampPalette(c("#d6c2a1", "#cbb493", "#bfa785", "#b49a77", 
                                  "#a98d6a", "#9e815d", "#937551", "#886946", 
                                  "#7e5e3c", "#735332", "#694929", "#5f3f21", 
                                  "#553518", "#4b2b11", "#41210a", "#371804"))(115)
#second_half <- oleron_palette[(half_length + 1):length(oleron_palette)]
#second_half = second_half[1:115]

# Reverse the second half
#reversed_second_half <- rev(second_half)

# Combine the first half with the reversed second half
custom_palette <- c(first_half, second_half)


# Probability Benthic by Season
# g = ggplot(data = resl) +
#   
#   # Add 200m contour
#   geom_contour(data = bf, 
#                aes(x=x, y=y, z=z),
#                breaks=c(-200),
#                linewidth=c(0.3),
#                colour="gray10", lty = 1)+
#   
#   # Add raster
#   geom_sf(data = resl, mapping = aes(fill = PropBenthic, color = PropBenthic), alpha = 0.75) +
#    
#   #scale_fill_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
#   #scale_color_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
#   scale_fill_gradientn(colors = custom_palette, 
#                        guide = guide_colorbar(theme = theme(legend.title.position = 'top',
#                                                             legend.key.width = unit(2.5, 'in')), 
#                                               title = 'Proportion Benthic / Demersal Diving', direction = 'horizontal')) +
#   scale_color_gradientn(colors = custom_palette, ) +
#   # Add depth contours
#   geom_contour(data = bf, 
#                aes(x=x, y=y, z=z),
#                breaks=c(-50),
#                linewidth=c(0.3),
#                colour="gray1", lty = 3)+
#   
#   
#   # Add wind energy areas
#   geom_sf(data = we, color = "orangered3", fill = NA, linewidth = 0.5) +
#   
#   # Add land
#   geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
#   
#   
#   coord_sf(xlim = c(-75, -63.5), ylim = c(38.75, 45)) +
#   
#   theme_bw() +
#   
#   guides(color = 'none') +
#   
#   xlab(NULL) +  # Remove x-axis label
#   ylab(NULL) +  # Remove y-axis label
#   ggtitle(NULL) +  # Remove plot title
#   
#   # remove title, axis text, legend, set text size, remove facet title strips 
#   theme(axis.text = element_blank(), 
#         text = element_text(size=18),
#         strip.background.x = element_blank(),
#         strip.text.x = element_blank(),
#         strip.placement = 'outside',
#         legend.direction = "horizontal") +
#   facet_wrap(~Season, ncol = 2, labeller = as_labeller(custom_labeller))
# 
# 
# g
# 
# 
# shift_legend3 <- function(p) {
#   pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
#     with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))
# 
#   if( length(pnls) == 0 ) stop( "No empty facets in the plot" )
# 
#   lemon::reposition_legend(p, "center", panel=names(pnls) )
# }
# 
# gshift = shift_legend3(g)
# 
# ggsave(filename =  '../plots/manuscript/Hg_Pre-Construction-ProbBenthic_HexagonalTessellation.png',
#        plot = gshift, device = 'png', dpi = 600, scale = 1, width = 10, height = 8)
# 
# 
# ggsave(filename =  '../plots/manuscript/Hg_Pre-Construction-ProbBenthic_HexagonalTessellation.svg',
#        plot = gshift, device = 'svg', dpi = 600, scale = 1, width = 10, height = 8)
# 


#################################################################################################################
############################################# Get Sediment Grain Sizes for 4th panel ############################
#################################################################################################################
sedi = raster::raster('C:/Users/Eleanor.heywood/Documents/ArcGIS/Projects/ManuscriptFigure1/SoftSediments.tif')

benthic = data
benthic$SediGrainSize = raster::extract(x = sedi, y = benthic)
hist(benthic$SediGrainSize, breaks = 50)

benthic = benthic[!is.na(benthic$SediGrainSize),]

benthic$SediClassification = NA


benthic$SediClassification[which(benthic$SediGrainSize > 0 & benthic$SediGrainSize <= 0.06)] = "Clay/Silt"
benthic$SediClassification[which(benthic$SediGrainSize > 0.06 & benthic$SediGrainSize <= 0.125)] = "Very Fine Sand"
benthic$SediClassification[which(benthic$SediGrainSize > 0.125 & benthic$SediGrainSize <= 0.25)] = "Fine Sand"
benthic$SediClassification[which(benthic$SediGrainSize > 0.25 & benthic$SediGrainSize <= 0.5)] = "Medium Sand"
benthic$SediClassification[which(benthic$SediGrainSize > 0.5 & benthic$SediGrainSize <= 1)] = "Coarse Sand"
benthic$SediClassification[which(benthic$SediGrainSize > 1 & benthic$SediGrainSize <= 2)] = "Very Coarse Sand"
benthic$SediClassification[which(benthic$SediGrainSize > 2)] = "Gravel/Granule"

benthic$SediClassification = factor(benthic$SediClassification, 
                                    levels = c('Clay/Silt', 'Very Fine Sand', 
                                               'Fine Sand', 'Medium Sand', 'Coarse Sand',
                                               'Very Coarse Sand', 'Gravel/Granule'), 
                                    labels = c('Cl/Si', 'VF Sand', 
                                               'F Sand', 'M Sand', 'C Sand',
                                               'VC Sand', 'Gravel'))

benthic$BenthDem = factor(as.character(benthic$BenthDem), levels = c("0", '1'), labels = c("Pelagic", 'Benthic/Demersal'))
e = ggplot(data = benthic) +
  geom_histogram(mapping = aes(x = SediClassification, fill = BenthDem, 
                               color = BenthDem, group = BenthDem), 
                 stat = "count", position = position_dodge(), alpha = 0.5) +
  scale_fill_manual(values = c("#d2dfe6", '#71644a')) +
  scale_color_manual(values = c("#4a6e76", '#71644a')) +
  scale_y_continuous(name = 'Number of Dives (K)', breaks = c(0, 5000,10000,15000,20000), labels = c('0','5', '10', '15', '20')) +
  ylab("Number of Dives")+
  theme_bw() +
  coord_cartesian(ylim = c(0,20500), expand = F) +
  theme(legend.position = 'none', axis.title.x = element_blank(), text = element_text(size = 14),
        axis.text.x = element_text(face = 'bold', angle = 0, vjust = 0.63))
e
# create tesselation panels for each Season
a = ggplot(data = resl[resl$Season == 'Winter',]) +
  
  # Add 200m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               linewidth=c(0.3),
               colour="gray20", lty = 1)+
  
  # Add raster
  geom_sf(data = resl[resl$Season == 'Winter',], mapping = aes(fill = PropBenthic, color = PropBenthic), alpha = 0.75) +
  
  #scale_fill_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  #scale_color_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  scale_fill_gradientn(colors = custom_palette, 
                       guide = guide_colorbar(theme = theme(legend.title.position = 'top',
                                                            legend.key.width = unit(2.5, 'in')), 
                                              title = 'Proportion Benthic / Demersal Diving', direction = 'horizontal')) +
  scale_color_gradientn(colors = custom_palette, ) +
  # Add depth contours
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               linewidth=c(0.3),
               colour="gray1", lty = 3)+

  # Add wind energy areas
  geom_sf(data = we, color = "black", fill = NA, linewidth = 0.65) +
  
  # Add land
  geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
  
  
  coord_sf(xlim = c(-75, -63.5), ylim = c(38.75, 45)) +
  
  theme_bw() +
  
  guides(color = 'none') +
  
  xlab(NULL) +  # Remove x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle(NULL) +  # Remove plot title
  
  # remove title, axis text, legend, set text size, remove facet title strips 
  theme(axis.text = element_blank(), 
        text = element_text(size=14),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = 'outside',
        legend.position =  'none')

b = ggplot(data = resl[resl$Season == 'Spring',]) +
  
  # Add 200m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               linewidth=c(0.3),
               colour="gray20", lty = 1)+
  
  # Add raster
  geom_sf(data = resl[resl$Season == 'Spring',], mapping = aes(fill = PropBenthic, color = PropBenthic), alpha = 0.75) +
  
  #scale_fill_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  #scale_color_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  scale_fill_gradientn(colors = custom_palette, 
                       guide = guide_colorbar(theme = theme(legend.title.position = 'top',
                                                            legend.key.width = unit(2.5, 'in')), 
                                              title = 'Proportion Benthic / Demersal Diving', direction = 'horizontal')) +
  scale_color_gradientn(colors = custom_palette, ) +
  # Add depth contours
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               linewidth=c(0.3),
               colour="gray1", lty = 3)+
  
  # Add wind energy areas
  geom_sf(data = we, color = "black", fill = NA, linewidth = 0.65) +
  
  # Add land
  geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
  
  
  coord_sf(xlim = c(-75, -63.5), ylim = c(38.75, 45)) +
  
  theme_bw() +
  
  guides(color = 'none') +
  
  xlab(NULL) +  # Remove x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle(NULL) +  # Remove plot title
  
  # remove title, axis text, legend, set text size, remove facet title strips 
  theme(axis.text = element_blank(), 
        text = element_text(size=14),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = 'outside',
        legend.position =  'none')

c = ggplot(data = resl[resl$Season == 'Summer',]) +
  
  # Add 200m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               linewidth=c(0.3),
               colour="gray20", lty = 1)+
  
  # Add raster
  geom_sf(data = resl[resl$Season == 'Summer',], mapping = aes(fill = PropBenthic, color = PropBenthic), alpha = 0.75) +
  
  #scale_fill_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  #scale_color_scico(palette = 'brocO', alpha = 0.85, direction = 1) +
  scale_fill_gradientn(colors = custom_palette, 
                       guide = guide_colorbar(theme = theme(legend.title.position = 'top',
                                                            legend.key.width = unit(2.5, 'in')), 
                                              title = 'Proportion Benthic / Demersal Diving', direction = 'horizontal')) +
  scale_color_gradientn(colors = custom_palette, ) +
  # Add depth contours
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               linewidth=c(0.3),
               colour="gray1", lty = 3)+
  
  # Add wind energy areas
  geom_sf(data = we, color = "black", fill = NA, linewidth = 0.65) +
  
  # Add land
  geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
  
  
  coord_sf(xlim = c(-75, -63.5), ylim = c(38.75, 45)) +
  
  theme_bw() +
  
  guides(color = 'none') +
  
  xlab(NULL) +  # Remove x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle(NULL) +  # Remove plot title
  
  # remove title, axis text, legend, set text size, remove facet title strips 
  theme(axis.text = element_blank(), 
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = 'outside',
        legend.position = "inside", legend.position.inside = c(0.7, 0.075), 
        legend.direction = "vertical", 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"),           # Reduce the size of the legend keys (squares/lines)
        legend.spacing.x = unit(0.5, 'cm'),             # Reduce spacing between legend items
        legend.spacing.y = unit(0.5, 'cm'),
        legend.title = element_blank())              # Reduce spacing between legend rows)

################################################# MS - sediment map

sedidf = as.data.frame(terra::rast(sedi), xy=T)

# Custom bins 
cust_bins = c(0, 0.06, 0.125, 0.25, 0.5, 1, 2, max(sedidf$SoftSediments))
sedidf_bin = sedidf %>% mutate(grainsize2 = cut(SoftSediments, breaks = cust_bins))
sedidf_bin$grainsize2 = factor(sedidf_bin$grainsize2, labels = c('Clay/Silt', 'Very Fine Sand',
                                                                 'Fine Sand', 'Medium Sand', 'Coarse Sand', 
                                                                 'Very Coarse Sand', 'Gravel'))
unique(sedidf_bin$grainsize2)


d = ggplot() +
  # Add land
  geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
  # Add raster
  geom_raster(data = sedidf_bin, mapping = aes(x=x, y=y, fill = grainsize2)) +
  scale_fill_manual(values = color_ramp <- colorRampPalette(c("#6b2c1a", "#9b421c", "#c3561c", "#e08214", 
                                                              "#fdb863", "#fee08b", "#fff3b0"))(7)) +
  # Add 200m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               linewidth=c(0.3),
               colour="gray10", lty = 1)+
  
  # Add depth contours
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               linewidth=c(0.3),
               colour="gray1", lty = 3)+
  
  # Add wind energy areas
  geom_sf(data = we, color = "black", fill = NA, linewidth = 0.65) +
  
  coord_sf(xlim = c(-75, -63.5), ylim = c(38.75, 45)) +
  theme_bw() +
  guides(color = 'none') +
  xlab(NULL) +  # Remove x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle(NULL) +  # Remove plot title
  # remove title, axis text, legend, set text size, remove facet title strips 
  theme(axis.text = element_blank(), 
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = 'outside', legend.title = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.675, 0.075), 
        legend.direction = "horizontal", 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"),           # Reduce the size of the legend keys (squares/lines)
        legend.spacing.x = unit(0.5, 'cm'),             # Reduce spacing between legend items
        legend.spacing.y = unit(0.5, 'cm')              # Reduce spacing between legend rows
  )

d

a = a + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 14), plot.tag.location ='panel')
b = b + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 14), plot.tag.location ='panel')
c = c + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 14), 
              legend.text = element_text(size = 10), plot.tag.location ='panel', 
              legend.title = element_blank(), legend.box.margin = margin(0,0,0,0), 
              legend.background = element_blank())
e = e + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 14), plot.tag.location ='panel')

d = d + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 14), legend.text = element_text(size = 9),
              plot.tag.location ='panel', legend.box.margin = margin(0,0,0,0), legend.background = element_blank(),
              legend.direction = 'vertical', legend.position.inside = c(0.87, 0.15))
d
library(patchwork)

layout = '
AB
CD
EE
'
mainfig7 = wrap_plots(D = d, C = c, B = b, A=a, E=e, design = layout, heights = c(1,1,0.75), widths = c(1,1)) +
  plot_annotation(tag_levels = 'a')

# mainfig7 = a + b + c + e + plot_layout(ncol = 2, nrow = 2) + 
#   plot_annotation(tag_levels = 'a')

mainfig7
ggsave(filename = '../plots/manuscript/Hg_Pre-Construction-ProbBenthic_HexagonalTessellation_withSediment.pdf', 
       plot = mainfig7, device = 'pdf', scale = 1.5, width = 170, height = 170, units = 'mm')

ggsave(filename = '../plots/manuscript/Hg_Pre-Construction-ProbBenthic_HexagonalTessellation_withSediment.png', 
       plot = mainfig7, device = 'png', scale = 1.5, width = 170, height = 170, units = 'mm')

############################################################################################################
######################################### SEDIMENT MAP #####################################################
############################################################################################################

ggsave(filename = '../plots/manuscript/Supplemental_SedimentGrainSize.png', 
       plot = d, device = 'png', scale = 1, width = 170, height = 130, dpi = 320, units = 'mm')

# WEB - sediment map
ggsave(filename = '../plots/manuscript/Supplemental_SedimentGrainSize.pdf', 
       plot = d, device = 'pdf', scale = 2, width = 1200, height = 1000, units = 'px')

#################################################################################################################
#################################################################################################################
################################################## SCRATCH ######################################################
#################################################################################################################
####################################################### Minimum Convex Hull #################################################
# # Create a MCP and then use it as the prediction grid for kriging
# 
# # Create MCP
# data.mcp = st_convex_hull(x = st_union(data))
# 
# plot(data.mcp)
# 
# # Create bounding box and grid
# grd = data.mcp %>% st_as_stars(dx = 5000) %>% st_crop(data.mcp)

############################################# KRIGING ########################################################################

# # Variogram - spatial predictions using geostatistical methods
# 
# v = variogram(PropWC~1, locations = data, data = data)
# 
# plot(v, plot.numbers = TRUE, xlab = "distance h [m]",
#      ylab = expression(gamma(h)),
#      xlim = c(0, 1.055 * max(v$dist)))
# 
# # Fit variogram model
# v.m = fit.variogram(v, vgm(1, model = 'Exp', range = 50000, nugget = 1))
# plot.new()
# plot(v, v.m, cutoff = 30000, pch =20, col ="red")
# 
# # Kriging
# rm(list = setdiff(ls(), c('v.m', 'v', 'data', 'grd', 'bathy_raster', 'bathy_stars', 'bathydf', 'landUTM')))
# 
# k = krige(formula = PropWC ~ 1, locations = data, newdata = grd, model = v.m, maxdist = 10000)
# 
# xlim = c(min(st_coordinates(data)[,1]) -1000, max(st_coordinates(data)[,1]) +1000)
# ylim = c(min(st_coordinates(data)[,2])-1000, max(st_coordinates(data)[,2])+1000)
# 
# library(viridis)
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_viridis(option = "E", name = 'Prop. Water Column') +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() +
#   theme(legend.position = 'bottom') 
#   
# 
# g
# 
# ggsave(filename = './plots/Hg-2019-2023-Continuous-DiveWaterColumn-Krige.jpeg', plot = g, device = 'jpeg', dpi = 300)
# 
# 
# class(tmp)
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_binned(type = 'viridis', breaks = c(0.05, 0.25, 0.50, 0.75, 0.95)) +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() 
# 
# g
# 
# 
# ggsave(filename = './plots/Hg-2019-2023-Binned-DiveWaterColumn-Krige.jpeg', plot = g, device = 'jpeg', dpi = 300)
# 
# 
# fulldatap = ggplot() +  geom_sf(data = data) +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim)

# # Jan - March
# data$Month = format(data$date, '%b')
# 
# 
# winter = data[data$Month %in% c('Jan', 'Feb'),]
# # Create MCP
# data.mcp = st_convex_hull(x = st_union(winter))
# 
# plot(data.mcp)
# 
# # Create bounding box and grid
# grd = data.mcp %>% st_as_stars(dx = 5000) %>% st_crop(data.mcp)
# 
# v = variogram(PropWC~1, locations = winter, data = winter)
# 
# plot(v, plot.numbers = TRUE, xlab = "distance h [m]",
#      ylab = expression(gamma(h)),
#      xlim = c(0, 1.055 * max(v$dist)))
# 
# # Fit variogram model
# v.m = fit.variogram(v, vgm(1, model = 'Exp', range = 50000, nugget = 2))
# plot.new()
# plot(v, v.m, cutoff = 30000, pch =20, col ="red")
# 
# # Kriging
# 
# k = krige(formula = PropWC ~ 1, locations = winter, newdata = grd, model = v.m, maxdist = 10000)
# 
# 
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_viridis(option = "E", name = 'Prop. Water Column') +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() +
#   theme(legend.position = 'bottom') 
# 
# 
# g
# 
# ggsave(filename = './plots/Hg-2019-2023-JANFEB-Continuous-DiveWaterColumn-Krige.jpeg', plot = g, device = 'jpeg', dpi = 300)
# 
# 
# 
# 
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_binned(type = 'viridis', breaks = seq(0.25, 0.95, 0.10), name = 'Prop. Water Column') +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() +
#   theme(legend.position = 'right') 
# 
# g
# 
# 
# 
# ##################### MARCH-May
# winter = data[data$Month %in% c('Mar', 'Apr', 'May'),]
# # Create MCP
# data.mcp = st_convex_hull(x = st_union(winter))
# 
# plot(data.mcp)
# 
# # Create bounding box and grid
# grd = data.mcp %>% st_as_stars(dx = 5000) %>% st_crop(data.mcp)
# 
# v = variogram(PropWC~1, locations = winter, data = winter)
# 
# plot(v, plot.numbers = TRUE, xlab = "distance h [m]",
#      ylab = expression(gamma(h)),
#      xlim = c(0, 1.055 * max(v$dist)))
# 
# # Fit variogram model
# v.m = fit.variogram(v, vgm(1, model = "Sph", range = 50000, nugget = 1))
# plot.new()
# plot(v, v.m, cutoff = 30000, pch =20, col ="red")
# 
# # Kriging
# 
# k = krige(formula = PropWC ~ 1, locations = winter, newdata = grd, model = v.m, maxdist = 10000)
# 
# 
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_viridis(option = "E", name = 'Prop. Water Column') +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() +
#   theme(legend.position = 'bottom') 
# 
# 
# g
# 
# ggsave(filename = './plots/Hg-2019-2023-JANFEB-Continuous-DiveWaterColumn-Krige.jpeg', plot = g, device = 'jpeg', dpi = 300)
# 
# 
# 
# 
# g = ggplot() +  
#   geom_stars(data = k, aes(fill = var1.pred, x = x, y = y), na.action = na.omit) +
#   scale_fill_binned(type = 'viridis', breaks = seq(0.25, 0.95, 0.10), name = 'Prop. Water Column') +
#   xlab(NULL) + ylab(NULL) +
#   geom_sf(data = landUTM) + 
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_bw() +
#   theme(legend.position = 'right') 
# 
# g
# 


###############################################################################################################################
###############################################################################################################################
####################################################### Kernel Density of Benthic Dives ####################################################
###############################################################################################################################
###############################################################################################################################

# library(splancs)
# library(marmap)
# library(ggplot2)
# df = data[data$DiveType == 'Benthic', ]
# df$t = as.numeric(format(df$date, "%j"))
# 
# dfpoints = df %>% mutate(x = st_coordinates(.)[,2], y = st_coordinates(.)[,1]) %>% st_drop_geometry()
# 
# dfpoints = as.points(dfpoints)
# 
# sealpts = coordinates(dfpoints)
# 
# boundingbox = bbox(sealpts)
# 
# xgr = seq(boundingbox[1,1], boundingbox[1,2], 10000) # values of x (lon) at which to compute the kernel
# ygr = seq(boundingbox[2,1], boundingbox[2,2], 10000) # values of y (lat) at which to compute the kernel
# zgr = seq(min(df$t) + 8, max(df$t) - 8, 15)
# b3d = kernel3d(pts = sealpts, times = df$t, xgr = xgr, ygr = ygr, zgr = zgr, hxy = 50*1000, hz = 8)
# 
# 
# # BATHY AND LAND
# 
# bathy = marmap::getNOAA.bathy(lon1 = min(data$lon) - 1, lon2 = max(data$lon) + 1, 
#                               lat1 = min(data$lat) - 1, lat2 = max(data$lat) + 1, resolution = 1)
# 
# bf = fortify.bathy(bathy)
# 
# bf = bf %>% st_as_sf(coords = c('x', 'y')) %>% st_set_crs(4326) %>% st_transform(projUTM) %>%
#   mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% st_drop_geometry()
# 
# 
# xlim = round(c(min(st_coordinates(df)[,1]), max(st_coordinates(df)[,1])), digits = 0)
# ylim = round(c(min(st_coordinates(df)[,2]), max(st_coordinates(df)[,2])), digits = 0)
# 
# 
# 
# library(raster)
# 
# 
# for (i in 1:length(zgr)) {
#   
#   m = list(x=xgr, y = ygr, z = b3d$v[,,i])
#   
#   r <- raster::raster(m, crs = CRS(projUTM))
#   
#   r1 = as(r, 'SpatRaster')
#   
#   r1 = terra::
#   
#   contours = amt::hr_isopleths(r1, levels = c(0.50, 0.95))
#   
#   g = ggplot()+
#     
#     # add 100m contour
#    # geom_contour(data = bf, 
#   #               aes(x=x, y=y, z=z),
#   #               breaks=c(-100),
#   #               size=c(0.3),
#   #               colour="grey")+
# 
#     # # add 250m contour
#     # geom_contour(data = bf, 
#     #              aes(x=x, y=y, z=z),
#     #              breaks=c(-250),
#     #              size=c(0.6),
#     #              colour="grey")+
#     # 
#     # # add coastline
#     geom_sf(data = landUTM, 
#                   fill = "darkgrey", color = NA) + 
#    
#     
#     # add male polygon
#     geom_sf(data = contours, mapping = aes(color = factor(level), fill = factor(level)), alpha = 0.25) +
#     scale_fill_manual(values = c('violet', 'purple4')) +
#     scale_color_manual(values = c('violet', 'purple4')) +
#     
#     # Add Boem wind areas
#     #geom_sf(data = wea, color = 'red', fill = "red", linewidth = 0.5, alpha = 0.5)+
#     
#     # configure projection and plot domain
#     coord_sf(xlim = xlim, ylim = ylim)+
#     
#     # formatting
#     ylab("")+xlab("")+
#     theme_bw() +
#     theme(legend.position = 'None')
#   g 
#   }
