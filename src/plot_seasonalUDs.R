# MANUSCRIPT FIGURES FOR SEASONAL UDs
# EIH
# June 6, 2024
rm(list=ls())
gc()
library(ggplot2)
library(dplyr)
library(sf)
library(marmap)
library(metR)

setwd("~/seal_telemetry/")

###################################### LOAD DATA ############################################################################################################################################################################
##############################################################################################################################################################################################################################
# read ud polygons
FUD = st_read('./data/L3/UD/Hg_PreConstruction_FemaleCollectiveUDs_AllSeasonPolygons.shp')
MUD = st_read('./data/L3/UD/Hg_PreConstruction_MaleCollectiveUDs_AllSeasonPolygons.shp')


# Load the land spatial polygons dataframe and subset based on relevant landmasses
land = rnaturalearth::ne_states()
land = land[land$admin %in% c("Canada", 'United States of America'),]
land = land[which(land$postal %in% c("ME", "NH", "VT", "RI", "CT", "MA", "PA", "NY", "DE", "NJ", "VA", "MD", "NC", 'OH', 'WV', 'SC', 'KY', 'TN', 'IL', 'WI', 'MN', 'IA', 'IN', 'MI') | land$region == "Eastern Canada"),]

## minor islands (islands less than 2 square km)
minor_islands = st_read("./READ-PSB-MoveSeals/ne_10m_minor_islands/ne_10m_minor_islands.shp") %>%
  mutate(name = featurecla) %>%
  dplyr::select(name, geometry)

## get the bounding box of the land and subset minor islands based on this bounding box
bounding_box <- st_bbox(land) %>% st_as_sfc()
minor_islands = st_intersection(minor_islands, bounding_box)

# Combine the land and minor islands
hires_land = bind_rows(land, minor_islands)
rm(minor_islands, land)

# BATHYMETRY
# GET BATHY FOR PLOTTING CONTOURS
xlim = c(-80, -55)
ylim = c(37, 50)
bathy = marmap::getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2], 
                              lat1 = ylim[1], lat2 = ylim[2], resolution = 4)
bf = fortify.bathy(bathy)

# WIND ENERGY AREAS
we = st_read("./data/shapefiles/BOEM_Wind_Leases_8_30_2023.shp")
we = st_intersection(we, bounding_box)
plot(we['STATE'])

we = st_union(we)

################################ PLOTTING ##############################################################
########################################################################################################
# Get female 95 and 50% ud
#FUD = FUD[FUD$PrcntVl %in% c('50', '95'), ]
FUD = FUD[FUD$Season %in% c("Winter", "Spring", 'Summer'), ]

FUD$Season = factor(FUD$Season, levels = c("Winter", "Spring", 'Summer'))

# Get male 95 and 50% ud
#MUD = MUD[MUD$PrcntVl %in% c('50', '95'), ]
MUD = MUD[MUD$Season %in% c("Winter", "Spring", 'Summer'), ]
MUD$Season <- factor(MUD$Season, levels = c("Winter", "Spring", 'Summer'))

# combine
AUD = bind_rows(MUD, FUD)

#### Plot both in facet grid by Season and Sex

# DEFINE COLORS
# males
dark_M = "#377eb8"
light_M = '#89b5e9'
# dark_M <- "#8B7355"  # Lighter Dark BurlyWood
# light_M <- "#CDAA7D"       # Lighter BurlyWood
# light_M <- "#FFDEAD" # Lighter BurlyWood

# females
dark_F = '#e41a1c'
light_F = '#f4a3a3'
# dark_F <- "#5A4EAB"
# light_F <- "#7A6FD6"
# light_F <- "#A3B8FF"


colorsM = c(dark_M, light_M)
colorsF = c(dark_F, light_F)
colorsM = data.frame(colors = colorsM, Sex = 'Male', PrcntVl = c('50', '95'))
colorsF = data.frame(colors = colorsF, Sex = 'Female', PrcntVl = c('50', '95'))

colors = bind_rows(colorsM, colorsF)
AUD$sex[AUD$Sex == 'Female'] = 'F'
AUD$sex[AUD$Sex == 'Male'] = 'M'

# get color mapping
AUD = left_join(AUD, colors, by = join_by(Sex, PrcntVl)) %>% 
  group_by(Season, Sex) %>% mutate(Label = paste(sex, ' = ', N, sep = '')) %>% 
  ungroup() %>%
  mutate(colorfactor = paste(Sex, PrcntVl))

# Create named vector for color mapping
color_mapping <- setNames(colors$colors, paste(colors$Sex, colors$PrcntVl))


# Plot with ggpattern for hatched patterns with colored lines and outlines
AUDplot <- ggplot() +
  
  # Add uds and color by sex and percent volume
  geom_sf(data = AUD, mapping = aes(fill = colorfactor, colour = colorfactor), alpha = 0.6) +

  scale_color_manual(values = color_mapping) +
  
  scale_fill_manual(values = color_mapping) +
  
  # Add depth contours
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               linewidth=c(0.3),
               colour="gray1", lty = 3)+
  
  geom_text_contour(data = bf, 
                    #nudge_y = -0.2,
                    aes(x = x, y = y, z = z),
                    label = '-50m',
                    breaks = c(-50),
                    min.size = 100, # Ensure a minimum length for labeling
                    check_overlap = TRUE,
                    size = 3,
                    color = "#696969",
                    fontface = "italic",
                    label.placer = label_placer_fraction(frac = 0.88, rot_adjuster = isoband::angle_identity())) + # Control the number of labels
  
  
  
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               linewidth=c(0.3),
               colour="gray10", lty = 1)+
  
  geom_text_contour(data = bf, 
                    nudge_y = -0.20,
                    aes(x = x, y = y, z = z),
                    label = '-200m',
                    breaks = c(-200),
                    min.size = 100, # Ensure a minimum length for labeling
                    check_overlap = TRUE,
                    size = 3,
                    color = "#696969",
                    fontface = "italic",
                    label.placer = label_placer_fraction(frac = 0.885, rot_adjuster = isoband::angle_identity())) + # Control the number of labels
  
  # Add wind energy areas
  geom_sf(data = we, color = "black", fill = NA, linewidth = 0.65) +
  
  
  geom_sf(data = hires_land, fill = 'gray95', color = 'gray75') +  # Set fill to NA for background polygons
  
  coord_sf(xlim = c(-75, -61), ylim = c(38.5, 45.1)) +
  
  theme_classic() +
  
  xlab(NULL) +  # Remove x-axis label
  ylab(NULL) +  # Remove y-axis label
  ggtitle(NULL) +  # Remove plot title
  
  # remove title, axis text, legend, set text size, remove facet title strips 
  theme(axis.text = element_blank(), 
        legend.position = 'none', 
        text = element_text(size=21),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 90, size = 18), 
        strip.placement = 'outside') +  # Control strip label appearance
  
  facet_grid(Season ~ Sex, switch = 'y') +
  
  geom_text(data = AUD, aes(label = Label), x = -64, y = 38.75, 
            size = 5, color = "black", hjust = 0, vjust = 0.5)

AUDplot


# SAVE AT MULTIPLE TYPES ACCORDIND TO MOVEMENT ECOLOGY JOURNAL
ggsave(filename = './plots/manuscript/Figure2.png', 
       plot = AUDplot, 
       device = 'png', dpi = 320, width = 170, height = 150, units = 'mm', scale = 1.25)

ggsave(filename = './plots/manuscript/Figure2.pdf', 
       plot = AUDplot, 
       device = 'pdf', width = 1000, height = 900, units = 'px', scale = 3)

##########################################################################################################################################

