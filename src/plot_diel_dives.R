# DIVES BY DIEL

############################# SETUP #############################################
#################################################################################
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

library(RchivalTag)

setwd('~/seal_telemetry/')

# Load Data
dives = read_csv(file = './data/L3/dive/Hg_2019-2023_DiveTypeClassification.csv')

# Subset to pre-construction
dives = dives[which(dives$date < as.POSIXct('2023-06-01', tz = 'UTC')),]

# The number of ptts per sex for dive data
tmp = dives[,c('ptt', 'sex')]
tmp = distinct(tmp)
table(tmp$sex)

shallowcutoff = -20

dives$DiveType = ifelse(dives$PWCProb_Benthic > 0.5, 'Benthic/Demersal', 'Pelagic')
dives$DiveType[which(dives$meanbathy_m > shallowcutoff)] = 'Shallow Water (<20m)'


dives$month = format(dives$date, '%b')
dives$month = factor(dives$month, 
                     levels = c("Jan", 'Feb', "Mar", 
                                "Apr", "May", "Jun", 
                                "Jul", "Aug", "Sep",
                                "Oct", "Nov", "Dec"))

dives$week = as.numeric(format(dives$date, "%W"))
dives$DOY = as.numeric(format(dives$date, "%j"))


dives = dives[dives$meanbathy_m <= -5, ]

#dives = dives[dives$DiveType != 'Shallow Water (<20m)', ]

# Get shape
# shape = read_csv('./data/L1/dive/Hg_2019-2023_BEHDiveRecords_QAQC.csv') %>% 
#   dplyr::mutate(diveID = SDPairID, DiveDur2 = DiveDur, Depth2 = -Depth, shape = shape_Dive) %>%
#   dplyr::select(diveID, DiveDur2, Depth2, shape)
# 
# dives = left_join(dives, shape, by = "diveID")
#all(dives$Depth == dives$Depth2) # TRUE - join is good

# 

dives$DiveType = factor(dives$DiveType, levels = c("Pelagic", 'Benthic/Demersal', 'Shallow Water (<20m)'))

# get lat lon
library(sf)
proj = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'
dives = st_as_sf(dives,coords = c('meanX_km', 'meanY_km')) %>% st_set_crs(st_crs(proj)) %>%
  st_transform(4326) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])


# Get daytime limits
lims = get_DayTimeLimits(pos = data.frame(datetime = dives$date, Lon = dives$lon, Lat = dives$lat))
classified = classify_DayTime(lims, twilight.set = 'naut')

divesinit = dives
# mutate to dives
dives = cbind(divesinit, classified)


# proportion dive type by hour 


dives$Season = NA
dives$Season[dives$month %in% c("Jan", "Feb")] = 'Winter'
dives$Season[dives$month %in% c("Mar", "Apr", 'May')] = 'Spring'
dives$Season[dives$month %in% c("Jun", "Jul", 'Aug')] = 'Summer'
dives$Season[dives$month %in% c("Sep", "Oct", 'Nov', 'Dec')] = 'Fall'

dives$Season = factor(dives$Season, levels = c('Winter', "Spring", "Summer", "Fall"))


# Get start and end periods of night
#nightperiods = as.data.frame(table(dives$daytime, dives$hour))
seasonaldielperiods = dives %>% st_drop_geometry() %>%
  mutate(jday = format(date, '%j', tz = 'EST'),
         HMS_sunrise = format(dives$sunrise, '%H:%M:%S', tz = "EST"), 
         HMS_sunset = format(dives$sunset, '%H:%M:%S', tz = "EST"),
         HMS_dawn = format(dives$dawn.naut, '%H:%M:%S', tz = "EST"),
         HMS_dusk = format(dives$dusk.naut, '%H:%M:%S', tz = "EST"),) 
# get decimal hour for sunrise
tmp = strsplit(seasonaldielperiods$HMS_sunrise, ':')
decHour = as.numeric(sapply(tmp, '[[', 1)) + as.numeric(sapply(tmp, '[[', 2))/60
seasonaldielperiods$sunrise_hour = decHour

# sunset
tmp = strsplit(seasonaldielperiods$HMS_sunset, ':')
decHour = as.numeric(sapply(tmp, '[[', 1)) + as.numeric(sapply(tmp, '[[', 2))/60
seasonaldielperiods$sunset_hour = decHour
# dawn
tmp = strsplit(seasonaldielperiods$HMS_dawn, ':')
decHour = as.numeric(sapply(tmp, '[[', 1)) + as.numeric(sapply(tmp, '[[', 2))/60
seasonaldielperiods$dawn_hour = decHour

tmp = strsplit(seasonaldielperiods$HMS_dusk, ':')
decHour = as.numeric(sapply(tmp, '[[', 1)) + as.numeric(sapply(tmp, '[[', 2))/60
seasonaldielperiods$dusk_hour = decHour

seasonaldielperiods = seasonaldielperiods %>%
  group_by(Season) %>% 
  summarise(Sunrise = mean(sunrise_hour),
            Sunset = mean(sunset_hour),
            Dawn = mean(dawn_hour),
            Dusk = mean(dusk_hour))


dives = left_join(dives, seasonaldielperiods, by = c('Season'))
dives$hour = as.numeric(format(dives$date, '%H', tz = 'EST'))

# Add shallow and deep dives to pelagic category
dives$DiveType = factor(dives$DiveType, levels = c('Shallow Water (<20m)', 'Pelagic', 'Benthic/Demersal'))

res = dives %>% st_drop_geometry() %>%
  group_by(Season, week, hour) %>%
  mutate(Ndives = n(), Sunrise = mean(Sunrise), Sunset = mean(Sunset), Dawn = mean(Dawn), Dusk = mean(Dusk)) %>% ungroup() %>%
  group_by(Season, week, hour, Ndives, DiveType, Sunrise, Sunset, Dawn, Dusk, .drop=F) %>%
  reframe(PropDiveType = n() / Ndives[1], meanDepth = mean(Depth, na.rm=T)) %>% ungroup() %>%
  filter(Season != 'Fall')


res$PropDiveType[is.na(res$PropDiveType)] = 0
#res$shape = factor(res$shape, levels = c('U', 'V', 'Square'))

library(patchwork)

res_minus_24 = res %>% mutate(hour = hour - 24)
res_plus_24 = res %>% mutate(hour = hour + 24)

  
res_wrapped = bind_rows(res, res_plus_24, res_minus_24)


labs = as.character(seq(0, 24, 2))
labs[13] = '0'
g2 = ggplot(subset(res_wrapped, Season == "Winter"), mapping = aes(x = hour, y = PropDiveType, group = DiveType, color = DiveType, fill = DiveType)) +
 
  # Nighttime
  geom_rect(aes(xmin = 0, xmax = Dawn, ymin = -Inf, ymax = Inf, alpha = 0.75), fill = '#1b2a49') +
  geom_rect(aes(xmin = Dusk, xmax = 24, ymin = -Inf, ymax = Inf), fill = '#1b2a49', alpha = 0.25) +
  
  # Crepuscular
  geom_rect(aes(xmin = Dawn, xmax = Sunrise, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_rect(aes(xmin = Sunset, xmax = Dusk, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_vline(mapping = aes(xintercept = Sunrise), color = 'orange', linewidth = 1.5) +
  geom_vline(mapping = aes(xintercept = Sunset), color = 'orange', linewidth = 1.5) +
  
  scale_x_continuous(
   breaks = seq(0, 24, by = 2),
   labels = labs) +
  #scale_x_continuous(limits = c(0,23), expand = c(0, 0), breaks = seq(0,22,2), labels = seq(0,22,2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_point(alpha = 0.5) +
  geom_smooth(alpha = 0.50, method = 'loess', span = 0.2) +
  scale_fill_manual(values = c("gray90", "#d2dfe6", '#71644a')) +
  scale_color_manual(values = c('gray90', "#d2dfe6", '#71644a')) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle(label = 'winter') +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0,24), expand = F)
g2


g3 = ggplot(data = subset(res_wrapped, Season == "Spring"), mapping = aes(x = hour, y = PropDiveType, color = DiveType, fill = DiveType)) +
  # Nighttime
  geom_rect(aes(xmin = 0, xmax = Dawn, ymin = -Inf, ymax = Inf, alpha = 0.75), fill = '#1b2a49') +
  geom_rect(aes(xmin = Dusk, xmax = 24, ymin = -Inf, ymax = Inf), fill = '#1b2a49', alpha = 0.25) +
  
  # Crepuscular
  geom_rect(aes(xmin = Dawn, xmax = Sunrise, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_rect(aes(xmin = Sunset, xmax = Dusk, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_vline(mapping = aes(xintercept = Sunrise), color = 'orange', linewidth = 1.5) +
  geom_vline(mapping = aes(xintercept = Sunset), color = 'orange', linewidth = 1.5) +
  scale_x_continuous(limits = c(0,2*pi), breaks = seq(0,2*pi, by = pi / 12)[2:25], labels = seq(1,24,1)) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 2),
    labels = labs) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_smooth(alpha = 0.50, method = 'loess', span = 0.2) +
  #geom_boxplot() +
  scale_fill_manual(values = c("gray90", "#d2dfe6", '#71644a')) +
  scale_color_manual(values = c('gray90', "#d2dfe6", '#71644a')) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle(label = 'spring') +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)
g3


# res$hour = factor(res$hour)
# g3 = ggplot(data = subset(res_wrapped, Season == "Spring"), mapping = aes(x = factor(hour), y = PropDiveType, color = DiveType, fill = DiveType)) +
#   # Nighttime
# #  geom_rect(aes(xmin = 0, xmax = Dawn, ymin = -Inf, ymax = Inf, alpha = 0.75), fill = '#1b2a49') +
# #  geom_rect(aes(xmin = Dusk, xmax = 23, ymin = -Inf, ymax = Inf), fill = '#1b2a49', alpha = 0.25) +
#   
#   # Crepuscular
#   #geom_rect(aes(xmin = Dawn, xmax = Sunrise, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
#   #geom_rect(aes(xmin = Sunset, xmax = Dusk, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
#   #geom_vline(mapping = aes(xintercept = Sunrise), color = 'orange', linewidth = 1.5) +
#   #geom_vline(mapping = aes(xintercept = Sunset), color = 'orange', linewidth = 1.5) +
#   #scale_x_continuous(limits = c(0,23), expand = c(0, 0), breaks = seq(0,24,2), labels = seq(0,24,2)) +
#   #scale_x_continuous(breaks = seq(0, 22, by = 2),
# #                     labels = c("0000","0200", "0400", "0600", "0800", "1000", "1200", 
# #                                "1400", "1600", "1800", "2000", "2200"), limits = c(0, 24)) +
#   #scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
#   #geom_smooth(alpha = 0.50, method = 'loess') +
#   geom_boxplot(position = position_dodge()) +
#   #geom_point() +
#   scale_fill_manual(values = c("gray90", "#aec8d6", '#8e7e5e')) +
#   scale_color_manual(values = c('gray90', "#aec8d6", '#8e7e5e')) +
#   ylab(NULL) +
#   xlab(NULL) +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   coord_cartesian(xlim = c(0,24))
# g3



g4 = ggplot(subset(res_wrapped, Season == "Summer"), mapping = aes(x = hour, y = PropDiveType, group = DiveType, color = DiveType, fill = DiveType)) +
  # Nighttime
  geom_rect(aes(xmin = 0, xmax = Dawn, ymin = -Inf, ymax = Inf, alpha = 0.75), fill = '#1b2a49') +
  geom_rect(aes(xmin = Dusk, xmax = 24, ymin = -Inf, ymax = Inf), fill = '#1b2a49', alpha = 0.25) +
  
  # Crepuscular
  geom_rect(aes(xmin = Dawn, xmax = Sunrise, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_rect(aes(xmin = Sunset, xmax = Dusk, ymin = -Inf, ymax = Inf), fill = '#3e5b96', color = '#3e5b96') +
  geom_vline(mapping = aes(xintercept = Sunrise), color = 'orange', linewidth = 1.5) +
  geom_vline(mapping = aes(xintercept = Sunset), color = 'orange', linewidth = 1.5) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 2),
    labels = labs) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_point(alpha = 0.5) +
  geom_smooth(mapping = aes(x = hour, y = PropDiveType, group = DiveType, color = DiveType, fill = DiveType),
              alpha = 0.50, method = 'loess', span = 0.3) +
  scale_fill_manual(values = c("gray90", "#d2dfe6", '#71644a')) +
  scale_color_manual(values = c('gray90', "#d2dfe6", '#71644a')) +
  ylab('Proportion') +
  xlab("Hour (EST)") +
  ggtitle(label = 'summer') +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0,24), expand = F)
g4


dives$daytime.long[dives$daytime.long %in% c('Dusk', 'Dawn')] = 'Crepuscular'
dives$DiveType = factor(dives$DiveType, levels = c('Shallow Water (<20m)', 'Pelagic', 'Benthic/Demersal'))

dives = dives[dives$Season != 'Fall', ]
g = ggplot(dives, mapping = aes(x = daytime, group = DiveType, fill = DiveType, color = DiveType)) +
  geom_histogram(stat = 'count', position = position_dodge(), alpha = 0.50) +
  scale_fill_manual(values = c("gray90", "#d2dfe6", '#71644a')) +
  scale_color_manual(values = c("gray80", "#4a6e76", '#71644a')) +
  theme_bw() +
  xlab(NULL) +
  ylab('# Dives') +
  theme(legend.position = 'none') +
  facet_wrap(~Season) +
  scale_y_continuous(limits = c(0,30000), expand = c(0,0)) +
  coord_cartesian(expand = F)
g

# Combine the plots into a 2x2 grid

# Standardize text size
library(patchwork)
g2 = g2 + theme(text = element_text(size = 16))
g3 = g3 + theme(text = element_text(size = 16))
g4 = g4 + theme(text = element_text(size = 16))
g = g + theme(text = element_text(size = 16))

combined_plot <- g2 + g3 + g4 + g +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'a')
combined_plot

ggsave(filename = './plots/manuscript/Figure5_dieldive_withpts.png', 
       plot = combined_plot, device = 'png', dpi = 300, 
       width = 170, height = 170, units = 'mm', scale = 1.5)

sampesizesummary = dives %>% st_drop_geometry() %>% group_by(Season, sex) %>% summarise(Nptt = length(unique(ptt)), Ndives = n())
t = flextable::flextable(sampesizesummary)
flextable::save_as_docx(t, path = './plots/manuscript/AnimalBiotelemetrySubmission/TableDiveSampleSizeSummary.docx')
