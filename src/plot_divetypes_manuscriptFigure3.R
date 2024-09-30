# EIH
# Revised 2024-09-19
# PLOTTING DIVE TYPES OVER TIME BY SEX AND ACROSS SEASONS
# MANUSCRIPT FIGURE 3
# SUPPLEMENTAL FIGURES

############################# SETUP #############################################
#################################################################################
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
setwd('~/seal_telemetry/')

# Load all imputed dive data with distance to shore and bathymetry
dives = read_csv(file = './data/L3/dive/Hg_2019-2023_DiveTypeClassification_withDistanceToShore.csv')

# Subset to pre-construction
dives = dives[which(dives$date < as.POSIXct('2023-06-01', tz = 'UTC')),]

# Get those dives in bathymetry >5m
dives = dives[dives$meanbathy_m <= -5, ]

# The number of ptts per sex for dive data
tmp = dives[,c('ptt', 'sex')]
tmp = distinct(tmp)
table(tmp$sex)

shallowcutoff = -20

dives$DiveType = ifelse(dives$PWCProb_Benthic > 0.5, 'Benthic/Demersal', 'Pelagic')
dives$DiveType[which(dives$meanbathy_m > shallowcutoff)] = 'Shallow Water (<20m)'
dives$DiveType = factor(dives$DiveType, levels = c("Pelagic", 'Benthic/Demersal', 'Shallow Water (<20m)'))



# ISO WEEK and Year
library(ISOweek)
dives$ISOweek = date2ISOweek(dives$date)
dives$ISOweek = substr(dives$ISOweek, 1, 8)
dives$week = as.numeric(format(dives$date, "%V"))
dives$year = as.numeric(format(dives$date, "%Y"))

restable = dives %>% group_by(week) %>%
  # number of ptts in each given ISO week
  mutate(Nptt = length(unique(ptt))) %>%
  ungroup() %>%
  # get the number of dives per week, ptt, andn divetype
  group_by(year, week, ISOweek, sex, ptt, Nptt, DiveType) %>% 
  count(DiveType, name = 'No_dives_perdivetype', .drop = F) 

# Get the total number of dives for each ptt and week
ndivesperweekptt = dives %>% 
  group_by(week, ptt) %>%
  count(ptt, name = 'Ndivesperweekptt', .drop = FALSE)

restable = left_join(restable, ndivesperweekptt, by = join_by(ptt, week))

restable$Proportion = restable$No_dives_perdivetype / restable$Ndivesperweekptt

restable = restable %>% arrange(week, year)
# get those weeks with 3 or more ptts representing the week
restable = restable[which(restable$Nptt >= 3), ]


# LABEL THE WEEK-YEAR DATES
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

w = paste(restable$ISOweek, 1, sep = '-')
restable$ISOdate = ISOweek2date(weekdate = w)
restable$day = as.numeric(format(restable$ISOdate, '%d'))
restable = restable %>% group_by(week) %>% mutate(mday = mode(day))
restable$lab = paste(format(restable$ISOdate, '%b'), restable$mday, sep = '-')


datelabs = restable %>% ungroup() %>% 
  filter(week %in% seq(3, 31, 4), day == mday) %>% dplyr::select(week, lab) %>% distinct()

restable$week = as.factor(restable$week)

restable = restable %>% arrange(week, ptt)
p = ggplot(data = restable, 
           mapping = aes(x = week, 
                         y = Proportion, 
                         fill = DiveType, 
                         color = DiveType)) +
  stat_boxplot(outliers = FALSE, alpha = 0.75, linewidth = 0.55, coef = 1) +
  scale_fill_manual(values = c("#d2dfe6", '#71644a', "gray90"), labels = levels(dives$DiveType)) +
  scale_color_manual(values = c("#4a6e76", '#71644a', "gray80"), labels = levels(dives$DiveType)) +
  scale_x_discrete(breaks = datelabs$week, 
                   labels= datelabs$lab) +
  

  
  ylab('Proportion') +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.89,0.89),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("segment", 
           x = seq(0,30,1), xend = seq(0,30,1), #adjust the length of ticks as you need
           y = -0.05, yend = -0.025) +
  annotate("segment",
           x = seq(1,31,4), xend = seq(1,31,4), #adjust the length of ticks as you need
           y = -0.078, yend = -0.025, linewidth = 1) +
  ggh4x::coord_axes_inside(ylim = c(0,1), clip = 'off') 

p


# PLOT MIXTURE MODEL 
x = dives$medianPWC[dives$medianPWC < 1.1]
wmm <- mixR::mixfit(x = x, ncomp = 2, family = 'weibull')
p3 <- plot(wmm, 
           legend.position = 'none', title = NULL,
           xlab = 'Proportion of Water Column', 
           trans = 0.75, breaks = 30)+ 
  scale_fill_manual(values = c("#d2dfe6", '#71644a')) +
  
  annotate(geom = 'segment', x = 0.874, xend = 0.874,
           y = 0,yend = 0.75, linewidth = 1.5, color = 'black', linetype = 1) +
  annotate(geom = 'text', x = 0.874-0.025, y = 0.9, label = '~0.87') +
  scale_x_continuous(breaks = seq(0, 1.1, 0.1)) +
  coord_cartesian(ylim = c(0,5.25), xlim = c(0,1.12), expand = F, clip = 'on')

p3

##################################### Add third panel of distance and bathymetry #######################################

trans <- 1  # Assuming you have a proper value for 'trans' to match the scales of 'distkm' and 'bathy'.

# get those weeks with 3 or more ptts
dives1 = dives[which(dives$week %in% unique(restable$week)),]

db = ggplot(data = dives1, aes(x=as.factor(week))) + 
  geom_boxplot(aes(y=distanttoshore/1000), color = 'grey40', fill = 'grey70',outliers = F, coef=1) +
  geom_boxplot(aes(y = meanbathy_m / trans), color = 'grey10', fill = 'grey30', outliers = F, coef=1) +  # Invert bathy values here
  scale_y_continuous(
    name = 'Dist. offshore(km)',
    sec.axis = sec_axis(~ .*trans, name = 'Bathymetry(m)')  # Invert back for the label
  ) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = datelabs$week, 
                   labels= datelabs$lab) +
  
  
  xlab('Week') +
  annotate("segment",
           x = seq(0,30,1), xend = seq(0,30,1), #adjust the length of ticks as you need
           y = -300, yend = -285) +
  annotate("segment",
           x = seq(1,31,4), xend = seq(1,31,4), #adjust the length of ticks as you need
           y = -312, yend = -285, linewidth = 1) +
  ggh4x::coord_axes_inside(ylim = c(-275,225)
                          , clip = 'off', labels_inside = F)

db

p = p + theme(axis.title.x = element_blank(), 
              axis.text.x = element_blank(), 
              text = element_text(size = 11),
              axis.text = element_text(size = 11),
              plot.margin = margin(0,0,0,0), 
              legend.position.inside = c(1,1))
db = db + theme(text = element_text(size = 11),
                axis.text = element_text(size = 11),
                plot.margin = margin(0,0,0,0))
p3 = p3 + theme(text = element_text(size = 11), 
                axis.text = element_text(size = 11),
                plot.margin = margin(0,0,0,0))



g = (p3) / 
  (p) /
  (db) + plot_annotation(tag_levels = 'a')
g

ggsave(filename = './plots/manuscript/Figure3.pdf', 
       plot = g, 
       device = 'pdf',width = 1200, height = 1200,units = 'px', scale = 2)


################################################################################################################
##################################### SCRATCH ##################################################################
################################################################################################################
###################################### PLOT PSEUDO PROFILES #######################

# dL = split(data6, data6$id)
# 
# 
# # Get Haulout Times for plotting
# haulout = read_csv("./L1/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT_V2_vmax2mps_updated_disthreshold.csv") %>% 
#   filter(is.na(TripID)) %>%
#   group_by(id, HauloutID) %>% dplyr::reframe(s = min(date), e = max(date)) %>% distinct()
# 
# 
# haulout = read_csv("./Hg-PTT253201_HAULOUT_TRIPID_ASSIGNMENT.csv") %>% 
#   filter(is.na(TripID)) %>%
#   group_by(id, HauloutID) %>% dplyr::reframe(s = min(datetime), e = max(datetime)) %>% distinct()
# 
# 
# 
# for (i in 1:length(dL)){
#   d = dL[[i]]
#   id = unique(d$id)
#   ho = haulout[which(haulout$id == id), ]
#   if (min(d$Depth) > -30){
#     next
#   }
#   pelagic = d[d$DiveType == "Pelagic", ]
#   benthic = d[d$DiveType == "Benthic", ]
#   shallow = d[d$DiveType == "Shallow", ]
#   
#   # Set y plot limits
#   if(min(d$bathydepth) < -200){
#     ylimits = c(-200,0)
#   }else{ylimits = c(min(d$bathydepth), 0)}
#   
#   g = ggplot() +
#     geom_rect(data = ho, aes(xmin=s, xmax=e, ymin=ylimits[1], ymax=0), alpha=.4, color = 'lightgray', fill = 'lightgray')+
#     geom_point(data = d, mapping = aes(date, bathydepth, color = 'Water Depth'), size = 1, shape = 8) +
#     geom_point(data = shallow, mapping = aes(date, Depth, color = 'Shallow Dive'), alpha = 0.5) +
#     
#     geom_point(data = pelagic, mapping = aes(date, Depth, color = 'Pelagic Dive')) +
#     geom_point(data = benthic, mapping = aes(date, Depth, color = 'Benthic Dive')) +
#     ylim(ylimits) +
#     scale_color_manual(values = c('Pelagic Dive' = 'royalblue1', 'Benthic Dive' = 'brown', 'Shallow Dive' = 'gray', 'Water Depth' = 'darkblue')) +
#     ylab('Depth (m)') +
#     xlab('Month') +
#     theme_bw()
#   
#   # Create a manual legend without replacing the scale
#   g <- g +
#     guides(
#       color = guide_legend(
#         title = unique(d$id),
#         override.aes = list(
#           color = c('brown',  'royalblue1', 'gray', 'darkblue'),
#           title = c('Benthic Dive', 'Pelagic Dive',  'Shallow Dive', 'Water Depth')
#         )
#       )
#     ) +
#     theme(legend.position = 'none') + # Remove the default legend
#     theme(legend.box = "horizontal", legend.justification = c(0, 0), legend.position = c(0.01, 0.01)) # Adjust legend position    
#   g
#   ggsave(g, filename = paste0('../plots/BenthicPelagicPlots/', unique(d$id), "_benthicpelagicdiveplot_withHOShading.jpg"), dpi = 300, width = 10, height = 7)
# }


