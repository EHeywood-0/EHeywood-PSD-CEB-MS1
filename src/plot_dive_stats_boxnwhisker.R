# Plot Supplemental Figs. 4 5 and 6
# Plot main manuscript Fig. 4
# EIH
# 2024-09-20

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
setwd('C:/Users/Eleanor.heywood/Documents/seal_telemetry/')

###############################################################################################################
##################################### SETUP & DATA LOAD #######################################################
###############################################################################################################

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
dives$doy = as.numeric(format(dives$date, '%j'))

###############################################################################################################
##################### PLOT SUPPLEMENTAL FIGURE 4 ##############################################################
###############################################################################################################
library(ggpattern)

g = ggplot() + 
  geom_histogram(data = dives[dives$DiveType!='Shallow Water (<20m)', ], mapping = aes(x = PWCProb_Benthic), 
                 bins = 30, color = 'gray25', fill = 'gray75') +
  theme_bw() +
  xlab('Probability Benthic / Demersal') +
  ylab('# Dives')
g
ggsave(filename = './plots/manuscript/Supplemental Figures/FigureS4.png', 
       plot = g, device = "png", dpi = 300, units = 'mm',scale = 0.8, width = 170, height = 140)


###############################################################################################################
###############################################################################################################

###############################################################################################################
##################################### SUPPLEMENTAL FIGURE 5 ###################################################
###############################################################################################################

# Summerize dive data by sex, ptt, and month
smonth = dives %>% group_by(month, sex, ptt) %>%
  summarise(`Mean Max Depth (m)` = round(mean(Depth, na.rm=T), 2),
            `Mean Duration (s)` = round(mean(DiveDur, na.rm=T), 2),
            `Mean IDI (s)` = round(mean(IDI, na.rm=T), 2)) %>%
  arrange(month, ptt) %>%
    group_by(month, sex) %>%
    mutate(SealsCount = n_distinct(ptt)) %>% ungroup()

# Get rid of months without sufficient data for visualizations
smonth = smonth[-which(smonth$month %in% c('Sep',"Oct", 'Nov', 'Dec')),]
smonth$month = factor(smonth$month, levels = c('Jan', 'Feb', "Mar", "Apr", "May", "Jun", "Jul", 'Aug'))


# Create a new variable for labeling with both M and F counts
smonth$label <- paste(smonth$sex, smonth$SealsCount, sep = '=')
smonth$sex = factor(smonth$sex, levels = c('F', 'M'))
labels = smonth %>% group_by(month, sex) %>% summarise(Label = unique(label)) %>%
  arrange(month, sex)

bluemale = '#377eb8'
redfemale = '#e41a1c'

gg1 <- ggplot() +
  geom_boxplot(smonth, mapping = aes(x = month, 
                                     y = `Mean Max Depth (m)`, 
                                     color = sex, fill = sex), 
               alpha = 0.25, outliers = F) +
  ylab('Mean max depth (m)') +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  geom_text(aes(x = labels$month, y = c(rep(c(-7,-2),8))), 
                label = labels$Label, 
                color = rep('gray30', 16), group = 'B', hjust = 0, nudge_x = -0.25) +
  theme_bw() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.1, 0.175),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'white', colour = 'gray20'),
        axis.title.x = element_blank()) 



gg2 <- ggplot() +
  geom_boxplot(smonth, mapping = aes(x = month, y = `Mean Duration (s)`, 
                                     color = sex, fill = sex), 
               alpha = 0.25, outliers = F) +
  ylab('Mean Duration (s)') +
  xlab('Month') +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_y_continuous(breaks = seq(50,320,20))+
  theme_bw() +
  theme(legend.position = 'none') 


gg3 <- ggplot() +
  geom_boxplot(smonth, mapping = aes(x = month, y = `Mean IDI (s)`, 
                                     color = sex, fill = sex), 
               alpha = 0.25, outliers = F) +
  ylab('Mean IDI (s)') +
  xlab('Month') +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_y_continuous(breaks = seq(50,750,50))+
  theme_bw() +
  theme(legend.position = 'none') 



# Combine using patchwork
comb = (gg1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 12))) /
        (gg2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 12))) /
        (gg3 + theme(axis.title.x = element_blank(), text = element_text(size = 12))) + 
  plot_annotation(tag_levels = 'a')

comb 
ggsave(filename = "./plots/manuscript/Supplemental Figures/FigureS5.png", plot = comb, 
       device = 'png', dpi = 300, width = 170, height = 170, units = 'mm',scale = 1.2)

###############################################################################################################
###############################################################################################################
###############################################################################################################

##################################################################################################################
########################################## SUPPLEMENTAL FIGURE 6 ##########################################################
##################################################################################################################
# DIVE TYPES BY SEX & SEASON
dives$month = format(dives$date, '%b')
dives$month = factor(dives$month, 
                     levels = c("Jan", 'Feb', "Mar", 
                                "Apr", "May", "Jun", 
                                "Jul", "Aug", "Sep",
                                "Oct", "Nov", "Dec"))


dives$season = NA
dives$season[which(dives$month %in% month.abb[c(1,2)])] = 'Winter'
dives$season[which(dives$month %in% month.abb[3:5])] = 'Spring'
dives$season[which(dives$month %in% month.abb[6:8])] = 'Summer'
dives$season[which(dives$month %in% month.abb[9:11])] = 'Fall'

dives = dives[!is.na(dives$season), ]
dives = dives[dives$season != 'Fall', ]

restable2 = dives %>%
  # get the number of dives per month per ptt
  group_by(season, month, sex, ptt) %>% 
  mutate(totaldives = n()) %>% ungroup() %>%
  group_by(season, month, sex, ptt, totaldives, DiveType) %>%
  count(DiveType, name = 'DTcount', .drop = F)

restable2 = restable2[!is.na(restable2$totaldives), ]

restable2$DTProp = restable2$DTcount / restable2$totaldives 

restable2$season = factor(restable2$season, levels = c('Winter', 'Spring', 'Summer'))

restable2$DiveType = as.character(restable2$DiveType)
unique(restable2$DiveType)

restable2$DiveType = factor(restable2$DiveType, levels = c('Shallow Water (<20m)', 'Benthic/Demersal', 'Pelagic'))



p4 = ggplot(data = restable2, mapping = aes(x = DiveType, y = DTProp, color = sex, fill = sex)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_color_manual(values = c(redfemale, bluemale), labels = c('Female', 'Male')) +
  scale_x_discrete(labels = c('SW', 'B/D', 'P')) +
  theme_bw() +
  coord_cartesian(clip = 'on') +
  facet_grid(~season) +
  ylab('Proportion of Total Dives') +
  xlab('') +
  theme(legend.title = element_blank(), 
        legend.position = 'inside',
        legend.position.inside = c(0.925,0.90),
        legend.background = element_rect(fill = 'white', colour = 'gray50'),
        text = element_text(size = 18),
        #axis.text.x = element_text(angle = 45, vjust = 0.55),
        strip.background = element_blank(),
        strip.placement = 'outside')   # Control strip label appearance



p4

ggsave(filename = './plots/manuscript/Supplemental Figures/FigureS6.png', 
       plot = p4, device = 'png', dpi = 300, width = 170, height = 85, units = 'mm',scale = 1.5)
ggsave(filename = './plots/manuscript/Supplemental Figures/FigureS6.pdf', 
       plot = p4, device = 'pdf', width = 1000, height = 650, units = 'px',scale = 3)


##################################################################################################################
##################################### MAIN MANUSCRIPT FIGURE 4 ###################################################
##################################################################################################################

Nppts = dives %>% group_by(week) %>% 
  summarise(Nptt = length(unique(ptt))) %>% filter(Nptt >= 3) 

dat = dives %>% filter(dives$week %in% Nppts$week) %>%
  arrange(week)
dat$week = factor(dat$week)
# LABEL THE WEEK-YEAR DATES
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

w = paste(dat$ISOweek, 1, sep = '-')
dat$ISOdate = ISOweek2date(weekdate = w)
dat$day = as.numeric(format(dat$ISOdate, '%d'))
dat = dat %>% group_by(week) %>% mutate(mday = mode(day))
dat$lab = paste(format(dat$ISOdate, '%b'), dat$mday, sep = '-')


datelabs = dat %>% ungroup() %>% 
  filter(week %in% seq(3, 31, 4), day == mday) %>% dplyr::select(week, lab) %>% distinct()

# Plot dive depth, duration and idi by dive type. 
gg1 <- ggplot() +
  geom_boxplot(dat, mapping = aes(x = week, 
                                     y = Depth, 
                                     fill = DiveType, colour = DiveType), 
               alpha = 0.75, outliers = F, coef = 1) + 

  ylab('Max depth (m)') +
  scale_fill_manual(values = c("#d2dfe6", '#71644a', "gray90", "gray90"), labels = levels(dat$DiveType)) +
  scale_color_manual(values = c("#4a6e76", '#71644a', "gray80", "gray90"), labels = levels(dat$DiveType)) +

  
  scale_x_discrete(breaks = sort(datelabs$week),
                  labels= datelabs$lab) +
  scale_y_continuous(breaks = seq(-180, 0, 20)) +
  

  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.15,0.2),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        #text = element_text(size = 18)
  ) +
  annotate("segment", 
           x = seq(0,30,1), xend = seq(0,30,1), #adjust the length of ticks as you need
           y = -189.5, yend = -185) +
  annotate("segment", 
           x = seq(1,31,4), xend = seq(1,31,4), #adjust the length of ticks as you need
           y = -195, yend = -185, linewidth = 1) +
  ggh4x::coord_axes_inside(ylim = c(-180,0), clip = 'off')


gg1


gg2 <- ggplot() +
  geom_boxplot(dat, mapping = aes(x = week, 
                                  y = DiveDur, 
                                  color = DiveType, fill = DiveType), 
               alpha = 0.75, outliers = F, coef = 1) +
  ylab('Duration (s)') +
  xlab('Week') +
  stat_boxplot(outliers = FALSE, alpha = 0.75) +
  scale_fill_manual(values = c("#d2dfe6", '#71644a', "gray90", "gray90"), labels = levels(dat$DiveType)) +
  scale_color_manual(values = c("#4a6e76", '#71644a', "gray80", "gray90"), labels = levels(dat$DiveType)) +
  
  
  
  scale_x_discrete(breaks = datelabs$week, 
                   labels= datelabs$lab) +
  
  scale_y_continuous(breaks = seq(60, 320, 20)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = 'none',
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        #text = element_text(size = 18)
  ) +
  annotate("segment", 
           x = seq(0,30,1), xend = seq(0,30,1), #adjust the length of ticks as you need
           y = 36, yend = 45) +
  annotate("segment", 
           x = seq(1,31,4), xend = seq(1,31,4), #adjust the length of ticks as you need
           y = 28, yend = 45, linewidth = 1) +
  ggh4x::coord_axes_inside(ylim = c(50,320), clip = 'off')


gg2

gg3 <- ggplot() +
  geom_boxplot(dat[dat$IDI < 1000,], mapping = aes(x = week, 
                                  y = IDI, 
                                  color = DiveType, fill = DiveType), 
               alpha = 0.75, outliers = F, coef = 1) +
  ylab('Inter-dive interval (s)') +
  stat_boxplot(outliers = FALSE, alpha = 0.5) +
  scale_fill_manual(values = c("#d2dfe6", '#71644a', "gray90", "gray90"), labels = levels(dat$DiveType)) +
  scale_color_manual(values = c("#4a6e76", '#71644a', "gray80", "gray90"), labels = levels(dat$DiveType)) +
  
  scale_x_discrete(breaks = datelabs$week, 
                   labels= datelabs$lab) +
  scale_y_continuous(breaks = seq(0,180,20)) +
  
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = 'none',
        legend.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
        #text = element_text(size = 18)
  ) +
  annotate("segment", 
           x = seq(0,30,1), xend = seq(0,30,1), #adjust the length of ticks as you need
           y = -8, yend = -4) +
  annotate("segment", 
           x = seq(1,31,4), xend = seq(1,31,4), #adjust the length of ticks as you need
           y = -12, yend = -4, linewidth = 1) +
  ggh4x::coord_axes_inside(ylim = c(0,170), clip = 'off')



gg3


gg1 = gg1 + theme(text = element_text(size = 14), plot.margin = margin(0,0,0,0))
gg2 = gg2 + theme(text = element_text(size = 14), plot.margin = margin(0,0,0,0))
gg3 = gg3 + theme(text = element_text(size = 14), plot.margin = margin(0,0,0,0))

comb = (gg1)  / 
  (gg2) / 
  (gg3) + plot_annotation(tag_level = 'a')
comb

ggsave(filename = "./plots/manuscript/Main Manuscript Figures/Figure4.png", plot = comb, 
       device = 'png', dpi = 320,
       width = 170, 
       height = 170,
       units = 'mm',
       scale = 1.2)


##################################################################################################################
##################################################################################################################
