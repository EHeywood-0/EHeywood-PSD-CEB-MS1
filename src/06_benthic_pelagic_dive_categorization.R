# CATEGORIZING DIVES AS BENTHIC OR PELAGIC
# EIH
# 2025-01-03

# Purpose: Read in data from the dive position imputation and rejection sampling
# Calculate the proportion of water column use
# Grouping by diveID, calculate mean lat, lon, and bathymetry - get se around bathymetry, and prop water column dive
# Categorize as benthic / pelagic
# Plot
# MAYBE extract other features to imputations and see if trend... sediment type... 

################################################################################################
################################################# SETUP ########################################
################################################################################################
rm(list = ls())
gc()
setwd('~/seal_telemetry')

# Data manipulation
library(readr)
library(dplyr)
library(tidyr)

# For querying bathymetry data: NOAA National Centers for Environmental Information (2022) ETOPO
# 2022 15 Arc-Second Global Relief Model. NOAA National Centers for Environmental Information. URL https://doi.org/10.25921/fd45-gt74
library(marmap)

# Spatial Libraries
library(sf)
library(stars)

# Plot
library(ggplot2)

################################################################################################
################################################# LOAD DATA ####################################
################################################################################################
# fnames = list.files(path = './data/L3/dive/imputeddivepos', 
#                     pattern = '*.csv$', all.files = TRUE, 
#                     recursive = TRUE, full.names = TRUE)
# 
# data = purrr::map_dfr(.x = fnames, .f = function(.x){
#   
#   tmp = read_csv(.x) %>% group_by(diveID) %>% 
#     mutate(Nimputes = n()) %>%
#     rename(bathy_m = bathydepth)
#   return(tmp)
#   
# })

data = read_csv("./data/L3/dive/Hg_2019-2023_DivePosEstimates_Imputations_BestEst.csv")

data = data %>%  
  mutate(totaldivespossible = n_dives_lost + n_dives_imputed) %>%
  mutate(percloss = n_dives_lost / totaldivespossible)

mean(data$percloss) # 4% average loss

# Subset to get those dive positions which have more than 50 imputations
data_ = data[data$Nimputes > 50, ]
totallost = data_ %>% group_by(id) %>% 
  summarise(AbsPercLoss = (totaldivespossible[1] - length(unique(diveID))) / totaldivespossible[1],
            NRetained = length(unique(diveID)), NPossible = totaldivespossible[1])
mean(totallost$AbsPercLoss)

################################################################################################
################################################# CALCULATE METRICS ############################
################################################################################################

dives = data_ %>% mutate(PercWaterCol = Depth / bathy_m) %>%
  
  group_by(id, diveID, date, Depth, DiveDur, IDI, Nimputes) %>% 
  summarise(meanbathy_m = mean(bathy_m),
            sebathy_m = sd(bathy_m) / Nimputes[1],
            
            meanPWC = mean(PercWaterCol),
            PW25 = quantile(PercWaterCol, 0.25),
            PW75 = quantile(PercWaterCol, 0.75),
            PW50 = quantile(PercWaterCol, 0.5),

            sePWC = sd(PercWaterCol) / Nimputes[1],
            # PWCLowerQuant = quantile(PercWaterCol, probs = 0.025),
            # PWCUpperQuant = quantile(PercWaterCol, probs = 0.975),
            
            meanX_km = mean(x),
            seX_km = sd(x) / Nimputes[1],
            meanY_km = mean(y),
            seY_km = sd(y) / Nimputes[1]) %>% 
  ungroup() %>%
  mutate(LowerCI = meanPWC - 1.96*sePWC, 
         UpperCI = meanPWC + 1.96*sePWC)
dives = dives[which(dives$PW50 < 1.1), ]
# USE PROBABILITY BENTHIC
hist(dives$PW50, breaks = 100, xlim = c(0,1.2))
hist(dives$PW25, breaks = 100, xlim = c(0,1.2))
hist(dives$PW75, breaks = 100, xlim = c(0,1.2))

hist(dives$PW50, breaks = 100, xlim = c(0,1.2))
abline(v = 0.875, col = 'red')

# BIMODAL - fit mixture model to determine the two underlying distributions

library(mixR)
x = dives$PW50

# xbeta = x
# xbeta[xbeta >1] = 1
# dives$x = dives$PW50
# dives$x[dives$x>1]=1
# oneset = dives[dives$ptt == dives$ptt[1],]
# betamm = betareg::betamix(formula = x ~ 1 | 1, data = oneset, k = 2)

# fit diff dist families
lnmm <- mixR::mixfit(x = x, ncomp = 2, family = 'lnorm')
wmm <- mixR::mixfit(x = x, ncomp = 2, family = 'weibull')
gmm <- mixR::mixfit(x = x, ncomp = 2, family = 'gamma')
nmm <- mixR::mixfit(x = x, ncomp = 2, family = 'normal')

# plot
p1 <- plot(nmm, title = 'Gaussian Mixture k=2', 
           legend.position = 'none', trans = 0.4) + 
  scale_fill_manual(values = c("cornflowerblue",'darkblue')) +
  annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(nmm$bic), 0)) +
  xlab('') +
  coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))
p2 <- plot(lnmm, title = 'Lognormal Mixture k=2', 
           legend.position = 'none', trans = 0.4)+ 
  scale_fill_manual(values = c("cornflowerblue",'darkblue'))+
  annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(lnmm$bic), 0)) +
  xlab('') +
  coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))

p3 <- plot(wmm, title = 'Weibull Mixture k=2', 
           legend.position = 'none', 
           xlab = '% Water column reached at max dive depth', trans = 0.4)+ 
  scale_fill_manual(values = c("cornflowerblue",'darkblue'))+
  annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(wmm$bic), 0)) +
  coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))

p4 <- plot(gmm, title = 'Gamma Mixture k=2', 
           legend.position = 'none', 
           xlab = '% Water column reached at max dive depth', trans = 0.4)+ 
  scale_fill_manual(values = c("cornflowerblue",'darkblue'))+
  annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(gmm$bic), 0)) +
  coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))

comb = gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2) 

ggsave(filename = './plots/manuscript/AnimalBiotelemetrySubmission/SupplementalMixtureModeling.png', 
       plot = comb, 
       device = 'png', dpi = 320, width = 170, height = 150, units = 'mm', scale = 1.25)

ggsave(filename = './plots/manuscript/AnimalBiotelemetrySubmission/SupplementalMixtureModeling.pdf', 
       plot = comb, 
       device = 'pdf', width = 1000, height = 900, units = 'px', scale = 3)

# # PLOT K= 3
# wmm3 <- mixR::mixfit(x = x, ncomp = 4, family = 'weibull')
# gmm3 <- mixR::mixfit(x = x, ncomp = 3, family = 'gamma')
# plot(wmm3)
# 
# p43 <- plot(gmm3, title = 'Gamma Mixture k=3', 
#            legend.position = 'none', 
#            xlab = '% Water column reached at max dive depth', trans = 0.4)+ 
#   scale_fill_manual(values = c("cornflowerblue",'red','darkblue'))+
#   annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(gmm3$bic), 0)) +
#   coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))
# p43
# p33 <- plot(wmm3, title = 'Weibull Mixture k=3', 
#             legend.position = 'none', 
#             xlab = '% Water column reached at max dive depth', trans = 0.4)+ 
#   scale_fill_manual(values = c("cornflowerblue",'red','darkblue'))+
#   annotate(geom = 'text', x = 0.2, y = 4,label = paste0('BIC = ', round(wmm3$bic), 0)) +
#   coord_cartesian(xlim = c(0,1.2), ylim = c(0,5))
# p33
# # plot
# 
# gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2) 


# Model selection
# select_gamma = select(x, ncomp = 1:3, family = 'gamma')
# plot(select_gamma)
# select_weibull = select(x, ncomp = 1:3, family = 'weibull')
# plot(select_weibull)


# SOLVE INTERSECTION
# Define the PDFs of the two Weibull distributions
weibull1 <- function(x, k1, lambda1) {
  (k1 / lambda1) * (x / lambda1)^(k1 - 1) * exp(-(x / lambda1)^k1)
}

weibull2 <- function(x, k2, lambda2) {
  (k2 / lambda2) * (x / lambda2)^(k2 - 1) * exp(-(x / lambda2)^k2)
}

# Define the parameters for the two Weibull distributions
k1 <- wmm$k[1]
lambda1 <- wmm$lambda[1]
k2 <- wmm$k[2]
lambda2 <- wmm$lambda[2]

# Define the function for the difference between the two PDFs
diff_weibull <- function(x) {
  weibull1(x, k1, lambda1) - weibull2(x, k2, lambda2)
}

# Plot the PDFs to visualize where they intersect
x_vals <- seq(0, 1.2, length.out = 1000)
y1 <- weibull1(x_vals, k1, lambda1)
y2 <- weibull2(x_vals, k2, lambda2)

plot(x_vals, y1, type = "l", col = "blue", ylim = c(0, max(y1, y2)), ylab = "Density", xlab = "x", main = "Intersection of Two Weibull Distributions")
lines(x_vals, y2, col = "red")
legend("topright", legend = c("Weibull 1", "Weibull 2"), col = c("blue", "red"), lty = 1)

# Find the intersection points using uniroot
# You might need to specify multiple intervals if there are multiple intersections
interval <- c(0.7,1)

intersection_points <- uniroot(diff_weibull, interval)$root

# Print the intersection points
print(intersection_points)

# Add the intersection points to the plot
points(intersection_points, weibull1(intersection_points, k1, lambda1), pch = 19, col = "green")

intersxn = round(intersection_points, digits = 2)

p3 <- plot(wmm, 
           legend.position = 'none', title = NULL,
           xlab = 'Max Dive Depth Percent of Water Column', trans = 0.4)+ 
  labs(subtitle = 'Weibull Mixture Model k=2')+
  scale_fill_manual(values = c("cornflowerblue",'darkblue')) +
  annotate(geom = 'segment', x = 0.87, xend = 0.87,
           y = 0,yend = 1, linewidth = 1, color = 'tomato', linetype = 1) +
  annotate(geom = 'text', x = intersection_points-0.05, y = 1.5, label = '~0.87') +
  coord_cartesian(xlim = c(0,1.1), ylim = c(0,5))
  
p3
#######################################################################################################
#######################################################################################################
#######################################################################################################

# USING ESTABLISHED DATA DRIVEN CUTOFF GET SOME DIVE CHARACTERISTICS
dives = data_ %>% mutate(PercWaterCol = Depth / bathy_m) %>%
  
  group_by(id, diveID, date, Depth, DiveDur, IDI, Nimputes) %>% 
  summarise(meanbathy_m = mean(bathy_m),
            sebathy_m = sd(bathy_m) / Nimputes[1],
            
            meanPWC = mean(PercWaterCol),
            medianPWC = median(PercWaterCol),
            PWC25 = quantile(PercWaterCol, 0.25),
            PWC75 = quantile(PercWaterCol, 0.75),
            sePWC = sd(PercWaterCol) / Nimputes[1],

            meanX_km = mean(x),
            seX_km = sd(x) / Nimputes[1],
            meanY_km = mean(y),
            seY_km = sd(y) / Nimputes[1],
            PWCProb_Benthic = length(which(PercWaterCol >= 0.87)) / n(),) %>% 
  ungroup() %>%
  mutate(LowerCI = meanPWC - 1.96*sePWC, 
         UpperCI = meanPWC + 1.96*sePWC)


# Anything that is over 1.10 perc water column is a very shallow water.. and due to tolerance set at 
#hist(dives$Depth[dives$meanPWC > 1.08])

#length(which(dives$meanPWC > 1.10))/nrow(dives) # ~700 dives are 10% or greater deeper than extracted bathymetry... 




dives$ptt = sapply(strsplit(dives$id, '-'), '[[',1)


meta = read.csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv", stringsAsFactors = FALSE)
meta = meta[,c('ptt', 'masskg', 'sex')]

data6 = merge(dives, meta, by = 'ptt')

write.csv(x = data6, file = './data/L3/dive/Hg_2019-2023_DiveTypeClassification.csv', row.names = FALSE)





