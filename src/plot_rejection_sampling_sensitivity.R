# Plot variance over N and mean diff of the x-coordinate
# EIH
# updated: 2024-09-30
rm(list = ls())
gc()
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

setwd("~/seal_telemetry/")

bootres = read_csv('./data/L3/dive/imputation_sensitivity_testing/bootstrapped_data.csv')

# each row in bootres represents the mean, sd, and se of x, y and bathymetry of N imputations
# for each diveID and each N (seq(10,200, 2)) there are 200 bootstrap replicates
# the following plots look at trends in those means and sds. 

# summarise the mean differences 
md_sum = bootres %>% group_by(N, diveID) %>%
  mutate(diffux = abs(mean(ux) - ux),
         diffuy = abs(mean(uy) - uy),
         diffub = abs(mean(ubathy) - ubathy)) %>% ungroup() %>%
  group_by(N) %>%
  summarise(diffUx = mean(diffux),
            diffUy = mean(diffuy),
            diffUb = mean(diffub),
            upperx = quantile(diffux, 0.75),
            lowerx = quantile(diffux, 0.25),
            uppery = quantile(diffuy, 0.75),
            lowery = quantile(diffuy, 0.25),
            upperb = quantile(diffub, 0.75),
            lowerb = quantile(diffub, 0.25))
            # iqrx = IQR(diffux, type = 0.95),
            # iqry = IQR(diffuy, type = 0.95),
            # iqrb = IQR(diffub, type = 0.95))

gx = ggplot(data = md_sum, mapping = aes(x = factor(N), y = diffUx)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowerx, ymax = upperx)) +
  labs(x = "", 
       y = expression("|" ~ bar(bar(x)) - bar(x)[i] ~ '|' ~ (m))) +
  scale_x_discrete(breaks = seq(10, 200, 10))+
  theme_minimal()

gx

gy = ggplot(data = md_sum, mapping = aes(x = factor(N), y = diffUy)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowery, ymax = uppery)) +
  labs(x = "", 
       y = expression("|" ~ bar(bar(y)) - bar(y)[i] ~ '|' ~ (m))) +
  scale_x_discrete(breaks = seq(10, 200, 10))+
  theme_minimal()

gy

gb = ggplot(data = md_sum, mapping = aes(x = factor(N), y = diffUb)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowerb, ymax = upperb)) +
  labs(x = "number of imputations", 
       y = expression("|" ~ bar(bar(bathy)) - bar(bathy)[i] ~ '|' ~ (m))) +
  scale_x_discrete(breaks = seq(10, 200, 10))+
  theme_minimal()

gb

# combine with a represenative caterpillar plot
sub = bootres[bootres$diveID == unique(bootres$diveID)[27], ]
sub$replicate = as.numeric(sub$replicate)

# Calculate average deviation from the mean
subN = sub %>% group_by(N) %>% 
  mutate(uxdiff = mean(ux) - ux, Ux = mean(ux)) %>% mutate(meanuxdiff = mean(uxdiff))

# example dive id
ymax = max(subN$Ux + subN$uxdiff) 
ymin = min(subN$Ux - subN$uxdiff) 

g <- ggplot(data = subN, mapping = aes(x = N, y = Ux, group = replicate)) +
  # each line is one of 200 bootstrap reps where the x coordinate represents the mean of N imputations
  geom_line(mapping = aes(x = N, y = ux), color = 'gray40', alpha = 0.1, linewidth = 1) +
  geom_line(linetype = 1, alpha = 1, linewidth = 2, color = 'black') +
  annotate(geom = 'text', x = 125, y = ymax-0.1, label = '200 bootstrapped datasets for one dive') +
  ylab(expression(bar(x) ~ coordinate ~ (km))) +
  xlab('number of imputations') +
  scale_x_continuous(breaks = seq(10,200,10)) +
  ylim(c(ymin,ymax)) +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(expand = F)
g

library(patchwork)
design = 'AB
          CD'


fig = wrap_plots(A = gx + theme(xlab(element_blank())), B = gy + theme(xlab(element_blank())), 
                 C = gb, D = g, 
                     design = design, 
                     heights = c(1,1), widths = c(1,1)) + plot_annotation(tag_levels = "a")
fig
ggsave(filename = './plots/manuscript/FigureS3.png',
       plot = fig, device = 'png', scale = 1.5, width = 170, height = 170, units = 'mm')
# for each dive id plot the difference from the mean of all replicates in N group.
samp = 1:length(unique(bootres$diveID))

for (i in samp){
  sub = bootres[bootres$diveID == unique(bootres$diveID)[i], ]
  sub$replicate = as.numeric(sub$replicate)
  
  # Calculate average deviation from the mean
  subN = sub %>% group_by(N) %>% 
    mutate(uxdiff = abs(mean(ux) - ux), Ux = mean(ux)) %>% mutate(meanuxdiff = mean(uxdiff))
  
  # MEAN DIFF
  ymax = max(subN$uxdiff)
  ymin = 0
  g <- ggplot(data = subN, mapping = aes(x = N, y = uxdiff, group = replicate)) +
    geom_line(linetype = 1, alpha = 0.1, linewidth = 1, color = 'gray30') +
    geom_line(data = subN, aes(x = N, y = meanuxdiff), color = "black", linewidth = 1, linetype = 1) +
    #scale_color_gradient2(midpoint = median(sub$replicate), low = "darkblue", mid = "green", high = "lightyellow") +
    annotate(geom = 'text', x = 125, y = ymax, label = paste0('200 Bootstrap Replicates for Dive ID: ', unique(sub$diveID))) +
    ylab(expression('|'~bar(x)[i] - bar(bar(x)) ~ '|' ~ (km))) +
    ylim(c(ymin,ymax)) +
    theme_bw() +
    theme(legend.position = 'none')
  g
  ggsave(filename = sprintf('./plots/rejection_sampling/difffrommean/RejectionSamplingMinNTest_xi-ux_centered0_%d.jpeg', i),
         plot = g, scale = 0.75, dpi = 300, device = 'jpeg')
  
  # MEAN +/- diffmean
  ymax = max(subN$ux) 
  ymin = min(subN$ux) 
  
  g <- ggplot(data = subN, mapping = aes(x = N, y = Ux, group = replicate)) +
    # each line is one of 200 bootstrap reps where the x coordinate represents the mean of N imputations
    geom_line(mapping = aes(x = N, y = ux), color = 'gray40', alpha = 0.1, linewidth = 1) +
    geom_line(linetype = 1, alpha = 1, linewidth = 2, color = 'black') +
    annotate(geom = 'text', x = 125, y = ymax, label = paste0('200 bootstrap replicates for dive ID: ', unique(sub$diveID))) +
    ylab(expression(bar(x) ~ coordinate ~ (km))) +
    xlab('number of imputations') +
    scale_x_continuous(breaks = seq(10,200,10)) +
    ylim(c(ymin,ymax)) +
    theme_bw() +
    theme(legend.position = 'none')
  g
  ggsave(filename = sprintf('./plots/rejection_sampling/difffrommean/RejectionSamplingMinNTest_xi-ux_centeredUX_%d.jpeg', i),
         plot = g, scale = 0.75, dpi = 300, device = 'jpeg')
  
  
} # end loop
