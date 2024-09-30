# PLOT UD MINIMUM SAMPLE SIZE 
setwd("~/seal_telemetry/")

bsize = 18

library(patchwork)
library(ggplot2)
library(SDLfilter)

# LOAD FEMALE
load('./data/L3/Hg2019-2023-OverlapProbability100PVC_BRB1kmres_Female.RData')

aF <- asymptote(overlapF, upper.degree = 5, estimator = 'glm', family = binomial)

fo = ggplot(data = aF$results, aes(x = x))+
  geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'gray70', linewidth = 1) + 
  geom_point(aes(y = y), size = 2) + 
  annotate("segment", x = aF$min.n, xend = aF$min.n, y = 0.3, yend = aF$h.asymptote*0.95, color = "gray10", linetype = 2,linewidth = 1) +
  geom_hline(yintercept = aF$h.asymptote*0.95, linetype = 2, color = 'gray40') +
 
  scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30)) +
  scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
  
  # Set plot attributes
  theme_light(base_size = bsize) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

fo

# LOAD MALE
load('./data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Male.RData')
aM <- asymptote(data = overlapM, upper.degree = 10, estimator = 'glm', family = binomial)
aM
mo = ggplot(data = aM$results, aes(x = x))+
  geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'gray70', linewidth = 1) + 
  
  geom_point(aes(y = y), size = 2) + 
  annotate("segment", x = aM$min.n, xend = aM$min.n, y = 0.3, yend = aM$h.asymptote*0.95, color = "gray10", linetype = 2,linewidth = 1) +
  geom_hline(yintercept = aM$h.asymptote*0.95, linetype = 2, color = 'gray40') +
  #geom_hline(yintercept = aM$h.asymptote, linetype = 3, color = 'black') +
  scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30), name = "N tracked") +
  scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
  theme_light(base_size = bsize)
mo

# LOAD ALL
load(file = './data/L3/HgPre-Construction_OverlapProbability100PVC_BRB1kmres-ALL.RData')

a  <- SDLfilter::asymptote(data = overlap, upper.degree = 10, estimator = 'glm', family = binomial)

g = ggplot(data = a$results, aes(x = x))+
  geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'gray70', linewidth = 1) + 
  
  geom_point(aes(y = y), size = 2) + 
  annotate("segment", x = a$min.n, xend = a$min.n, y = 0.3, yend = a$h.asymptote*0.95, color = "gray10", linetype = 2,linewidth = 1) +
  geom_hline(yintercept = a$h.asymptote*0.95, linetype = 2, color = 'gray40') +
  scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30)) +
  scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
  theme_light(base_size = bsize) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) 
g

# Combined plot of all, female, and male
combined_p = g/fo/mo +
  plot_annotation(tag_levels = 'a')
combined_p
ggsave(filename = "./plots/manuscript/Supplemental Figures/FigureS3.png", 
       plot = combined_p, device = 'png', dpi = 320, 
       width = 170, height = 170, units = 'mm',scale = 1.25)

ggsave(filename = "./plots/manuscript/Supplemental Figures/FigureS3.pdf", 
       plot = combined_p, device = 'pdf', 
       width = 1200, height = 1200, units = 'px', scale = 2)

########################################## SCRATCH SEASONAL ##################################################
# Seasonal break outs failed to estimate asymptote
# # Winter
# load('./data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Winter.RData')
# 
# a  <- SDLfilter::asymptote(data = Woverlap, upper.degree = 5, estimator = 'glm', family = binomial)
# 
# w = ggplot(data = a$results, aes(x = x))+
#   geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'darkgrey', linewidth = 1) + 
#   
#   geom_point(aes(y = y), size = 2) + 
#   annotate("segment", x = a$min.n, xend = a$min.n, y = 0.3, yend = a$h.asymptote*0.95, color = "darkgreen", linetype = 1,linewidth = 1) +
#   geom_hline(yintercept = a$h.asymptote*0.95, linetype = 2, color = 'darkgray') +
#   scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30)) +
#   scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
#   theme_light(base_size = bsize) +
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank()) 
# w
# 
# # Spring
# load('./data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Spring.RData')
# 
# a  <- SDLfilter::asymptote(data = Woverlap, upper.degree = 150, estimator = 'glm', family = binomial)
# 
# sp = ggplot(data = a$results, aes(x = x))+
#   geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'darkgrey', linewidth = 1) + 
#   
#   geom_point(aes(y = y), size = 2) + 
#   annotate("segment", x = a$min.n, xend = a$min.n, y = 0.3, yend = a$h.asymptote*0.95, color = "darkgreen", linetype = 1,linewidth = 1) +
#   geom_hline(yintercept = a$h.asymptote*0.95, linetype = 2, color = 'darkgray') +
#   scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30)) +
#   scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
#   theme_light(base_size = bsize) +
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank()) 
# sp
# 
# # Spring
# load('./data/L3/Hg_Pre-Construction-OverlapProbability100PVC_BRB1kmres_Summer.RData')
# 
# a  <- SDLfilter::asymptote(data = Woverlap, upper.degree = 50, estimator = 'glm', family = binomial)
# 
# su = ggplot(data = a$results, aes(x = x))+
#   geom_errorbar(aes(ymin = y_lwr, ymax = y_upr), width = 0.4, colour = 'darkgrey', linewidth = 1) + 
#   
#   geom_point(aes(y = y), size = 2) + 
#   annotate("segment", x = a$min.n, xend = a$min.n, y = 0.3, yend = a$h.asymptote*0.95, color = "darkgreen", linetype = 1,linewidth = 1) +
#   geom_hline(yintercept = a$h.asymptote*0.95, linetype = 2, color = 'darkgray') +
#   scale_x_continuous(breaks = seq(3, 30, 2), limits = c(2,30)) +
#   scale_y_continuous(limits = c(0.3,1), name = "Overlap probability") +
#   theme_light(base_size = bsize) +
#   theme(axis.text.x = element_blank(), axis.title.x = element_blank()) 
# su


