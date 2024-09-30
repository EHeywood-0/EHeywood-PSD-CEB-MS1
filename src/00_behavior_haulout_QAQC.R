# HAULOUT BEHAVIOR CSV PRE-PROCESSING
# EIH
# inputs: wildlifecomputers_readdata_fxns.R; *HaulOut.csvs
# updated: 2024-09-27
rm(list = ls())
gc()
setwd('~/seal_telemetry/')

library(readr)
library(dplyr)
library(tidyr)
ho_dur_max = 24*7

source('./READ-PSB-MoveSeals/src/fxns_wildlifecomputers_readdata.R')

# Read in meta data on deployments
meta = read_csv(file = "./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv")
colnames(meta) = tolower(colnames(meta))

# Get rid of bad ptts 177509 (no locations), 240183 (rehab), 240186 (tag malfunction), 
keeps = setdiff(unique(meta$ptt), c(177509, 240186, 240183))
meta = meta[meta$ptt %in% keeps, ]

files = list.files(path = "./", pattern = '*HaulOut.csv', full.names = T, recursive = T)
ho = WC_read_haulout(filenames = files)

# Filter out data occurring outside deployment
ho = left_join(ho, meta, by = join_by(ptt == ptt)) %>% filter(startmin >= cutoffstart, endmax <= cutoffend)

# Per discussion on Oct 13, 2023 - remove those rows where start min and max do not equal or end min and max do not equal
ho2 = ho[which(ho$startmin == ho$startmax),]
ho3 = ho2[which(ho2$endmax == ho2$endmin), ]

# get rid of negative durations and ones that don't make any biological sense (haulout periods for > 7 days)
ho3$duration_hour = as.numeric(difftime(ho3$endmax, ho3$startmax, units = 'hours'))
haulouts = ho3[ho3$duration_hour > 0 & ho3$duration_hour <= ho_dur_max, ]

write.csv(x = haulouts, file = './data/L1/haulout/Hg_2019-2023_Beh_Haulout_Events.csv', row.names = FALSE)
