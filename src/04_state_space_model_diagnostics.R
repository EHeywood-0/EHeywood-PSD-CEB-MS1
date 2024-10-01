# MODEL DIAGNOSTICS
# EIH
# 2024-05-17

rm(list = ls())
gc()
setwd('~/seal_telemetry/')
library(parallel)
library(dplyr)
library(tidyr)
source("./READ-PSB-MoveSeals/src/fxns_helper.R")

####################################### DIAGNOSTICS ##########################################
##############################################################################################


### SSMCRW - diags for 2 hour regularized estimated positions
# computationally expensive due to one-step ahead residuals
# 524 trip segmets
load(file = './data/L2/SSM/Hg-2019-2023-SSM-ModelObjects.RData')

# THIS TAKES FOREVER>>>> 524 trip segments... 
mclapply(X = 1:nrow(segfits), FUN = function(x){
  plot_ssm_diags(segfits[x,], writedir = './plots/ssm_diags/', model = 'crw')
})




  
load(file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_modelobjects_refits.RData')

mclapply(X = 1:nrow(segrefits), FUN = function(x){
  plot_ssm_diags(segrefits[x, ], writedir = './plots/ssm_validation/ssm_dive_diags/', model = 'rw')
})

load(file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_modelobjects_fits.RData')
# plts = list.files("./plots/ssm_validation/ssm_dive_diags")
# plts = trimws(plts, which = 'both')
# tripnames = sapply(strsplit(plts, ' '), '[[', 1)
# plts
# 
# plts[duplicated(tripnames)]
# missing = setdiff(fits$id, tripnames)
# which(fits$id %in% missing)
mclapply(X = 1:nrow(fits), FUN = function(x){
  plot_ssm_diags(fits[x, ], writedir = './plots/ssm_validation/ssm_dive_diags/', model = 'crw')
})


#######################################################################################################
#######################################################################################################
############################################# SCRATCH #################################################
#######################################################################################################
#######################################################################################################

evalimages = list.files('./plots/ssm_dive_diags/', recursive = T)
ptt = trimws(sapply(strsplit(evalimages, '-'), '[[', 1))
TripID = paste(ptt, sapply(strsplit(evalimages, '-'), '[[',2), sep = '-')

eval_df = data.frame(ptt = ptt, TripID = TripID, Image = evalimages,
                     Keep = NA, Notes = NA)

write_csv(x = eval_df, file = './data/L2/SSM/MPM_Eval_Worksheet.csv')
