# EIH 
# estimate dive positions with the ssmcrw
# updated: 2024-09-27
# inputs: fxn_ssm_plot_diagnostics.R
#       : Hg_2019-2023_BEHDiveRecords_QAQC.csv
#       : HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv
rm(list = ls())
gc()
setwd('~/seal_telemetry/')

library(sf)
library(dplyr)
library(aniMotum)

# source ssm diagnostic plot function
source('./READ-PSB-MoveSeals/src/fxn_ssm_plot_diagnostics.R')


# load data (locs and dive data)
dives = read.csv("./data/L1/dive/Hg_2019-2023_BEHDiveRecords_QAQC.csv", stringsAsFactors = F)

dives$start_Dive = as.POSIXct(dives$start_Dive, format = "%Y-%m-%d %H:%M:%S",tz = 'UTC')
dives$end_Dive = as.POSIXct(dives$end_Dive, format = "%Y-%m-%d %H:%M:%S",tz = 'UTC')

# get list of ptts with behavior data
ptts = unique(dives$ptt)


# Read in Locations with haulout and trip id assigned
df = read.csv("./data/L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv", stringsAsFactors = F) 
df$date = as.POSIXct(df$datetime, format = "%Y-%m-%dT%H:%M:%SZ",tz = 'UTC')
df$Year = format(df$date, "%Y")
df = df[df$id %in% ptts, ]

# subset to at sea behavior
atsea = df[!df$final_haulout,c("Year","id", "TripID", "SegID", "date", "lc", "lon", "lat", 
                                 "error.radius", "smaj", "smin", "eor")]


length(unique(atsea$id)) # 34

# subset to ptt segment ids with >50 locs for modelling
seglengths = as.data.frame(table(atsea$SegID))
qualsegs = unique(seglengths$Var1[which(seglengths$Freq >= 50)])
anilocs = atsea[atsea$SegID %in% qualsegs,]

# get the start and end of segments
grouped = split(anilocs, anilocs$SegID)

tmp = lapply(X = grouped, FUN = function(segdf){
  
  segdf$SegStart = min(segdf$date)
  segdf$SegEnd = max(segdf$date)
  return(segdf)
  
})

anilocs = do.call(rbind,tmp)

# calculate the time differences between locations for each trip
grouped = split(anilocs, anilocs$TripID)

tmp = lapply(X = grouped, FUN = function(segdf){
  
  segdf$locTimeDiff= c(0, difftime(time1 = segdf$date[2:nrow(segdf)],
                                   time2 = segdf$date[1:(nrow(segdf)-1)], units = 'hours'))
  return(segdf)
  
})

anilocs = do.call(rbind,tmp)

# clean workspace
rm(list = setdiff(ls(), c('anilocs', 'dives', 'diagplots')))

######################################### DIVE TIME MODELING ###############################################

# 1. fit a state space model to the locations from each segment
# 2. predict the locations using a vector of the mid-point of all dives within modellable segments 

##################################### step 1: time-match dives to trip segments ######################################

# for matching to regularized modeled Segments
# calculate the temporal midpoint if each dive as representative of likely bottom time
starttime = dives$start_Dive
endtime = dives$end_Dive
mid_dive_dt = as.POSIXct((as.numeric(starttime) + as.numeric(endtime)) / 2, origin = '1970-01-01', tz = 'UTC')

tsvec = data.frame(id = dives$ptt, diveID = dives$SDPairID, mid_dive_dt = mid_dive_dt)
tsvec$id = as.character(tsvec$id)
tsvec$diveID = as.character(tsvec$diveID)

dives$mid_dive_dt = mid_dive_dt

# split by SegID 
grouped = split(anilocs, anilocs$SegID)

# Iterate over each segment id in seglocsL
tmp = lapply(grouped, FUN = function(df){
  
  # subset the dives dataframe based on the unique ptt id in segdf
  divesdf = tsvec[tsvec$id == as.character(unique(df$id)),]
  
  # get the dives in this segment by filtering the dives between the segment start and end times
  dives_in_seg = divesdf[which(divesdf$mid_dive_dt <= df$SegEnd[1] & divesdf$mid_dive_dt >= df$SegStart[1]),]
  
  # if there are dives in this segment, assign segID to dives_in_seg df and store the dives_in_seg df in our storage list
  if (nrow(dives_in_seg) > 0) {
    dives_in_seg$TripID = NA
    dives_in_seg$TripID = unique(df$TripID)
    dives_in_seg$SegID = NA
    dives_in_seg$SegID = unique(df$SegID)
  }
  
  return(dives_in_seg)
  
})
  
# bind the rows of our divesinseg storage list ' tmp'
divesinseg = do.call(rbind, tmp)
rm(list = setdiff(ls(), c('anilocs', 'divesinseg', 'dives')))
head(divesinseg)
rownames(divesinseg) = NULL
rownames(anilocs) = NULL
##################################################################################################

####### step 2: fit a model using each segid and use a vector of dive times instead of a regularized time step 
divesinseg = divesinseg[order(divesinseg$SegID, divesinseg$mid_dive_dt),]

locs_for_fit = anilocs[,c('SegID', 'date', 'lc', 'lon', 'lat', 'error.radius', 'smaj', 'smin', 'eor')]
locs_for_fit = locs_for_fit[order(locs_for_fit$SegID, locs_for_fit$date),]

colnames(locs_for_fit)[1] = 'id'


tvec_df = data.frame(id = divesinseg$SegID, date = divesinseg$mid_dive_dt)

tvec_df = tvec_df %>% arrange(id, date)

locs_subset = locs_for_fit[which(locs_for_fit$id %in% unique(divesinseg$SegID)),] # 55069

locs_subset = locs_subset[order(locs_subset$id, locs_subset$date), ] # 54458

############### fit ssm crw models
set.seed(28)
vmax = 3.5 # Gallon et al. 2007
fits = fit_ssm(locs_subset, vmax = vmax, model = 'crw', time.step = tvec_df) 

############### look at convergence issues and standard error estimation issues
segfits = fits
s = as.data.frame(summary(segfits)$Stattab)

if (any(s$converged == FALSE)){ # FALSE they all converged
nonconverged = segfits$id[!segfits$converged] # empty
}

# look at standard error estimations that failed to be estimated
idsforrefit = sapply(X = 1:nrow(segfits), FUN = function(x){
  tmppar = as.data.frame(segfits$ssm[[x]]$par)
  par_mles = as.data.frame(segfits$ssm[[x]]$tmb$env$last.par.best)
  
  if (any(is.na(tmppar$`Std. Error`)) | any(is.na(par_mles))){
    idforrefit = segfits$id[x]
  }else{idforrefit = NA}
  
  return(idforrefit)
  
})

idsforrefit = idsforrefit[!is.na(idsforrefit)]

if (exists('nonconverged')){
idsforrefit = unique(c(nonconverged, idsforrefit))
}else{
  idsforrefit = idsforrefit
}
# refit these if there are idsforrefit, fit without estimating psi
if (length(idsforrefit) > 0){
  refitanilocs = locs_subset[locs_subset$id %in% idsforrefit, ]
  tv = tvec_df[tvec_df$id %in% idsforrefit, ]
  segrefits = fit_ssm(x = refitanilocs, 
                      vmax = vmax, 
                      model = 'crw', 
                      time.step = tv, 
                      map = list(psi = factor(NA)))
  

  nose = sapply(X = 1:nrow(segrefits), FUN = function(x){
    tmppar = as.data.frame(segrefits$ssm[[x]]$par)
    par_mles = as.data.frame(segrefits$ssm[[x]]$tmb$env$last.par.best)
    
    if (any(is.na(tmppar$`Std. Error`)) | any(is.na(par_mles))){
      idforrefit = segrefits$id[x]
    }else{idforrefit = NA}
    
    return(idforrefit)
    
  })
  nose = nose[!is.na(nose)]
}


save(object = fits, file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_modelobjects_fits.RData')

save(object = segrefits, file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_modelobjects_refits.RData')

####################################################################################################################################
# for dive location imputation with rejection sampling - SAVE OUT RANDOM AND FIXED EFFECT PARAMETERS 

parres = lapply(X = fits$ssm, FUN = function(x){
  
  curssm = x
  predicted = x$predicted
  
  par_mles <- x$tmb$env$last.par.best
  
  # random effects
  re <- x$tmb$env$random
  par_re <- par_mles[re]
  
  h_re <- x$tmb$env$spHess(par_mles, random = TRUE) ## conditional prec. of u | theta
  
  returnlist = list(ssm = curssm, predicted = predicted, par_mles = par_mles, par_re = par_re, h_re = h_re)
  
  return(returnlist)
  
  
})

names(parres) = fits$id

parresrefits = lapply(X = segrefits$ssm, FUN = function(x){
  
  curssm = x
  predicted = x$predicted
  
  par_mles <- x$tmb$env$last.par.best
  
  # random effects
  re <- x$tmb$env$random
  par_re <- par_mles[re]
  
  h_re <- x$tmb$env$spHess(par_mles, random = TRUE) ## conditional prec. of u | theta
  
  returnlist = list(ssm = curssm, predicted = predicted, par_mles = par_mles, par_re = par_re, h_re = h_re)
  
  return(returnlist)
  
  
})

names(parresrefits) = segrefits$id

# Get rid of the refit ids from parres 
parres = parres[-c(which(names(parres) %in% idsforrefit))]

# Add in the refitted ids
parres = c(parres, parresrefits)

# Save out list of params
save(object = parres, file = './data/L2/SSM/Hg-2019-2023-SSM-DivePosition_paramobjects.RData')
####################################################################################################################################
################################################ get predictions ###################################################################
preds = grab(fits, what = 'predicted')
predsrefit = grab(segrefits, what = 'predicted')
preds[preds$id %in% idsforrefit, ] = predsrefit

# clean workspace
rm(list = setdiff(ls(), c('segrefits','fits','divesinseg', 'dives', 'preds', 'locs_for_fit', 'idsforrefit')))


# merge back with dive info
dive_preds = preds
dive_preds$SegID = dive_preds$id
dive_preds$id = sapply(strsplit(dive_preds$SegID, '-'), "[", 1)

# get dive id back in the mix
dis = divesinseg[,c('SegID', 'diveID', 'mid_dive_dt')]
colnames(dis)[3] = 'date'
dive_preds_ = left_join(x = dive_preds, y = dis, by = join_by(SegID, date))

dives$ptt = as.character(dives$ptt)
dives$diveID = as.character(dives$SDPairID)
modeled_dives = left_join(dive_preds_, dives, by = join_by(id==ptt, diveID))

# remove some unnecessary columns
colstoeliminate = which(colnames(modeled_dives) %in% c("percentdecoded", "passes", "percentargosloc", "msgperpass", "ds", "mininterval", 
                    'tagmodel', 'notes', 'xmitstart', 'xmitend', "totxmitdays",'datastart', 'dataend'))
modeled_dives = modeled_dives[, -colstoeliminate]

modeled_dives$TripID = paste(modeled_dives$id, sapply(strsplit(x = modeled_dives$SegID, '-'), '[[', 2), sep = '-')

# reorder for neatness
modeled_dives = modeled_dives[,c(1, 28,38, 15:17, 2, 37, 3:14, 18:25)]

# write out dive positions
write.csv(modeled_dives, file = "./data/L2/SSM/Hg-2019-2023-SSM-DivePositions.csv", 
          row.names = FALSE)
# end