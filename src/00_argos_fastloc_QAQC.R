# EIH
# updated 2024-09-27
# Location and FAST GPS Pre-PROCESSING
# inputs: # wildlife computers portal downloads (L0) Argos and Fastloc data files 
#       : the output of 0000_process_wc_summary_files.R for the data cutoff dates

# setup
rm(list = ls())
gc()
library(readr)
library(tidyr)
library(dplyr)

getOption(x = 'digits.secs', default = 6)
setwd('~/seal_telemetry')

# source some helper functions
source("./READ-PSB-MoveSeals/src/fxns_wildlifecomputers_readdata.R")

# Will subset all data >= Deployment Date and <= Final Transmission

# Load Meta Data
meta = read_csv(file = "./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv")


# Get rid of bad ptts 177509 (no locations), 240183 (rehab), 240186 (tag malfunction), 
keeps = setdiff(unique(meta$ptt), c(177509, 240186, 240183))
meta = meta[meta$ptt %in% keeps, ]


# list fastloc files

# List fastgps files
files = list.files("./data/L0/", pattern = '*-1-FastGPS.csv', full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
filedf = data.frame(fname = sapply(strsplit(files, '/'), '[[', 6), ptt = sapply(strsplit(files, '/'), '[[', 5))
# Filter by keeps
filedf = filedf[filedf$ptt %in% as.character(keeps), ]

fastlocs = WC_read_fastloc(files = filedf$fname, datadir = './data/L0/')
fastlocs$ptt = as.numeric(fastlocs$ptt)

# List all location data files
files = list.files("./data/L0/", pattern = '*-Locations.csv', recursive = TRUE, full.names = TRUE)

ptts = as.numeric(sapply(strsplit(x = files, split = '/'), "[[", 5))
length(unique(ptts))
fnames = sapply(strsplit(x = files, split = '/'), "[[", 6)
filedf = data.frame(ptt = ptts, fnames = fnames, fullnames = files)
filedf = filedf[as.numeric(filedf$ptt) %in% keeps, ]

# If there is both a locations file and a [single digit]-Locations.csv file... keep the 1-Locations.csv file
fastlocFiles = filedf[grepl(pattern = '-\\d-Locations.csv', x = filedf$fnames), ]
argosonly = filedf[-c(which(filedf$ptt %in% fastlocFiles$ptt)), ]

fs = rbind(fastlocFiles, argosonly)
fs = fs[fs$fnames != '142351-Locations.csv', ]

# LOAD 61 TAGS - 2 failed to transmit
locs = WC_read_locs(files = fs$fullnames) # 142115
# Checks
length(unique(locs$ptt))
which(!fs$ptt %in% unique(locs$ptt))

############################################################################################
# MERGE BY TIME AND PTT TO GET FASTLOC ANCILLARY DATA
# Time match to merge in sat and GPS error data
storageL = list()
fastloc_ptts = fastlocFiles$ptt
for (i in 1:length(fastloc_ptts)){
  id = fastloc_ptts[i]
  tm1 = locs[locs$ptt == id, ]
  tm1 = tm1[order(tm1$datetime),]
  tm2 = fastlocs[fastlocs$ptt == id, ]
  tm2 = tm2[order(tm2$datetime), ]
  tm2$RowIndex = 1:nrow(tm2)
  tm1$timematch_idx = findInterval(x = tm1$datetime, vec = c(-Inf, head(tm2$datetime, -1)) + c(0, diff(tm2$datetime)/2))
  tm1$timematch_idx[tm1$type == 'Argos'] = NA
  
  colnames(tm2)[9] = 'fastloc_datetime'
  tm2 = tm2[,c('ptt', 'fastloc_datetime','hauled.out', 'satellites', 'residual', 'time.error', 'RowIndex')]
  
  merger = left_join(tm1, tm2, by = join_by(ptt, timematch_idx == RowIndex))
  
  storageL[[i]] = merger
  
}

merged_fltags = do.call(rbind, storageL)

argosdata = locs[!locs$ptt %in% fastloc_ptts, ]

# Recombine the merged fastloc data and the argos only tags
locs_ = bind_rows(argosdata, merged_fltags)

##### ELIMINATE PRE DEPLOYMENT DATA AND POST-DEPLOYMENT DATA 

# merge with meta to get the transmission cutoff dates 

locs2 = left_join(locs_, meta, by = 'ptt') %>% 

  filter(datetime >= cutoffstart, datetime <= cutoffend)%>%
  rename(lat=latitude,lon=longitude,smaj=error.semi.major.axis, smin=error.semi.minor.axis,
         lc=quality, eor= error.ellipse.orientation,id=ptt)


# eliminated 1379 locations left on tags pre-deployment or post-deployment transmissions that are not on animal 
nrow(locs) - nrow(locs2)

# write out data
write_csv(x = locs2, file = './data/L1/locs/Hg_2019-2023_CombinedDeploymentLocs.csv')

# clean workspace
rm(list = setdiff(ls(), 'locs2'))

# Eliminate FastGPS positions that have both bad sats (fewer than 6) and high residuals (>= 30) (Dujon et al)
badgps = which(locs2$residual > 30 & locs2$satellites < 6)

if(length(badgps) > 0 ){
  locs3 = locs2[-c(badgps), ]
}else{locs3 = locs2}


#Impossible positions no matter what LC class
locs4 = locs3[which(locs3$lon > -77 & locs3$lon <= -54 & locs3$lat >= 35.0),] 

# Calculate percentage location classes
table(locs4$lc, useNA = 'ifany')
round(table(locs2$lc, useNA = 'ifany') / nrow(locs4) * 100)

# Remove Zs
locs4 = locs4[locs4$lc != 'Z', ]
table(locs4$lc, useNA = 'ifany')
round(table(locs4$lc, useNA = 'ifany') / nrow(locs4) * 100)

# Get rid of NA
locs4 = locs4[!is.na(locs4$lc), ]

# remove duplicates
locs4 = locs4 %>% arrange(id, datetime, desc(error.radius))

dupes = which(duplicated(locs4[,c('id', 'datetime')]))

locs5 = locs4[-dupes,]
table(locs5$lc, useNA = 'ifany')
length(unique(locs5$id))

# Remove those with comments indicating erroneous location
locs6 = locs5[is.na(locs5$comment), ]

table(locs6$type)

# GET NUMBERS FOR PRE WINDFARM CONTSTRUCTION PERIOD
nrow(locs6[locs6$datetime < as.POSIXct(x = '2023-06-01 00:00:00', tz = 'UTC'),])
table(locs6$type[locs6$datetime < as.POSIXct(x = '2023-06-01 00:00:00', tz = 'UTC')])
length(unique(locs6$id))
# write out filtered data
write.csv(x = locs6, file = './data/L1/locs/Hg_2019-2023_Prefiltered_Locs.csv', row.names = F)

# SCRATCH
# A PRIORI SPEED FILTER 

# APRIORI SPEED FILTER FROM ARGOS FILTER 
# data = locs3
# 
# ll = split(data, data$id)
# 
# resL = lapply(X = ll, FUN = function(x, vmax = 10){
#   
#   Gs = x[x$lc == 'G', ]
#   tmp = x[x$lc != 'G',]
#   tmp = tmp[order(tmp$date), ]
#   filt = argosfilter::sdafilter(lat = tmp$lat, 
#                                 lon = tmp$lon, 
#                                 dtime = tmp$datetime, 
#                                 lc = tmp$lc, 
#                                 vmax = vmax, ang = -1)
#   
#   dat = tmp[which(filt %in% c('not', 'end_location')), ]
#   print('n removed:')
#   print(nrow(tmp) - nrow(dat))
#   fin = rbind(dat, Gs)
#   fin = fin[order(fin$datetime), ]
#   
#   return(fin)
#   rm(tmp, dat, fin)
#   
# })

#locs_filt = bind_rows(resL)
#rm(ll, data)
