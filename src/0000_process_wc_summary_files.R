# Process wildlife computers transmission summary files to get some useful metadata 
# EIH
# updated: 2024-09-27
# create date cutoffs for data based on start and end of data transmission records 
# inputs: wildlife computer summary files, '*-Summary.csv' files, '*.All.csv files, and location files
#         Tag_Deployments_21719_summaries.csv
# performs some cleaning and joining with meta data to produce cleaned summary and meta data records
# WILDLIFE COMPUTER DEFINITIONS:
# https://static.wildlifecomputers.com/2019/05/10152339/Spreadsheet-File-Descriptions.pdf

# But Briefly, I will describe relevant columns for subsetting DEPLOYMENT records 
# EarliestXmitTime & LatestXmitTime: the time of the first and last received Argos transmission (should be the initial and final record in ALL)
######## This be inaccurate if the tag was turned on for testing prior to deployment. 
# As such, to fix this issue, we will subset ALL.csv files to all Argos records day of or after deployment date and occurring in the same deploy year. Then take the min and max of this

# EarliestDataTime: The time stamp of the first data point, including status messages. This should be close to the deployment date.
# This can be error prone as well as is the case with Robs tags where he does not clear data. To remedy, I will pull the first location from the locations file that occurs during or after deployment day
# LatestDataTime: The time stamp of the last data point, excluding status messages
# This can be error prone due to redeployments in subsequent years etc. As such I will replace these errors with the final location timestamp occurring in the year of the target animal deployment.
rm(list = ls())
gc()
setwd("~/seal_telemetry/")
library(readr)
library(tidyr)
library(dplyr)
library(lubridate)


# load Summary files (N=63)
files = list.files(path = './data/L0/', pattern = '-Summary.csv', recursive = TRUE, full.names = TRUE, all.files = TRUE)

trx = lapply(X = files, FUN = function(x){
  tmp = read.csv(x)
  tmp$DeployID = as.character(tmp$DeployID)
  tmp$EarliestDataTime = as.POSIXct(tmp$EarliestDataTime, format = "%H:%M:%S %d-%b-%Y", tz = 'UTC')
  tmp$EarliestXmitTime = as.POSIXct(tmp$EarliestXmitTime, format = "%H:%M:%S %d-%b-%Y", tz = 'UTC')
  tmp$LatestDataTime = as.POSIXct(tmp$LatestDataTime, format = "%H:%M:%S %d-%b-%Y", tz = 'UTC')
  tmp$LatestXmitTime = as.POSIXct(tmp$LatestXmitTime, format = "%H:%M:%S %d-%b-%Y", tz = 'UTC')
  
  tmp = tmp[,colSums(is.na(tmp)) == 0]
  return(tmp)
  
})

# bind into dataframe & select relevant columns
Xmits = dplyr::bind_rows(trx) %>% dplyr::select(-c(SW, DeployDate, ReleaseDate, ReleaseType))

# read in raw meta data
meta = read_csv("./data/meta/Tag_Deployments_21719_summaries.csv") %>% dplyr::select(PTT, deploydate)
meta$deploydate = as.POSIXct(meta$deploydate, format = '%m/%d/%Y', tz = 'UTC')

# join Xmits from Summary files with meta to get deploydate
Xmits = left_join(Xmits, meta, by = join_by(Ptt == PTT))


# LOAD ALL FILES (Contains all attempted Argos transmits) - this is where the XmitTimes come from (the time of the first and last received argos messages)
files = list.files(path = './data/L0/', pattern = '-All.csv', recursive = TRUE, full.names = TRUE, all.files = TRUE)
all = lapply(X = files, FUN = function(x){
  tmp = read_csv(x, show_col_types = FALSE) %>% dplyr::select(`Platform ID No.`, `Msg Date`)
  tmp$`Msg Date` = as.POSIXct(tmp$`Msg Date`, format = '%m/%d/%Y %H:%M:%S', tz = 'UTC')
  return(tmp)
  
})

#
all = bind_rows(all) %>% left_join(., meta, by = join_by(`Platform ID No.` == PTT))

# Filter all messages so that they occur within the year of deployment and occur after or eqaual to calendar day deployment
firstlast = all %>% filter(as.Date(`Msg Date`) >= as.Date(deploydate), year(deploydate) == year(`Msg Date`)) %>% group_by(`Platform ID No.`) %>% summarise(FirstAllTime = min(`Msg Date`), LastAllTime = max(`Msg Date`))


# Bind in the All First and Last transmit dates
Xmits = left_join(Xmits, firstlast, by = join_by(Ptt == `Platform ID No.`))

# Create new transmit start columns, preserving the original where it is correct
Xmits$XmitStart = if_else(as.Date(Xmits$EarliestXmitTime) >= as.Date(Xmits$deploydate), Xmits$EarliestXmitTime, Xmits$FirstAllTime)
Xmits$XmitEnd = if_else(year(Xmits$LatestXmitTime) == year(Xmits$deploydate), Xmits$LatestXmitTime, Xmits$LastAllTime)
Xmits$totxmitdays = round(as.numeric(difftime(Xmits$XmitEnd, Xmits$XmitStart, units = 'days')), digits = 1)


# WC defines this as the first and last data point (earliest includes status messages, latest does not include status messages)
# We will define it as the first and last timestamp of the location data when the EarliestDataTime is incorrect for our deployment

Xmits$DataStart = if_else(as.Date(Xmits$EarliestDataTime) >= as.Date(Xmits$deploydate), Xmits$EarliestDataTime, NA)
Xmits$DataEnd = if_else(year(Xmits$LatestDataTime) == year(Xmits$deploydate), Xmits$LatestDataTime, NA)


ptts = unique(Xmits$Ptt[is.na(Xmits$DataStart) | is.na(Xmits$DataEnd)])

files = list.files(path = './data/L0/', recursive = T, full.names = T, pattern = '-Locations.csv')

grps = sapply(X = ptts, FUN = function(x){grepl(pattern = x, x = files)})

indices = which(rowSums(grps) == 1)

files = files[indices]

filedf = data.frame(ptt = sapply(strsplit(files, '/'), '[[', 5), file = files)
filedf = filedf[order(filedf$ptt, filedf$file),]
floc = filedf[grepl(pattern = '-1-Locations.csv', x = filedf$file),]
argos = filedf[!filedf$ptt %in% floc$ptt, ]

filedfmain = bind_rows(floc, argos)

locs = lapply(X = filedfmain$file, FUN = function(x){
  tmp = read_csv(x) %>% dplyr::select(Ptt, Date) 
  tmp$Date = as.POSIXct(tmp$Date, format = "%H:%M:%OS %d-%b-%Y", tz = 'UTC')
  return(tmp)
})

locs = bind_rows(locs)
locs = left_join(locs, meta, by = join_by(Ptt == PTT))

# Get the missing datastart times as the first and last location of on deployment
datacalcs = locs %>% filter(as.Date(Date) >= as.Date(deploydate), year(Date) == year(deploydate)) %>% 
  group_by(Ptt, deploydate) %>% summarize(DataStart = min(Date), DataEnd = max(Date)) %>% arrange(Ptt)

baddatastartptts = Xmits$Ptt[is.na(Xmits$DataStart)]
baddataendptts = Xmits$Ptt[is.na(Xmits$DataEnd)]


# Merge df1 and df2 based on the "Ptt" column
Xmits_ <- merge(Xmits, datacalcs[,c('Ptt', 'DataStart', 'DataEnd')], by = "Ptt", suffixes = c("_df1", "_df2"), all.x = TRUE)

# Replace NA values in df1$DataStart with values from df2$DataStart where applicable
Xmits_$DataStart_df1[is.na(Xmits_$DataStart_df1)] <- Xmits_$DataStart_df2[is.na(Xmits_$DataStart_df1)]
Xmits_$DataEnd_df1[is.na(Xmits_$DataEnd_df1)] <- Xmits_$DataEnd_df2[is.na(Xmits_$DataEnd_df1)]

Xmits_ = Xmits_ %>% rename(DataStart = DataStart_df1, DataEnd = DataEnd_df1) %>% 
  dplyr::select(Ptt, PercentDecoded, Passes, PercentArgosLoc, MsgPerPass, DS, MinInterval,deploydate, XmitStart, XmitEnd, totxmitdays, DataStart, DataEnd)
Xmits_$totdatadays = round(as.numeric(difftime(Xmits_$DataEnd, Xmits_$DataStart, units = 'days')), digits = 1)

kmmeta = read_csv("./data/meta/deployment summary_20192023.csv") %>%
  dplyr::select(PTT, tagmodel, sex, deployloc, masskg, lengthcm, girthcm, notes)
kmmeta$tagmodel[kmmeta$PTT == 176859] = 'SPLASH10'

# FIX TAG MODEL
unique(kmmeta$tagmodel)


# THESE SHOULD BE FASTLOC ENABLED (THEY ALL HAD FASTLOC FILES, N=23)
# 224161 224162 224163 224164 224165 224166 224167 224168 235616 235617 235618 235619
# 240184 240185 142351 225840 225841 225842 237647 237648 237649 237651 237652

# 240186 did not transmit fastloc data but is fastloc capable


deffastloc = c(224161, 224162, 224163, 224164, 224165, 224166, 224167, 224168, 235616, 235617, 235618, 235619,
               240184, 240185, 142351, 225840, 225841, 225842, 237647, 237648, 237649, 237651, 237652, 240186)

kmmeta$Fastloc = FALSE
kmmeta$Fastloc[kmmeta$PTT %in% deffastloc] = TRUE

# ALL DIFF SPLASH TAG MODELS AS WRITTEN IN META
# [1] "SPLASH10-297"   "SPLASH10"      
# [6] "SPLASH tag/GPS" "SPLASH10-296F" "SPLASH297" "SPLASH10F-297A"
# [11] "SPLASH10-351F"  "SPLASH Fastloc" "SPLASH10F-393A" "SPLASH10F-296A" "SPLASH10F"
kmmeta$`Tag Model New` = NA
kmmeta$`Tag Model New`[kmmeta$tagmodel %in% c('SPOT6', 'SPOT tag', 'SPOT293', 'SPOT293A')] = 'SPOT-293'


##### NO FASTLOC
kmmeta$`Tag Model New`[!kmmeta$Fastloc & kmmeta$tagmodel %in% c('SPLASH10-297', 'SPLASH10', 'SPLASH10297', 'SPLASH297')] = 'SPLASH10-297'

# SPLASH10-296F - 206722, 206723, 206724
kmmeta$`Tag Model New`[!kmmeta$Fastloc & kmmeta$tagmodel == 'SPLASH10-296F'] = 'SPLASH10-296'

# SPLASH10-351 config - does not have fastloc
kmmeta$`Tag Model New`[!kmmeta$Fastloc & kmmeta$tagmodel %in% c("SPLASH10-351F")] = 'SPLASH10-351'

##### FASTLOC
# FASTLOC-297 config - we are going to assume that general SPLASH10F or GPS are SPLASH10-F-297 configs as this was the most purchased configuration
kmmeta$`Tag Model New`[kmmeta$Fastloc & kmmeta$tagmodel %in% c('SPLASH10F-297A', 'SPLASH10F', 'SPLASH tag/GPS', 'SPLASH Fastloc')] = 'SPLASH10-F-297'

# FASTLOC-296 config
kmmeta$`Tag Model New`[kmmeta$Fastloc & kmmeta$tagmodel == 'SPLASH10F-296A'] = 'SPLASH10-F-296'

# SPLASH10F-393A
kmmeta$`Tag Model New`[kmmeta$Fastloc & kmmeta$tagmodel == 'SPLASH10F-393A'] = 'SPLASH10-F-393'


Xmits2 = full_join(Xmits_, kmmeta, by = join_by(Ptt == PTT))

colnames(Xmits2) = tolower(colnames(Xmits2))

# DEFINE DATACUTOFF PERIODS
Xmits2$cutoffstart = Xmits2$datastart
Xmits2$cutoffend = Xmits2$dataend

# RESOLVE PTT 195526 - That is the animal that died in the gillnet. I estimated the date it died to be around 6/16 based on the dive records. The tag was brought up on a boat and continued pinging back to the dock and back to the fishers house which was how I was able to find it and get it back. 
# It says transmissions went to 7/27, so there was a full month of locational data collected when the tag wasn't on the animal.
Xmits2$cutoffend[Xmits2$ptt == 195526] = as.POSIXct(x = '2020-06-16 00:00:00', tz = 'UTC')

# Calculate dry period as the difference between deployment and first transmit
Xmits2$dryprd_days = round(as.numeric(difftime(Xmits2$xmitstart, Xmits2$deploydate, units = 'days')), digits = 2)

# write out cleaned summary and metadata file
write_csv(x = Xmits2, file = './data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv')

