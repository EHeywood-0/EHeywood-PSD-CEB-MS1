# Determine Haulout Periods from Percent Dry Timelines
# EIH
# inputs: wildlifecomputers_readdata_fxns.R; 
#         find_haulout_fxn.R;
#         Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv; 
#         *.Histos.csv
# updated: 2024-09-27
# method reference: Irani et al. 2024 (in press)
# INFO
# Within the Histos.csv spreadsheet, the Percent Timeline data have a HistType
# = “Percent.” 
# This timeline records the percentage of time spent above a threshold depth within 24 1-hour blocks per day.
# In other words, there are 24 bins presented in the Histos.csv file, each bin represents a 1-hour period of
# the day. Within each bin is the percentage of time the tag was “dry” or above
# a user-specified depth. When programming your tag, you can choose the
# “Low-Resolution” option, which would round the percentages to the nearest 10%, or the “High-Resolution” option, which
# would round the percentages to the nearest 1%.
# These Percent Timelines are often referred to as Percent Dry Timelines if the wet/dry sensor is used instead of a threshold
# depth. This is usually done to identify haul-out periods or surface behavior for air- breathing animals.
rm(list = ls())
gc()
setwd("~/seal_telemetry/")
library(dplyr)
library(readr)
library(tidyr)

source("./READ-PSB-MoveSeals/src/fxns_wildlifecomputers_readdata.R")
source("./READ-PSB-MoveSeals/src/fxn_find_haulout.R")
# Set max haulout duratin to 7 days
ho_dur_max = 24*7
# Read in meta data on deployments
meta = read_csv(file = "./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv")
colnames(meta) = tolower(colnames(meta))

# Get rid of bad ptts 177509 (no locations), 240183 (rehab), 240186 (tag malfunction), 
keeps = setdiff(unique(meta$ptt), c(177509, 240186, 240183))
meta = meta[meta$ptt %in% keeps, ]

# List all Histos files
files = list.files(path = "./data/L0/", pattern = 'Histos.csv', full.names = TRUE, recursive = TRUE)
filesdf = data.frame(ptt = as.numeric(sapply(strsplit(files, '/'), '[[', 5)), fname = sapply(strsplit(files, '/'), '[[', 6), fullfname = files)
filesdf = filesdf[filesdf$ptt %in% keeps, ]

data_wide = WC_read_histos(files = filesdf$fullfname)
length(unique(data_wide$ptt)) == 61

data_wide$ptt = as.numeric(data_wide$ptt)
data_wide = data_wide[order(data_wide$ptt, data_wide$datetime),]

# Merge deployment summary data with main and make sure we're only using data post deployment dates
data_wide = left_join(data_wide, meta[,c('ptt', 'cutoffstart', 'cutoffend')], "ptt") %>% 
  filter(datetime >= cutoffstart, datetime <= cutoffend) %>%
  arrange(ptt, datetime)

# Duplicates in data wide 
dupes = duplicated(x = data.frame(x = data_wide$ptt, y = as.Date(data_wide$datetime)))
which(dupes)
dupeidx = sort(c(which(dupes), which(dupes)-1))

# remove dupes
data_wide = data_wide[-which(dupes), ]

length(unique(data_wide$ptt)) == 61
#unique(data_wide$ptt)[which(!unique(data_wide$ptt) %in% unique(data_wide_$ptt))]
# PTT 194405 did not transmit beyond January and the only percent timelines transmitted occurred prior to deployment
############################################# HAULOUT DEFINITION USING PERCENT TIMELINES #################################################
##########################################################################################################################################
# DETERMINE HAULOUT 
# (1) ≥ 50%, i.e., the head of the animal was dry for at least 30 minutes within that hour (similar to Tucker et al. in press). 
###### This threshold allows us to eliminate the contribution of percent dry values that occurred due to short surface intervals, 
###### when the seal is at the surface of the water to breathe or rest. Considering the bimodal distribution of hourly percent dry values, 
# the first mode of this distribution (~15%, range: 0-40%: Figure 4.B) corresponds to these surface times. 
# Here, values <50% were assigned a haul-out proportion of 0%.

#(2) Adjacent to an hourly bin ≥ 95%, i.e., the seal was hauled out for almost the entire hour (≥ 57 minutes). 
# We include these periods as the probability that the hour-long or multi-hour haul-out event began/ended in the hour itself is much lower than if it were to have begun or ended in an adjacent hour. 
# As such, the percentage dry of these ‘’tail’’ bins are used as the estimated amount of time 
# spent hauled out entering or leaving an hour-long or multi-hour haul-out event.

# Convert from wide to long form dataframe
data = pivot_longer(data = data_wide, cols = starts_with('bin'), names_to = 'bin', values_to = 'PercentTimeDry') %>% arrange(ptt, datetime)



# Sometimes no data is NA, sometimes it is 0, make it all zero
data$PercentTimeDry[is.na(data$PercentTimeDry)] = 0
data$PercentTimeDry = as.numeric(data$PercentTimeDry)
# Method one: ≥ 50%, i.e., the head of the animal was dry for at least 30 minutes within that hour

# Apply function
data = data %>% arrange(ptt, datetime)
dataL = split(data, data$ptt)

res1 = bind_rows(lapply(dataL, FUN = function(x){
  tmp = find_haulout(percentdry = x$PercentTimeDry, timevec = x$datetime)
  tmp$id = unique(x$ptt)
  rownames(tmp) = NULL
  return(tmp)
}))

res1 = res1 %>% arrange(id, Date)
res1$ConID = paste(res1$id, res1$ConID, sep = '-')
dupes = duplicated(x = res1[,c('id', 'Date')])
which(dupes)

# Combine the list of data frames into a single data frame

df <- res1 %>% group_by(id, ConID) %>% arrange(Date) %>%
  mutate(HauloutID = paste(ConID, cumsum(if_else(Haulout==TRUE, 0,1)), sep = "-"),
         AtSeaID = paste(ConID, cumsum(if_else(Haulout == FALSE, 0,1)), sep = '-')) %>%
  arrange(id, Date)

df$HauloutID[!df$Haulout] <- NA
df$AtSeaID[df$Haulout] <- NA

df$PercentDry = df$PercentDry/100

# Get haulout start and end times + at sea start and end times
df = df %>% group_by(id, HauloutID) %>% arrange(Date) %>% 
  mutate(
    st = if_else((PercentDry[1] >= 0.95), Date[1], Date[1]+(3600-PercentDry[1]*3600)), 
    en = if_else(n() == 1, Date[n()] + 3599, Date[n()]+ PercentDry[n()]*3600)) %>% 
  ungroup() %>% group_by(id, ConID) %>% arrange(Date) %>%
  
  # get at-sea start and end times... 
  
 mutate(stas = case_when(!Haulout & lag(Haulout) ~lag(en),
                         !Haulout[1] ~ Date[1]),
        enas = case_when(!Haulout & lead(Haulout) ~ lead(st),
                        !Haulout[n()] ~ Date[n()] + 3600)) %>%
  arrange(id, Date)

# Make st and en times associated with at sea times NA
df$st[!df$Haulout] = NA
df$en[!df$Haulout] = NA

# get just haulout
hotmp = df %>% filter(Haulout) %>% dplyr::select(-c(stas, enas))

# get just atsea
astmp = df %>% filter(!Haulout) %>% group_by(id, AtSeaID) %>% arrange(Date) %>%
  mutate(st = stas[1], en = enas[n()]) %>% dplyr::select(-c(stas, enas)) %>% arrange(id, Date)

# bind these again
df = bind_rows(hotmp, astmp) %>% arrange(id, Date)

#################################### get haulout events ###############
percentdry_atsea_times = df %>% filter(!Haulout) %>%
  group_by(id, ConID, AtSeaID) %>%
  reframe(EventID = paste('AtSea-', AtSeaID, sep = ''), Start = st[1], End = en[1]) %>% 
  dplyr::select(id, ConID, EventID, Start, End) %>% distinct()


percentdry_haulout_times = df %>% filter(Haulout) %>%
  group_by(id, ConID, HauloutID) %>%
  reframe(EventID = paste('Haulout-', HauloutID, sep = ''), Start = st[1], End = en[1]) %>% 
  dplyr::select(id, ConID, EventID, Start, End) %>% distinct()




hauloutbouts = bind_rows(percentdry_atsea_times, percentdry_haulout_times) %>% arrange(id, Start)

hauloutbouts$What = sapply(strsplit(hauloutbouts$EventID, '-'), '[[', 1)

hauloutbouts$duration_hour = as.numeric(difftime(hauloutbouts$End, hauloutbouts$Start, unit = 'hours'))

keepindx = which(hauloutbouts$duration_hour > 0 & hauloutbouts$duration_hour <= ho_dur_max & hauloutbouts$What == 'Haulout')
keepindx = c(keepindx, which(hauloutbouts$What == 'AtSea'))

hauloutbouts_ = hauloutbouts[keepindx, ]

hauloutbouts_ = hauloutbouts_[order(hauloutbouts_$id, hauloutbouts_$Start), ]

write_csv(hauloutbouts_, file = "./data/L1/haulout/Hg_2019-2023_PercentDry_HauloutAtSea_Events.csv")  
