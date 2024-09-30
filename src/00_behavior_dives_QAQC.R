# First level of dive behavior QAQC and processing
# EIH
# updated:2024-09-27
# inputs: wildlifecomputers_readdata_fxns.R
# 1. Read raw Behavior records for all tags 2019-2023
# 2. Assign a row / event id for each ptt beh dataframe, preserving the order of the data
# 3. cleans dive data
rm(list = ls())
gc()

# Set working directory to raw data on server
setwd('~/seal_telemetry/')
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)

source('./READ-PSB-MoveSeals/src/fxns_wildlifecomputers_readdata.R')
########## DETAILS OF BEHAVIOR CSV

#- DeployID Friendly (and should be unique) name given to the tag by the user. If no specific friendly name is given, this is the PTT ID.

#- PTT Argos Platform Transmitter Terminal identifier, which is a unique number identifying your instrument.

#- DepthSensor Indicates the resolution of the depth sensor. Typical values are 0.1m, 0.5, and 1.0m. 
#             The value in this column is blank if the gain was unavailable - in this case a 0.5-meter resolution is assumed.

#- Source Indicates where the data came from, or how it was generated.

#- Instr Wildlife Computers instrument family.

#- Count Total number of times this particular data item was received, verified, and successfully
# decoded.

#- Start Timestamp when the particular behavior item began, to the nearest 30 seconds.

#- End Timestamp when the particular behavior item ended.

#- What Kind of element - Surface, Haulout, Dive, Message, or Unrecognized.

#- Number The number of times the behavior element described by the following Shape, Depth,
# and Duration columns occurred back-to-back.

#- Shape If the behavior is a dive, this is the shape of the dive determined by the tag (see below).

#- DepthMin The dive Depth is between DepthMin and DepthMax meters.

#- DepthMax The dive Depth is between DepthMin and DepthMax meters.

#- DurationMin The duration of the element is between DurationMin and DurationMax seconds.

#- DurationMax The duration of the element is between DurationMin and DurationMax seconds.

################################## THESE ARE RELICS I BELIEVE
#- Shallow In between dives, this is the spent ABOVE the threshold to determine start and end of a
# dive.

#- Deep In between dives, this is the time spent BELOW the threshold to determine start and
# end of a time and ABOVE the depth required to qualify a dive.

#################################################################################################
# NOTES ABOUT ADL and POTENTIAL DIVE DATA ALIASING DUE TO 10Sec sampling interval for our juvenile gray seals.

# FOR ALL DIVES: plot IDI vs Duration - if 100% of these don't make sense... likely some depth aliasing. Meaning
# if the animal isn't spending time at the surfacing recovering from a 30 min dive, it likely isn't real.
# Vertical ascent / decent rate can be > 1m/s thus the animal could surface
# (5 m to sea level, quick breath, and go back down in <10 seconds) resulting in dive aliasing (two different dives getting concatenated)

########################### DEFINITIONS and BACKGROUND ######################################################
# ADL: the point during a dive when lactic acid begins to increase above resting levels, i.e. the threshold where the animal goes from exclusively oxidative metabolism to including anaerobic metabolism (Costa et al., 2023)

#Beyond this ???threshold,??? a substantial increase in the post-dive surface interval
# is required to clear the lactate acid to restore it to pre-dive levels. 
# Unfortunately, the ADL has only been directly measured in a few species because of logistical challenges: 
# Weddell seals (Kooyman et al., 1980)
# Baikal seals (Pusa sibirica) (Ponganis et al., 1997a)
# California sea lions (Ponganis et al., 1997c), 
# emperor penguins (Aptenodytes forsteri) (Ponganis et al., 1997b), 
# beluga whales (Delphinapterus leucas) (Shaffer et al., 1996) and bottlenose dolphins (Williams et al., 1999).

# When animals dive aerobically (below the ADL), the post-dive surface interval is short relative to the dive duration.
# However, when they exceed the ADL, the post-dive surface interval 
# (or recovery time) becomes a more significant proportion of the dive cycle. 
# Thus, an animal spends more time at the surface recovering from a dive than underwater in a dive, 
# reducing its overall efficiency. 
# This led to the concept that animals should dive aerobically for routine sustained diving, 
# reserving anaerobic dives for rare or extreme events (Kooyman et al., 1980).


# Noren et al. 2005 
# -CADL (calculated/estimated by dividing usable oxygen stores (measured) by the estimated rate of O2 consumption) 
# reports in the literature:
# Aerobic dive limit for western Atlantic adult gray seals: 11.9 +/- 0.9 minutes
# cADL Yearling: 6.1 +/- 0.3 minutes
# 24 Days Post Weaning: 3.8 +/- 0.1 minutes

# Personal Comms.
# err on the side of a time-dependant threshold informed by 2024 higher resolution data. 
# longest 2024 dive is 8 mins 35 seconds - 
# anything beyond 10 mins and definitely 12 mins, should be treated with caution. 
# as a field ecological contextualization of diving physiology 
# is so poor that we're just not there yet in terms of getting a more physiologically-informed 
# basis for what might be an ADL cut-off would be, or what that actually means for the animal during repeated diving. 

# Beck et al. 2003 dive behavior sable island grays

# 90 - 92% of dives < 8 minutes duration
# Maximum durations were 20.3 minutes (M) and 22.0 min (F)

# 95-97% of dives < 120m depth
# Max dive depth was 412 for M and 354 F

#Jessopp et al. 2013 
# Max dive depth in eastern Atlantic grays: 455m

##############################################################################################################
# Load Deployment Metadata
meta = read_csv(file = "./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv")
colnames(meta) = tolower(colnames(meta))

# Get rid of bad ptts 177509 (no locations), 240183 (rehab), 240186 (tag malfunction), 
# 177509 was Argos Only and failed to transmit locations
# 240183 was rehab
# 240186 transmitted some data but suffered tag malfunction after less than one day
keeps = setdiff(unique(meta$ptt), c(177509, 240186, 240183))
meta = meta[meta$ptt %in% keeps, ]

# Load BEH Data (N=35)
files = list.files(path = './', pattern = '-Behavior.csv', full.names = T, recursive = T)
# read in all files using custom WC_read_beh script
beh = WC_read_beh(files)

# Remove rehab animal (240186), malfunctioning tag (240183), and no locs (177509)
beh = beh[beh$ptt %in% keeps, ]
length(unique(beh$ptt)) #34

table(meta$sex[meta$ptt %in% unique(beh$ptt)])
# 17M and 17F

############################# Calculate a Continuous ID as data may be continuous between messages from the same PTT
# Filter for rows where 'what' is "Message"
filtered <- subset(beh, what == "Message")

# Sort the filtered data frame by 'ptt' and 'roworder'
sorted <- filtered[order(filtered$ptt, filtered$roworder), ]

# Group by 'ptt'
grouped <- split(sorted, sorted$ptt)

# Calculate time differences within each 'ptt'
conidL <- lapply(grouped, function(df) {
  deltat = difftime(df$start[2:nrow(df)], df$end[1:nrow(df)-1], units = 'secs')
  deltat = c(0, deltat)
  # Create 'conid'
  conid = paste(df$ptt, cumsum(ifelse(deltat > 60, 1, 0)), sep = '--')
  
})

sorted$conid = unlist(conidL)

sorted = sorted[,c('mid', "conid")]

beh2 = merge(beh, sorted, by = "mid", all.x = TRUE)

beh2 = beh2[ ,c('ptt', 'mid', 'conid', 'roworder', 'start', 'end', 'what', 'shape', 'depth', 'duration')]

beh2 = beh2[order(beh2$ptt, beh2$roworder), ]


######################################### PAIR DIVES AND IDIs#################################################################
##############################################################################################################################
# create dive / surface pair IDS based on conid (continuous behavior stretches) 
# this is not representative of ALL collected dives or surface intervals, 
# as ones which cannot be paired are eliminated
idi = beh2[beh2$what != 'Message', ]

grouped = split(idi, idi$conid)

processed = lapply(grouped, FUN = function(df){
  
  # Within a continuous behavior stretch, check to make sure dive / surface intervals are every other
  # If they are not, throw out one, as two dives occurring back to back likely bad data
  keep = c(ifelse(df$what[1:nrow(df)-1] != df$what[2:nrow(df)], TRUE, FALSE), TRUE)
  
  # Filter tmp by keep (boolean)
  tmp = df[keep, ]
  
  # if the continuous event grouping starts with a surface, eliminate it
  if (tmp$what[1] == "Surface") {
    tmp <- tmp[2:nrow(tmp), ]
  }
  # if it ends with a dive, eliminate it
  if (tmp$what[nrow(tmp)] == "Dive") {
    tmp <- tmp[1:nrow(tmp)-1, ]
  }
  
  return(tmp)
  
})

# Combine into a single data frame
filtered <- dplyr::bind_rows(processed)
print("Number of dives or surface intervals filtered out due to a lack of a corresponding event (dive-surface interval matches")
nrow(idi) - nrow(filtered) # 137 unmatched . . . do we want to eliminate these??? 

idi$eventid = paste(idi$ptt, idi$roworder, sep = '--')
filtered$eventid = paste(filtered$ptt, filtered$roworder, sep = '--')

eventids_elim = setdiff(unique(idi$eventid), unique(filtered$eventid))

datanomatch = idi[idi$eventid %in% eventids_elim, ]
table(datanomatch$what)

# Eliminated 50 dives and 62 surfacing periods with no match, 41 of the 50 dives are greater than 15 minutes which is highly unlikely, if not impossible
length(which(datanomatch$duration >= 15*60))

# Group by 'conid' and 'what'
grouped_conid_what <- split(filtered, list(filtered$conid, filtered$what))

# Add 'SDPairID' column
grouped_conid_what <- lapply(grouped_conid_what, function(df) {
  df$SDPairID <- paste(df$conid, seq_len(nrow(df)), sep = "-")
  df
})

res = dplyr::bind_rows(grouped_conid_what)

idi0 = idi
idi = res[order(res$ptt, res$roworder), ]

rm(list=setdiff(ls(), c("idi", 'meta')))


##################################################### CONVERT FROM LONG TO WIDE DATA FORMAT ############################################################
########################################################################################################################################################

# Make new reshape function that is much faster than the base R stats version
new_reshape <- function(...){
  interaction <- function(x, drop) do.call(paste0, x)
  environment(reshape) <- environment()
  reshape(...)
}


# Get dive / surface start times as a wide dataframe
stlong <- idi[, c("SDPairID", "what", 'depth', 'duration', 'start', 'end')]

idiW <- new_reshape(data = stlong, 
                    direction = 'wide',
                    idvar = "SDPairID", 
                    timevar = "what",
                    sep = "_")

IDI = idiW[,c('SDPairID', 'duration_Dive', 'duration_Surface', 'depth_Dive','start_Dive', 'end_Dive', 'start_Surface', 'end_Surface')]
colnames(IDI)[2:4] = c('DiveDur', 'IDI', 'Depth')
rownames(IDI) = NULL


IDI$ptt = as.numeric(sapply(strsplit(IDI$SDPairID, split = '--'), '[[', 1))

write.csv(x = IDI, 
          file = '~/seal_telemetry/data/L1/dive/Hg_2019-2023-PairedDivesIDIs.csv', 
          row.names = FALSE)


####################################### ELIMINATE DIVES ##########################################################
# First eliminate by the data cutoffdates
IDI = left_join(IDI, meta, by = 'ptt') %>% filter(start_Dive >= cutoffstart, end_Surface <= cutoffend)

# Eliminate dives > 15 minutes following Beck et al., 2003 and our own unpublished 2014 data
# Eliminate dives deeper than 455 meters (Jessopp et al. 2013) - max recorded depth for the species
idisub = IDI[which(IDI$DiveDur/60 < 15), ]
idisub = idisub[which(idisub$Depth <= 455),]
idisub = idisub[which(idisub$Depth >= 5), ]


# 806 dives with durations over 15 minutes...
length(which(IDI$DiveDur/60 >= 15))
nrow(IDI) - nrow(idisub) #eliminated 6479 dives
########################## EXPLORATORY PLOTS ############################

length(which(idisub$DiveDur/60 > 8)) / nrow(idisub) * 100

# Zoom in on axis
g = ggplot(data = idisub, mapping = aes(x = DiveDur/60, y = IDI/60)) +
  geom_point() +
  xlab('Dive Dur (min)') +
  ylab('IDI (min)')+
  ylim(c(0, 2.5)) +
  xlim(c(0,15)) +
  theme_bw()
g


##########################################################################################################################################
beh = idisub
beh$duration = beh$DiveDur/60

# 791 dives greater than 8 minutes (but less than 15) or ~4% (data 2019-2023) 
length(which(beh$duration > 8))  / nrow(beh) * 100

# Number of PTTs for which we have recorded dive behavior
length(unique(beh$ptt)) # 34 ptts with data

# Visualize the distribution of dive depths (meters) and duration (minutes)
hist(beh$Depth, breaks = 50, xlim = c(0,475),xlab = 'Maximum Dive Depth (m)', main = NULL)
print(paste('Max Dive Depth Recorded: ', max(beh$Depth), 'm; for a duration of ', round(beh$duration[beh$Depth == max(beh$Depth)]), ' min.'))
hist(beh$duration, breaks = 50, xlim = c(0,15),xlab = 'Maximum Dive Duration (minutes)', main = NULL)
print(paste('Max Dive Duration Recorded: ', max(beh$duration), 'min; at a depth of:', 
            round(beh$Depth[beh$duration == max(beh$duration)]), ' meters, and a surface interval of: ',
            round(beh$IDI[beh$duration == max(beh$duration)]/60)))


write.csv(x = beh, file = './data/L1/dive/Hg_2019-2023_BEHDiveRecords_QAQC.csv', row.names = FALSE)
