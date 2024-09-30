# EIH
# Summarise the following for 2023, NEFSC Data for the BOEM REPORT
## For each PTT: Number of Trips, Average Trip Distance, Average Trip Duration, Percentage of trips inside/outside Vineyard lease area
rm(list = ls())

setwd("~/seal_telemetry/")
library(readr)
library(tidyr)
library(dplyr)
library(sf)


###################################################################################################################
########################################## LOAD AND SUBSET DATA ###################################################
###################################################################################################################
################get timing information from full set, get spatial info from modellable trips only #################
###################################################################################################################

# READ IN SSM REGULARIZED COMPLETE TRIPS, REROUTED AROUND LAND
data = read_csv("./data/L2/SSM/Hg-2019-2023-SSM-Tracks-CompleteTrips-Reroute.csv")
data$ptt = sapply(strsplit(data$id, '-'), '[[', 1)
data$TripID = paste(data$ptt, sapply(strsplit(data$id, '-'), '[[', 2), sep = '-')
data$SegID = data$id
data = data %>% dplyr::select(ptt, TripID, SegID, date, lon, lat) %>% arrange(ptt, date)

# Subset to our study time period 2019 - May 31, 2023. Get trips that end prior to May 31
data = data %>% group_by(TripID) %>% mutate(TripStart = min(date), TripEnd = max(date)) %>%
  ungroup() %>% arrange(ptt, date)

# 384 Qualifying Trips - to report spatial characteristics
qual_trips = unique(data$TripID[which(data$TripEnd < as.POSIXct(x = '2023-06-01 00:00:00', tz = "UTC"))])

# Subset data by these qualifying trip ids
data_ = data[data$TripID %in% qual_trips, ]
length(unique(data_$TripID))

# read in meta data on tag deployments
meta = read_csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv") %>% 
  dplyr::select(-c(percentdecoded, passes, percentargosloc, msgperpass, ds, mininterval, xmitstart, xmitend, datastart, dataend))

# merge
data_ = merge(data_, meta, by = 'ptt')
spatialtripdata = data_ %>% arrange(ptt, date)
rm(data, data_, qual_trips)

################################################ TEMPORAL TRIP DATA ####################################################################
# Get trip freq from full trip file (Total Trips per PTT / TotDatDays)
# SUBSET THOSE TRIPS THAT ARE BOOKENDED - 

# Load data, filter out post construction and ptts not in spatial (we won't report on these)
# Mutate a colum 'totlocdays' - we will use this to calculate the Trip Frequency within the study period (total # trips / total days of locations)
df = read_csv("./data/L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv") %>%
  filter(id %in% unique(spatialtripdata$ptt), datetime < as.POSIXct('2024-06-01', tz = 'UTC')) %>% arrange(id, datetime) %>% 
  group_by(id) %>% 
  mutate(totlocdays = as.numeric(difftime(max(datetime), min(datetime), units = 'days'))) %>% 
  ungroup()
  

## Get qualifying trips (must be bookended by haulout periods)
df$eventID = NA
df$eventID[is.na(df$HauloutID)] = paste("T-", df$TripID[is.na(df$HauloutID)], sep="")
df$eventID[is.na(df$TripID)] = paste("H-", df$HauloutID[is.na(df$TripID)], sep = "")

hoL = split(df, df$id)
qual_trip_ids = c()
for (i in 1:length(hoL)){
  tmp = hoL[[i]]
  tmp = tmp[order(tmp$datetime), ]
  events = unique(tmp$eventID)
  ht = sapply(strsplit(events, "-"), "[", 1)
  ts = which(ht == "T")
  hs = which(ht=="H")
  # if the last trip idx is greater than the last haulout idx, remove it
  if (ts[length(ts)] > hs[length(hs)]){
    quals = ts[1:length(ts)-1]
  }else{quals = ts}
  
  qual_ids = events[quals]
  qual_ids = gsub("T-", "", qual_ids)
  qual_trip_ids = c(qual_trip_ids, qual_ids)
  rm(qual_ids, quals, hs, ts, ht, events, tmp)
}

# Filter to qualifying trips and get temporal info
temporaltripdata = df %>% filter(TripID %in% qual_trip_ids) %>%
  group_by(id, TripID) %>% 
  mutate(TripDur = round(as.numeric(difftime(max(datetime), min(datetime), units = 'days')), 2),
         TripStart = min(datetime),
         TripEnd = max(datetime),
         TripISOWeek = lubridate::isoweek(min(datetime)),
         TripMonth = format(min(datetime), '%b')) %>% ungroup() %>%
  left_join(meta, by = join_by(id ==ptt)) %>%
  group_by(id) %>% 
  mutate(TotalTrips = length(unique(TripID)), TripFreq = length(unique(TripID))/totlocdays) %>% 
  dplyr::select(id, TripID, TotalTrips, TripFreq, totlocdays,TripDur, TripStart, TripEnd, TripISOWeek, TripMonth) %>% distinct() %>%
  dplyr::filter(TripDur >= 1)


rm(list = setdiff(ls(), c('spatialtripdata', 'temporaltripdata', 'meta')))

########################################################################################################################################
########################################## SPATIAL TRIP METRICS ########################################################################
########################################################################################################################################
########################################### DISTANCE TO SHORE ##########################################################################
library(geosphere)
library(rnaturalearth)

# convert to sf object
dfsf <- spatialtripdata %>% st_as_sf(coords = c('lon','lat')) %>% 
  st_set_crs(4326)

coastline = ne_coastline(scale = 'medium')

# Define bounding box for the US East Coast and Canada
east_coast_bbox <- st_as_sfc(st_bbox(c(xmin = -76, ymin = 37, xmax = -59, ymax = 46), crs = 4326))

# Subset coastline data to include only the US East Coast and Canada
east_coast <- st_intersection(x = coastline, east_coast_bbox) %>%
  st_union() 
plot(east_coast)

# Convert the sf MULTILINESTRING object to a SpatialLines object
east_coast_sp <- as(east_coast, "Spatial")

# use dist2Line from geosphere - only works for WGS84 to get the distance to shore of each location in a trip
dist <- geosphere::dist2Line(p = st_coordinates(dfsf), line = east_coast_sp)

#combine initial data with distance to coastline
dist_ = as.data.frame(dist)
dist_ = dist_ %>% rename(LANDLON = lon, LANDLAT = lat, distanttoshore = distance) %>% dplyr::select(-c(ID))

# Bind distance information to original dataset
spatialtripdata = cbind(spatialtripdata, dist_)

#plot
library(ggplot2)
p <- ggplot() + 
  geom_sf(data=east_coast) +
  geom_sf(data=dfsf) +
  geom_segment(data=spatialtripdata,aes(x=lon,y=lat,xend=LANDLON,yend=LANDLAT)) 

p
########################################################################################################################################
######################################### GET TRIP LENGTHS (in km) #####################################################################
# Cast to sf line string object
data_sf = st_as_sf(spatialtripdata, coords = c("lon", 'lat')) %>% st_set_crs(4326) %>% 
  group_by(TripID) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

# Get trip length in kilometers
data_sf$TripLength_km = as.numeric(st_length(data_sf) / 1000)



########################################################################################################################################
################################# Calculate TRIP LENGTHS AND PERCENTAGES OF TRIPS OCCURRING INSIDE WEAS ################################
########################################################################################################################################

# READ IN BOEM LEASE ARES
s = st_read(dsn = "./data/shapefiles/BOEM_Wind_Leases_8_30_2023.shp")

# DEFINE A BOUNDING BOX
dfbbox = st_bbox(data_sf)
# INTERSECT TO GET WEAs OCCURRING ONLY WITHIN THE SPATIAL BOUNDS OF THE SEALS
WEA_NE = st_intersection(s, y = st_as_sfc(dfbbox))
plot(WEA_NE['COMPANY'])
WEA_union = st_union(WEA_NE)
# Take trip lines data_sf and intersect with the wind energy areas
data_int = st_intersection(x = data_sf, y = WEA_NE) 

data_int$len_in_WEA = as.numeric(st_length(data_int) / 1000)


########################################################################################################################################
###################################### CALCULATE TIME IN WEA AND PERCENT TIME IN WEA ###################################################
########################################################################################################################################
data_pts = st_as_sf(spatialtripdata, coords = c("lon", 'lat')) %>% st_set_crs(4326)

# Take trip lines data_sf and intersect with the wind energy areas
timeinwea = st_intersection(x = data_pts, y = WEA_union) 
timeinwea = timeinwea %>% mutate(inwea = TRUE) %>% dplyr::select(SegID, date, inwea) %>% st_drop_geometry()


data_pts = left_join(data_pts, timeinwea, by = join_by(SegID, date)) %>% st_drop_geometry()
data_pts$inwea[is.na(data_pts$inwea)] = FALSE
data_pts$inweaID = paste(data_pts$SegID, cumsum(ifelse(data_pts$inwea, 0, 1)), sep = '-')
data_pts$inweaID[!data_pts$inwea] = NA
data_pts = data_pts %>% group_by(TripID, inweaID) %>% summarise(TimeinWeaSecs = as.numeric(difftime(max(date), min(date), units = 'secs')))
data_pts$TimeinWeaSecs[is.na(data_pts$inweaID)] = NA


# Total Time by Trip
timeinwea = data_pts %>% group_by(TripID) %>% summarise(TimeinWeaDay = sum(TimeinWeaSecs, na.rm = T)/60/60/24)
########################################################################################################################################
###################################### AGGREGATE AND SUMMARIZE ###################################################
########################################################################################################################################


data_int_s = data_int %>%  group_by(TripID) %>% # Here you need to insert all the columns from your shapes
  dplyr::summarize(TripLen_in_WEA_km = sum(len_in_WEA)) %>%
  st_drop_geometry()

# Merge with timeinwea data
weadat = left_join(timeinwea, data_int_s, by = 'TripID')

# merge back with other tripdata and get the percentage in WEA

aggtripdata = left_join(data_sf, weadat, 'TripID')

aggtripdata$`Percent in WEA` = aggtripdata$TripLen_in_WEA_km / aggtripdata$TripLength_km * 100

aggtripdata$InWEA = ifelse(is.na(aggtripdata$TripLen_in_WEA_km), FALSE, TRUE)

# MERGE THIS AGG TRIPDATA TO spatialtripdata
spatialtripdata_ = left_join(spatialtripdata, aggtripdata, by = 'TripID') %>%
  group_by(TripID) %>% 
  mutate(MeanDisttoShore = mean(distanttoshore)/1000, 
         MaxDisttoShore = max(distanttoshore)/1000, 
         MinDisttoShore = min(distanttoshore)/1000,
         SDDisttoShore = sd(distanttoshore)/1000,
         PercTimeinWea = TimeinWeaDay / (as.numeric(difftime(TripEnd, TripStart, units = 'days'))) * 100) %>% ungroup() %>%
  dplyr::select(ptt, TripID, TripStart, TripEnd, 
                TripLength_km, MeanDisttoShore, MaxDisttoShore, MinDisttoShore, SDDisttoShore,
                InWEA,TripLen_in_WEA_km, `Percent in WEA`, TimeinWeaDay, PercTimeinWea, 
                sex, masskg, lengthcm, girthcm) %>% distinct() %>%
  arrange(TripID)

spatialtripdata_$TimeinWeaDay[!spatialtripdata_$InWEA] = NA
spatialtripdata_$PercTimeinWea[!spatialtripdata_$InWEA] = NA

write_csv(x = spatialtripdata_, file = './data/L3/Hg-PreConstruction-SpatialTripMetrics_CompleteModelableTrips.csv')


###################################################################################################################################
################################################ MERGE SPATIAL AND TEMPORAL #######################################################
###################################################################################################################################

test = left_join(temporaltripdata, spatialtripdata_ %>% dplyr::select(-c(TripStart, TripEnd, ptt, sex, masskg, lengthcm, girthcm)), by = 'TripID')

fullsummary = left_join(test, meta %>% dplyr::select(ptt, tagmodel, sex, masskg, lengthcm, girthcm), by = join_by(id == ptt))

write_csv(x = fullsummary, file = './data/L3/Hg-PreConstruction-SpatialTemporalTripMetrics.csv')
write_csv(x = temporaltripdata, file = './data/L3/Hg-PreConstruction-TemporalTripMetrics.csv')

################################################# FINAL SUMMARY ###################################################################
pttsummtripmetrix = fullsummary %>%
  group_by(id, sex, masskg, lengthcm, TotalTrips, TripFreq, totlocdays) %>% 
  summarise(`Mean Trip Duration (days)` = round(mean(TripDur, na.rm=T)),
            `Max Trip Duration (days)` = round(max(TripDur)),
            `Min Trip Duration (days)` = round(min(TripDur)),
            `SD Trip Duration (days)` = round(sd(TripDur, na.rm=T)),
            `Mean Trip Length (km)`= round(mean(TripLength_km, na.rm=T)),
            `SD Trip Length (km)` = round(sd(TripLength_km, na.rm=T)),
            `No. Trips in WEA`= sum(InWEA),
            `Mean Trip Length in WEA (km)` = round(mean(TripLen_in_WEA_km, na.rm = T)),
            `SD Trip Length in WEA (km)` = round(sd(TripLen_in_WEA_km, na.rm=T)),
            `Mean Percent Time in WEA (km)` = round(mean(PercTimeinWea, na.rm=T)),
            `SD Percent Time in WEA (km)` = round(sd(PercTimeinWea, na.rm=T)),
            `Mean Max Distance to Shore (km)` = round(mean(MaxDisttoShore, na.rm=T)),
            `SD Max Distance to Shore (km)` = round(sd(MaxDisttoShore, na.rm=T)),
            `Mean Distance to Shore (km)`= round(mean(MeanDisttoShore, na.rm=T)),
            `SD Mean Distance to Shore (km)` = round(sd(MeanDisttoShore, na.rm=T)),
            `Trip Months` = paste(unique(TripMonth), collapse = ','))



write_csv(pttsummtripmetrix, "./data/L3/Hg-Pre-Construction-AllTripMetrix_byPTT.csv")

################################################### PLOTTING #######################################################################

# PLOTS BY WEEK
t = fullsummary %>% dplyr::select(id, sex, TripID, TripLength_km, TripDur, TripISOWeek, TripMonth, TripStart)
t = t[which(t$TripMonth != 'Dec'),]
t$TripMonth = factor(t$TripMonth, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
# Likely need to only include completed trips in here (bookended by haulout event)

# Create some labels
# Calculate the number of seals transmitting in each month
tlabs <- t %>%
  group_by(TripMonth) %>%
  mutate(SealsCount = n_distinct(id), TripCount = n_distinct(TripID)) %>% ungroup()

# Create a new variable for labeling with both Male and Female counts
tlabs$label <- paste('N: ', tlabs$SealsCount, '\n','Trips: ', tlabs$TripCount)

labels = tlabs %>% group_by(TripMonth) %>% summarise(Label = unique(label))
names(t)
g = ggplot() +
  geom_boxplot(data = t, mapping = aes(x = TripMonth, y = TripLength_km, fill = sex)) + 
  ylab('Trip Length (km)') +
  xlab('Month') +
  geom_text(aes(x = labels$TripMonth, y = rep(2700, nrow(labels)), label = labels$Label), 
            position = position_dodge(0.9), group = 'B') +
  theme_bw()
g
ggsave(filename = "./plots/manuscript/2024Manuscript-TripDistanceBoxPlot.jpg", 
       plot = g, dpi = 300, width = 10, height = 8, scale = 1)

g = ggplot() +
  geom_boxplot(data = t, mapping = aes(x = TripMonth, y = TripDur, fill = sex), outliers = F) +
  ylab('Trip Duration (day)') +
  xlab('Month') +
  geom_text(aes(x = labels$TripMonth, y = rep(35, nrow(labels)), label = labels$Label), 
            position = position_dodge(0.9), group = 'B') +
  theme_bw() +
  theme(legend.position = 'NONE') 
g
ggsave(filename = "./plots/manuscript/2024Manuscript-TripDurationBoxPlot.jpg", plot = g, dpi = 300,
       width = 10, height = 8, scale = 1)


########### PLOTTING SCRATCH
dat = fullsummary
#dat = dat[dat$TotalTrips>1, ]
#dat = dat[dat$TripDur > 1, ]

# AGAINST MASS
plot(dat$masskg, dat$TripFreq)
plot(dat$masskg, dat$TripDur)
plot(dat$masskg, dat$TripLength_km)
plot(dat$masskg, dat$MinDisttoShore)
plot(dat$masskg, dat$MaxDisttoShore)

# AGAINST LENGTH
plot(dat$lengthcm, dat$TripFreq)
plot(dat$lengthcm, dat$TripDur)
plot(dat$lengthcm, dat$TripLength_km)
plot(dat$lengthcm, dat$MinDisttoShore)
plot(dat$lengthcm, dat$MaxDisttoShore)
plot(dat$lengthcm, dat$MeanDisttoShore)

# AGAINST Girth
plot(dat$girthcm, dat$TripFreq)
plot(dat$girthcm, dat$TripDur)
plot(dat$girthcm, dat$TripLength_km)
plot(dat$girthcm, dat$MinDisttoShore)
plot(dat$girthcm, dat$MaxDisttoShore)
plot(dat$girthcm, dat$MeanDisttoShore)

# AGAINST WEEK
#plot(dat$TripISOWeek, dat$TripFreq)
plot(dat$TripISOWeek, dat$TripDur)
plot(dat$TripISOWeek, dat$TripLength_km)
plot(dat$TripISOWeek, dat$MinDisttoShore)
plot(dat$TripISOWeek, dat$MaxDisttoShore)
plot(dat$TripISOWeek, dat$MeanDisttoShore)

# AGAINST SEX
#plot(dat$TripISOWeek, dat$TripFreq)
boxplot(TripFreq ~ sex, dat)
boxplot(TripDur ~ sex, dat)
boxplot(TripLength_km ~ sex, dat)
boxplot(MeanDisttoShore ~ sex, dat)
boxplot(MinDisttoShore ~ sex, dat)
boxplot(MaxDisttoShore ~ sex, dat)


# LOOKS LIKE THERE IS AN AFFECT OF LENGTH / GIRTH / MASS / SEX on TripFreq

