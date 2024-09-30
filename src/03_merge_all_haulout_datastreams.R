# combine all haul out data streams, merge overlapping events
# EIH
# inputs: Hg_2019-2023_PercentDry_HauloutAtSea_Events.csv;
#         Hg_2019-2023_Beh_Haulout_Events.csv;
#         HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv;
#         Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv
# updated: 2024-09-27
rm(list = ls())
gc()
library(readr)
library(dplyr)

# Set working dir
setwd("~/seal_telemetry/")
#load all haulout files

dat1<-read_csv("./data/L1/haulout/Hg_2019-2023_PercentDry_HauloutAtSea_Events.csv") %>% 
  filter(What == 'Haulout') %>% 
  mutate(ptt = id, s = Start, e = End, type = 'percent dry') %>% 
  dplyr::select(ptt, type, s, e) %>% arrange(ptt, s)


dat2<-read_csv("./data/L1/haulout/Hg_2019-2023_Beh_Haulout_Events.csv") %>% 
  mutate(s = startmin, e = endmin, type = 'behavior') %>% 
  dplyr::select(ptt, type, s, e) %>% arrange(ptt, s)


# Spatial (haulouts from tracking data)
dat3 = read_csv("./data/L1/locs/HG_2019-2023_HAULOUT_TRIPID_ASSIGNMENT.csv") %>%
  filter(!is.na(HauloutID)) %>% mutate(type = 'spatial', ptt = id) %>% 
  group_by(ptt, type, HauloutID) %>% reframe(s = min(datetime), e = max(datetime)) %>%
  dplyr::select(-c(HauloutID))
  
all = rbind(dat1, dat2, dat3) %>% arrange(ptt, s, e)

# split
allL = split(all, all$ptt)

res = lapply(X = allL, FUN = function(x, nlags = 1){
  x = x[order(x$s, x$e),]
 
  if (nrow(x) < 4){
    print('not enough events to look at lags') 
    return(x)
    }
  
  if (nrow(x) <= nlags & nrow(x) >= 4){
    
    nlags = nrow(x) - 1
    ncx = ncol(x)
    for(i in 1:nlags){
      dt = c(rep(0, i), difftime(x$s[(i+1):nrow(x)], x$e[1:(nrow(x)-i)], units = 'mins'))
      which(dt < 0)
      x[,(ncx+i)] = dt
      colnames(x)[ncx+i] = paste('Lag', i)
   
       
    }
    return(x)
    } # end second if
  
  if (nrow(x) > nlags){
    ncx = ncol(x)
    for(i in 1:nlags){
      dt = c(rep(0, i), difftime(x$s[(i+1):nrow(x)], x$e[1:(nrow(x)-i)], units = 'mins'))
      which(dt < 0)
      x[,(ncx+i)] = dt
      colnames(x)[ncx+i] = paste('Lag', i, sep = '')
      }    
    return(x)  
  }
})

res = bind_rows(res)
res = res[order(res$ptt, res$s),]

##########################
# LOOK AT LAGS - take the the index of each negative first lag - 1 (END)

find_first_n_consecutive_negatives <- function(idx, s, e) {
  
  
  eidx = e[idx-1] # the value from which to compare everything
  n = length(e)
  
  if (idx == n){
    count = 0
    return(count)
    }
  
  lags = s[c((idx+1):n)] - eidx
  
  # Create a logical vector indicating if elements are negative
  is_negative <- lags < 0
  
  # Initialize counter and storage for the result
  count <- 0
  
  # Iterate through the vector
  for(i in 1:length(is_negative)){
    
    if (is_negative[i]){
      count = count+1
    }else{break}
    
  }
  
  return(count)
}

res$Lag1[is.na(res$Lag1)] = 0

lagsL = split(res, res$ptt)

res2L = lapply(X = lagsL, function(x){

  repeat{
    # Order by start time
    x = x[order(x$s), ]
    x$oid = 0
    negidx = which(x$Lag1 < 0)
    
    if (length(negidx) == 0) {
      return(x)
      break
    }else{
      for (idx in negidx) {
        num = find_first_n_consecutive_negatives(idx = idx, s = x$s, e = x$e)
        x$oid[c((idx-1):(idx+num))] = paste(x$ptt[1], idx, sep = '-')
        } # end for loop
      
      # Compute new start and end times
    
      overlaps = x[x$oid != 0, ]
      overlaps = overlaps[order(overlaps$s), ]
      overlaps = overlaps %>%
        group_by(oid) %>%
        mutate(s = min(s), e = max(e), type = 'Beh + PercDry') %>%
        ungroup() %>%
        dplyr::select(ptt, type, s, e) %>%
        distinct()
    
    none = x[x$oid == 0, ] %>%
      dplyr::select(ptt, type, s, e)
    
    comb = bind_rows(overlaps, none) %>%
      arrange(s)
    
    comb$Lag1 = c(0, difftime(comb$s[2:nrow(comb)], comb$e[1:(nrow(comb)-1)], units = 'mins'))
    
    x = comb # Update x with comb for the next iteration
    } # end else
  } # end repeat
})


# Bind test store
res2 = bind_rows(res2L)


finalHO = res2 %>% group_by(ptt) %>% arrange(s) %>%
  mutate(startho = s, endho = e, HID = paste(ptt, row_number(), sep = '-'), ho_dur_hours = round(as.numeric(difftime(e, s, units = 'hours')), digits = 2)) %>%
  dplyr::select(ptt, HID, type, startho, endho, ho_dur_hours)

meta = read_csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv") %>%
  dplyr::select(ptt, sex, masskg, lengthcm, girthcm, tagmodel, deploydate, xmitstart, xmitend, totxmitdays, datastart, dataend, totdatadays, cutoffstart, cutoffend, dryprd_days)

finalHO = left_join(finalHO, meta, by = 'ptt') %>% arrange(ptt, startho)

# TAG 195526 which was bycaught and died and has data cutoff of 06/16 has one percent dry haulout with an end date of 06/17
# We will eliminate this
finalHO = finalHO[-c(which(finalHO$ptt == 195526 & finalHO$endho > finalHO$cutoffend)),]

# get rid of - duration haulouts
testing = finalHO %>% group_by(ptt) %>% reframe(MeanHODur = mean(ho_dur_hours))
testing1 = finalHO %>% filter(ho_dur_hours > 0) %>% group_by(ptt) %>% reframe(MeanHODur = mean(ho_dur_hours))

# this is due to the one loc spatial haulouts.
finalHO = finalHO[finalHO$ho_dur_hours > 0, ]

write_csv(x = finalHO, file = "./data/L1/haulout/Hg_2019-2023_ALL_HauloutMethods_Merged.csv")
