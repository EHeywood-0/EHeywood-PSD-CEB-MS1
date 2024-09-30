######################################### IDENTIFY HAULOUT BOUTS FROM PERCENT DRY DATA ###################################################

#' Identify haulout bouts from the percent dry data from the Histos output of Wildlife Computers Spot or Splash tags
#' # DETERMINE HAULOUT 
#' (1) ≥ 50%, i.e., the head of the animal was dry for at least 30 minutes within that hour (similar to Tucker et al. in press). 
#'##### This threshold allows us to eliminate the contribution of percent dry values that occurred due to short surface intervals, 
#'##### when the seal is at the surface of the water to breathe or rest. Considering the bimodal distribution of hourly percent dry values, 
#' the first mode of this distribution (~15%, range: 0-40%: Figure 4.B) corresponds to these surface times. 
#' Here, values <50% were assigned a haul-out proportion of 0%.

#'(2) Adjacent to an hourly bin ≥ 95%, i.e., the seal was hauled out for almost the entire hour (≥ 57 minutes). 
#' We include these periods as the probability that the hour-long or multi-hour haul-out event began/ended in the hour itself is much lower than if it were to have begun or ended in an adjacent hour. 
#' As such, the percentage dry of these ‘’tail’’ bins are used as the estimated amount of time 
#' spent hauled out entering or leaving an hour-long or multi-hour haul-out event.

#' @param percentdry An object of class "numeric". A vector of percentages from Histos representing percent of time in a given hour that the tag was dry
#' @param timevec An object of class "POSIXCT". A time vector of START_TIMES, of the same length as percent dry
#' @param threshold1 An object of class numeric. A vector of size one indicating the percent threshold dry above which that hour period is considered hauled out
#' @param threshold2 An object of class numeric. A vector of size one indicating the percent threshold dry for which adjacent percent dry readings will be considered hauled out if they are adjacent in time
#' @return Returns an object of class "BOOLEAN". A vector of size percentdry indicating which observations are considered hauled out (TRUE) or not hauled out (FALSE)
#' @examples
#' # Add some code illustrating how to use the function

find_haulout = function(percentdry, timevec, threshold1 = 50, threshold2 = 95){
  require(lubridate)  
  require(dplyr)
  # ERROR TRAPPING
  if(is.numeric(percentdry) == FALSE && length(percentdry) < 1){
    stop("The input provided was not a numeric vector!")
  }
  if(length(timevec) != length(percentdry)){
    stop("The input timevec and percentdry vectors are not of equal length")
  }
  if(is.POSIXct(timevec) == FALSE){
    stop("The provided timevec is not of class POSIXct")
  }
  
  # Enumerate timevec
  tv = data.frame(tv = timevec) %>% group_by(tv) %>% summarise(n = n())
  
  ### MAKE THIS A FOR LOOP
  dateL = list()
  for (i in 1:nrow(tv)){
    date = seq.POSIXt(from = as.POSIXct(strptime(as.Date(tv$tv[i]), format = "%Y-%m-%d", tz = "UTC")), by = "60 min", length.out = tv$n[i])
    date = as.data.frame(date)
    dateL[[i]] = date
    
  }
  
  date = bind_rows(dateL)
  
  df = data.frame(Start = timevec, Date = date$date, PercentDry = percentdry)
  
  # Solve the issue where the start time is not 00:00, slicing only the relevant hours in that day bin
  df$Hour = hour(df$Start)
  df$n_slice = ifelse(df$Hour == 0, 24, 24-df$Hour)
  dfL = split(df, df$Start)
  for (d in 1:length(dfL)){
    tmp = dfL[[d]]
    nslice = unique(tmp$n_slice)
    tokeep = ifelse(tmp$n_slice < 24, c(rep(FALSE, nrow(tmp) - nslice), rep(TRUE, nslice)), rep(TRUE, nrow(tmp)))
    dfL[[d]] = tmp[tokeep, ]
  }
  df = do.call(rbind, dfL)
  
  nobs = nrow(df)
  
  df$Haulout = ifelse(df$PercentDry >= threshold1, TRUE, FALSE)
  
  df$timediff = as.numeric(c(0,difftime(df$Date[2:nobs], df$Date[1:nobs-1], unit = 'hours')))
  
  df$ConID = cumsum(ifelse(df$timediff == 1, 0, 1))
  
  dfL = split(df, df$ConID)
  storage = lapply(dfL, FUN = function(x){
    idx95 = which(x$PercentDry >= threshold2)
    idx_haulout = unique(c(idx95, idx95+1, idx95-1)) 
    idx_haulout = idx_haulout[idx_haulout %in% seq(1,nrow(x),1)]
    x$Haulout[idx_haulout] = TRUE
    return(x)
  })
  res = bind_rows(storage)
  return(res)
  
} # End of function