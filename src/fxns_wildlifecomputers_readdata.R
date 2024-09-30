# Wildlife Computers data reading functions
# EIH
# 2024-03-27

# Description: A series of functions to read in and process vectors of filenames direct from WC portal downloaded data
# Recommend: WC Portal downloads should be unzipped with an archive utility that has checksum functionality (e.g.7-zip) 
# so that corrupted files from download are captured

# Once established a clean dir of portal downloads, this family of functions takes a list of filenames, reads in and does some basic processing / data cleaning

#################################################### LOAD LOCATION DATA #############################################################
#' Load Argos/FastGPS Wildlife Computer Location Data
#' 
#' @param files a vector of full path file names to be loaded
#' @param datetime_format A character string. The default for the formatting Wildlife Computers date times is "%H:%M:%OS %d-%b-%Y" but others can be specified by POSIX standard
#' @param tz A character string indicating the time zone (default is Universal Time Coordinated or UTC)
#' 
#' 
#' 
#' @return Data frame with positional data
WC_read_locs <- function(files, datetime_format = '%H:%M:%OS %d-%b-%Y', tz = 'UTC') {
  
  
  dataL = lapply(X = files, FUN = function(x){
    tmp = read.csv(x, header = TRUE, stringsAsFactors = FALSE, na.strings = "")
    colnames(tmp) = tolower(colnames(tmp))
    tmp$quality[tmp$type == 'FastGPS'] = 'G'
    tmp$datetime = as.POSIXct(strptime(x = tmp$date, format = datetime_format, tz = tz))
    if (anyNA(tmp$datetime)) {
      
      warning('NAs introduced converting to POSIXct - check date time formatting assumptions')
      
    }
    names(tmp)
    cols_keep = setdiff(colnames(tmp), c("count", "gpe.msd", "gpe.u", "offset", "offset.orientation"))
    tmp = tmp[,cols_keep]
    return(tmp)
  })
  
  data = do.call(rbind, dataL)
  
  return(data)
}

#' Load FastGPS Wildlife Computer Fastloc Data and subset columns for quality control
#' 
#' @param files a vector of file names to be loaded
#' @param datetime_format A character string. The default for the formatting Wildlife Computers date times is "%H:%M:%OS %d-%b-%Y" but others can be specified by POSIX standard
#' @param tz A character string indicating the time zone (default is Universal Time Coordinated or UTC)
#' 
#' 
#' 
#' @return Data frame with positional data
WC_read_fastloc <- function(files, datadir,datetime_format = '%H:%M:%OS %d-%b-%Y', tz = 'UTC') {
  
  require(dplyr)
  
  dataL = lapply(X = files, FUN = function(x){
    fullpath = list.files(path = datadir, recursive = T, pattern = x, full.names = T)
    tmp = read.csv(fullpath, header = TRUE, stringsAsFactors = FALSE, na.strings = "")
    print(x)
    print(fullpath)
    
    if (!all(is.na(tmp[1,]))){
      colnames(tmp) = tolower(colnames(tmp))
      tmp = tmp[,c('name', 'locnumber', 'hauled.out', 'satellites', 'residual', 'time.error', 'day', 'time')]
      tmp$datetime = as.POSIXct(strptime(paste(tmp$time, tmp$day, sep = ' '), format = datetime_format, tz = 'UTC'))
      tmp$ptt = sapply(strsplit(x, '-'), '[[', 1)
      tmp$ptt = as.character(tmp$ptt)
      
      return(tmp)
    }
    
  })
  
  data = do.call(rbind, dataL)
  
  return(data)
}





#' Load Wildlife Computer Histos Data and subset specific user defined histtypes
#' 
#' @param files a vector of full path file names to be loaded
#' @param HistType a character string indicating which WC histos type is desired (Percent, TAD, TAT, DiveDepth, DiveDuration)
#' @param datetime_format A character string. The default for the formatting Wildlife Computers date times is "%H:%M:%OS %d-%b-%Y" but others can be specified by POSIX standard
#' @param tz A character string indicating the time zone (default is Universal Time Coordinated or UTC)
#' 
#' 
#' 
#' @return Data frame with histogram data

WC_read_histos <- function(files, HistType = c('1Percent', 'Percent'),datetime_format = '%H:%M:%S %d-%b-%Y', tz = 'UTC') {
  
  colClasses=c(rep('character', 87))
  
  dataL = lapply(X = files, FUN = function(x){
    tmp = read.csv(x, header = TRUE, stringsAsFactors = FALSE, na.strings = "", colClasses = colClasses)
    colnames(tmp) = tolower(colnames(tmp))
    tmp = tmp[tmp$histtype %in% HistType, ]
    tmp$datetime = as.POSIXct(strptime(x = tmp$date, format = datetime_format, tz = tz))
    if (anyNA(tmp$datetime)) {
      
      warning('NAs introduced converting to POSIXct - check date time formatting assumptions')
      
    }
    names(tmp)
    cols_keep = c("ptt", "datetime", "time.offset", paste('bin', 1:24, sep = ''))
    tmp = tmp[,cols_keep]
    return(tmp)
  })
  
  data = do.call(rbind, dataL)
  
  return(data)
}

#' Load Wildlife Computer Behavior Data
#' 
#' @param files a vector of full path file names to be loaded
#' @param datetime_format A character string. The default for the formatting Wildlife Computers date times is "%H:%M:%OS %d-%b-%Y" but others can be specified by POSIX standard
#' @param tz A character string indicating the time zone (default is Universal Time Coordinated or UTC)
#' 
#' 
#' 
#' @return Data frame with positional data



WC_read_beh = function(files, datetime_format = '%H:%M:%S %d-%b-%Y', tz = 'UTC'){
  
  dataL = lapply(X = files, FUN = function(x){
    tmp = read.csv(x, header = T, stringsAsFactors = F)
    tmp$RowOrder = 1:nrow(tmp) # preserve the order of the events
    
    colnames(tmp) = tolower(colnames(tmp))
    
    tmp$start = as.POSIXct(strptime(x = tmp$start, format = datetime_format, tz = tz))
    tmp$end = as.POSIXct(strptime(x = tmp$end, format = datetime_format, tz = tz))
    
    # get a message id
    tmp$mid = paste(tmp$ptt, cumsum(ifelse(tmp$what == 'Message', 1, 0)), sep = "_")
    
    # calculate the average between the min and max estimates as the true value is likely between these (WC) 
    tmp$depth = apply(X = data.frame(min = tmp$depthmin, max = tmp$depthmax), MARGIN = 1, FUN = mean)
    tmp$duration = apply(X = data.frame(min = tmp$durationmin, max = tmp$durationmax), MARGIN = 1, FUN = mean)
    
    tmp = tmp[,c('ptt', 'mid', 'roworder','start', 'end', 'what', 'shape', 'depth', 'duration')]
    
    return(tmp)
  })
  # Bind rows of all behavior data
  beh = do.call(rbind, dataL)
  return(beh)
} # end function



#################### HAULOUT CSVs ############################################################
#' Load Argos/FastGPS Wildlife Computer Haulout CSV Data
#' 
#' @param files a vector of full path file names to be loaded
#' @param datetime_format A character string. The default for the formatting Wildlife Computers date times is "%H:%M:%OS %d-%b-%Y" but others can be specified by POSIX standard
#' @param tz A character string indicating the time zone (default is Universal Time Coordinated or UTC)
#' 
#' @return Data frame with positional data
WC_read_haulout = function(filenames, datetime_format = '%H:%M:%S %d-%b-%Y', tz = 'UTC') {
  require(dplyr)
  dataL = lapply(X = filenames, FUN = function(x){
    tmp = read.csv(x, stringsAsFactors = F, na.strings = '')
    
    if(nrow(tmp)>=1){
      colnames(tmp) = tolower(colnames(tmp))
      tmp$startmax = as.POSIXct(strptime(tmp$startmax, format = datetime_format, tz = tz))
      tmp$startmin = as.POSIXct(strptime(tmp$startmin, format = datetime_format, tz = tz))
      tmp$endmax = as.POSIXct(strptime(tmp$endmax, format = datetime_format, tz = tz))
      tmp$endmin = as.POSIXct(strptime(tmp$endmin, format = datetime_format, tz = tz))
      
      tmp = tmp[,c('ptt', 'startmin', 'startmax', 'endmin', 'endmax', 'mindurationminutes', 'maxdurationminutes')]
      return(tmp)
    }else{next}
  })
  
  ho = bind_rows(dataL)
  return(ho)
  
}


