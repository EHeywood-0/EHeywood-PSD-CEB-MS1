# Helper Functions
# EIH
# updated: 2024-09-30

# Load functions for later - from Theo Michelot
#' Split track at gaps
#' 
#' @param data Data frame with (at least) columns for "ID" and "time"
#' @param max_gap Longest allowed gap, in minutes (track will be split at longer gaps)
#' @param shortest_track Shortest track to keep after splitting, in minutes. Shorter
#' tracks will be removed from the output data set.
#' 
#' @return Data frame with identical structure as input, where ID column
#' has been replaced by new ID for split tracks. Old ID still accessible as
#' ID_old column
split_at_gap <- function(data, max_gap = 60, shortest_track = 0) {
  # Number of tracks
  n_tracks <- length(unique(data$ID))
  
  # Save old ID and reinitialise ID column
  data$ID_old <- data$ID
  data$ID <- character(nrow(data))
  
  # Loop over tracks (i.e., over IDs)
  for(i_track in 1:n_tracks) {
    # Indices for this track
    ind_this_track <- which(data$ID_old == unique(data$ID_old)[i_track])
    track_length <- length(ind_this_track)
    
    # Time intervals in min
    dtimes <- difftime(data$time[ind_this_track[-1]], 
                       data$time[ind_this_track[-track_length]],
                       units = "mins")
    
    # Indices of gaps longer than max_gap
    ind_gap <- c(0, which(dtimes > max_gap), track_length)
    
    # Create new ID based on split track
    subtrack_ID <- rep(1:(length(ind_gap) - 1), diff(ind_gap))
    data$ID[ind_this_track] <- paste0(data$ID_old[ind_this_track], "-", subtrack_ID)
  }
  
  # Only keep sub-tracks longer than some duration
  track_lengths <- sapply(unique(data$ID), function(id) {
    ind <- which(data$ID == id)
    difftime(data$time[ind[length(ind)]], data$time[ind[1]], units = "min")
  })
  ID_keep <- names(track_lengths)[which(track_lengths >= shortest_track)]
  data <- subset(data, ID %in% ID_keep)
  
  return(data)
}




#' Output Diagnostic Plots for aniMotum::ssm_fit()
#' plots fitted and predicted locations with model error
#' calculates and plots one-step-ahead residuals, qq and acf plots
#' 
#' @param ssm_df ssm_df output from fit_ssm
#' @param writedir the directory where you want the plots written
#' @param model can be one of 'crw' 'rw' or 'mp'
#' 
plot_ssm_diags = function(ssm_df, writedir, model = 'crw') {
  require(aniMotum)
  require(patchwork)
  require(ggplot2)
  fit = ssm_df
  

  id = fit$id  
  
  f = plot(fit, what = "fitted")
  f = f[[1]] + ggtitle(paste(id, 'fitted', sep = ' - '))
  p = plot(fit, what = "predicted") 
  p = p[[1]] + ggtitle('predicted')
  p2 = plot(fit, what = "predicted", type = 2, alpha = 0.1)
  p2 = p2[[1]] + ggtitle('predicted')
  # OSAR RESIDUALS AND PLOTS
  res = osar(fit)
  r1 = plot(res, type = 'ts')
  r2 = plot(res, type = 'qq')
  r3 = plot(res, type = 'acf')
  
  if(model %in% c('crw', 'rw')){
    
    layout = 'AB
              CD
              EG'
    diagfig = wrap_plots(A = f, B = p, C = p2, D = r1, E = r2, G = r3, 
                         design = layout, 
                         heights = c(1,1,1), widths = c(1,1))
    diagfig
    
    ggsave(filename = paste(writedir,id, '-', model,'-SSM_Diagnostics.jpeg'), 
           plot = diagfig,
           device = 'jpeg',
           dpi = 300,
           width = 10, 
           height = 9)

    }
  
  if(model == 'mp'){
    
    p3 = plot(fit, what = 'predicted', type = 3, normalise = TRUE)
    p4 = map(fit, what = 'predicted', normalise = TRUE, silent = TRUE) + theme_bw() + theme(legend.position = 'none') 
    
    p2 = p2[[1]] + theme(title = element_blank())
    comb =  (p3[[1]] | p4) /
      (r1 | r3) /
      (r2 | p2)
    
    ggsave(filename = paste(writedir,id, '-', model,'-SSM_Diagnostics.jpeg', sep=''), 
           plot = comb,
           device = 'jpeg',
           dpi = 300,
           width = 10, 
           height = 10)
     
  }

 
  
}





#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
######################################################## SCRATCH ################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

## FROM: https://fawda123.github.io/WRTDStidal/reference/goodfit.html
## The goodness of fit measure for quantile regression is estimated as 1 minus the ratio between 
## the sum of absolute deviations in the fully parameterized models and the sum of absolute 
## deviations in the null (non-conditional) quantile model. 
## The values are useful for comparisons between quantile models, 
## but they are not comparable to standard coefficients of determination. 
## The latter is based on the variance of squared deviations, 
# whereas goodness of fit values for quantile regression are based on absolute deviations. 
# Goodness of fit values will always be smaller than R2 values.

#' #' Quantile regression goodness of fit
#' #'
#' #' Calculate quantile regression goodness of fit using residuals and non-conditional residuals
#' #' 
#' #' @param resid numeric vector of residuals from the conditional quantile model
#' #' @param resid_nl numeric vector of residuals from the non-conditional (null) quantile model
#' #' @param tau numeric value from zero to one for the estimated quantile
#' #' 
#' #' @export
#' #' 
#' #' @details The goodness of fit measure for quantile regression is estimated as 1 minus the ratio between the sum of absolute deviations in the fully parameterized models and the sum of absolute deviations in the null (non-conditional) quantile model.  The values are useful for comparisons between quantile models, but they are not comparable to standard coefficients of determination. The latter is based on the variance of squared deviations, whereas goodness of fit values for quantile regression are based on absolute deviations.  Goodness of fit values will always be smaller than R2 values. 
#' #' 
#' #' @return A numeric value from 0 to 1 indicating goodness of fit
#' #' 
#' #' @seealso \code{\link{wrtdsrsd}} for residuals
#' #' 
#' #' @references Koenker, R., Machado, J.A.F. 1999. Goodness of fit and related inference processes for quantile regression. Journal of the American Statistical Association. 94(448):1296-1310.
#' #' 
#' #' @examples
#' #' 
#' #' library(quantreg)
#' #' 
#' #' ## random variables
#' #' x <- runif(100, 0, 10)
#' #' y <- x + rnorm(100)
#' #' 
#' #' ## quantile model
#' #' mod <- rq(y ~ x, tau = 0.5)
#' #' res <- resid(mod)
#' #' 
#' #' ## non-conditional quantile model
#' #' mod_nl <- rq(y ~ 1, tau = 0.5)
#' #' rsd_nl <- resid(mod_nl)
#' #' 
#' #' goodfit(res, rsd_nl, 0.5)
#' #' 
#' #' ## r2 of mean model for comparison
#' #' mod_lm <- lm(y ~ x)
#' #' 
#' #' summary(mod_lm)$r.squared
#' goodfit <- function(resid, resid_nl, tau){
#'   
#'   # minimum sum of deviations
#'   V1 <- resid * (tau - (resid < 0))
#'   V1 <- sum(V1, na.rm = T) 
#'   
#'   # null sum of deviations
#'   V0 <- resid_nl * (tau - (resid_nl < 0))
#'   V0 <- sum(V0, na.rm = T) 
#'   
#'   # explained deviance
#'   out <- 1 - V1/V0
#'   
#'   # exceptions for output
#'   if(any(c(Inf, -Inf) %in% out)) out <- NA
#'   if(V1 > V0) out <- NA
#'   
#'   return(out)
#'   
#' }



#' #' Error catching function for animotum::fit_ssm
#' #' 
#' #' @param aniLocs dataframe of positions formatted for state-space modelling in the animotum R package
#' #' @param model the type of movement model to fit - either correlated random walk(crw) or random walk (rw)
#' #' @param ts Shortest track to keep after splitting, in minutes. Shorter
#' #' tracks will be removed from the output data set.
#' #' 
#' #' @return Data frame with identical structure as input, where ID column
#' #' has been replaced by new ID for split tracks. Old ID still accessible as
#' #' ID_old column
#' try_ssm <- function(aniLocs, model = "crw", ts = 2, vmax = 5) {
#'   
#'   # Inputs
#'   if (!(model %in% c('rw', 'crw')))
#'     stop('model must be one of "crw" or "rw"')
#'   if (!is.numeric(ts))
#'     stop('time step ts must be an integer')
#'   
#'   # Error handling
#'   
#'   # 'tryCatch()' will return the last evaluated expression
#'   # in case the "try" part was completed successfully
#'   fit = tryCatch(
#'     {
#'       # The return value of `fit_ssm()` is the actual value
#'       # that will be returned in case there is no condition
#'       # (e.g. warning or error).
#'       fit_ssm(x = aniLocs,
#'               spdf = TRUE, # pre-filtering on
#'               vmax = vmax, # max speed 5m/s
#'               model = model, #
#'               time.step = ts,
#'               control = ssm_control(verbose = 1)
#'       )
#'     },
#'     
#'     error = function(e) {
#'       message(aniLocs$id[1], "\nerror:\n", e)
#'       message("Here's the original error message:")
#'       message(conditionMessage(e))
#'       # Choose a return value in case of error
#'       'FAILED'
#'     },
#'     
#'     warning = function(w) {
#'       message(paste("SegID caused a warning:", aniLocs$id[1]))
#'       message("Here's the original warning message:")
#'       message(aniLocs$id[1], "\nwarning:\n", w)
#'       # Choose a return value in case of warning
#'       'WARNING'
#'     },
#'     
#'     finally = {
#'       message(paste("Processed SegID:", aniLocs$id[1]))
#'     }
#'   )
#'   return(fit)
#' }
# SUMMARY FUNCTION
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)


# summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
#                       conf.interval=.95, .drop=TRUE) {
#   library(plyr)
#   
#   # New version of length which can handle NA's: if na.rm==T, don't count them
#   length2 <- function (x, na.rm=FALSE) {
#     if (na.rm) sum(!is.na(x))
#     else       length(x)
#   }
#   
#   # This does the summary. For each group's data frame, return a vector with
#   # N, mean, and sd
#   datac <- ddply(data, groupvars, .drop=.drop,
#                  .fun = function(xx, col) {
#                    c(N    = length2(xx[[col]], na.rm=na.rm),
#                      mean = mean   (xx[[col]], na.rm=na.rm),
#                      sd   = sd     (xx[[col]], na.rm=na.rm)
#                    )
#                  },
#                  measurevar
#   )
#   
#   # Rename the "mean" column    
#   datac <- rename(datac, c("mean" = measurevar))
#   
#   datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
#   
#   # Confidence interval multiplier for standard error
#   # Calculate t-statistic for confidence interval: 
#   # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#   ciMult <- qt(conf.interval/2 + .5, datac$N-1)
#   datac$ci <- datac$se * ciMult
#   
#   return(datac)
# }

