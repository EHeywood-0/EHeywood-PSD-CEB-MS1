# Rejection Sampling Functions
# EIH
# 2024-09-30
####################################################################################################################################################
# FUNCTIONS FOR REJECTION SAMPLING OF SSM MODELED DIVE POSITIONS

# From Josh Hatch
# https://github.com/jmhatch-NOAA/aniMotum/blob/master/R/sample_path.R

##' @title rmvnorm_prec
##'
##' @description Sample from MVN parameterized using a sparse precision matrix
##'
##' @details Modified from VAST
##' 
##' @param n Number of observations
##' @param mu Mean vector
##' @param Omega Sparse precision matrix
##' 
##' @keywords internal
##' 
rmvnorm_prec <- function(n = 1, mu, Omega) {
  z = matrix(rnorm(length(mu) * n), ncol = n)
  L = Matrix::Cholesky(Omega, super = TRUE)
  z = Matrix::solve(L, z, system = "Lt") # z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") # z = Pt %*% z
  z = as.vector(z)
  return(mu + z)
}

##' @title sample_paths
##'
##' @description Sample predicted locations, centered on estimated values
##'
##' @param x a \code{ssm} fit object with class `ssm_df`
##' @param n number of samples to draw for each path
##` 
##' @return a `sf` object containing the sampled paths from a \code{ssm} fit object
##' 
##' @export
##' 
##' @examples
##' # sample 10 paths from an estimated track(s)
##' samp_paths <- sample_paths(x = fit, n = 10)
##' 
##' # From Josh Hatch
# https://github.com/jmhatch-NOAA/aniMotum/blob/master/R/sample_path.R
sample_paths <- function(x, n = 1) {
  
  # ensure n > 0
  stopifnot("n must be greater than 0" = (n > 0))
  
  # ensure x is a ssm fit object
  stopifnot("x must be an `ssm_df` fit object" = inherits(x, "ssm_df"))
  
  # ensure user supplies predictions
  stopifnot("sample_path only works on predicted locations, need to set time.step in `fit_ssm`" = all(sapply(X = x$ssm, FUN = function(X) !is.null(X$predicted), simplify = TRUE)))
  
  # loop through tracks
  samp_tracks_sf = purrr::map_dfr(.x = seq_along(x$ssm), .f = function(.x) {
    
    # mles from ssm fit object
    par_mles <- x$ssm[[.x]]$tmb$env$last.par.best
    
    # random effects
    re <- x$ssm[[.x]]$tmb$env$random
    par_re <- par_mles[re]
    h_re <- x$ssm[[.x]]$tmb$env$spHess(par_mles, random = TRUE) ## conditional prec. of u | theta
    
    # CRS
    prj <- sf::st_crs(x$ssm[[.x]]$predicted)
    
    # sample paths
    samp_paths <- purrr::map_dfr(.x = 1:n, .f = function(.y) { 
      
      # sample
      samp_fit <- rmvnorm_prec(n = 1, mu = par_re, Omega = h_re)
      samp_path <- samp_fit[which(names(par_re) == 'mu')]
      
      # re-order based on aniMotum ordering of random effects
      samp_locs <- cbind(samp_path[seq(1, length(samp_path), by = 2)],
                         samp_path[seq(2, length(samp_path), by = 2)])
      samp_locs <- as.data.frame(samp_locs, row.names = 1:nrow(samp_locs))
      
      # extract predicted locs
      samp_locs <- samp_locs[!x$ssm[[.x]]$isd, ]
      colnames(samp_locs) <- c('sample_x', 'sample_y')
      
      # output
      samp_locs |> dplyr::mutate(id = x$ssm[[.x]]$predicted$id,
                                 date = x$ssm[[.x]]$predicted$date,
                                 sample = .y)
      
    })
    
    # coerce to sf object
    sf::st_as_sf(samp_paths, coords = c('sample_x', 'sample_y'), remove = FALSE) |>
      sf::st_set_crs(prj) 
    
  })
  
  # return
  return(samp_tracks_sf)
  
}


# function to extract bathymetry to imputed position
######
##' @title impute_bathy
##'
##' @description help function to sample predicted locations, centered on estimated values
##'
##' @param curssm a \code{ssm} fit object with class `ssm`
##' @param iter the current iteration of 'n' which is the total number of imputations
##' @param prj the projection of the predicted values in x
##' @param par_re
##' @param h_re
##' @param bathymetry a bathymetric raster with spatial coverage of the predicted locations
##' @param dives a dataframe of dives with SegID, diveID, and dive depths / durations
##' @return a `sf` object containing the sampled paths from a \code{ssm} fit object with corresponding bathymetry values
##'
##' @export
##'
##' @examples
##'
impute_bathy = function(curssm, iter, prj, par_re, h_re, bathymetry, dives){
  stopifnot('curssm must be an aniMotum ssm fit object with class `ssm`' = inherits(curssm, 'ssm'))
  #stopifnot('there are no ids in curssm that match ids in dives - check input data' = !curssm$predicted$id[1] %in% unique(dives$SegID))
  
  # sample
  samp_fit <- rmvnorm_prec(n = 1, mu = par_re, Omega = h_re)
  samp_path <- samp_fit[which(names(par_re) == 'mu')]
  
  # re-order based on aniMotum ordering of random effects
  samp_locs <- cbind(samp_path[seq(1, length(samp_path), by = 2)],
                     samp_path[seq(2, length(samp_path), by = 2)])
  samp_locs <- as.data.frame(samp_locs, row.names = 1:nrow(samp_locs))
  
  # extract predicted locs
  samp_locs <- samp_locs[!curssm$isd, ]
  colnames(samp_locs) <- c('sample_x', 'sample_y')
  
  # output
  samp_locs<- samp_locs |> dplyr::mutate(id = curssm$predicted$id,
                                         date = curssm$predicted$date,
                                         sample = iter)
  
  
  # coerce to sf object with WGS84 coordinates
  samp_locs = sf::st_as_sf(samp_locs, coords = c('sample_x', 'sample_y'), remove = FALSE) |>
    sf::st_set_crs(prj) %>% st_transform(4326)
  
  # Get bathy
  bathy_depths = st_extract(bathymetry, samp_locs)
  
  # Bind depth data to data
  samp_locs$bathydepth = bathy_depths$layer
  
  # Merge with depth information contained in 'dives'
  divesub = dives[dives$SegID == unique(samp_locs$id), ]
  samp_locs = left_join(samp_locs, divesub, by = join_by(id == SegID, date))
  
}

######
##' @title sample_path_reject
##'
##' @description Sample predicted locations, centered on estimated values
##'
##' @param x a list of \code{ssm} fit objects with class `ssm`
##' @param n number of samples to draw for each path
##' @param tolerance in meters - allowable error in depth differences between bathy and dive depth where bathy is shallower than dives
##' @param minresample an integer indicating the number of dives that will be thrown out if stableiter is exceeded
##' @param stableiter an integer indicating the maximum number of resample attempts for stable resample lengths
##' 
##' @return a `sf` object containing the sampled paths from a \code{ssm} fit object
##'
##' @export
##'
##' @examples
##' # impute 10 paths from an estimated SSM model, allowing 10 meters error where bathy is shallower than dive depth
##' # and stop, and discard remaining dives if the number to resample (where bathy is shallower) is < 5 and there has been 100 stable iterations of the while loop
##' imputed_paths = sample_path_reject(x = x, bathymetry = bathymetry, dives = dives, n = 10, tolerance = 10, stableiter = 100, minresample = 5)

##'
sample_path_reject <- function(x, bathymetry, dives, n = 1, tolerance = 2, stableiter = 200, progress = FALSE) {
  require(purrr)
  # ensure n > 0
  stopifnot("n must be greater than 0" = (n > 0))
  
  # ensure the contents of list x are a ssm fit objects
  stopifnot("x must be a list containing an 'ssm' fit object" = inherits(x$ssm, "ssm"))
  
  # ensure user supplies predictions
  stopifnot("sample_path only works on predicted locations, need to set time.step in `fit_ssm`" = !is.null(x$predicted))
  
  # ensure that dive depths are negative, representing depth below sea level (0)
  if(all(dives$Depth > 0)){
    dives$Depth = -dives$Depth
  }
  
  ####################################### BEGIN REJECTION SAMPLING ##################################################################
  # mles from ssm fit object
  par_mles <- x$par_mles
  
  # random effects
  par_re <- x$par_re
  
  h_re <- x$h_re ## conditional prec. of u | theta
  
  # CRS
  prj <- sf::st_crs(x$predicted)
  
  # Initialize progress handler
  if (progress){
    p <- progressr::progressor(along = 1:n)
  }
  # impute n paths, iterating over 1:n in another map_dfr()
  samp_paths <- purrr::map_dfr(.x = 1:n, .f = function(.i) {   
    #message('Current n: ', .i)
    
    if (progress){
    p(sprintf("Processing n = %d", .i))
    }
    nsuccess = 0
    last_resamp_length = NULL
    stable_iterations = 0
    
    #Take initial sample using mini helper function impute_bathy
    samp_locs = impute_bathy(curssm = x$ssm, iter = .i, prj, par_re, h_re, bathymetry, dives)
    
    # for this iter ('n') we want one success for each thing, thus we will iterate until we have all sensible imputations
    while(nsuccess < 1) {
      
      # If any of the dive data is deeper than bathy, reject this one (implement some tolerance here?)
      if(any((samp_locs$bathydepth - samp_locs$Depth) > tolerance | is.na(samp_locs$bathydepth))){
        
        resamp = which((samp_locs$bathydepth - samp_locs$Depth) > tolerance | is.na(samp_locs$bathydepth))
        #message('Length Current Resample:', length(resamp))               
        
        # Resample with impute_bathy function
        resamp_locs = impute_bathy(curssm = x$ssm, iter = .i, prj, par_re, h_re, bathymetry, dives)
        
        # Check if length(resamp) is stable
        if (!is.null(last_resamp_length) && length(resamp) == last_resamp_length) {
          stable_iterations = stable_iterations + 1
          #message("Current Stable Iterations:", stable_iterations)
        }else{stable_iterations = 0}
        
        last_resamp_length = length(resamp)  # Update last_resamp_length
        
        # If stable for 500 iterations and length(resamp) < 5, remove those indices and exit
        if (stable_iterations >= stableiter) {
          samp_locs <- samp_locs[-resamp, ]
          return(samp_locs)
        }
        
        # Replace values at resamp indices
        samp_locs[resamp,] = resamp_locs[resamp,]
        
      }else{
        # Advance nsuccess
        nsuccess = nsuccess + 1
      }
      
    } # end while
    # return samp locs
    
    return(samp_locs)
    
  }) # end internal purrr::map_dfr()
  
  
  # return
  return(samp_paths)
  
} # end function



