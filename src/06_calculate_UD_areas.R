# Calculate UD areas 
library(raster)
#library(terra)
library(dplyr)
setwd('~/seal_telemetry/')


# DEFINE VOLUME CONTOUR FUNCTION
# ---- roxygen documentation ----
#' @title Volume contour from Raster
#'
#' @description
#'   Compute a percent volume contour polygon from a raster UD.
#' @details
#'   The volras function is a simpler version of the getvolumeUD function from the package \code{adehabitatHR} developed by C. Calenge. It allows the output to be a 'raster looking' polygon (i.e., the cells that are within the UD) or a simplified (smoothed) polygon.
#'   
#' @param x a \code{RasterLayer}
#' @param percent a percent value to get the volume contour, e.g., 95. Note: This is a simple function and only accepts one value at a time.
#' @param simplify (logical; default = TRUE) whether or not to simplify the output home range polygon using \code{gSimplify} from \code{rgeos} with a tolerance value of 1.5 times the spatial resolution of the UD. 
#' 
#' @return
#'   A \code{SpatialPolygonsDataFrame}.
#'
#' @seealso fbtgUD, rspUD, tgkde
#' @examples
#' data(m3)
#' ud <- tgkde(m3,disfun='inv',method='vanderWatt')
#' raster::plot(ud)
#' hr <- volras(ud,95)
#' sp::plot(hr,add=TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----

volras <- function(x,percent=95, simplify=TRUE, returnwhat = 'both'){
  
  
  require(smoothr)
  require(sf)
  require(raster)
  
  x[is.na(x)] <- 0
  pfs <- raster::crs(x)
  
  ## standardize it so that the total volume is 1 over the area
  v <- as.vector(values(x))
  index<-1:length(v)
  vord<-v[order(v, decreasing=TRUE)]
  vsu<-cumsum(vord)
  
  cont <- which(vsu > (percent/100)*max(vsu))[1]
  cont.lev <- vord[cont]
  
  #Get all the cells above the cont.lev and make 1
  m <- c(0, cont.lev, 0, cont.lev, max(x[]), 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  xr <- raster::reclassify(x, rclmat)
  
  #Convert to polygon and simplify if desired
  if (returnwhat %in% c('polygon', 'both')){
    hr <- st_as_sf(raster::rasterToPolygons(xr,
                                            fun=function(x){x==1},
                                            n = 4, dissolve=TRUE, 
                                            digits = 4))
    
    if(simplify){
      hr <- smoothr::smooth(x = hr, method = 'ksmooth', smoothness = 5)
    }
    
    if(returnwhat == 'polygon'){
      return(hr)
    }
    
  }
  
  if (returnwhat %in% c('area', 'both')){
    
    # sum reclassified raster to get the number of cells in contour
    numcells = as.numeric(raster::cellStats(x = xr, stat = 'sum', na.rm=T))
    resolution = res(xr)
    resolutionsqkm = round(resolution[1] * resolution[2] / 1000000, digits = 0)
    areasqkm = numcells * resolutionsqkm
    
    if(returnwhat == 'area'){
      return(areasqkm)
    }
    
  }
  
  if(returnwhat == 'both'){
    return(list(hr, areasqkm))
  }
  
}


seasons = c('Winter', "Spring", 'Summer', 'Fall')

# Define which volume contours
perc = c(95, 50)


# GET AREAS FOR FEMALES
for (season in 1:length(seasons)){
  
  FStackR = raster::stack(x = paste0('./data/L3/UD/UD_rasters/Hg_PreConstruction_FemaleUDs_', seasons[season],'.grd'))
  tmpL = list()
  
  # Calculate volume rasters for each layer in the stack to get individual 
  for (i in 1:nlayers(FStackR)){
    print(seasons[season])
    print(paste(i, ' of ', nlayers(FStackR)))
    curast = FStackR[[i]]
    
    res = lapply(X = 1:length(perc), FUN = function(x) {
      tmp = volras(x = curast, percent = perc[x], simplify = TRUE, returnwhat = 'area')
      return(tmp)
    })
    
    asqkm = unlist(res)
    tmp = data.frame(area_sqkm = as.numeric(asqkm))
    tmp$PercentVolume = as.character(perc)
    tmp$Sex = 'Female'
    tmp$Season = seasons[season]
    tmp$tripid = gsub(pattern = 'X', replacement = '', x = names(FStackR)[i])
    tmp$tripid = gsub(pattern = '\\.', replacement = '-', x = tmp$tripid)
    
    tmp$ptt = sapply(strsplit(tmp$tripid, '-'), '[[', 1)
    tmp = tmp[,c(6, 5, 4, 3, 2, 1)]
    tmp
    tmpL[[i]] = tmp
    }
  
  indPolys = bind_rows(tmpL)
  
  if(season == 1){
    areas = indPolys
  }else{
    areas = bind_rows(areas, indPolys)
  }
}

# Female UD areas
Fareas = areas

###############################################################################################
###############################################################################################
###############################################################################################
# NOW FOR MALES 
# GET AREAS FOR FEMALES
for (season in 1:length(seasons)){
  
  StackR = raster::stack(x = paste0('./data/L3/UD/UD_rasters/Hg_PreConstruction_MaleUDs_', seasons[season],'.grd'))
  
  tmpL = vector('list', nlayers(StackR))
  
  # Calculate volume rasters for each layer in the stack to get individual 
  for (i in 1:nlayers(StackR)){
    print(seasons[season])
    print(paste(i, ' of ', nlayers(StackR)))
    curast = StackR[[i]]
    
    res = lapply(X = 1:length(perc), FUN = function(x) {
      tmp = volras(x = curast, percent = perc[x], simplify = TRUE, returnwhat = 'area')
      return(tmp)
    })
    
    asqkm = unlist(res)
    tmp = data.frame(area_sqkm = as.numeric(asqkm))
    tmp$PercentVolume = as.character(perc)
    tmp$Sex = 'Male'
    tmp$Season = seasons[season]
    tmp$tripid = gsub(pattern = 'X', replacement = '', x = names(StackR)[i])
    tmp$tripid = gsub(pattern = '\\.', replacement = '-', x = tmp$tripid)
    
    tmp$ptt = sapply(strsplit(tmp$tripid, '-'), '[[', 1)
    tmp = tmp[,c(6, 5, 4, 3, 2, 1)]
    tmp
    tmpL[[i]] = tmp
  }
  
  indPolys = bind_rows(tmpL)
  
  if(season == 1){
    areas = indPolys
  }else{
    areas = bind_rows(areas, indPolys)
  }
}

Mareas = areas


# Together 
HR_AREAS = bind_rows(Mareas, Fareas)

library(readr)
write_csv(x = HR_AREAS,file =  './data/L4/Hg_PreConstruction_Individual_HOMERANGEAREAS.csv')



# meta
meta = readr::read_csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv") %>% dplyr::select(ptt, lengthcm, masskg, girthcm, cutoffstart, cutoffend)
meta$ptt = as.character(meta$ptt)
HR_AREAS = left_join(HR_AREAS, meta, 'ptt')

library(ggplot2)
HR_AREAS$Season = factor(HR_AREAS$Season, levels = c('Winter', 'Spring', 'Summer'))
ggplot(HR_AREAS, mapping = aes(x = Season, y = area_sqkm, fill = Sex)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~PercentVolume, scales = 'free')

# 
# library(ggpubr)
# library(rstatix)
# 
# HR_AREAS %>%
#   group_by(Season, Sex, PercentVolume) %>%
#   get_summary_stats(area_sqkm, type = "mean_sd")
# 
# ggqqplot(HR_AREAS[HR_AREAS$PercentVolume == '95',], "area_sqkm", ggtheme = theme_bw()) +
#   facet_grid(Season ~ Sex)
# ggqqplot(HR_AREAS[HR_AREAS$PercentVolume == '50',], "area_sqkm", ggtheme = theme_bw()) +
#   facet_grid(Season ~ Sex)
# 
# 
# # Two-way mixed ANOVA test
# library(geepack)
# library(patchwork)
# library(ggsignif)
# 
# HR_AREAS$Year = factor(format(HR_AREAS$cutoffstart, '%Y'))
# HR_AREAS$BCI = HR_AREAS$masskg / (HR_AREAS$lengthcm * HR_AREAS$girthcm)
# 
# # Fit a linear regression model
# bcimod <- lm(masskg ~ lengthcm + girthcm, data = HR_AREAS)
# 
# # Get the residuals (RBCI)
# HR_AREAS$RBCI <- resid(bcimod)
# 
# a50 = HR_AREAS[HR_AREAS$PercentVolume == '50',]
# a95 = HR_AREAS[HR_AREAS$PercentVolume == '95',]
# 
# a50 = a50 %>% arrange(ptt, tripid)
# a95 = a95 %>% arrange(ptt, tripid)
# 
# # Fit GEE with a rank-based method for non-parametric analysis
# model <- geeglm(area_sqkm ~ Season + girthcm, 
#                 id = ptt, 
#                 data = a95, 
#                 family = Gamma(link = 'identity'), 
#                 corstr = 'ar1')
# 
# # Summarize the results
# summary(model)
# plot(model)
# QIC(model)
# # Basic boxplot of area_sqkm by sex and season
# p <- ggplot(a95, aes(x = Sex, y = area_sqkm)) +
#   geom_boxplot(outliers = T, alpha = 0.5) +
#   #scale_fill_manual(values = c('slateblue', 'burlywood2'))+
#   scale_y_continuous(breaks = c(25000, 50000, 75000, 100000, 125000, 150000), labels = c('25', '50', '75','100', '125', '150')) +
#   ylab(expression("95% UD Area (thousands km"^2*")")) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         legend.position = 'None', text = element_text(size = 16))
# p
# 
# # Obtain estimated marginal means for the interaction of Sex and Season
# emmeans_combined <- emmeans(model, ~ Sex | Season)
# 
# # Perform pairwise comparisons
# pairwise_combined <- pairs(emmeans_combined)
# pairwise_combined_summary <- summary(pairwise_combined)
# 
# # Display the pairwise comparisons
# print(pairwise_combined_summary)
# 
# library(ggsignif)
# 
# # Determine which comparisons you want to annotate, based on the results of pairwise_combined_summary
# # For example, you may want to compare Male vs. Female within each Season
# 
# p1 <- p + geom_signif(
#   y_position = c(max(a95$area_sqkm) * 0.95
#                  #max(a95$area_sqkm) * 0.75, 
#                  #max(a95$area_sqkm) * 0.75
#                  ), 
#   xmin = 1,
#   xmax = 2,
#   annotations = c(ifelse(pairwise_combined_summary$p.value[1] < 0.05, "***", "ns")))
# 
# print(p1)
# 
# p2 = ggplot(a95, aes(x = Season, y = area_sqkm)) +
#   geom_boxplot(outliers = T, alpha = 0.5) +
#   scale_y_continuous(breaks = c(25000, 50000, 75000, 100000, 125000, 150000), labels = c('25', '50', '75','100', '125', '150')) +
# #  labs(y = "95% UD (thousands sq. km)") +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         legend.position = 'None', text = element_text(size = 16))
# 
# p2
# # Obtain estimated marginal means for the interaction of Sex and Season
# emmeans_combined <- emmeans(model, ~ Season | Sex)
# 
# # Perform pairwise comparisons
# pairwise_combined <- pairs(emmeans_combined)
# pairwise_combined_summary <- summary(pairwise_combined)
# 
# # Display the pairwise comparisons
# print(pairwise_combined_summary)
# 
# p2 = p2 + geom_signif(
#   y_position = c(max(a95$area_sqkm) * 0.875), 
#   xmin = c(1),
#   xmax = c(2),
#   annotations = c('*'), colour = 'black') +
#   geom_signif(y_position = c(max(a95$area_sqkm) * 0.95), 
#               xmin = c(1),
#               xmax = c(3),
#               annotations = c('*'), colour = 'black')
# 
# p2
# p95 = p1 | p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# p95
# #####################################################################################################################################
# # 50% UD
# 
# # Fit GEE with a rank-based method for non-parametric analysis
# model <- geeglm(area_sqkm ~ Season + girthcm, 
#                 id = ptt, 
#                 data = a50, 
#                 family = Gamma(link = 'identity'), 
#                 corstr = 'ar1')
# 
# # Summarize the results
# summary(model)
# plot(model)
# QIC(model)
# # Basic boxplot of area_sqkm by sex and season
# p <- ggplot(a50, aes(x = Sex, y = area_sqkm)) +
#   geom_boxplot(outliers = T, alpha = 0.5) +
#   #scale_fill_manual(values = c('slateblue', 'burlywood2'))+
#   scale_y_continuous(breaks = c(5000, 10000, 15000, 20000, 25000, 30000), labels = c('5', '10', '15','20', '25', '30')) +
#   ylab(expression("50% UD Area (thousands km"^2*")")) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         legend.position = 'None', text = element_text(size = 16))
# p
# 
# # Obtain estimated marginal means for the interaction of Sex and Season
# emmeans_combined <- emmeans(model, ~ Sex | Season)
# 
# # Perform pairwise comparisons
# pairwise_combined <- pairs(emmeans_combined)
# pairwise_combined_summary <- summary(pairwise_combined)
# 
# # Display the pairwise comparisons
# print(pairwise_combined_summary)
# 
# 
# # Determine which comparisons you want to annotate, based on the results of pairwise_combined_summary
# # For example, you may want to compare Male vs. Female within each Season
# 
# p1 <- p + geom_signif(
#   y_position = c(max(a50$area_sqkm) * 0.95), 
#   xmin = 1,
#   xmax = 2,
#   annotations = c(ifelse(pairwise_combined_summary$p.value[1] < 0.05, "***", "ns")))
# 
# print(p1)
# 
# p2 = ggplot(a50, aes(x = Season, y = area_sqkm)) +
#   geom_boxplot(outliers = T, alpha = 0.5) +
#   scale_y_continuous(breaks = c(5000, 10000, 15000, 20000, 25000, 30000), labels = c('5', '10', '15','20', '25', '30')) +
#   #  labs(y = "95% UD (thousands sq. km)") +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         legend.position = 'None', text = element_text(size = 16))
# 
# p2
# # Obtain estimated marginal means for the interaction of Sex and Season
# emmeans_combined <- emmeans(model, ~ Season | Sex)
# 
# # Perform pairwise comparisons
# pairwise_combined <- pairs(emmeans_combined)
# pairwise_combined_summary <- summary(pairwise_combined)
# 
# # Display the pairwise comparisons
# print(pairwise_combined_summary)
# 
# p2 = p2 + geom_signif(
#   y_position = c(max(a50$area_sqkm) * 0.95), 
#   xmin = c(1),
#   xmax = c(2),
#   annotations = c('*'), colour = 'black')
# p2
# p50 = p1 | p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# p50
# 
# 
# comb = p50 /
#   p95
# comb
# 
# ggsave(filename = './plots/manuscript/Hg_PreConstruction_GEEGLM_UD_PLOT.png', 
#        plot = comb, device = 'png', dpi = 300, width = 8, height = 8, scale = 1)

# WRITE OUT MODEL COMPARISON TABLE 
