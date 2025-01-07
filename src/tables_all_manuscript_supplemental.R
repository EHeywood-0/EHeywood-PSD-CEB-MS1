# MS SUPPLEMENTAL TABLES
# EIH
# Updated 2024-09-20

# SETUP
setwd('~/seal_telemetry/')
library(readr)
library(dplyr)
library(flextable)
library(officer)

# SET FLEXTABLE DEFAULTS
set_flextable_defaults(font.size = 10, font.family = 'Arial',
                       font.color = 'black',
                       border.color = 'black',
                       theme_fun = 'theme_booktabs', 
                       padding.bottom = 2, padding.top = 2, 
                       padding.left = 1, padding.right = 1,
                       line_spacing = 1, 
                       text.align = 'left', 
                       table_align = 'center')


sect_properties <- prop_section(
  page_size = page_size(orient = "landscape",
                        width = 8.3, height = 11.7),
  type = "continuous",
  page_margins = page_mar()
)

##########################################################################################################################################
#################### SUPPLEMENTAL TABLE 1: DEPLOYMENT CHARACTERISTICS ####################################################################
##########################################################################################################################################
# LOAD META DATA - # remove repeat row of 177309 
meta = read_csv(file = './data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv') %>%
  distinct(ptt, .keep_all = T)
meta$ptt = as.character(meta$ptt)
meta$totdatadays[meta$ptt == '177509'] = 0
meta$deploydate[meta$ptt == '177509'] = as.POSIXct('2019-01-26', tz = 'UTC')


ST1 = meta[,c('ptt', 'sex', 'deploydate', 'totdatadays','deployloc', 'tag model new', 'masskg', 'lengthcm', 'girthcm')]

ST1$deploydate = format(ST1$deploydate, '%Y-%m-%d')

nrow(ST1)
table(ST1$deployloc)

table(ST1$sex)

colnames(ST1) = c("ptt", 'sex', 'date', 'duration (d)','location', 'model', 'mass (kg)', 'length (cm)', 'girth (cm)')


# Add in info about types of data
ST1$`data type` = 'Argos Only'
ST1$`data type`[which(grepl(pattern = 'SPLASH', x = ST1$model))] = '+Time-Depth'
ST1$`data type`[which(grepl(pattern = '-F', x = ST1$model))] = '+Time-Depth +Fastloc'

# Select final columns
ST1 = ST1 %>% 
  dplyr::select(ptt, model, `data type`, date, location, `duration (d)`, sex, `mass (kg)`, `length (cm)`, `girth (cm)`)
ST1 = arrange(ST1, date)

table(ST1$`data type`)

ST1$sex <- ifelse(ST1$sex == "M", "\u2642", "\u2640")

ft = flextable::flextable(ST1)

# Add superscripts where a deployment failed (1) or was too short (2)
ft = compose(ft, j = 'ptt', i = which(ST1$ptt %in% c('177509', '240186')),
             value = as_paragraph(as_chunk(c('177509', '240186')), as_chunk("a", props = fp_text(vertical.align =  "superscript"))))
ft = compose(ft, j = 'ptt', i = which(ST1$ptt %in% c('194405', '194406')),
             value = as_paragraph(as_chunk( c('194405', '194406')), as_chunk("b", props = fp_text(vertical.align = "superscript"))))

# Add a footer with references for superscripts
ft <- add_footer_lines(ft, values = c(
  "a: Tag failure (N=2)",
  "b: No modeled data due to short deployment durations (N=2)"
))

# Set overarching headers
ft = add_header_row(x=ft, values = c('tag', 'deployment', 'morphometrics'), 
                    colwidths = c(3, 3, 4), top = T)
# Manually add vertical borders to partition the header sections

ft = vline(x = ft, i = c(1,2), part = 'header', j = c(3, 6), border = fp_border(width = 2))
# Center the header text
ft <- align(ft, part = "header", i = c(1,2),align = "center")
ft <- align(ft, part = 'body',align = "center")

ft = autofit(ft)
ft
save_as_docx(ft, path = "./plots/manuscript/TableS1.docx", pr_section = sect_properties)

##################################################################################################################################
################################### SUPPLEMENTAL TABLE 2: TRIP STATISTICS ########################################################
##################################################################################################################################
rm(list = setdiff(ls(), 'sect_properties'))

# READ IN THE CALCULATED PROBABILITY OF EACH TRIP BEING IN A WEA, ORGANIZED INTO TRIPID, SEASONS AND SEX 
# where a trip is classified as occurring in the season which contains most data points for a month in that season
tm = read_csv('./data/L4/Hg-PreConstruction-TripUDs-ProbWEA.csv') %>% dplyr::select(-c(Sex, Season))

# read in UD areas
udareas = read_csv('./data/L4/Hg_PreConstruction_Individual_HOMERANGEAREAS.csv') %>%
  mutate(TripID = tripid) %>% dplyr::select(-c(tripid, Sex, ptt))
udareas = tidyr::pivot_wider(data = udareas, names_from = PercentVolume, values_from = area_sqkm, names_prefix = 'UD_Area_sqkm_')

# READ IN SUMMARIZED TRIP DATA - both spatial and temporal metrics
dat = read_csv("~/seal_telemetry/data/L3/Hg-PreConstruction-SpatialTemporalTripMetrics.csv") %>%
  dplyr::select(-c(TripMonth)) 

dat = left_join(tm, dat, by = 'TripID') %>%  
  left_join(udareas, by = 'TripID')

# READ IN META DATA AND GET PTT DEPLOYDATE AND TOTAL DATA DAYS
meta = read_csv('~/seal_telemetry/data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv') %>% 
  dplyr::select(ptt, deploydate, totdatadays)

dat = left_join(dat, meta, by = join_by(id == ptt))
rm(udareas, tm) 


# trips data by ptt
dat$ptt = as.character(dat$id)
dat$year = format(dat$deploydate, '%Y')

ST2 = dat %>%
  group_by(ptt, sex) %>% arrange(year, sex) %>%
  reframe(n = format(length(unique(TripID)), nsmall = 0),
          `duration` = paste0(format(round(median(TripDur), digits = 1), nsmall = 1), 
                              " (", 
                              format(round(IQR(TripDur), digits = 1), nsmall = 1), ")"),
          `length` = paste0(format(round(median(TripLength_km, na.rm = T), digits = 1), nsmall = 1), 
                            " (", 
                            format(round(IQR(TripLength_km, na.rm = T), digits = 1), nsmall = 1), ")"),
          `dist. offshore` = paste0(format(round(median(MeanDisttoShore, na.rm = T), digits = 1), nsmall = 1), 
                                    " (", 
                                    format(round(IQR(MeanDisttoShore, na.rm=T), digits = 1), nsmall = 1), ")"),
          `core area` = paste0(format(round(median(UD_Area_sqkm_50, na.rm = T), digits = 1), nsmall = 1), 
                               " (", 
                               format(round(IQR(UD_Area_sqkm_50, na.rm = T), digits = 1), nsmall = 1), ")"),
          `home range area` = paste0(format(round(median(UD_Area_sqkm_95, na.rm = T), digits = 1), nsmall = 1), 
                                     " (", 
                                     format(round(IQR(UD_Area_sqkm_95, na.rm = T), digits = 1), nsmall = 1), ")"),
          `n intersect` = format(length(unique(TripID[Probability_in_Wea > 0])), nsmall = 0),
          `probability` = paste0(format(round(median(Probability_in_Wea, na.rm = T), digits = 2), nsmall = 2), 
                                 " (", 
                                 format(round(IQR(Probability_in_Wea, na.rm=T), digits = 2), nsmall = 2), ")")
  )


# where probability is very small make (< 0.01)
ST2$`probability`[which(grepl(pattern = '0.00', x = ST2$`probability`) & ST2$n > 1)] = gsub(pattern = '0.00', 
                                                                                            replacement = '<0.01', 
                                                                                            x = ST2$`probability`[which(grepl(pattern = '0.00', x = ST2$`probability`) & ST2$n > 1)])

ST2$`probability`[which(grepl(pattern = '\\(0.00\\)', x = ST2$`probability`) & ST2$n == 1)] = gsub(pattern = '\\(0.00\\)', 
                                                                                            replacement = '\\(NA\\)', 
                                                                                            x = ST2$`probability`[which(grepl(pattern = '\\(0.00\\)', x = ST2$`probability`) & ST2$n == 1)])
ST2$`probability`[which(grepl(pattern = '0.00 \\(NA\\)', x = ST2$`probability`) & ST2$n == 1)] = gsub(pattern = '0.00', 
                                                                                                   replacement = '<0.01', 
                                                                                                   x = ST2$`probability`[which(grepl(pattern = '0.00 \\(NA\\)', x = ST2$`probability`) & ST2$n == 1)])
ST2[] = apply(X = ST2, MARGIN = 2, FUN = function(x){
  gsub(pattern = '\\(0.0\\)', replacement = '\\(NA\\)', x = x)
})

# Get sex symbols
ST2$sex <- ifelse(ST2$sex == "M", "\u2642", "\u2640")

# Adjust colnames
colnames(ST2) = c('ptt', 'sex', 
                  'n', 'duration (d)', 'length (km)', 'dist. offshore (km)', # trips
                 'core range (km\u00B2)', 'home range (km\u00B2)', #UD
                 'n intersect', 'probability')



ft = flextable::flextable(ST2)
# Add superscripts to indicate median and iqr
# Apply superscripts
ft = compose(ft, i = 1, j = c(4:8), part = 'header',
             value = as_paragraph(as_chunk(names(ST2)[c(4:8)]), as_chunk(c(rep('a', 5)), props = fp_text(vertical.align =  "superscript"))))
ft = compose(ft, i = 1, j = c(10), part = 'header',
             value = as_paragraph(as_chunk(names(ST2)[c(10)]), as_chunk(c(rep('b', 1)), props = fp_text(vertical.align =  "superscript"))))
# Add a footer with references for superscripts
ft <- add_footer_lines(ft, values = c(
  "a: median(± interquartile range)",
  "b: mean(± standard deviation)"
))


# Set overarching headers
ft = add_header_row(x=ft, values = c('', 'trips', 'utilization distributions', 'in WEA'), 
                    colwidths = c(2, 4, 2, 2), top = T)
# Manually add vertical borders to partition the header sections

ft = vline(x = ft, i = c(1,2), part = 'header', j = c(2, 6, 8), border = fp_border(width = 2))
# Center the header text
ft <- align(ft, part = "header", i = c(1,2),align = "center")
ft <- align(ft, part = 'body',align = "center")

ft = autofit(ft)
ft

# Save as word docx
save_as_docx(ft, path = "./plots/manuscript/TableS2.docx", pr_section = sect_properties)


##################################################################################################################################
################################### SUPPLEMENTAL TABLE 3: HAULOUT STATISTICS #####################################################
##################################################################################################################################
rm(list = setdiff(ls(), 'sect_properties'))

ho = read_csv('./data/L1/haulout/Hg_2019-2023_ALL_HauloutMethods_Merged.csv') %>% 
  filter(ho_dur_hours > 0, ptt != 194405, ptt!=194406) 
ho$ptt = as.character(ho$ptt)
ho$year = format(ho$deploydate, '%Y')

# create table of % of total time spent hauled out, average haulout duration, SD, min and max
hotable = ho %>% arrange(year, sex) %>% 
  group_by(ptt, sex) %>% 
  summarise(n = format(length(unique(HID)), nsmall = 0), 
            `tag dur` = format(totdatadays[1], nsmall = 1),
            `total (d)` = format(round(sum(ho_dur_hours)/24, digits = 1), nsmall = 1),
            `duration (h)` = paste0(format(round(median(ho_dur_hours), digits = 0), nsmall = 0), 
                                    " (", 
                                    format(round(IQR(ho_dur_hours), digits = 0), nsmall = 0), ")"),
            `max duration (h)` = format(round(max(ho_dur_hours), digits = 0), nsmall = 0)
            )
hotable$`% time` = format(round(as.numeric(hotable$`total (d)`)/as.numeric(hotable$`tag dur`)*100, 1), nsmall = 1)


# Make nice table
ST3 = hotable %>% 
  dplyr::select(ptt, sex,  
                n, `total (d)`, `% time`, 
                `duration (h)`, `max duration (h)`)

# Get sex symbols
ST3$sex <- ifelse(ST3$sex == "M", "\u2642", "\u2640")


ft = flextable::flextable(ST3)

# Add superscripts to indicate median and iqr
ft = compose(ft,i = 1, j = 6, part = 'header',
             value = as_paragraph(as_chunk(names(ST3)[6]), as_chunk("1", props = fp_text(vertical.align =  "superscript"))))

# Add a footer with references for superscripts
ft <- add_footer_lines(ft, 
                       values = "1: median(± interquartile range)")
# Set overarching headers
ft = add_header_row(x=ft, values = c('', 'haulouts'), 
                    colwidths = c(2, 5), top = T)
# Manually add vertical borders to partition the header sections

ft = vline(x = ft, i = c(1,2), part = 'header', j = 2, border = fp_border(width = 2))
# Center the header text
ft <- align(ft, part = "header", i = c(1,2),align = "center")
ft <- align(ft, part = 'body',align = "center")

ft = autofit(ft)
ft
save_as_docx(ft, path = "~/seal_telemetry/plots/manuscript/TableS3.docx", pr_section = sect_properties)


######### SCRATCH 
#With adult data, create separate relationships for males and females 
#this is not necessary for the pups because they are similar size/weight

#all animals:
# vi_hgm <- lm(masskg ~ lengthcm * I(girthcm^2), dat = hotable)
# summary(vi_hgm) 
# hotable$bmi = residuals(vi_hgm)


# m = lm(`Percent Time` ~ sex, data = hotable[hotable$`Tag Duration (d)`> 30,])
# summary(m)
# ggplot(hotable[hotable$`Tag Duration (d)` > 30,], mapping = aes(x = bmi, y = `Percent Time`)) +
#   geom_point() +
#   geom_smooth(method = 'lm')

