# MAIN MANUSCRIPT TABLES
# EIH
# Updated 2024-09-20

# SETUP
setwd('~/seal_telemetry/')
library(readr)
library(dplyr)
library(tidyr)
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

#############################################################################################################################
#################### TABLE 1: DEPLOYMENT CHARACTERISTICS ####################################################################
#############################################################################################################################
# LOAD META DATA - # remove repeat row of 177309 
meta = read_csv(file = './data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv') %>%
  distinct(ptt, .keep_all = T)
meta$ptt = as.character(meta$ptt)
meta$totdatadays[meta$ptt == '177509'] = 0
meta$deploydate[meta$ptt == '177509'] = as.POSIXct('2019-01-26', tz = 'UTC')

meta$year = format(meta$deploydate, '%Y')

meta = meta[,c('year', 'ptt', 'sex', 'totdatadays', 'masskg', 'lengthcm', 'girthcm')]

nrow(meta)

table(meta$sex)

#colnames(meta) = c('year', 'ptt', 'sex', 'tag duration (d)', 'mass (kg)', 'length (cm)', 'girth (cm)')
meta$sex <- ifelse(meta$sex == "M", "\u2642", "\u2640")

# Summarize T1 by year and by sex - excluding those that failed to transmit
meta$year = factor(meta$year)
meta$sex = factor(meta$sex)

# rounding number
rn = 1
ns = 1
T1 = meta %>% group_by(year, sex, .drop = F) %>%
  summarise(n = n(),
            `mass (kg)` = paste0(format(round(mean(masskg, na.rm = T), digits = rn), nsmall = ns), 
                                 ' (', 
                                 format(round(sd(masskg), digits = rn), nsmall=ns), ')'),
            `length (cm)` = paste0(format(round(mean(lengthcm, na.rm = T), digits = rn), nsmall = ns), 
                                   ' (',
                                   format(round(sd(lengthcm), digits = rn), nsmall = ns), ')'),
            `girth (cm)` = paste0(format(round(mean(girthcm, na.rm = T), digits = rn), nsmall = ns), 
                                  ' (',
                                  format(round(sd(girthcm), rn), nsmall = ns), ')'),
            `tag duration (d)` = paste0(format(round(mean(totdatadays[!ptt %in% c("240186", "177509")], na.rm = T), digits = rn), nsmall = ns), 
                                        ' (',
                                        format(round(sd(totdatadays[!ptt %in% c("240186", "177509")]), digits = rn), nsmall = ns), ')'))


# Calculate totals (mean across sexes) for each variable
totalrow = meta %>% group_by(sex, .drop = F) %>%
  summarise(year = 'total',
            n = n(),
            `mass (kg)` = paste0(format(round(mean(masskg, na.rm = T), digits = rn), nsmall = ns), 
                                 ' (', 
                                 format(round(sd(masskg), digits = rn), nsmall=ns), ')'),
            `length (cm)` = paste0(format(round(mean(lengthcm, na.rm = T), digits = rn), nsmall = ns), 
                                   ' (',
                                   format(round(sd(lengthcm), digits = rn), nsmall = ns), ')'),
            `girth (cm)` = paste0(format(round(mean(girthcm, na.rm = T), digits = rn), nsmall = ns), 
                                  ' (',
                                  format(round(sd(girthcm), rn), nsmall = ns), ')'),
            `tag duration (d)` = paste0(format(round(mean(totdatadays[!ptt %in% c("240186", "177509")], na.rm = T), digits = rn), nsmall = ns), 
                                        ' (',
                                        format(round(sd(totdatadays[!ptt %in% c("240186", "177509")]), digits = rn), nsmall = ns), ')'))

# Append the total row to the summarized dataframe

T1 <- bind_rows(T1, totalrow)

T1[] <- lapply(T1, function(x) gsub("\\(NA\\)", "", x))

T1[] <- lapply(T1, function(x) gsub("NaN ", "NA",  x))



ft = flextable::flextable(T1)

ft = compose(ft, j = 3, i = 1,
             value = as_paragraph(as_chunk(c("5")), as_chunk("b", props = fp_text(vertical.align =  "superscript"))))
ft = compose(ft, j = 3, i = 10,
             value = as_paragraph(as_chunk(c("10")), as_chunk("b", props = fp_text(vertical.align =  "superscript"))))


ft = compose(ft, i = 1, j = c(4:7), part = 'header',
             value = as_paragraph(as_chunk(names(T1)[c(4:7)]), as_chunk(c(rep('a', 4)), props = fp_text(vertical.align =  "superscript"))))

# Add a footer with references for superscripts
ft <- add_footer_lines(ft, values = c(
  "a: mean(± standard deviation) where n>1",
  "b: tag failed to transmit and was removed for calculating mean tag duration (total = 1F, 1M)"
))


# Merge the cells in the "year" column where necessary to nest the rows by sex
ft <- merge_v(ft, j = "year")




# Manually add vertical borders to partition the header sections
ft <- theme_booktabs(ft)
ft = autofit(ft)

ft = hline(x = ft, i = c(2,4,6,8, 10), border = fp_border(width = 2))
# Center the header text
ft <- align(ft, part = "header", i = 1,align = "center")
ft <- align(ft, part = 'body',align = "center")
ft <- bold(ft, part = "header", i = 1, bold = TRUE)
ft <- bold(ft, part = "body", i = c(11,12), j = c(1,3:7),bold = TRUE)


ft

save_as_docx(ft, path = "./plots/manuscript/Table1.docx", pr_section = sect_properties)

#######################################################################################################################
#################### TABLE 2: TRIP CHARACTERISTICS ####################################################################
#######################################################################################################################
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

# SEASONAL TRACKING DATA DURS FOR HAULOUT PERCENT TIME CALCS
# READ IN HAULOUT DATA - ptts 194405 and 194406 did not produce any trip data so eliminate
ho = read_csv('./data/L1/haulout/Hg_2019-2023_ALL_HauloutMethods_Merged.csv') %>% 
  filter(ho_dur_hours > 0, ptt != 194405, ptt!=194406) %>% mutate(Month = format(startho, '%b')) %>% filter(Month != 'Dec') %>%
  mutate(Season = case_when(Month %in% c('Jan','Feb') ~ 'Winter',
                            Month %in% c('Mar', 'Apr', 'May') ~ 'Spring',
                            Month %in% c('Jun', 'Jul', 'Aug') ~ 'Summer',
                            Month %in% c('Sep', 'Oct', 'Nov') ~ 'Fall'))

ho$Year = format(ho$cutoffstart, '%Y')

ho = ho %>% group_by(Season, ptt) %>% 
  mutate(SeasonStart = case_when(Season == 'Winter' ~ as.POSIXct(paste0(Year, '-01-01 00:00:00'), tz = 'UTC'),
                                 Season == 'Spring' ~ as.POSIXct(paste0(Year, '-03-01 00:00:00'), tz = 'UTC'),
                                 Season == 'Summer' ~ as.POSIXct(paste0(Year, '-06-01 00:00:00'), tz = 'UTC'),
                                 Season == 'Fall' ~   as.POSIXct(paste0(Year, '-09-01 00:00:00'), tz = 'UTC')),
         SeasonEnd = case_when(Season == 'Winter' ~ as.POSIXct(paste0(Year, '-03-01 00:00:00'), tz = 'UTC')-1,
                               Season == 'Spring' ~ as.POSIXct(paste0(Year, '-06-01 00:00:00'), tz = 'UTC')-1,
                               Season == 'Summer' ~ as.POSIXct(paste0(Year, '-09-01 00:00:00'), tz = 'UTC')-1,
                               Season == 'Fall' ~   as.POSIXct(paste0(Year, '-11-01 00:00:00'), tz = 'UTC')-1))


seasonaltrackingdays = ho %>%
  group_by(Season, ptt) %>% 
  mutate(SeasonalTrackingDuration = case_when(cutoffend >= SeasonEnd & cutoffstart >= SeasonStart ~ difftime(SeasonEnd, cutoffstart, 'days'),
                                              cutoffend >= SeasonEnd & cutoffstart < SeasonStart ~ difftime(SeasonEnd, SeasonStart, 'days'),
                                              cutoffend < SeasonEnd & cutoffstart >= SeasonStart ~ difftime(cutoffend, cutoffstart, 'days'),
                                              cutoffend < SeasonEnd & cutoffstart < SeasonStart ~ difftime(cutoffend, SeasonStart, 'days'))) %>% 
  dplyr::select(ptt, Season, SeasonalTrackingDuration) %>% distinct()

seasonaltrackingdays$SeasonalTrackingDuration = as.numeric(seasonaltrackingdays$SeasonalTrackingDuration)


ho =  left_join(ho, seasonaltrackingdays, by = join_by(Season, ptt))




ho_sum = ho %>%
  group_by(sex, ptt, Season) %>%
  mutate(TotTimeHO = sum(ho_dur_hours)/24, 
         PercTimeHO = (sum(ho_dur_hours)/24) / SeasonalTrackingDuration[1]*100) %>% ungroup() %>% 
  group_by(Season, sex) %>%
  summarise(`Median Haulout Dur` = format(round(median(ho_dur_hours), digits = 1), nsmall = 1), 
            `Mean % Time HO` = format(round(mean(PercTimeHO), 1), nsmall = 1),
            iqrhodur = format(round(IQR(ho_dur_hours), digits = 1), nsmall = 1), 
            sdperctimeho = format(round(sd(PercTimeHO), digits = 1), nsmall = 1))
rm(meta)


# FILTER OUT NON-SSM TRIPS
gsum = dat %>% 
  group_by(Season, sex) %>%
  summarise(N = length(unique(id)),
            `Total Trips` = length(unique(TripID)),
            `Trip Duration Days` = paste0(format(round(median(TripDur), digits = 1), nsmall = 1), 
                                          " (", 
                                          format(round(IQR(TripDur), digits = 1), nsmall = 1), ")"),
            `Trip Length (km)` = paste0(format(round(median(TripLength_km, na.rm = T), digits = 1), nsmall = 1), 
                                        " (", 
                                        format(round(IQR(TripLength_km, na.rm = T), 1), nsmall = 1), ")"),
            `Trip 95% UD Area (sq.km)` = paste0(format(round(median(UD_Area_sqkm_95, na.rm = T), digits = 0), nsmall = 0), 
                                                " (", 
                                                format(round(IQR(UD_Area_sqkm_95, na.rm = T), digits = 0), nsmall = 0), ")"),
            `Trip 50% UD Area (sq.km)` = paste0(format(round(median(UD_Area_sqkm_50, na.rm = T), digits = 0), nsmall = 0), 
                                                " (", 
                                                format(round(IQR(UD_Area_sqkm_50, na.rm = T), digits = 0), nsmall = 0), ")"),
            
            `Distance to Shore (km)` = paste0(format(round(median(MeanDisttoShore, na.rm = T), digits = 1), nsmall = 1), 
                                              " (", 
                                              format(round(IQR(MeanDisttoShore, na.rm=T), digits = 1), nsmall = 1), ")"),
            `N Trips in WEA` = paste0(format(round(length(unique(TripID[Probability_in_Wea > 0])), digits = 0)), 
                                      ' (',
                                      format(round(length(unique(TripID[Probability_in_Wea > 0])) / length(unique(TripID))*100, digits = 0)),
                                      "%)"),
            `Prob. of Trip in WEA` = paste0(format(round(median(Probability_in_Wea, na.rm = T), digits = 2), nsmall = 2), 
                                            " (", 
                                            format(round(IQR(Probability_in_Wea, na.rm=T), digits = 2), nsmall = 2), ")")
  ) %>% arrange(Season, sex)

# merge with ho summary
hosum = tibble(Season = ho_sum$Season, 
               sex = ho_sum$sex, 
               `Median Haulout Dur` = paste(ho_sum$`Median Haulout Dur`, " (",ho_sum$iqrhodur, ")", sep = ""),
               `Mean % Time HO` = paste(ho_sum$`Mean % Time HO`, " (",ho_sum$sdperctimeho, ")", sep=""))
gsum = left_join(gsum, hosum, by = c('Season', 'sex'))


tot = dat %>% ungroup() %>% group_by(sex) %>%
  summarise(N = length(unique(id)),
            `Total Trips` = length(unique(TripID)),
            `Trip Duration Days` = paste0(format(round(median(TripDur), digits = 1), nsmall = 1), 
                                          " (", 
                                          format(round(IQR(TripDur), digits = 1), nsmall = 1), ")"),
            `Trip Length (km)` = paste0(format(round(median(TripLength_km, na.rm = T), digits = 1), nsmall = 1), 
                                        " (", 
                                        format(round(IQR(TripLength_km, na.rm = T), 1), nsmall = 1), ")"),
            `Trip 95% UD Area (sq.km)` = paste0(format(round(median(UD_Area_sqkm_95, na.rm = T), digits = 0), nsmall = 0), 
                                                " (", 
                                                format(round(IQR(UD_Area_sqkm_95, na.rm = T), digits = 0), nsmall = 0), ")"),
            `Trip 50% UD Area (sq.km)` = paste0(format(round(median(UD_Area_sqkm_50, na.rm = T), digits = 0), nsmall = 0), 
                                                " (", 
                                                format(round(IQR(UD_Area_sqkm_50, na.rm = T), digits = 0), nsmall = 0), ")"),
            
            `Distance to Shore (km)` = paste0(format(round(median(MeanDisttoShore, na.rm = T), digits = 1), nsmall = 1), 
                                              " (", 
                                              format(round(IQR(MeanDisttoShore, na.rm=T), digits = 1), nsmall = 1), ")"),
            `N Trips in WEA` = paste0(format(round(length(unique(TripID[Probability_in_Wea > 0])), digits = 0)), 
                                      ' (',
                                      format(round(length(unique(TripID[Probability_in_Wea > 0])) / length(unique(TripID))*100, digits = 0)),
                                      "%)"),,
            `Prob. of Trip in WEA` = paste0(format(round(median(Probability_in_Wea, na.rm = T), digits = 2), nsmall = 2), 
                                            " (", 
                                            format(round(IQR(Probability_in_Wea, na.rm=T), digits = 2), nsmall = 2), ")")
            
  )

totho = ho %>% group_by(sex, ptt) %>%
  mutate(PercTimeHO = (sum(ho_dur_hours)/24) / totdatadays[1]*100) %>% 
  ungroup() %>% group_by(sex) %>%
  summarise(`Median Haulout Dur` = paste(format(round(median(ho_dur_hours),digits = 1), nsmall = 1), 
                                         " (",
                                         format(round(IQR(ho_dur_hours),digits = 1), nsmall = 1), ")", sep = ""),
            `Mean % Time HO` = paste(format(round(mean(PercTimeHO),digits = 1), nsmall = 1), 
                                     " (",
                                     format(round(sd(PercTimeHO),digits = 1), nsmall = 1), ")", sep=""))

tot = left_join(tot, totho, 'sex')

gsum = bind_rows(gsum, tot)
gsum$Season[9:10] = 'total'
gsum$Season = factor(gsum$Season, levels = c('Winter', 'Spring', 'Summer', 'Fall', 'total'), labels = c('winter', 'spring', 'summer', 'fall', 'total'))
gsum = arrange(gsum, Season, sex)

# where probability is very small make (< 0.01)
gsum$`Prob. of Trip in WEA`[gsum$`Prob. of Trip in WEA` %in% c('0.00 (0.00)')] = '<0.01'
gsum$`Prob. of Trip in WEA`[grepl(pattern = '0.00', x = gsum$`Prob. of Trip in WEA`)] = gsub(pattern = '0.00', replacement = '<0.01', x = gsum$`Prob. of Trip in WEA`[grepl(pattern = '0.00', x = gsum$`Prob. of Trip in WEA`)])


# RENAME AND ARRANGE COLUMNS
T2 = gsum %>% dplyr::select(Season, sex, N, `Total Trips`, `Trip Duration Days`, 
                              `Trip Length (km)`, `Distance to Shore (km)`,
                              `Trip 50% UD Area (sq.km)`, `Trip 95% UD Area (sq.km)`,`N Trips in WEA`, `Prob. of Trip in WEA`, `Median Haulout Dur`, `Mean % Time HO`)

colnames(T2) = c('season', 'sex', 'N', 
                 'n', 'duration (d)', 'length (km)', 'dist. offshore (km)', # trips
                 'core range (km\u00B2)', 'home range (km\u00B2)', #UD
                 'n trips', 'probability', # wea
                 'duration (h)', '% time') # haulout
T2$sex <- ifelse(T2$sex == "M", "\u2642", "\u2640")


T2$sex = factor(T2$sex)

# Make flex table
ft = flextable::flextable(T2)

# Apply superscripts
ft = compose(ft, i = 1, j = c(5:9, 11:12), part = 'header',
             value = as_paragraph(as_chunk(names(T2)[c(5:9, 11:12)]), as_chunk(c(rep('a', 7)), props = fp_text(vertical.align =  "superscript"))))
ft = compose(ft, i = 1, j = c(13), part = 'header',
             value = as_paragraph(as_chunk(names(T2)[c(13)]), as_chunk(c(rep('b', 1)), props = fp_text(vertical.align =  "superscript"))))
# Add a footer with references for superscripts
ft <- add_footer_lines(ft, values = c(
  "a: median(± interquartile range)",
  "b: mean(± standard deviation)"
))

# Merge the cells in the "season" column where necessary to nest the rows by sex
ft <- merge_v(ft, j = "season")

# Set overarching headers
ft = add_header_row(x=ft, values = c('', 'trips', 'utilization distributions', 'in WEA', 'haulout'), 
                    colwidths = c(3, 4, 2, 2, 2), top = T)



# Manually add vertical borders to partition the header sections
ft <- theme_booktabs(ft)
# Manually add vertical borders to partition the header sections
ft = vline(x = ft, i = c(1,2), part = 'header', j = c(3, 7, 9, 11), border = fp_border(width = 2))
ft = hline(x = ft, i = c(2,4,6,8), border = fp_border(width = 1))



# Center the header text
ft <- align(ft, part = "header", i = c(1,2),align = "center")
ft <- align(ft, part = 'body',align = "center")
ft <- bold(ft, part = "header", i = c(1, 2), bold = TRUE)
ft <- bold(ft, part = "body", i = c(9,10), j = c(1),bold = TRUE)


ft = autofit(ft)



ft

save_as_docx(ft, path = "./plots/manuscript/Table2.docx", pr_section = sect_properties)


