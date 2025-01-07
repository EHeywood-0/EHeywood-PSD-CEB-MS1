library(readr)
library(dplyr)
library(tidyr)
setwd("~/seal_telemetry")

path = paste(getwd(), '/data/L3/dive/imputeddivepos/', sep = '')

fnames = list.files(path = path,
                    pattern = '*.csv$', all.files = TRUE,
                    recursive = TRUE, full.names = TRUE)

data = purrr::map_dfr(.x = fnames, .f = function(.x){

  tmp = read_csv(.x) %>% group_by(diveID) %>%
    mutate(Nimputes = n()) %>%
    rename(bathy_m = bathydepth)
  return(tmp)

})

data = data %>%
  mutate(totaldivespossible = n_dives_lost + n_dives_imputed) %>%
  mutate(percloss = n_dives_lost / totaldivespossible)


datmerge = data[data$bathy_m < 0, ]
View(head(datmerge))
write_csv(x = datmerge, file = './data/L3/dive/Hg_2019-2023_DivePosEstimates_Imputations_BestEst.csv')
