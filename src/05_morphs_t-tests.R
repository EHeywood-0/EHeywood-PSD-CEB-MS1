# MORPHOMETRICS SUMMARIES FOR INTRO MANUSCRIPT

library(readr)
library(dplyr)
library(tidyr)
library(effsize)

sealdat<-read_csv("./data/meta/Hg2019-2023_WC_Tag_summaryFiles+MetaData.csv") 
#sealdat = sealdat[!sealdat$notes %in% c("tag malfunction", "never signaled"),]
sealdat$year = format(sealdat$deploydate, '%Y')

# TOTALS
mean(sealdat$masskg)
range(sealdat$masskg)

tagsum<-sealdat %>%
  group_by(year, sex)%>%
  summarize(
    count = n(),
    meanwt = round(mean(masskg), digits = 1),
    sdwt = round(sd(masskg), digits = 1),
    meanlen = round(mean(lengthcm), digits = 1),
    sdlen = round(sd(lengthcm), digits = 1),
    meangirth = round(mean(girthcm), digits = 1),
    sdgir = round(sd(girthcm), digits = 1),
    meantagdur = round(mean(totdatadays), digits = 1),
    sdtagdur = round(sd(totdatadays), digits = 1)
  )

morphsumyearsex<-as.data.frame(tagsum)

print(morphsumyearsex)

totalsline<-sealdat %>%
  group_by(sex)%>%
  summarize(
    allmass=mean(masskg),
    allmasssd=sd(masskg),
    alllen=mean(lengthcm),
    alllensd=sd(lengthcm),
    allgir=mean(girthcm),
    allgirsd=sd(girthcm),
    alldays=mean(totdatadays, na.rm = T),
    alldyssd=sd(totdatadays, na.rm = T)
  )

print(totalsline)

## TTEST MASS
mmass = sealdat$masskg[sealdat$sex == 'M']
fmass = sealdat$masskg[sealdat$sex == 'F']
# Conduct the independent t-test
res <- t.test(x = mmass, y = fmass, paired = FALSE, var.equal = TRUE)
print(res)

cohen.d(mmass, fmass)

## TTEST LENGHT
m = sealdat$lengthcm[sealdat$sex == 'M']
f = sealdat$lengthcm[sealdat$sex == 'F']
# Conduct the independent t-test
res <- t.test(x = m, y = f, paired = FALSE, var.equal = TRUE)
print(res)

cohen.d(m, f)


## TTEST GIRTH
m = sealdat$girthcm[sealdat$sex == 'M']
f = sealdat$girthcm[sealdat$sex == 'F']
# Conduct the independent t-test
res <- t.test(x = m, y = f, paired = FALSE, var.equal = TRUE)
print(res)

cohen.d(m, f)


