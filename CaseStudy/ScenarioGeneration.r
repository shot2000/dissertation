library (JGR)
JGR ()
update.packages()

rm(list=ls(all=TRUE)) 
setwd("H:/01-Dissertation/R_Spatial_Analysis_2015/")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c('DeducerText','DeducerSpatial','DeducerPlugInScaling','DeducerExtras',"Deducer","JGR","ggplot2", "plyr", "reshape2", "RColorBrewer", "scales", "grid", "rJava", "RSQLite","sqldf","XLConnect","foreign")
ipak(packages)
### need to delete all extra "-"s in the file at column 26 before this step
HHhotspot <- read.dbf('H:/01-Dissertation/R_Spatial_Analysis_2015/HHhotspot.dbf')
head (HHhotspot)
##creating the increased density scenario

minGI<- min(HHhotspot$GiZScore)

summary (HHhotspot$GiZScore)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-7.0560 -1.2320  0.5782 -0.2811  1.6680  4.2700 

HHhotspot<- mutate (HHhotspot,GiABS1= GiZScore - minGI)
sumGI1<- sum(HHhotspot$GiABS1)
sumHH<- sum(HHhotspot$HH2010)
HHhotspot<- mutate (HHhotspot,GiFra=GiABS1/sumGI1)
HHhotspot<- mutate (HHhotspot,newHH_INCDensity=GiFra*sumHH)
#getting density
HHhotspot<- mutate (HHhotspot,HHdensity=HH2010/ACRES)
HHhotspot<- mutate (HHhotspot,HHdensityINC_den=newHH_INCDensity/ACRES)

ave (HHhotspot$HHdensity)
2.368521
ave (HHhotspot$HHdensityINC_den)
2.909786
###  0.2285244673785877 increase of HH density 



minGI<- min(HHhotspot$GiZScore)

summary (HHhotspot$GiZScore)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-7.0560 -1.2320  0.5782 -0.2811  1.6680  4.2700 

##creating the decreased density scenario
HHhotspot<- mutate (HHhotspot,GiABS2= GiZScore - minGI-1.232)
sumGI2<- sum(HHhotspot$GiABS2)
sumHH<- sum(HHhotspot$HH2010)
HHhotspot<- mutate (HHhotspot,GiFra=GiABS2/sumGI2)
HHhotspot<- mutate (HHhotspot,newHH_DECDensity=GiFra*sumHH)
#getting density
HHhotspot<- mutate (HHhotspot,HHdensityDEC_den=newHH_DECDensity/ACRES)

ave (HHhotspot$HHdensity)
2.368521
ave (HHhotspot$HHdensityDEC_den)
1.865315
###  Decrease of HH density 
###  0.2124559

##Export the data 
write.dbf(HHhotspot,'H:/01-Dissertation/R_Spatial_Analysis_2015/CaseStudy/HHScenarioData.dbf')



# 1. Denser urban development (shorter trip length, higher densities, more accessible transit)
# 2. Sparse suburban style development (longer trip length, lower densities, less accessible transit)
