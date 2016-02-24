###########

http://rstudio-pubs-static.s3.amazonaws.com/5032_4139cd3628e043e79ab2713f541a7463.html

##############
Reduce Colinearity
##  A VIF value >= 10 indicates high collinearity and inflated standard errors
https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/

## Clear the workspace
  
  rm(list=ls())

library(car)          ## variable colinerity check 
library(maps)         ## Projections
library(maptools)     ## Data management
library(sp)           ## Data management
library(spdep)        ## Spatial autocorrelation
library(gstat)        ## Geostatistics
library(splancs)      ## Kernel Density
library(spatstat)     ## Geostatistics
library(pgirmess)     ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(spgwr)        ## GWR
library(RANN)  
library(rgdal)  
library (foreign)
###3load the data
setwd("E:/01-Dissertation/R_Spatial_Analysis_2015")
## load the shapefile
taz<-readOGR(".","TAZ_CarbEM4Regression0514") 

taz <- read.table ("TAZ_CARBEM.csv", header=TRUE, sep=",")
#save(taz,file="Datasets.RData")
load("Datasets.RData")

##Change data type to double
taz[, 6]<-as.double(taz[, 6])
taz[,26]<-as.double(taz[,26])
taz[,23]<-as.double(taz[,23])
taz[,37]<-as.double(taz[,37])


data <- taz
names(data)

# Change Column names
names(data )[6]<-paste("AT")
names(data )[8]<-paste("TOTAL_HH")
names(data )[9]<-paste("AVGWK") 
names(data )[10]<-paste("AVGPER")
names(data )[11]<-paste("AVGAUTO")
names(data )[19]<-paste("EMP_L") 
names(data )[20]<-paste("EMP_M") 
names(data )[21]<-paste("EMP_H") 
names(data )[23]<-paste("POP") 
# read the spatial weight file

tazgal <- read.gal("TAZ_CarbEM4Regression0514.gal",override.id=TRUE)
attributes (tazgal)

# Creating a spatial weights list from a gal file
tazw <- nb2listw(tazgal)

###build the dataframe

n <- c("-50%", "-25%", "Baseline", "+25%", "+50%")
a <- c( "POP" , "TOTAL_HH" , "TOTAL_EMPL" , "POP_DENSIT" , "EMP_DENSIT" , "TOTAL_AUTO" , "EMP_M" , "AVGWK" , "AVGAUTO" , "Avg_CarbEM" , "Avg_TRIPSP")
df<- (n,a)
Master <- as.data.frame (matrix (nrow=11, ncol=5))
names(Master)[1:5]<-paste(n)
rownames(Master)[1:11]<-paste(a)

########
## Spatial Durbin Model
########

ev <- eigenw(similar.listw(tazw))
mod.SDEm <- errorsarlm(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP,  data,  tazw,method="eigen", control=list(pre_eig=ev),tol.solve=1.384e-14)

summary(mod.SDEm, correlation=TRUE)

####
## sensitivity test will scape ACRES and AT since they are not numerical
## each variable was given four scenarios for testing -50%, -25%, 25% and 50%
####################################################################################################

dataT_POP<- data #copy the data frame (so we don't mess up the original)

# Change the pop
dataT_POP_Minus50<- data
dataT_POP_Minus50$POP<- data$POP*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_Minus50, listw=tazw))   
Master[1,1] <- sum (new.pred$fit)   

dataT_POP_Minus25<- data              
dataT_POP_Minus25$POP<- data$POP*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_Minus25, listw=tazw))   
Master[1,2] <- sum (new.pred$fit) 


dataT_POP_Plus25<- data
dataT_POP_Plus25$POP<- data$POP*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_Plus25, listw=tazw))   
Master[1,4] <- sum (new.pred$fit)  

dataT_POP_Plus50<- data
dataT_POP_Plus50$POP<- data$POP*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_Plus50, listw=tazw))   
Master[1,5] <- sum (new.pred$fit) 

# Change the HH
dataT_HH_Minus50<- data
dataT_HH_Minus50$TOTAL_HH<- data$TOTAL_HH*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_HH_Minus50, listw=tazw))   
Master[2,1] <- sum (new.pred$fit)   

dataT_HH_Minus25<- data              
dataT_HH_Minus25$TOTAL_HH<- data$TOTAL_HH*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_HH_Minus25, listw=tazw))   
Master[2,2] <- sum (new.pred$fit) 


dataT_HH_Plus25<- data
dataT_HH_Plus25$TOTAL_HH<- data$TOTAL_HH*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_HH_Plus25, listw=tazw))   
Master[2,4] <- sum (new.pred$fit)  

dataT_HH_Plus50<- data
dataT_HH_Plus50$TOTAL_HH<- data$TOTAL_HH*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_HH_Plus50, listw=tazw))   
Master[2,5] <- sum (new.pred$fit) 

# Change the EMP
dataT_EMP_Minus50<- data
dataT_EMP_Minus50$TOTAL_EMPL<- data$TOTAL_EMPL*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_Minus50, listw=tazw))   
Master[3,1] <- sum (new.pred$fit)   

dataT_EMP_Minus25<- data              
dataT_EMP_Minus25$TOTAL_EMPL<- data$TOTAL_EMPL*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_Minus25, listw=tazw))   
Master[3,2] <- sum (new.pred$fit) 


dataT_EMP_Plus25<- data
dataT_EMP_Plus25$TOTAL_EMPL<- data$TOTAL_EMPL*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_Plus25, listw=tazw))   
Master[3,4] <- sum (new.pred$fit)  

dataT_EMP_Plus50<- data
dataT_EMP_Plus50$TOTAL_EMPL<- data$TOTAL_EMPL*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_Plus50, listw=tazw))   
Master[3,5] <- sum (new.pred$fit) 

# Change the POP Density
dataT_POP_DEN_Minus50<- data
dataT_POP_DEN_Minus50$POP_DENSIT<- data$POP_DENSIT*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_DEN_Minus50, listw=tazw))   
Master[4,1] <- sum (new.pred$fit)   

dataT_POP_DEN_Minus25<- data              
dataT_POP_DEN_Minus25$POP_DENSIT<- data$POP_DENSIT*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_DEN_Minus25, listw=tazw))   
Master[4,2] <- sum (new.pred$fit) 


dataT_POP_DEN_Plus25<- data
dataT_POP_DEN_Plus25$POP_DENSIT<- data$POP_DENSIT*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_DEN_Plus25, listw=tazw))   
Master[4,4] <- sum (new.pred$fit)  

dataT_POP_DEN_Plus50<- data
dataT_POP_DEN_Plus50$POP_DENSIT<- data$POP_DENSIT*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_POP_DEN_Plus50, listw=tazw))   
Master[4,5] <- sum (new.pred$fit) 

# Change the EMP Density
dataT_EMP_DEN_Minus50<- data
dataT_EMP_DEN_Minus50$EMP_DENSIT<- data$EMP_DENSIT*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_DEN_Minus50, listw=tazw))   
Master[5,1] <- sum (new.pred$fit)   

dataT_EMP_DEN_Minus25<- data              
dataT_EMP_DEN_Minus25$EMP_DENSIT<- data$EMP_DENSIT*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_DEN_Minus25, listw=tazw))   
Master[5,2] <- sum (new.pred$fit) 

dataT_EMP_DEN_Plus25<- data
dataT_EMP_DEN_Plus25$EMP_DENSIT<- data$EMP_DENSIT*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_DEN_Plus25, listw=tazw))   
Master[5,4] <- sum (new.pred$fit)  

dataT_EMP_DEN_Plus50<- data
dataT_EMP_DEN_Plus50$EMP_DENSIT<- data$EMP_DENSIT*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_DEN_Plus50, listw=tazw))   
Master[5,5] <- sum (new.pred$fit) 

# Change the Total Auto
dataT_TAUTO_DEN_Minus50<- data
dataT_TAUTO_DEN_Minus50$TOTAL_AUTO<- data$TOTAL_AUTO*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_TAUTO_DEN_Minus50, listw=tazw))   
Master[6,1] <- sum (new.pred$fit)   

dataT_TAUTO_DEN_Minus25<- data              
dataT_TAUTO_DEN_Minus25$TOTAL_AUTO<- data$TOTAL_AUTO*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_TAUTO_DEN_Minus25, listw=tazw))   
Master[6,2] <- sum (new.pred$fit) 

dataT_TAUTO_DEN_Plus25<- data
dataT_TAUTO_DEN_Plus25$TOTAL_AUTO<- data$TOTAL_AUTO*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_TAUTO_DEN_Plus25, listw=tazw))   
Master[6,4] <- sum (new.pred$fit)  

dataT_TAUTO_DEN_Plus50<- data
dataT_TAUTO_DEN_Plus50$TOTAL_AUTO<- data$TOTAL_AUTO*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_TAUTO_DEN_Plus50, listw=tazw))   
Master[6,5] <- sum (new.pred$fit) 

# Change the EMP M
dataT_EMP_M_Minus50<- data
dataT_EMP_M_Minus50$EMP_M<- data$EMP_M*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_M_Minus50, listw=tazw))   
Master[7,1] <- sum (new.pred$fit)   

dataT_EMP_M_Minus25<- data              
dataT_EMP_M_Minus25$EMP_M<- data$EMP_M*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_M_Minus25, listw=tazw))   
Master[7,2] <- sum (new.pred$fit) 

dataT_EMP_M_Plus25<- data
dataT_EMP_M_Plus25$EMP_M<- data$EMP_M*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_M_Plus25, listw=tazw))   
Master[7,4] <- sum (new.pred$fit)  

dataT_EMP_M_Plus50<- data
dataT_EMP_M_Plus50$EMP_M<- data$EMP_M*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP_M_Plus50, listw=tazw))   
Master[7,5] <- sum (new.pred$fit) 

# Change the Average Worker
dataT_AVGWK_Minus50<- data
dataT_AVGWK_Minus50$AVGWK<- data$AVGWK*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGWK_Minus50, listw=tazw))   
Master[8,1] <- sum (new.pred$fit)   

dataT_AVGWK_Minus25<- data              
dataT_AVGWK_Minus25$AVGWK<- data$AVGWK*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGWK_Minus25, listw=tazw))   
Master[8,2] <- sum (new.pred$fit) 

dataT_AVGWK_Plus25<- data
dataT_AVGWK_Plus25$AVGWK<- data$AVGWK*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGWK_Plus25, listw=tazw))   
Master[8,4] <- sum (new.pred$fit)  

dataT_AVGWK_Plus50<- data
dataT_AVGWK_Plus50$AVGWK<- data$AVGWK*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGWK_Plus50, listw=tazw))   
Master[8,5] <- sum (new.pred$fit) 

# Change the Average Auto
dataT_AVGAUTO_Minus50<- data
dataT_AVGAUTO_Minus50$AVGAUTO<- data$AVGAUTO*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGAUTO_Minus50, listw=tazw))   
Master[9,1] <- sum (new.pred$fit)   

dataT_AVGAUTO_Minus25<- data              
dataT_AVGAUTO_Minus25$AVGAUTO<- data$AVGAUTO*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGAUTO_Minus25, listw=tazw))   
Master[9,2] <- sum (new.pred$fit) 

dataT_AVGAUTO_Plus25<- data
dataT_AVGAUTO_Plus25$AVGAUTO<- data$AVGAUTO*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGAUTO_Plus25, listw=tazw))   
Master[9,4] <- sum (new.pred$fit)  

dataT_AVGAUTO_Plus50<- data
dataT_AVGAUTO_Plus50$AVGAUTO<- data$AVGAUTO*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_AVGAUTO_Plus50, listw=tazw))   
Master[9,5] <- sum (new.pred$fit) 

# Change the Average Carbon
dataT_Avg_CarbEM_Minus50<- data
dataT_Avg_CarbEM_Minus50$Avg_CarbEM<- data$Avg_CarbEM*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_CarbEM_Minus50, listw=tazw))   
Master[10,1] <- sum (new.pred$fit)   

dataT_Avg_CarbEM_Minus25<- data              
dataT_Avg_CarbEM_Minus25$Avg_CarbEM<- data$Avg_CarbEM*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_CarbEM_Minus25, listw=tazw))   
Master[10,2] <- sum (new.pred$fit) 

dataT_Avg_CarbEM_Plus25<- data
dataT_Avg_CarbEM_Plus25$Avg_CarbEM<- data$Avg_CarbEM*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_CarbEM_Plus25, listw=tazw))   
Master[10,4] <- sum (new.pred$fit)  

dataT_Avg_CarbEM_Plus50<- data
dataT_Avg_CarbEM_Plus50$Avg_CarbEM<- data$Avg_CarbEM*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_CarbEM_Plus50, listw=tazw))   
Master[10,5] <- sum (new.pred$fit) 

# Change the Average SPD
dataT_Avg_TRIPSP_Minus50<- data
dataT_Avg_TRIPSP_Minus50$Avg_TRIPSP<- data$Avg_TRIPSP*(1-0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_TRIPSP_Minus50, listw=tazw))   
Master[11,1] <- sum (new.pred$fit)   

dataT_Avg_TRIPSP_Minus25<- data              
dataT_Avg_TRIPSP_Minus25$Avg_TRIPSP<- data$Avg_TRIPSP*(1-0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_TRIPSP_Minus25, listw=tazw))   
Master[11,2] <- sum (new.pred$fit) 

dataT_Avg_TRIPSP_Plus25<- data
dataT_Avg_TRIPSP_Plus25$Avg_TRIPSP<- data$Avg_TRIPSP*(1+0.25)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_TRIPSP_Plus25, listw=tazw))   
Master[11,4] <- sum (new.pred$fit)  

dataT_Avg_TRIPSP_Plus50<- data
dataT_Avg_TRIPSP_Plus50$Avg_TRIPSP<- data$Avg_TRIPSP*(1+0.5)
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_Avg_TRIPSP_Plus50, listw=tazw))   
Master[11,5] <- sum (new.pred$fit) 

Master[,3]<- sum (data$CarEMTAZ)

library(xlsx) #load the package
write.xlsx(x=Master, file = "data.xlsx",sheetName = "TestSheet", row.names = TRUE,col.names = TRUE)