## for sensitivity
dat2$SPEND[dat2$COUNTRY=="PL"] <- dat2$SPEND[dat2$COUNTRY=="PL"]*2

## Clear the workspace

rm(list=ls())


## Install packages

install.packages("maps")
install.packages("maptools")
install.packages("sp")
install.packages("spdep")
install.packages("gstat")
install.packages("splancs")
install.packages("spatstat")
install.packages("lattice")
install.packages("pgirmess")
install.packages("RColorBrewer")
install.packages("classInt")
install.packages("spgwr")
install.packages("rgdal")

## Load spatial packages

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

setwd("G:/01-Dissertation/R_Spatial_Analysis_2015")
## load the shapefile
taz<-readOGR(".","TAZ_CarbEM4Regression0514")  
#save(taz,file="Datasets.RData")
load("Datasets.RData")

data <- taz
names(data)

# CHange Column names
names(data )[6]<-paste("AT")
names(data )[8]<-paste("TOTAL_HH")
names(data )[9]<-paste("AVGWK") 
names(data )[10]<-paste("AVGPER")
names(data )[11]<-paste("AVGAUTO")
names(data )[19]<-paste("EMP_L") 
names(data )[20]<-paste("EMP_M") 
names(data )[21]<-paste("EMP_H") 
names(data )[23]<-paste("POP") 

###check for colinearity


########
## Linear Model
########

mod.lm2 <- lm(CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data=data)

mod.lm <- lm(CarEMTAZ ~ ACRES + AT  + TOTAL_HH + POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP, data=data)


##VIF > 10 is highly colinear

vif (mod.lm2)
ACRES         AT        POP   TOTAL_HH TOTAL_EMPL POP_DENSIT EMP_DENSIT 
1.574900   2.691958  24.899501  22.078868   1.476185   1.854916   1.622575 
TOTAL_AUTO      EMP_M      AVGWK    AVGAUTO Avg_CarbEM Avg_TRIPSP 
14.317840   1.112205   2.574727   4.949845   3.496820   4.039569 

vif (mod.lm)

####droped colinear variables


summary(mod.lm)

##########
Call:
  lm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
       POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
       Avg_TRIPSP, data = data)

Residuals:
  Min       1Q   Median       3Q      Max 
-1.60216 -0.28434 -0.03248  0.21790  2.77159 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept) -5.388e-01  1.034e-01  -5.213 2.46e-07 ***
  ACRES        6.821e-05  4.711e-05   1.448 0.148083    
AT           1.986e-01  4.875e-02   4.074 5.17e-05 ***
  TOTAL_EMPL  -5.288e-05  2.609e-05  -2.027 0.043085 *  
  TOTAL_HH     2.178e-03  6.533e-05  33.347  < 2e-16 ***
  POP_DENSIT  -2.787e-02  5.029e-03  -5.543 4.27e-08 ***
  EMP_DENSIT   9.039e-04  2.929e-04   3.086 0.002110 ** 
  EMP_M        2.510e-01  7.441e-02   3.374 0.000783 ***
  AVGWK        2.471e-01  8.789e-02   2.811 0.005076 ** 
  AVGAUTO      7.397e-03  6.862e-02   0.108 0.914193    
Avg_CarbEM  -6.423e+01  2.211e+01  -2.905 0.003791 ** 
  Avg_TRIPSP   1.979e-02  3.749e-03   5.279 1.75e-07 ***
  ---
  Signif. codes:  0 ?**?0.001 ?*?0.01 ??0.05 ??0.1 ??1

Residual standard error: 0.5269 on 681 degrees of freedom
Multiple R-squared:  0.7776,  Adjusted R-squared:  0.774 
F-statistic: 216.4 on 11 and 681 DF,  p-value: < 2.2e-16

###########
AIC(mod.lm)
[1] 1092.482
##getting BIC
bic=AIC(mod.lm,k = log(length(data)))
bic
[1] 1151.515
##Extract Log-Likelihood from an lm Object
logLik(mod.lm)
'log Lik.' -533.2409 (df=13)

## Plot residuals

res <- mod.lm$residuals
png(filename="plots/Residuals from OLS Model.png", width = 800, height = 600)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)

classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-3.5,-1.5,-0.5,0.5,1.5,3.5), rtimes = 1)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, main="Residuals from OLS Model", pretty=T, border="grey")
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from OLS Model",ncol=5)

dev.off()

############ Residual Autocorrelation

# read the spatial weight file
tazgal <- read.gal("TAZ_CarbEM4Regression0514.gal",override.id=TRUE)
attributes (tazgal)
# Creating a spatial weights list from a gal file
tazw <- nb2listw(tazgal)
moran.test(res, listw=tazw, zero.policy=T)

# the results
#Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = 1.2099, p-value = 0.1132
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
0.026397503      -0.001445087       0.000529529 
	
########
## SAR Model (WARNING: This takes a while to run)
########

mod.sar <- lagsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                      POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                      Avg_TRIPSP, data = data, listw=tazw , zero.policy=T, tol.solve=2e-14)
summary(mod.sar)

##results

Call:lagsarlm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                Avg_TRIPSP, data = data, listw = tazw, zero.policy = T, tol.solve = 2e-14)

Residuals:
  Min        1Q    Median        3Q       Max 
-1.604924 -0.286649 -0.025675  0.210998  2.804901 

Type: lag 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -4.9812e-01  1.0328e-01 -4.8228 1.416e-06
ACRES        4.2356e-05  4.6832e-05  0.9044 0.3657685
AT           2.0750e-01  4.8133e-02  4.3109 1.626e-05
TOTAL_EMPL  -5.8549e-05  2.5799e-05 -2.2695 0.0232398
TOTAL_HH     2.2117e-03  6.6080e-05 33.4705 < 2.2e-16
POP_DENSIT  -2.8329e-02  4.9641e-03 -5.7067 1.152e-08
EMP_DENSIT   8.3464e-04  2.9130e-04  2.8652 0.0041670
EMP_M        2.4627e-01  7.3405e-02  3.3550 0.0007937
AVGWK        2.6212e-01  8.7068e-02  3.0106 0.0026077
AVGAUTO      4.0065e-02  6.8673e-02  0.5834 0.5596147
Avg_CarbEM  -6.4586e+01  2.1808e+01 -2.9616 0.0030604
Avg_TRIPSP   2.0154e-02  3.6985e-03  5.4492 5.059e-08

Rho: -0.085439, LR test value: 6.0644, p-value: 0.013794
Asymptotic standard error: 0.0349
z-value: -2.4481, p-value: 0.014361
Wald statistic: 5.9932, p-value: 0.014361

Log likelihood: -530.2087 for lag model
ML residual variance (sigma squared): 0.2701, (sigma: 0.51971)
Number of observations: 693 
Number of parameters estimated: 14 
AIC: 1088.4, (AIC for lm: 1092.5)
LM test for residual autocorrelation
test value: 8.5313, p-value: 0.0034908

###########residual

res <- mod.sar$residuals
png(filename="plots/Residuals from Spatial Auto-Regressive Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Spatial Auto-Regressive Model",ncol=5)

dev.off()

## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)
#Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = 2.5224, p-value = 0.005829
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
0.0565954544     -0.0014450867      0.0005294812 

###########
AIC(mod.sar)
[1] 1088.417

bic=AIC(mod.sar,k = log(length(data)))
bic
[1] 1151.992

logLik(mod.sar)
'log Lik.' -530.2087 (df=14)


########
## SEM Model (WARNING: This takes a while to run)
########


mod.sem <- errorsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                        POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                        Avg_TRIPSP, data = data, listw = tazw, zero.policy=T, tol.solve=1e-15)


summary(mod.sem)

Call:errorsarlm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                  POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                  Avg_TRIPSP, data = data, listw = tazw, zero.policy = T, tol.solve = 1e-15)

Residuals:
  Min        1Q    Median        3Q       Max 
-1.623820 -0.282586 -0.024815  0.215013  2.771203 

Type: error 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -5.4202e-01  1.0554e-01 -5.1358 2.810e-07
ACRES        8.9334e-05  4.8223e-05  1.8525 0.0639518
AT           1.9249e-01  5.0028e-02  3.8476 0.0001193
TOTAL_EMPL  -5.1489e-05  2.6063e-05 -1.9756 0.0482018
TOTAL_HH     2.1897e-03  6.5412e-05 33.4757 < 2.2e-16
POP_DENSIT  -2.8027e-02  5.0869e-03 -5.5097 3.594e-08
EMP_DENSIT   8.6429e-04  2.9763e-04  2.9039 0.0036855
EMP_M        2.5134e-01  7.4339e-02  3.3810 0.0007222
AVGWK        2.5527e-01  8.8026e-02  2.8999 0.0037323
AVGAUTO      1.4246e-03  6.9162e-02  0.0206 0.9835658
Avg_CarbEM  -6.4301e+01  2.1871e+01 -2.9400 0.0032820
Avg_TRIPSP   1.9545e-02  3.7187e-03  5.2559 1.473e-07

Lambda: 0.075007, LR test value: 1.379, p-value: 0.24027
Asymptotic standard error: 0.060765
z-value: 1.2344, p-value: 0.21706
Wald statistic: 1.5237, p-value: 0.21706

Log likelihood: -532.5514 for error model
ML residual variance (sigma squared): 0.272, (sigma: 0.52154)
Number of observations: 693 
Number of parameters estimated: 14 
AIC: 1093.1, (AIC for lm: 1092.5)

#################################
res <- mod.sem$residuals

classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
cols <- findColours(classes_fx,pal)
png(filename="plots/Residuals from Spatial Error Model.png", width = 800, height = 600)
par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Spatial Error Model",ncol=5)

dev.off()

## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)
# Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = 0.0352, p-value = 0.486
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
-0.0006347203     -0.0014450867      0.0005295174 

###########
AIC(mod.sem)
[1] 1093.103

bic=AIC(mod.sem,k = log(length(data)))
bic
[1] 1156.677

logLik(mod.sem)
'log Lik.' -532.5514 (df=14)

########
## SDM Model 
########


mod.sdm <- lagsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                      POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                      Avg_TRIPSP, data = data, listw = tazw, zero.policy=T, type="mixed", tol.solve=1.47932e-15)
summary(mod.sdm)

Call:lagsarlm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                Avg_TRIPSP, data = data, listw = tazw, type = "mixed", zero.policy = T, 
              tol.solve = 1.47932e-15)

Residuals:
  Min        1Q    Median        3Q       Max 
-1.466746 -0.284153 -0.028162  0.228628  3.025738 

Type: mixed 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept)    -4.0864e-01  1.7091e-01 -2.3909 0.0168049
ACRES           2.1560e-04  5.8955e-05  3.6571 0.0002551
AT              3.4062e-01  8.6818e-02  3.9234 8.730e-05
TOTAL_EMPL     -6.0398e-05  2.6613e-05 -2.2695 0.0232384
TOTAL_HH        2.2104e-03  6.8218e-05 32.4015 < 2.2e-16
POP_DENSIT     -3.0544e-02  5.6927e-03 -5.3655 8.074e-08
EMP_DENSIT      5.8882e-04  3.5813e-04  1.6441 0.1001452
EMP_M           2.0042e-01  7.4702e-02  2.6829 0.0072981
AVGWK           3.1376e-01  9.0577e-02  3.4640 0.0005322
AVGAUTO         7.2833e-02  7.8350e-02  0.9296 0.3525843
Avg_CarbEM     -7.1337e+01  2.1386e+01 -3.3357 0.0008509
Avg_TRIPSP      1.9671e-02  3.6391e-03  5.4055 6.463e-08
lag.ACRES      -3.8645e-04  7.9765e-05 -4.8449 1.267e-06
lag.AT         -7.7019e-02  1.1388e-01 -0.6763 0.4988327
lag.TOTAL_EMPL -4.7395e-05  5.0153e-05 -0.9450 0.3446636
lag.TOTAL_HH   -4.2401e-04  1.8164e-04 -2.3344 0.0195765
lag.POP_DENSIT  1.2737e-02  1.0032e-02  1.2697 0.2042069
lag.EMP_DENSIT  7.1609e-04  5.6102e-04  1.2764 0.2018168
lag.EMP_M      -3.4103e-02  1.5380e-01 -0.2217 0.8245161
lag.AVGWK      -2.4517e-01  1.7841e-01 -1.3742 0.1693732
lag.AVGAUTO     6.2541e-02  1.3683e-01  0.4571 0.6476094
lag.Avg_CarbEM -3.5291e+01  5.0533e+01 -0.6984 0.4849405
lag.Avg_TRIPSP  1.2023e-02  8.3791e-03  1.4349 0.1513209

Rho: 0.021338, LR test value: 0.12855, p-value: 0.71994
Asymptotic standard error: 0.061175
z-value: 0.3488, p-value: 0.72724
Wald statistic: 0.12166, p-value: 0.72724

Log likelihood: -507.6612 for mixed model
ML residual variance (sigma squared): 0.25338, (sigma: 0.50337)
Number of observations: 693 
Number of parameters estimated: 25 
AIC: 1065.3, (AIC for lm: 1063.5)
LM test for residual autocorrelation
test value: 5.939, p-value: 0.014809

#############################################

res <- mod.sdm$residuals
png(filename="plots/Residuals from Spatial Durbin Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-4.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Spatial Durbin Model",ncol=5)

dev.off()
## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)
# 	Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = -0.3182, p-value = 0.6248
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
-0.0087638893     -0.0014450867      0.0005291067 

###########
AIC(mod.sdm)
[1] 1065.322

bic=AIC(mod.sdm,k = log(length(data)))
bic
[1] 1178.848


logLik(mod.sdm)
'log Lik.' -507.6612 (df=25)


########
## Kelejian-Prucha Model 
########

##mod.kpm <- GMerrorsar(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP, data,  tazw,verbose=TRUE)
mod.kpm <- gstsls(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                    POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                    Avg_TRIPSP, data,  tazw,verbose=TRUE)


summary(mod.kpm, correlation = FALSE, Hausman=FALSE)

##results

Call:gstsls(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
              POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
              Avg_TRIPSP, data = data, listw = tazw, verbose = TRUE)

Residuals:
  Min        1Q    Median        3Q       Max 
-1.640256 -0.280245 -0.038826  0.213125  2.802804 

Type: GM SARAR estimator
Coefficients: (GM standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
Rho_Wy      -1.1896e-01  4.3162e-02 -2.7561 0.0058498
(Intercept) -4.7917e-01  1.1300e-01 -4.2405 2.231e-05
ACRES        8.5420e-05  5.1917e-05  1.6453 0.0999054
AT           1.9892e-01  5.2817e-02  3.7662 0.0001658
TOTAL_EMPL  -5.5997e-05  2.6352e-05 -2.1250 0.0335883
TOTAL_HH     2.2290e-03  6.6770e-05 33.3830 < 2.2e-16
POP_DENSIT  -2.8531e-02  5.2183e-03 -5.4674 4.567e-08
EMP_DENSIT   7.1598e-04  3.0943e-04  2.3139 0.0206738
EMP_M        2.4240e-01  7.4995e-02  3.2322 0.0012283
AVGWK        2.8616e-01  8.9290e-02  3.2048 0.0013514
AVGAUTO      3.7707e-02  7.2372e-02  0.5210 0.6023568
Avg_CarbEM  -6.4410e+01  2.1751e+01 -2.9612 0.0030641
Avg_TRIPSP   1.9536e-02  3.7140e-03  5.2601 1.440e-07

Lambda: 0.17902
Residual variance (sigma squared): 0.27024, (sigma: 0.51984)
GM argmin sigma squared: 0.26729
Number of observations: 693 
Number of parameters estimated: 15 


###########residual Plot

res <- mod.kpm$residuals
png(filename="plots/Residuals from Kelejian-Prucha Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Kelejian-Prucha Model",ncol=5)

dev.off()

## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)

#Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = 0.048, p-value = 0.4808
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
-0.0003400694     -0.0014450867      0.0005294290 



########
## Spatial Durbin Error Model 
########

ev <- eigenw(similar.listw(tazw))
mod.SDEm <- errorsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                         POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                         Avg_TRIPSP, data,  tazw,method="eigen", control=list(pre_eig=ev),tol.solve=1.384e-14)

summary(mod.SDEm, correlation=TRUE)

### Results


Call:errorsarlm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                  POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                  Avg_TRIPSP, data = data, listw = tazw, method = "eigen", 
                tol.solve = 1.384e-14, control = list(pre_eig = ev))

Residuals:
  Min        1Q    Median        3Q       Max 
-1.623820 -0.282586 -0.024815  0.215013  2.771203 

Type: error 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -5.4202e-01  1.0554e-01 -5.1358 2.810e-07
ACRES        8.9334e-05  4.8223e-05  1.8525 0.0639518
AT           1.9249e-01  5.0028e-02  3.8476 0.0001193
TOTAL_EMPL  -5.1489e-05  2.6063e-05 -1.9756 0.0482018
TOTAL_HH     2.1897e-03  6.5412e-05 33.4757 < 2.2e-16
POP_DENSIT  -2.8027e-02  5.0869e-03 -5.5097 3.594e-08
EMP_DENSIT   8.6429e-04  2.9763e-04  2.9039 0.0036855
EMP_M        2.5134e-01  7.4339e-02  3.3810 0.0007222
AVGWK        2.5527e-01  8.8026e-02  2.8999 0.0037323
AVGAUTO      1.4246e-03  6.9162e-02  0.0206 0.9835658
Avg_CarbEM  -6.4301e+01  2.1871e+01 -2.9400 0.0032820
Avg_TRIPSP   1.9545e-02  3.7187e-03  5.2559 1.473e-07

Lambda: 0.075007, LR test value: 1.379, p-value: 0.24027
Asymptotic standard error: 0.060765
z-value: 1.2344, p-value: 0.21706
Wald statistic: 1.5237, p-value: 0.21706

Log likelihood: -532.5514 for error model
ML residual variance (sigma squared): 0.272, (sigma: 0.52154)
Number of observations: 693 
Number of parameters estimated: 14 
AIC: 1093.1, (AIC for lm: 1092.5)

Correlation of coefficients 
sigma lambda (Intercept) ACRES AT    TOTAL_EMPL TOTAL_HH POP_DENSIT
lambda      -0.03                                                              
(Intercept)  0.00  0.00                                                        
ACRES        0.00  0.00   0.07                                                 
AT           0.00  0.00  -0.31       -0.29                                     
TOTAL_EMPL   0.00  0.00  -0.32        0.09 -0.14                               
TOTAL_HH     0.00  0.00   0.02       -0.14  0.00  0.04                         
POP_DENSIT   0.00  0.00  -0.36        0.22  0.15  0.06      -0.43              
EMP_DENSIT   0.00  0.00  -0.16       -0.11  0.30 -0.40       0.11     0.04     
EMP_M        0.00  0.00  -0.42        0.02  0.11 -0.03      -0.11    -0.02     
AVGWK        0.00  0.00  -0.29        0.05  0.05  0.13      -0.07    -0.08     
AVGAUTO      0.00  0.00  -0.17       -0.08 -0.44  0.06       0.00     0.17     
Avg_CarbEM   0.00  0.00  -0.02       -0.10  0.01  0.01       0.06     0.01     
Avg_TRIPSP   0.00  0.00   0.01        0.00  0.00  0.05      -0.26    -0.06     
EMP_DENSIT EMP_M AVGWK AVGAUTO Avg_CarbEM
lambda                                               
(Intercept)                                          
ACRES                                                
AT                                                   
TOTAL_EMPL                                           
TOTAL_HH                                             
POP_DENSIT                                           
EMP_DENSIT                                           
EMP_M       -0.10                                    
AVGWK        0.00       0.01                         
AVGAUTO     -0.01      -0.01 -0.62                   
Avg_CarbEM  -0.03      -0.01  0.01  0.01             
Avg_TRIPSP   0.02      -0.04 -0.05 -0.06   -0.79     


###########residual

res <- mod.SDEm$residuals
png(filename="plots/Residuals from Spatial Durbin Error Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Spatial Durbin Error Model ",ncol=5)

dev.off()

## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)

#Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = 0.0352, p-value = 0.486
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
-0.0006347203     -0.0014450867      0.0005295174 


########
AIC(mod.SDEm)
[1] 1093.103

bic=AIC(mod.SDEm,k = log(length(data)))
bic
[1] 1156.677

logLik(mod.SDEm)
'log Lik.' -532.5514 (df=14)

########

## Manski Model 
########

mod.mam <- sacsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                      POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                      Avg_TRIPSP,  data,  tazw, type="sacmixed")
summary(mod.mam, correlation=TRUE)

###results
Call:sacsarlm(formula = CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                Avg_TRIPSP, data = data, listw = tazw, type = "sacmixed")

Residuals:
  Min        1Q    Median        3Q       Max 
-1.426094 -0.275260 -0.023359  0.228974  2.951609 

Type: sacmixed 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept)    -2.7369e-01  1.5500e-01 -1.7658 0.0774339
ACRES           2.4911e-04  6.3526e-05  3.9214 8.804e-05
AT              3.5582e-01  8.9861e-02  3.9597 7.506e-05
TOTAL_EMPL     -5.3670e-05  2.7459e-05 -1.9546 0.0506348
TOTAL_HH        2.2306e-03  7.0608e-05 31.5918 < 2.2e-16
POP_DENSIT     -3.1451e-02  5.8862e-03 -5.3432 9.134e-08
EMP_DENSIT      4.3618e-04  3.7842e-04  1.1526 0.2490563
EMP_M           1.9674e-01  7.7183e-02  2.5490 0.0108043
AVGWK           3.0823e-01  9.3551e-02  3.2947 0.0009851
AVGAUTO         7.3902e-02  8.1098e-02  0.9113 0.3621525
Avg_CarbEM     -7.0975e+01  2.1448e+01 -3.3092 0.0009358
Avg_TRIPSP      1.9301e-02  3.6875e-03  5.2341 1.658e-07
lag.ACRES      -3.9709e-04  7.7920e-05 -5.0962 3.466e-07
lag.AT         -1.6150e-01  1.1692e-01 -1.3813 0.1671889
lag.TOTAL_EMPL -3.4811e-05  4.7441e-05 -0.7338 0.4630831
lag.TOTAL_HH   -1.0165e-03  3.5345e-04 -2.8760 0.0040275
lag.POP_DENSIT  2.0380e-02  9.8507e-03  2.0689 0.0385523
lag.EMP_DENSIT  5.4544e-04  5.5989e-04  0.9742 0.3299636
lag.EMP_M      -9.8952e-02  1.4670e-01 -0.6745 0.4999876
lag.AVGWK      -2.4264e-01  1.6522e-01 -1.4686 0.1419529
lag.AVGAUTO     1.7511e-02  1.2888e-01  0.1359 0.8919203
lag.Avg_CarbEM -1.1202e+01  4.7304e+01 -0.2368 0.8128071
lag.Avg_TRIPSP  4.5336e-03  8.5463e-03  0.5305 0.5957850

Rho: 0.31018
Asymptotic standard error: 0.1662
z-value: 1.8663, p-value: 0.062
Lambda: -0.34871
Asymptotic standard error: 0.20558
z-value: -1.6962, p-value: 0.089845

LR test value: 53.728, p-value: 6.7467e-07

Log likelihood: -506.3768 for sacmixed model
ML residual variance (sigma squared): 0.24285, (sigma: 0.4928)
Number of observations: 693 
Number of parameters estimated: 26 
AIC: 1064.8, (AIC for lm: 1092.5)

Correlation of coefficients 
sigma rho   lambda (Intercept) ACRES AT    TOTAL_EMPL
rho            -0.62                                                
lambda          0.63 -0.95                                          
(Intercept)    -0.32  0.51 -0.49                                    
ACRES          -0.13  0.21 -0.20   0.09                             
AT              0.03 -0.04  0.04  -0.01       -0.21                 
TOTAL_EMPL     -0.05  0.08 -0.07  -0.02        0.01 -0.01           
TOTAL_HH       -0.09  0.14 -0.13   0.05       -0.23  0.00  0.04     
POP_DENSIT      0.02 -0.03  0.03  -0.01        0.24 -0.01  0.09     
EMP_DENSIT      0.05 -0.08  0.08   0.00        0.00  0.01 -0.46     
EMP_M           0.00  0.00  0.00  -0.01       -0.02  0.08 -0.02     
AVGWK          -0.02  0.04 -0.04   0.06        0.04  0.03  0.12     
AVGAUTO         0.02 -0.04  0.04   0.03       -0.11 -0.14  0.04     
Avg_CarbEM     -0.03  0.05 -0.04   0.05       -0.05 -0.03 -0.01     
Avg_TRIPSP      0.05 -0.08  0.07  -0.06       -0.05 -0.02  0.04     
lag.ACRES      -0.04  0.07 -0.06  -0.01       -0.74  0.18  0.04     
lag.AT          0.20 -0.32  0.31  -0.24        0.10 -0.80 -0.03     
lag.TOTAL_EMPL -0.18  0.29 -0.27  -0.12        0.10 -0.01 -0.41     
lag.TOTAL_HH    0.59 -0.94  0.89  -0.50       -0.14  0.03 -0.06     
lag.POP_DENSIT -0.24  0.38 -0.36   0.00       -0.08  0.00 -0.05     
lag.EMP_DENSIT  0.20 -0.32  0.31  -0.28       -0.05 -0.02  0.27     
lag.EMP_M       0.17 -0.27  0.25  -0.48       -0.07 -0.03 -0.01     
lag.AVGWK       0.09 -0.14  0.13  -0.29       -0.06  0.02 -0.07     
lag.AVGAUTO     0.08 -0.13  0.12  -0.14        0.02  0.10 -0.04     
lag.Avg_CarbEM -0.18  0.29 -0.27   0.16        0.05 -0.01  0.08     
lag.Avg_TRIPSP  0.29 -0.47  0.45  -0.25       -0.03  0.01 -0.05     
TOTAL_HH POP_DENSIT EMP_DENSIT EMP_M AVGWK AVGAUTO
rho                                                              
lambda                                                           
(Intercept)                                                      
ACRES                                                            
AT                                                               
TOTAL_EMPL                                                       
TOTAL_HH                                                         
POP_DENSIT     -0.45                                             
EMP_DENSIT     -0.03     0.02                                    
EMP_M          -0.14     0.03      -0.02                         
AVGWK          -0.03    -0.07       0.05      -0.01              
AVGAUTO         0.08     0.07       0.00      -0.05 -0.44        
Avg_CarbEM      0.03     0.04       0.05       0.00  0.02  0.02  
Avg_TRIPSP     -0.21    -0.08      -0.03      -0.04 -0.04 -0.05  
lag.ACRES       0.22    -0.19      -0.03       0.02 -0.03  0.06  
lag.AT         -0.04     0.01       0.04      -0.08 -0.03  0.13  
lag.TOTAL_EMPL  0.04    -0.03       0.20       0.00 -0.06 -0.04  
lag.TOTAL_HH   -0.29     0.12       0.08       0.03 -0.03  0.02  
lag.POP_DENSIT  0.31    -0.63      -0.01      -0.03  0.04 -0.09  
lag.EMP_DENSIT -0.01     0.00      -0.63      -0.01 -0.07  0.00  
lag.EMP_M       0.04    -0.01      -0.02      -0.47 -0.02  0.04  
lag.AVGWK      -0.02     0.03      -0.05       0.02 -0.55  0.21  
lag.AVGAUTO    -0.07    -0.04       0.02       0.04  0.23 -0.62  
lag.Avg_CarbEM  0.05    -0.06      -0.13      -0.01  0.03 -0.01  
lag.Avg_TRIPSP -0.03     0.09       0.12       0.01 -0.04  0.04  
Avg_CarbEM Avg_TRIPSP lag.ACRES lag.AT lag.TOTAL_EMPL
rho                                                                 
lambda                                                              
(Intercept)                                                         
ACRES                                                               
AT                                                                  
TOTAL_EMPL                                                          
TOTAL_HH                                                            
POP_DENSIT                                                          
EMP_DENSIT                                                          
EMP_M                                                               
AVGWK                                                               
AVGAUTO                                                             
Avg_CarbEM                                                          
Avg_TRIPSP     -0.79                                                
lag.ACRES       0.03       0.03                                     
lag.AT          0.00       0.06      -0.29                          
lag.TOTAL_EMPL  0.04      -0.02       0.11     -0.19                
lag.TOTAL_HH   -0.05       0.11      -0.10      0.33  -0.22         
lag.POP_DENSIT -0.04       0.04       0.19     -0.09   0.09         
lag.EMP_DENSIT -0.11       0.10      -0.12      0.28  -0.36         
lag.EMP_M       0.00       0.01       0.08      0.12  -0.12         
lag.AVGWK      -0.01       0.00       0.05      0.04   0.02         
lag.AVGAUTO    -0.03       0.06      -0.04     -0.21   0.06         
lag.Avg_CarbEM -0.23       0.19      -0.04     -0.08   0.11         
lag.Avg_TRIPSP  0.18      -0.26      -0.08      0.18  -0.11         
lag.TOTAL_HH lag.POP_DENSIT lag.EMP_DENSIT lag.EMP_M
rho                                                                
lambda                                                             
(Intercept)                                                        
ACRES                                                              
AT                                                                 
TOTAL_EMPL                                                         
TOTAL_HH                                                           
POP_DENSIT                                                         
EMP_DENSIT                                                         
EMP_M                                                              
AVGWK                                                              
AVGAUTO                                                            
Avg_CarbEM                                                         
Avg_TRIPSP                                                         
lag.ACRES                                                          
lag.AT                                                             
lag.TOTAL_EMPL                                                     
lag.TOTAL_HH                                                       
lag.POP_DENSIT -0.49                                               
lag.EMP_DENSIT  0.38        -0.09                                  
lag.EMP_M       0.20        -0.14          -0.04                   
lag.AVGWK       0.11        -0.14           0.04           0.14    
lag.AVGAUTO     0.11         0.18           0.04          -0.03    
lag.Avg_CarbEM -0.23         0.07          -0.12          -0.03    
lag.Avg_TRIPSP  0.36        -0.19           0.17           0.06    
lag.AVGWK lag.AVGAUTO lag.Avg_CarbEM
rho                                                
lambda                                             
(Intercept)                                        
ACRES                                              
AT                                                 
TOTAL_EMPL                                         
TOTAL_HH                                           
POP_DENSIT                                         
EMP_DENSIT                                         
EMP_M                                              
AVGWK                                              
AVGAUTO                                            
Avg_CarbEM                                         
Avg_TRIPSP                                         
lag.ACRES                                          
lag.AT                                             
lag.TOTAL_EMPL                                     
lag.TOTAL_HH                                       
lag.POP_DENSIT                                     
lag.EMP_DENSIT                                     
lag.EMP_M                                          
lag.AVGWK                                          
lag.AVGAUTO    -0.53                               
lag.Avg_CarbEM -0.04     -0.09                     
lag.Avg_TRIPSP  0.02      0.03       -0.81


###########residual

res <- mod.mam$residuals
png(filename="plots/Residuals from Manski Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)
cols <- findColours(classes_fx,pal)

par(mar=rep(0,4))
plot(data,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1.5,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from Manski Model ",ncol=5)


dev.off()

## Residual Autocorrelation

moran.test(res, listw=tazw, zero.policy=T)
#Moran's I test under randomisation

data:  res  
weights: tazw  

Moran I statistic standard deviate = -0.167, p-value = 0.5663
alternative hypothesis: greater
sample estimates:
  Moran I statistic       Expectation          Variance 
-0.0052871289     -0.0014450867      0.0005291471 

##############
AIC(mod.mam)
[1] 1064.754

bic=AIC(mod.mam,k = log(length(data)))
bic
[1] 1182.82

logLik(mod.mam)
'log Lik.' -506.3768 (df=26)

