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

setwd("H:/01-Dissertation/R_Spatial_Analysis_2015")
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
########
## Linear Model
########

mod.lm <- lm(CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data=data)
summary(mod.lm)

##########
Call:
lm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data = data)

Residuals:
     Min       1Q   Median       3Q      Max 
-2.05135 -0.25754 -0.03712  0.21129  2.66938 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept) -1.606e-01  1.083e-01  -1.482 0.138675    
ACRES        7.610e-05  4.484e-05   1.697 0.090160 .  
AT           1.868e-01  4.710e-02   3.965 8.10e-05 ***
POP          4.565e-04  1.006e-04   4.535 6.80e-06 ***
TOTAL_HH     4.199e-04  2.174e-04   1.932 0.053826 .  
TOTAL_EMPL  -3.776e-05  2.484e-05  -1.520 0.128934    
POP_DENSIT  -2.128e-02  4.870e-03  -4.369 1.44e-05 ***
EMP_DENSIT   4.214e-04  2.838e-04   1.485 0.138110    
TOTAL_AUTO   4.961e-04  1.040e-04   4.772 2.24e-06 ***
EMP_M        2.396e-01  7.069e-02   3.390 0.000740 ***
AVGWK        1.583e-01  8.502e-02   1.862 0.062971 .  
AVGAUTO     -1.940e-01  7.480e-02  -2.594 0.009688 ** 
Avg_CarbEM  -7.543e+01  2.106e+01  -3.582 0.000365 ***
Avg_TRIPSP   1.989e-02  3.558e-03   5.589 3.32e-08 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Residual standard error: 0.5001 on 679 degrees of freedom
Multiple R-squared: 0.8002,	Adjusted R-squared: 0.7964 
F-statistic: 209.2 on 13 and 679 DF,  p-value: < 2.2e-16 

###########
AIC(mod.lm)
[1] 1022.084
##getting BIC
bic=AIC(mod.lm,k = log(length(data)))
##Extract Log-Likelihood from an lm Object
logLik(mod.lm)
'log Lik.' -496.0418 (df=15)
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

## Residual Autocorrelation
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
 
Moran I statistic standard deviate = 0.11, p-value = 0.4562
alternative hypothesis: greater 
sample estimates:
Moran I statistic       Expectation          Variance 
     0.0010859147     -0.0014450867      0.0005292909 
	
########
## SAR Model (WARNING: This takes a while to run)
########

mod.sar <- lagsarlm(CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL +
+ POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO +
+ Avg_CarbEM + Avg_TRIPSP, data = data, listw=tazw , zero.policy=T, tol.solve=2e-14)
summary(mod.sar)

##results

Call:lagsarlm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, zero.policy = T, 
    tol.solve = 2e-14)

Residuals:
      Min        1Q    Median        3Q       Max 
-2.031904 -0.250650 -0.034679  0.206391  2.695469 

Type: lag 
Coefficients: (asymptotic standard errors) 
               Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -1.3751e-01  1.0741e-01 -1.2803 0.2004344
ACRES        5.6604e-05  4.4711e-05  1.2660 0.2055115
AT           1.9277e-01  4.6517e-02  4.1440 3.414e-05
POP          4.5876e-04  9.9445e-05  4.6132 3.966e-06
TOTAL_HH     4.6802e-04  2.1592e-04  2.1676 0.0301884
TOTAL_EMPL  -4.2244e-05  2.4621e-05 -1.7158 0.0861948
POP_DENSIT  -2.1807e-02  4.8179e-03 -4.5263 6.004e-06
EMP_DENSIT   3.7842e-04  2.8109e-04  1.3462 0.1782265
TOTAL_AUTO   4.7792e-04  1.0364e-04  4.6112 4.004e-06
EMP_M        2.3586e-01  6.9770e-02  3.3805 0.0007236
AVGWK        1.6982e-01  8.4205e-02  2.0168 0.0437169
AVGAUTO     -1.6245e-01  7.5533e-02 -2.1508 0.0314955
Avg_CarbEM  -7.5381e+01  2.0780e+01 -3.6276 0.0002861
Avg_TRIPSP   2.0156e-02  3.5120e-03  5.7392 9.514e-09

Rho: -0.065259, LR test value: 3.8462, p-value: 0.049858
Asymptotic standard error: 0.03343
    z-value: -1.9521, p-value: 0.050927
Wald statistic: 3.8107, p-value: 0.050927

Log likelihood: -494.1186 for lag model
ML residual variance (sigma squared): 0.24351, (sigma: 0.49347)
Number of observations: 693 
Number of parameters estimated: 16 
AIC: 1020.2, (AIC for lm: 1022.1)
LM test for residual autocorrelation
test value: 1.31, p-value: 0.2524

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
 
Moran I statistic standard deviate = 1.0418, p-value = 0.1488
alternative hypothesis: greater 
sample estimates:
Moran I statistic       Expectation          Variance 
     0.0225214042     -0.0014450867      0.0005292397 

###########
AIC(mod.sar)
[1] 1020.237

bic=AIC(mod.sar,k = log(length(data)))
bic
[1] 1092.894

logLik(mod.sar)
'log Lik.' -494.1186 (df=16)
########
## SEM Model (WARNING: This takes a while to run)
########


mod.sem <- errorsarlm(CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, zero.policy=T, tol.solve=1e-15)
summary(mod.sem)

Call:errorsarlm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + 
    TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + 
    AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, 
    zero.policy = T, tol.solve = 1e-15)

Residuals:
     Min       1Q   Median       3Q      Max 
-2.05138 -0.25715 -0.03673  0.21124  2.66950 

Type: error 
Coefficients: (asymptotic standard errors) 
               Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -1.6075e-01  1.0734e-01 -1.4976 0.1342316
ACRES        7.6827e-05  4.4454e-05  1.7282 0.0839488
AT           1.8654e-01  4.6688e-02  3.9955 6.455e-05
POP          4.5731e-04  9.9709e-05  4.5864 4.509e-06
TOTAL_HH     4.1992e-04  2.1537e-04  1.9498 0.0512003
TOTAL_EMPL  -3.7769e-05  2.4594e-05 -1.5357 0.1246102
POP_DENSIT  -2.1299e-02  4.8254e-03 -4.4139 1.015e-05
EMP_DENSIT   4.2090e-04  2.8127e-04  1.4965 0.1345346
TOTAL_AUTO   4.9510e-04  1.0299e-04  4.8074 1.529e-06
EMP_M        2.3961e-01  6.9998e-02  3.4230 0.0006192
AVGWK        1.5857e-01  8.4201e-02  1.8832 0.0596732
AVGAUTO     -1.9404e-01  7.4093e-02 -2.6189 0.0088227
Avg_CarbEM  -7.5409e+01  2.0845e+01 -3.6176 0.0002974
Avg_TRIPSP   1.9876e-02  3.5225e-03  5.6424 1.677e-08

Lambda: 0.0032987, LR test value: 0.0024828, p-value: 0.96026
Asymptotic standard error: 0.062279
    z-value: 0.052967, p-value: 0.95776
Wald statistic: 0.0028055, p-value: 0.95776

Log likelihood: -496.0405 for error model
ML residual variance (sigma squared): 0.24505, (sigma: 0.49502)
Number of observations: 693 
Number of parameters estimated: 16 
AIC: 1024.1, (AIC for lm: 1022.1)

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
 
Moran I statistic standard deviate = 0.0639, p-value = 0.4745
alternative hypothesis: greater 
sample estimates:
Moran I statistic       Expectation          Variance 
     2.511625e-05     -1.445087e-03      5.292904e-04 

###########
AIC(mod.sem)
[1] 1024.081

bic=AIC(mod.sem,k = log(length(data)))
bic
[1] 1096.738

logLik(mod.sem)
'log Lik.' -496.0405 (df=16)
########
## SDM Model (WARNING: This takes a while to run)
########


mod.sdm <- lagsarlm(CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, zero.policy=T, type="mixed", tol.solve=1.47932e-15)
summary(mod.sdm)

Call:lagsarlm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
    POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
    Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, type = "mixed", 
    zero.policy = T, tol.solve = 1.47932e-15)

Residuals:
      Min        1Q    Median        3Q       Max 
-1.954091 -0.263275 -0.047956  0.217257  2.794910 

Type: mixed 
Coefficients: (asymptotic standard errors) 
                  Estimate  Std. Error z value  Pr(>|z|)
(Intercept)    -6.3735e-02  1.7691e-01 -0.3603 0.7186379
ACRES           2.1679e-04  5.6788e-05  3.8176 0.0001347
AT              3.0049e-01  8.3073e-02  3.6172 0.0002978
POP             5.5010e-04  1.0738e-04  5.1227 3.011e-07
TOTAL_HH        6.0516e-04  2.3397e-04  2.5865 0.0096956
TOTAL_EMPL     -5.0795e-05  2.5444e-05 -1.9964 0.0458939
POP_DENSIT     -2.5502e-02  5.5321e-03 -4.6099 4.029e-06
EMP_DENSIT      3.7456e-04  3.4305e-04  1.0918 0.2749051
TOTAL_AUTO      2.5572e-04  1.1208e-04  2.2816 0.0225122
EMP_M           1.8658e-01  7.1481e-02  2.6102 0.0090493
AVGWK           1.9110e-01  8.8355e-02  2.1629 0.0305500
AVGAUTO        -8.0334e-02  8.7730e-02 -0.9157 0.3598300
Avg_CarbEM     -7.5783e+01  2.0472e+01 -3.7019 0.0002140
Avg_TRIPSP      1.9839e-02  3.4686e-03  5.7195 1.068e-08
lag.ACRES      -3.3558e-04  7.7041e-05 -4.3558 1.326e-05
lag.AT         -2.8099e-02  1.1030e-01 -0.2547 0.7989266
lag.POP        -4.0606e-04  1.8026e-04 -2.2527 0.0242809
lag.TOTAL_HH   -1.7116e-04  3.9631e-04 -0.4319 0.6658255
lag.TOTAL_EMPL -1.5352e-05  4.9389e-05 -0.3109 0.7559136
lag.POP_DENSIT  1.5029e-02  9.7010e-03  1.5492 0.1213218
lag.EMP_DENSIT  2.5217e-04  5.6187e-04  0.4488 0.6535716
lag.TOTAL_AUTO  6.1609e-04  2.0102e-04  3.0649 0.0021775
lag.EMP_M      -9.2180e-03  1.4680e-01 -0.0628 0.9499302
lag.AVGWK      -1.2454e-01  1.7357e-01 -0.7175 0.4730560
lag.AVGAUTO    -1.0240e-01  1.5746e-01 -0.6503 0.5154900
lag.Avg_CarbEM -4.7620e+01  4.8346e+01 -0.9850 0.3246325
lag.Avg_TRIPSP  1.0390e-02  8.0088e-03  1.2973 0.1945310

Rho: -0.048167, LR test value: 0.60202, p-value: 0.43781
Asymptotic standard error: 0.062502
    z-value: -0.77065, p-value: 0.44091
Wald statistic: 0.59391, p-value: 0.44091

Log likelihood: -474.5257 for mixed model
ML residual variance (sigma squared): 0.2302, (sigma: 0.47979)
Number of observations: 693 
Number of parameters estimated: 29 
AIC: 1007.1, (AIC for lm: 1005.7)
LM test for residual autocorrelation
test value: 10.52, p-value: 0.0011811

#############################################

res <- mod.sdm$residuals
png(filename="plots/Residuals from Spatial Durbin Model.png", width = 800, height = 600)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=c(-2.5,-1.5,-0.5,0.5,1.5,2.5), rtimes = 1)
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
 
Moran I statistic standard deviate = -0.422, p-value = 0.6635
alternative hypothesis: greater 
sample estimates:
Moran I statistic       Expectation          Variance 
    -0.0111517507     -0.0014450867      0.0005290678 


###########
AIC(mod.sdm)
[1] 1024.081

bic=AIC(mod.sdm,k = log(length(data)))
bic
[1] 1138.741

logLik(mod.sdm)
'log Lik.' -474.5257 (df=16)


########
## Kelejian-Prucha Model 
########

##mod.kpm <- GMerrorsar(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP,  
                      data,  tazw,verbose=TRUE)
mod.kpm <- gstsls(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP,  
                      data,  tazw,verbose=TRUE)


summary(mod.kpm, correlation = FALSE, Hausman=FALSE)

##results

Call:gstsls(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
              POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
              Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, verbose = TRUE)

Residuals:
  Min        1Q    Median        3Q       Max 
-2.031648 -0.257512 -0.040764  0.206467  2.693383 

Type: GM SARAR estimator
Coefficients: (GM standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
Rho_Wy      -6.0503e-02  3.9722e-02 -1.5232 0.1277177
(Intercept) -1.4244e-01  1.1185e-01 -1.2735 0.2028571
ACRES        7.3671e-05  4.7562e-05  1.5490 0.1213907
AT           1.8828e-01  4.8592e-02  3.8746 0.0001068
POP          4.7496e-04  1.0199e-04  4.6568 3.211e-06
TOTAL_HH     4.6491e-04  2.2205e-04  2.0937 0.0362885
TOTAL_EMPL  -4.1916e-05  2.5070e-05 -1.6720 0.0945227
POP_DENSIT  -2.2093e-02  4.9560e-03 -4.4577 8.283e-06
EMP_DENSIT   3.7225e-04  2.9029e-04  1.2823 0.1997330
TOTAL_AUTO   4.5638e-04  1.0577e-04  4.3146 1.599e-05
EMP_M        2.3543e-01  7.1105e-02  3.3111 0.0009294
AVGWK        1.7369e-01  8.5950e-02  2.0209 0.0432940
AVGAUTO     -1.6451e-01  7.8168e-02 -2.1045 0.0353312
Avg_CarbEM  -7.4727e+01  2.0979e+01 -3.5619 0.0003681
Avg_TRIPSP   1.9902e-02  3.5564e-03  5.5962 2.190e-08

Lambda: 0.068083
Residual variance (sigma squared): 0.24841, (sigma: 0.49841)
GM argmin sigma squared: 0.24391
Number of observations: 693 
Number of parameters estimated: 17 


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

Moran I statistic standard deviate = -0.0133, p-value = 0.5053
alternative hypothesis: greater
sample estimates:
Moran I statistic       Expectation          Variance 
     -0.001751955      -0.001445087       0.000529229 



########
## Spatial Durbin Error Model 
########

ev <- eigenw(similar.listw(tazw))
mod.SDEm <- errorsarlm(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP,  data,  tazw,method="eigen", control=list(pre_eig=ev),tol.solve=1.384e-14)

summary(mod.SDEm, correlation=TRUE)

### Results

Call:errorsarlm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + 
                  TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + 
                  AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, 
                method = "eigen", tol.solve = 1.384e-14, control = list(pre_eig = ev))

Residuals:
  Min       1Q   Median       3Q      Max 
-2.05138 -0.25715 -0.03673  0.21124  2.66950 

Type: error 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept) -1.6075e-01  1.0734e-01 -1.4976 0.1342316
ACRES        7.6827e-05  4.4454e-05  1.7282 0.0839488
AT           1.8654e-01  4.6688e-02  3.9955 6.455e-05
POP          4.5731e-04  9.9709e-05  4.5864 4.509e-06
TOTAL_HH     4.1992e-04  2.1537e-04  1.9498 0.0512003
TOTAL_EMPL  -3.7769e-05  2.4594e-05 -1.5357 0.1246102
POP_DENSIT  -2.1299e-02  4.8254e-03 -4.4139 1.015e-05
EMP_DENSIT   4.2090e-04  2.8127e-04  1.4965 0.1345346
TOTAL_AUTO   4.9510e-04  1.0299e-04  4.8074 1.529e-06
EMP_M        2.3961e-01  6.9998e-02  3.4230 0.0006192
AVGWK        1.5857e-01  8.4201e-02  1.8832 0.0596732
AVGAUTO     -1.9404e-01  7.4093e-02 -2.6189 0.0088227
Avg_CarbEM  -7.5409e+01  2.0845e+01 -3.6176 0.0002974
Avg_TRIPSP   1.9876e-02  3.5225e-03  5.6424 1.677e-08

Lambda: 0.0032987, LR test value: 0.0024828, p-value: 0.96026
Asymptotic standard error: 0.062279
z-value: 0.052967, p-value: 0.95776
Wald statistic: 0.0028055, p-value: 0.95776

Log likelihood: -496.0405 for error model
ML residual variance (sigma squared): 0.24505, (sigma: 0.49502)
Number of observations: 693 
Number of parameters estimated: 16 
AIC: 1024.1, (AIC for lm: 1022.1)

Correlation of coefficients 
sigma lambda (Intercept) ACRES AT    POP   TOTAL_HH TOTAL_EMPL POP_DENSIT EMP_DENSIT TOTAL_AUTO
lambda       0.00                                                                                          
(Intercept)  0.00  0.00                                                                                    
ACRES        0.00  0.00   0.05                                                                             
AT           0.00  0.00  -0.25       -0.30                                                                 
POP          0.00  0.00   0.08        0.07 -0.17                                                           
TOTAL_HH     0.00  0.00  -0.33       -0.07  0.08 -0.70                                                     
TOTAL_EMPL   0.00  0.00  -0.27        0.11 -0.15  0.06 -0.06                                               
POP_DENSIT   0.00  0.00  -0.24        0.20  0.16 -0.03 -0.23     0.06                                      
EMP_DENSIT   0.00  0.00  -0.22       -0.12  0.30 -0.05  0.20    -0.40       0.00                           
TOTAL_AUTO   0.00  0.00   0.34       -0.05  0.14 -0.44 -0.28     0.01       0.19      -0.15                
EMP_M        0.00  0.00  -0.39        0.02  0.11 -0.04 -0.01    -0.03      -0.02      -0.11       0.02     
AVGWK        0.00  0.00  -0.28        0.03  0.08 -0.19  0.13     0.11      -0.07       0.01       0.06     
AVGAUTO      0.00  0.00  -0.30       -0.04 -0.45  0.17  0.18     0.05       0.06       0.07      -0.49     
Avg_CarbEM   0.00  0.00  -0.04       -0.10  0.00  0.02  0.06     0.01      -0.01      -0.02      -0.08     
Avg_TRIPSP   0.00  0.00   0.01        0.00  0.00 -0.01 -0.07     0.05      -0.05       0.02       0.01     
EMP_M AVGWK AVGAUTO Avg_CarbEM
lambda                                    
(Intercept)                               
ACRES                                     
AT                                        
POP                                       
TOTAL_HH                                  
TOTAL_EMPL                                
POP_DENSIT                                
EMP_DENSIT                                
TOTAL_AUTO                                
EMP_M                                     
AVGWK        0.02                         
AVGAUTO     -0.01 -0.56                   
Avg_CarbEM  -0.01  0.01  0.05             
Avg_TRIPSP  -0.04 -0.05 -0.06   -0.79    

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

Moran I statistic standard deviate = 0.0639, p-value = 0.4745
alternative hypothesis: greater
sample estimates:
Moran I statistic       Expectation          Variance 
2.512062e-05     -1.445087e-03      5.292904e-04 


########
AIC(mod.SDEm)
[1] 1024.081

bic=AIC(mod.SDEm,k = log(length(data)))
bic
[1] 1096.738

logLik(mod.SDEm)
'log Lik.' -496.0405 (df=16)

########
## Manski Model 
########

mod.mam <- sacsarlm(CarEMTAZ ~  ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + Avg_TRIPSP,  data,  tazw, type="sacmixed")
summary(mod.mam, correlation=TRUE)

###results
Call:sacsarlm(formula = CarEMTAZ ~ ACRES + AT + POP + TOTAL_HH + TOTAL_EMPL + 
                POP_DENSIT + EMP_DENSIT + TOTAL_AUTO + EMP_M + AVGWK + AVGAUTO + 
                Avg_CarbEM + Avg_TRIPSP, data = data, listw = tazw, type = "sacmixed")

Residuals:
  Min        1Q    Median        3Q       Max 
-1.938546 -0.265238 -0.050193  0.209100  2.758267 

Type: sacmixed 
Coefficients: (asymptotic standard errors) 
Estimate  Std. Error z value  Pr(>|z|)
(Intercept)    -3.1665e-02  1.4991e-01 -0.2112 0.8327148
ACRES           2.4125e-04  6.1591e-05  3.9169 8.970e-05
AT              3.1126e-01  8.5876e-02  3.6245 0.0002895
POP             5.7653e-04  1.1256e-04  5.1220 3.023e-07
TOTAL_HH        6.3272e-04  2.4113e-04  2.6239 0.0086922
TOTAL_EMPL     -4.9611e-05  2.6061e-05 -1.9037 0.0569541
POP_DENSIT     -2.6776e-02  5.7182e-03 -4.6826 2.832e-06
EMP_DENSIT      3.2166e-04  3.5927e-04  0.8953 0.3706238
TOTAL_AUTO      2.0821e-04  1.1780e-04  1.7675 0.0771479
EMP_M           1.8258e-01  7.3569e-02  2.4818 0.0130724
AVGWK           1.8440e-01  9.1021e-02  2.0259 0.0427734
AVGAUTO        -6.2193e-02  9.0389e-02 -0.6881 0.4914155
Avg_CarbEM     -7.4597e+01  2.0602e+01 -3.6208 0.0002937
Avg_TRIPSP      1.9578e-02  3.5209e-03  5.5606 2.689e-08
lag.ACRES      -3.5144e-04  7.5495e-05 -4.6552 3.237e-06
lag.AT         -8.9435e-02  1.1967e-01 -0.7474 0.4548472
lag.POP        -5.0723e-04  1.8663e-04 -2.7179 0.0065706
lag.TOTAL_HH   -3.1228e-04  3.8535e-04 -0.8104 0.4177248
lag.TOTAL_EMPL -6.9386e-06  4.7345e-05 -0.1466 0.8834844
lag.POP_DENSIT  2.0215e-02  9.6881e-03  2.0865 0.0369315
lag.EMP_DENSIT  1.6161e-04  5.5173e-04  0.2929 0.7695875
lag.TOTAL_AUTO  5.3039e-04  2.4024e-04  2.2078 0.0272590
lag.EMP_M      -5.7492e-02  1.4681e-01 -0.3916 0.6953399
lag.AVGWK      -1.0398e-01  1.6413e-01 -0.6335 0.5263834
lag.AVGAUTO    -1.0840e-01  1.5434e-01 -0.7024 0.4824597
lag.Avg_CarbEM -2.5863e+01  4.9251e+01 -0.5251 0.5995018
lag.Avg_TRIPSP  4.6759e-03  9.1173e-03  0.5129 0.6080466

Rho: 0.18309
Asymptotic standard error: 0.23527
z-value: 0.77824, p-value: 0.43643
Lambda: -0.2717
Asymptotic standard error: 0.26903
z-value: -1.0099, p-value: 0.31253

LR test value: 45.217, p-value: 7.0756e-05

Log likelihood: -473.4331 for sacmixed model
ML residual variance (sigma squared): 0.22533, (sigma: 0.47469)
Number of observations: 693 
Number of parameters estimated: 30 
AIC: 1006.9, (AIC for lm: 1022.1)

Correlation of coefficients 
sigma rho   lambda (Intercept) ACRES AT    POP   TOTAL_HH TOTAL_EMPL POP_DENSIT EMP_DENSIT
rho            -0.59                                                                                     
lambda          0.60 -0.97                                                                               
(Intercept)    -0.11  0.18 -0.18                                                                         
ACRES          -0.16  0.27 -0.26   0.06                                                                  
AT              0.03 -0.05  0.05   0.00       -0.20                                                      
POP            -0.10  0.17 -0.16   0.04        0.08 -0.06                                                
TOTAL_HH        0.00  0.01 -0.01  -0.01       -0.02  0.08 -0.67                                          
TOTAL_EMPL     -0.02  0.03 -0.03  -0.06        0.01 -0.01  0.01  0.00                                    
POP_DENSIT      0.05 -0.08  0.08  -0.01        0.19 -0.02 -0.03 -0.24     0.08                           
EMP_DENSIT      0.04 -0.06  0.06   0.04        0.01  0.02 -0.03  0.07    -0.46       0.01                
TOTAL_AUTO      0.11 -0.19  0.18  -0.03       -0.15 -0.01 -0.47 -0.29     0.00       0.19      -0.06     
EMP_M           0.00 -0.01  0.01  -0.03       -0.02  0.08 -0.07 -0.02    -0.02       0.04      -0.02     
AVGWK          -0.02  0.03 -0.03   0.04        0.03  0.05 -0.20  0.17     0.11      -0.08       0.05     
AVGAUTO        -0.01  0.01 -0.01   0.07       -0.02 -0.10  0.13  0.28     0.03      -0.04       0.04     
Avg_CarbEM     -0.05  0.08 -0.08   0.04       -0.03 -0.03  0.05  0.04    -0.01       0.02       0.06     
Avg_TRIPSP      0.06 -0.09  0.09  -0.04       -0.05 -0.02 -0.02 -0.06     0.05      -0.07      -0.03     
lag.ACRES       0.02 -0.03  0.03  -0.03       -0.71  0.17 -0.01  0.00     0.04      -0.15      -0.03     
lag.AT          0.26 -0.43  0.42  -0.15        0.02 -0.72 -0.05 -0.04    -0.03       0.05       0.04     
lag.POP         0.24 -0.41  0.40  -0.07       -0.11  0.03 -0.60  0.35     0.01       0.04       0.02     
lag.TOTAL_HH    0.15 -0.25  0.24  -0.28       -0.10 -0.02  0.32 -0.55    -0.03       0.16      -0.02     
lag.TOTAL_EMPL -0.14  0.24 -0.24  -0.15        0.13 -0.02  0.06  0.01    -0.36      -0.04       0.19     
lag.POP_DENSIT -0.22  0.37 -0.36  -0.08       -0.03  0.01  0.06  0.13    -0.06      -0.62       0.01     
lag.EMP_DENSIT  0.12 -0.21  0.20  -0.28       -0.07 -0.02 -0.04 -0.05     0.26       0.03      -0.61     
lag.TOTAL_AUTO  0.37 -0.62  0.60   0.13       -0.08  0.05  0.12  0.11    -0.02      -0.02       0.09     
lag.EMP_M       0.21 -0.35  0.34  -0.41       -0.10 -0.02 -0.02  0.01    -0.01       0.00      -0.02     
lag.AVGWK       0.08 -0.13  0.13  -0.24       -0.06  0.02  0.10 -0.09    -0.06       0.04      -0.06     
lag.AVGAUTO    -0.14  0.24 -0.23  -0.18        0.04  0.04 -0.04 -0.17    -0.02       0.00      -0.05     
lag.Avg_CarbEM -0.26  0.44 -0.42   0.08        0.10 -0.02  0.07 -0.01     0.07      -0.09      -0.13     
lag.Avg_TRIPSP  0.35 -0.59  0.57  -0.13       -0.10  0.03 -0.13  0.04    -0.03       0.11       0.11     
TOTAL_AUTO EMP_M AVGWK AVGAUTO Avg_CarbEM Avg_TRIPSP lag.ACRES lag.AT lag.POP lag.TOTAL_HH
rho                                                                                                      
lambda                                                                                                   
(Intercept)                                                                                              
ACRES                                                                                                    
AT                                                                                                       
POP                                                                                                      
TOTAL_HH                                                                                                 
TOTAL_EMPL                                                                                               
POP_DENSIT                                                                                               
EMP_DENSIT                                                                                               
TOTAL_AUTO                                                                                               
EMP_M           0.07                                                                                     
AVGWK           0.04       0.00                                                                          
AVGAUTO        -0.50      -0.07 -0.37                                                                    
Avg_CarbEM     -0.10      -0.01  0.02  0.06                                                              
Avg_TRIPSP      0.04      -0.04 -0.04 -0.06   -0.79                                                      
lag.ACRES       0.10       0.03 -0.03  0.00    0.00       0.03                                           
lag.AT          0.11      -0.07 -0.03  0.05   -0.02       0.07      -0.25                                
lag.POP         0.34       0.07  0.10 -0.08   -0.07       0.05       0.10      0.06                      
lag.TOTAL_HH    0.20       0.00 -0.10 -0.18   -0.03       0.05      -0.08      0.23  -0.49               
lag.TOTAL_EMPL -0.07       0.00 -0.06  0.00    0.05      -0.02       0.11     -0.22   0.03   -0.22       
lag.POP_DENSIT -0.14      -0.04  0.04  0.00   -0.02       0.03       0.15     -0.14  -0.14   -0.32       
lag.EMP_DENSIT  0.12      -0.01 -0.07 -0.08   -0.12       0.10      -0.11      0.27   0.02    0.35       
lag.TOTAL_AUTO -0.31      -0.04 -0.05  0.24    0.00       0.05      -0.02      0.30  -0.08   -0.06       
lag.EMP_M       0.03      -0.42 -0.03  0.03   -0.01       0.02       0.10      0.18   0.10    0.06       
lag.AVGWK      -0.03       0.01 -0.51  0.17    0.00       0.00       0.05      0.08  -0.13    0.14       
lag.AVGAUTO     0.25       0.06  0.17 -0.58   -0.03       0.02       0.00     -0.32   0.09    0.13       
lag.Avg_CarbEM -0.07      -0.01  0.03  0.00   -0.15       0.12      -0.07     -0.18  -0.17   -0.07       
lag.Avg_TRIPSP  0.11       0.01 -0.03  0.01    0.11      -0.17      -0.03      0.27   0.25    0.07       
lag.TOTAL_EMPL lag.POP_DENSIT lag.EMP_DENSIT lag.TOTAL_AUTO lag.EMP_M lag.AVGWK lag.AVGAUTO
rho                                                                                                       
lambda                                                                                                    
(Intercept)                                                                                               
ACRES                                                                                                     
AT                                                                                                        
POP                                                                                                       
TOTAL_HH                                                                                                  
TOTAL_EMPL                                                                                                
POP_DENSIT                                                                                                
EMP_DENSIT                                                                                                
TOTAL_AUTO                                                                                                
EMP_M                                                                                                     
AVGWK                                                                                                     
AVGAUTO                                                                                                   
Avg_CarbEM                                                                                                
Avg_TRIPSP                                                                                                
lag.ACRES                                                                                                 
lag.AT                                                                                                    
lag.POP                                                                                                   
lag.TOTAL_HH                                                                                              
lag.TOTAL_EMPL                                                                                            
lag.POP_DENSIT  0.10                                                                                      
lag.EMP_DENSIT -0.38          -0.10                                                                       
lag.TOTAL_AUTO -0.08          -0.12          -0.05                                                        
lag.EMP_M      -0.13          -0.17          -0.05           0.24                                         
lag.AVGWK      -0.01          -0.14           0.03           0.17           0.15                          
lag.AVGAUTO     0.09           0.19           0.07          -0.57          -0.15     -0.52                
lag.Avg_CarbEM  0.13           0.12          -0.12          -0.30          -0.11     -0.07      0.09      
lag.Avg_TRIPSP -0.12          -0.23           0.15           0.36           0.15      0.04     -0.16      
lag.Avg_CarbEM
rho                          
lambda                       
(Intercept)                  
ACRES                        
AT                           
POP                          
TOTAL_HH                     
TOTAL_EMPL                   
POP_DENSIT                   
EMP_DENSIT                   
TOTAL_AUTO                   
EMP_M                        
AVGWK                        
AVGAUTO                      
Avg_CarbEM                   
Avg_TRIPSP                   
lag.ACRES                    
lag.AT                       
lag.POP                      
lag.TOTAL_HH                 
lag.TOTAL_EMPL               
lag.POP_DENSIT               
lag.EMP_DENSIT               
lag.TOTAL_AUTO               
lag.EMP_M                    
lag.AVGWK                    
lag.AVGAUTO                  
lag.Avg_CarbEM               
lag.Avg_TRIPSP -0.83      

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

Moran I statistic standard deviate = -0.2155, p-value = 0.5853
alternative hypothesis: greater
sample estimates:
Moran I statistic       Expectation          Variance 
-0.0064017950     -0.0014450867      0.0005290743 

##############
AIC(mod.mam)
[1] 1006.866

bic=AIC(mod.mam,k = log(length(data)))
bic
[1] 1143.097

logLik(mod.mam)
'log Lik.' -473.4331 (df=30)

