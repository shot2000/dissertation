Interpreting Spatial Lag Models

Load libraries and data

library(classInt)
## Warning: package 'classInt' was built under R version 2.15.3
## Loading required package: class
## Loading required package: e1071
## Warning: package 'e1071' was built under R version 2.15.3
library(spdep)
## Warning: package 'spdep' was built under R version 2.15.3
## Loading required package: sp
## Loading required package: boot
## Loading required package: Matrix
## Loading required package: lattice
## Attaching package: 'lattice'
## The following object(s) are masked from 'package:boot':
## 
## melanoma
## Loading required package: MASS
## Loading required package: nlme
## Loading required package: maptools
## Warning: package 'maptools' was built under R version 2.15.3
## Loading required package: foreign
## Loading required package: grid
## Checking rgeos availability: FALSE Note: when rgeos is not available,
## polygon geometry computations in maptools depend on gpclib, which has a
## restricted licence. It is disabled by default; to enable gpclib, type
## gpclibPermit()
## Loading required package: deldir
## deldir 0.0-21
## Loading required package: coda
## Warning: package 'coda' was built under R version 2.15.3
## Loading required package: splines
library(RColorBrewer)
library(gstat)
## Warning: package 'gstat' was built under R version 2.15.3
load("C:/Users/molly/Documents/1_Coursework/courses by topic/Geography/Methods in Geography/Data/soco.rda")
Build weights matrix

soco_nbq <- poly2nb(soco)  #Queen's neighborhood, row standarized
soco_nbq_w <- nb2listw(soco_nbq)
Examine variogram of outcome of interest, percent children living in povery (PPOV)

plot(variogram(soco$PPOV ~ 1, locations = coordinates(soco), data = soco, cloud = F), 
    type = "b", pch = 16, main = "Variogram of PPOV")
plot of chunk unnamed-chunk-3

Examine Moran's I Build weights matrix

moran.mc(soco$PPOV, soco_nbq_w, nsim = 999)
## 
##  Monte-Carlo simulation of Moran's I
## 
## data:  soco$PPOV 
## weights: soco_nbq_w  
## number of simulations + 1: 1000 
##  
## statistic = 0.5893, observed rank = 1000, p-value = 0.001
## alternative hypothesis: greater
Variogram and Moran's I indicate there is spatial autocorrelation in PPOV.

Build a spatial lag model

soco_LAG <- lagsarlm(PPOV ~ PFHH + PUNEM + PBLK + P65UP, data = soco, soco_nbq_w)
summary(soco_LAG)
## 
## Call:lagsarlm(formula = PPOV ~ PFHH + PUNEM + PBLK + P65UP, data = soco, 
##     listw = soco_nbq_w)
## 
## Residuals:
##        Min         1Q     Median         3Q        Max 
## -0.2457653 -0.0284360 -0.0028762  0.0262169  0.2374894 
## 
## Type: lag 
## Coefficients: (asymptotic standard errors) 
##              Estimate Std. Error  z value  Pr(>|z|)
## (Intercept) -0.100260   0.007375 -13.5946 < 2.2e-16
## PFHH         0.429404   0.040246  10.6695 < 2.2e-16
## PUNEM        1.354637   0.065959  20.5374 < 2.2e-16
## PBLK        -0.069046   0.015335  -4.5025 6.716e-06
## P65UP        0.291192   0.035210   8.2701 2.220e-16
## 
## Rho: 0.5172, LR test value: 491.5, p-value: < 2.22e-16
## Asymptotic standard error: 0.02118
##     z-value: 24.42, p-value: < 2.22e-16
## Wald statistic: 596.3, p-value: < 2.22e-16
## 
## Log likelihood: 2239 for lag model
## ML residual variance (sigma squared): 0.002192, (sigma: 0.04682)
## Number of observations: 1387 
## Number of parameters estimated: 7 
## AIC: -4464, (AIC for lm: -3974)
## LM test for residual autocorrelation
## test value: 4.85, p-value: 0.027651
Results indicate that predictors and the lag are significant. the z-test is based on the standard error of Moran's I (0.51719/0.02118=24.42). The LR test is a test of the model with and without the spatial lag. The reported LR suggests that the addition of the lag is an improvement, it is formally equivalent to:

lm1 <- lm(PPOV ~ PFHH + PUNEM + PBLK + P65UP, data = soco)
anova(soco_LAG, lm1)
##          Model df   AIC logLik Test L.Ratio p-value
## soco_LAG     1  7 -4464   2239    1                
## lm1          2  6 -3974   1993    2     491       0
What would happen to PPOV in the southeast if the unemployment rate rose from 6% to 75% in Jefferson County, Alabama.

soco.new <- soco  #copy the data frame (so we don't mess up the original)

# Change the unemployment rate
soco.new@data[soco.new@data$CNTY_ST == "Jefferson County  AL", "PUNEM"] <- 0.75

# The original predicted values
orig.pred <- as.data.frame(predict(soco_LAG))

# The predicted values with the new unemployment rate in Alabama
new.pred <- as.data.frame(predict(soco_LAG, newdata = soco.new, listw = soco_nbq_w))

# the difference between the predicted values
jCoAL_effect <- new.pred$fit - orig.pred$fit
el <- data.frame(name = soco$CNTY_ST, diff_in_pred_PPOV = jCoAL_effect)
soco.new$jee <- el$diff_in_pred_PPOV

# sort counties by absolute value of the change in predicted PPOV
el <- el[rev(order(abs(el$diff_in_pred_PPOV))), ]
el[1:10, ]  #show the top 10 counties
##                       name diff_in_pred_PPOV
## 37    Jefferson County  AL           0.99079
## 4          Bibb County  AL           0.11637
## 58    St. Clair County  AL           0.11165
## 5        Blount County  AL           0.10881
## 64       Walker County  AL           0.10738
## 59       Shelby County  AL           0.10516
## 63   Tuscaloosa County  AL           0.09641
## 1050    El Paso County  TX          -0.08629
## 1197     Sutton County  TX          -0.08108
## 437         Lee County  KY          -0.07285
We observe that counties as far away as Kentucky experienced substantial change in their predicted poverty rate. The changes in PPOV are both positive and negative.

Map changes in predicted poverty

breaks <- c(min(soco.new$jee), -0.03, 0.03, max(soco.new$jee))
labels <- c("Negative effect (greater than .03)", "No Effect (effect -.03 to .03)", 
    "Positive effect (greater than .03)")


np <- findInterval(soco.new$jee, breaks)
colors <- c("red", "yellow", "blue")

# Draw Map
plot(soco.new, col = colors[np])
mtext("Effects of a change in Jefferson County, AL (set PUNEM = .75)\n on predicted values in a spatial lag model", 
    side = 3, line = 1)
legend("topleft", legend = labels, fill = colors, bty = "n")
plot of chunk unnamed-chunk-8

Map the magnitude of the changes caused by altering the PPOV in Jefferson County, AL.

pal5 <- brewer.pal(6, "Spectral")
cats5 <- classIntervals(soco.new$jee, n = 5, style = "jenks")
colors5 <- findColours(cats5, pal5)
plot(soco.new, col = colors5)
legend("topleft", legend = round(cats5$brks, 2), fill = pal5, bty = "n")
mtext("Effects of a change in Jefferson County, AL (set PUNEM = .75)\n on predicted values in a spatial lag model", 
    side = 3, line = 1)
plot of chunk unnamed-chunk-9

Use the impacts() function to understand the direct (local), indirect(spill-over), and total effect of a unit change in each of the predictor variables. The changes reported by impacts are the global average impact:

impacts(soco_LAG, listw = soco_nbq_w)
## Impact measures (lag, exact):
##         Direct Indirect   Total
## PFHH   0.45695  0.43243  0.8894
## PUNEM  1.44154  1.36419  2.8057
## PBLK  -0.07348 -0.06953 -0.1430
## P65UP  0.30987  0.29325  0.6031
The output from impacts says that a 100% increase in unemployment leads to a 280% increase in childhood poverty rate (since 1 unit here = 100%)