data(boston) 
gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) 
  +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT) + TOWN, 
  data=boston.c, nb2listw(boston.soi), method="eigen", tol.solve=1e-12) 
  
  boston<- boston.c				
boston$CRIM <- boston.c$CRIM/4
  
p2 <-as.data.frame(predict(gp2, newdata=boston, listw=nb2listw(boston.soi)))


ev <- eigenw(similar.listw(tazw))
mod.SDEm <- errorsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                         POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                         Avg_TRIPSP, data,  tazw,method="eigen", control=list(pre_eig=ev),tol.solve=1.384e-14)
						 
mod.mam <- sacsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                      POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                      Avg_TRIPSP,  data,  tazw, listw2 = NULL, type="sac",method = "eigen", quiet = NULL, zero.policy = NULL, tol.solve = 1e-10,llprof=NULL, interval1=NULL, interval2=NULL, trs1=NULL, trs2=NULL,control = list())


mod.mam <- sacsarlm(CarEMTAZ ~ ACRES + AT + TOTAL_EMPL + TOTAL_HH + 
                      POP_DENSIT + EMP_DENSIT + EMP_M + AVGWK + AVGAUTO + Avg_CarbEM + 
                      Avg_TRIPSP,  data,  tazw, type="sacmixed")
						 
dataT_EMP<- data #copy the data frame (so we don't mess up the original)

# Change the unemployment rate
dataT_EMP$TOTAL_EMPL <- data$TOTAL_EMPL/2

# The original predicted values
orig.pred <- as.data.frame(predict(mod.mam))

# The predicted values with the new Total emp
new.pred <- as.data.frame(predict(mod.SDEm, newdata = dataT_EMP, listw=tazw))

new.pred <- as.data.frame(predict(mod.mam, newdata = dataT_EMP, listw=tazw))








#stepwise VIF function with preallocated vectors
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
   
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
   
  #get initial vif value for all comparisons of variables
  vif_init<-vector('list', length = ncol(in_frame))
  names(vif_init) <- names(in_frame)
  for(val in names(in_frame)){
      form_in<-formula(paste(val,' ~ .'))
      vif_init[[val]]<-VIF(lm(form_in,data=in_frame,...))
      }
  vif_max<-max(unlist(vif_init))
   
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
        prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
        cat('\n')
        cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
        }
    return(names(in_frame))
    }
  else{
  
    in_dat<-in_frame
  
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
       
      vif_vals<-vector('list', length = ncol(in_dat))
      names(vif_vals) <- names(in_dat)
  
      for(val in names(in_dat)){
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat,...))
        vif_vals[[val]]<-vif_add
        }
 
      max_row<-which.max(vif_vals)[1]
  
      vif_max<-vif_vals[max_row]
  
      if(vif_max<thresh) break
       
      if(trace==T){ #print output of each iteration
        vif_vals <- do.call('rbind', vif_vals)
        vif_vals
        prmatrix(vif_vals,collab='vif',rowlab=row.names(vif_vals),quote=F)
        cat('\n')
        cat('removed: ', names(vif_max),unlist(vif_max),'\n\n')
        flush.console()
        }
  
      in_dat<-in_dat[,!names(in_dat) %in% names(vif_max)]
  
      }
  
    return(names(in_dat))
     
    }
   
  }