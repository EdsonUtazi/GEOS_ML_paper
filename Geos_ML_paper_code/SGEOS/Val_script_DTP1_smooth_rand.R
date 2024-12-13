#Implemented on a HPC but can be run locally

#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(rgdal)
library(ggplot2)

filePathData <- "/.../"
filePathData1 <- "/.../Shapefiles/"
filePathData2 <- "/.../DTP1/"

set.seed(500) #set seed - doesn't work for INLA functions

#Source CRPS code
source(paste0(filePathData, "CRPS_func.R"))

#Load the processed input data
vaxdata <- read.csv(paste0(filePathData,"NG_Vax_DTP1.csv"), header=TRUE) #Vaccination data
vaxcov  <- read.csv(paste0(filePathData,"Covariates_merged_selected_transformed.csv"), header = TRUE) #Covariate data

#Merge both data sets
data.merge <- merge(vaxcov, vaxdata, by = "DHSCLUST")
data.merge <- na.omit(data.merge)

#Split the merged data 
vaxcov.all  <- data.merge[,c(3,6:18,4,5)]   
vaxdata.all <- data.merge[,20:31]

#Convert Urban-rural to 0's and 1's 
Urban_rural <- as.factor(vaxcov.all$URBAN_RURA)
Urban_rural <- as.numeric(Urban_rural)
Urban_rural[Urban_rural==1] <- 0
Urban_rural[Urban_rural==2] <- 1
vaxcov.all$URBAN_RURA <- Urban_rural

#Delete clusters where TotChild is <=1 
zero.clust <- which(is.na(vaxdata.all$Age12_23_Tot)|vaxdata.all$Age12_23_Tot<=1)
if (length(zero.clust)>0){
  vaxdata.all <- vaxdata.all[-zero.clust,]
  vaxcov.all  <- vaxcov.all[-zero.clust,]
}

#Numbers vaccinated and numbers sampled
Numvacc.all    <- vaxdata.all$Age12_23_Vax
weights.all    <- vaxdata.all$Age12_23_Tot

#Observation coordinates
coords.all    <- cbind(vaxcov.all$LONGNUM,vaxcov.all$LATNUM)

#Define logit and inverse logit function
inverse_logit <- function(aa){
  out = exp(aa)/(1+exp(aa))
  return(out)
}

logit <- function(aa){
  out = log(aa/(1-aa))
  return(out)
}

#Back-transforming the covariates
  vaxcov.all$l_n2_Mean     <- exp(vaxcov.all$l_n2_Mean) - 0.05
  vaxcov.all$l_n28_Mean     <- exp(vaxcov.all$l_n28_Mean) - 0.2     
  vaxcov.all$l_n51_Mean	  <- exp(vaxcov.all$l_n51_Mean) - 0.05 #_gp
  vaxcov.all$l_n58_Mean     <- exp(vaxcov.all$l_n58_Mean) - 0.05 #_gp
  vaxcov.all$logit_wealth.prop     <- inverse_logit(vaxcov.all$logit_wealth.prop)   #_gp
  vaxcov.all$logit_wealth.prop[vaxcov.all$logit_wealth.prop==0.99] <- 1; vaxcov.all$logit_wealth.prop[vaxcov.all$logit_wealth.prop==0.01] <- 0
  vaxcov.all$logit_educ.prop     <- inverse_logit(vaxcov.all$logit_educ.prop) #_gp
  vaxcov.all$logit_educ.prop[vaxcov.all$logit_educ.prop==0.99] <- 1; vaxcov.all$logit_educ.prop[vaxcov.all$logit_educ.prop==0.01] <- 0
  vaxcov.all$logit_sba.prop     <- inverse_logit(vaxcov.all$logit_sba.prop) #_gp
  vaxcov.all$logit_sba.prop[vaxcov.all$logit_sba.prop==0.99] <- 1; vaxcov.all$logit_sba.prop[vaxcov.all$logit_sba.prop==0.01] <- 0
  vaxcov.all$logit_phone_internet.prop     <- inverse_logit(vaxcov.all$logit_phone_internet.prop) #_gp
  vaxcov.all$logit_phone_internet.prop[vaxcov.all$logit_phone_internet.prop==0.99] <- 1; vaxcov.all$logit_phone_internet.prop[vaxcov.all$logit_phone_internet.prop==0.01] <- 0


#Grouping the covariates into bins for the smooth functions
gg <- 10 #Number of groupd tuned to avoid numerical errors 
vaxcov.all$n12_Mean_gp <- inla.group(vaxcov.all$n12_Mean, n = 5, method = "quantile")# , method = "quantile")
vaxcov.all$l_n51_Mean_gp <- inla.group(vaxcov.all$l_n51_Mean, n = 5, method = "quantile") #10 chosen because range of values is small
vaxcov.all$l_n28_Mean_gp <- inla.group(vaxcov.all$l_n28_Mean, n = 5, method = "quantile")
vaxcov.all$n69_Mean_gp <- inla.group(vaxcov.all$n69_Mean, n = gg, method = "quantile")
vaxcov.all$logit_wealth.prop_gp <- inla.group(vaxcov.all$logit_wealth.prop, n = 10, method = "quantile")
vaxcov.all$logit_educ.prop_gp <- inla.group(vaxcov.all$logit_educ.prop, n = 10, method = "quantile")
vaxcov.all$health_card_doc.prop_gp <- inla.group(vaxcov.all$health_card_doc.prop, n = gg, method = "quantile")
vaxcov.all$logit_sba.prop_gp <- inla.group(vaxcov.all$logit_sba.prop, n = 10, method = "quantile")
vaxcov.all$l_n2_Mean_gp <- inla.group(vaxcov.all$l_n2_Mean, n = 10, method = "quantile")


#######################################################
#Start k-fold cross-validation loop
cv <- 10   #i.e. 10-fold cross validation
lim <- floor(nrow(coords.all)/cv)

#Create a random arrangement of the observations
srand <- sample(1:nrow(coords.all),nrow(coords.all), replace=FALSE)

val.out <- matrix(0, cv, 7) #Output matrix for val metrics
obspred.out <- data.frame() #To store CV predictions

for (kk in 1:cv){
  if (kk < cv) {qq <- (((kk-1)*lim)+1):(lim*kk); samp.c <- srand[qq]}
  if (kk == cv) {qq <- (((kk-1)*lim)+1):nrow(coords.all); samp.c <- srand[qq]}
  
  #Validation data
  coords.nc 	<- coords.all[samp.c,]
  Numvacc.nc	<- Numvacc.all[samp.c]
  weights.nc	<- weights.all[samp.c]
  vaxcov.nc	<- vaxcov.all[samp.c,]
  yp.nc=np.nc=rep(NA, length(Numvacc.nc))
  
  #Training data
  coords  <- coords.all[-samp.c,]
  Numvacc <- Numvacc.all[-samp.c] 
  weights <- weights.all[-samp.c]
  vaxcov  <- vaxcov.all[-samp.c,]
  
  #Covariates for model-fitting
  xp1     <- vaxcov$URBAN_RURA
  xp2     <- vaxcov$logit_phone_internet.prop
  xp3	  <- vaxcov$n26_Mean
  xp4     <- vaxcov$l_n58_Mean
  xp5     <- vaxcov$media.prop
  xp6	  <- vaxcov$n12_Mean_gp
  xp7	  <- vaxcov$l_n51_Mean_gp
  xp8     <- vaxcov$l_n28_Mean_gp
  xp9     <- vaxcov$n69_Mean_gp
  xp10     <- vaxcov$logit_wealth.prop_gp
  xp11     <- vaxcov$logit_educ.prop_gp
  xp12     <- vaxcov$health_card_doc.prop_gp
  xp13     <- vaxcov$logit_sba.prop_gp
  xp14     <- vaxcov$l_n2_Mean_gp
  
  #Covariates for validation
  xv1     <- vaxcov.nc$URBAN_RURA
  xv2    <- vaxcov.nc$logit_phone_internet.prop
  xv3	  <- vaxcov.nc$n26_Mean
  xv4     <- vaxcov.nc$l_n58_Mean
  xv5     <- vaxcov.nc$media.prop
  xv6	  <- vaxcov.nc$n12_Mean_gp
  xv7	  <- vaxcov.nc$l_n51_Mean_gp
  xv8     <- vaxcov.nc$l_n28_Mean_gp
  xv9     <- vaxcov.nc$n69_Mean_gp
  xv10     <- vaxcov.nc$logit_wealth.prop_gp
  xv11     <- vaxcov.nc$logit_educ.prop_gp
  xv12     <- vaxcov.nc$health_card_doc.prop_gp
  xv13     <- vaxcov.nc$logit_sba.prop_gp
  xv14     <- vaxcov.nc$l_n2_Mean_gp

  #Adm0 shapefile
  shp_ng  <- readShapePoly(paste0(filePathData1,"gadm36_NGA_0.shp"))
  
  #meshfit: fine triangulated mesh
  shp_df <- fortify(shp_ng)
  shp.bnd <- cbind(shp_df$long, shp_df$lat)
  c.bnd <- as.matrix(shp.bnd)	
  meshfit <- inla.mesh.2d(loc=coords,loc.domain=c.bnd, max.edge=c(0.05, 0.6),cutoff=0.1)		
  
  #For priors
  nu <- 1 #Matern smoothness parameter, redundant here as it implies alpha=2
  alpha <- 2
  
  r0 <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction (i.e. 0.05*(ymax-ymin))
  
  #Matern SPDE model object using inla.pcmatern
  spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0, 0.01), prior.sigma=c(3, 0.01))
  
  # Training data - covariates and stack object
  X0 <- model.matrix(~ -1 + factor(xp1) + xp2 + xp3 + xp4 + xp5 + xp6 + xp7 + xp8 +
                       xp9 + xp10 + xp11 + xp12 + xp13 + xp14)
  Xobs1 <-  as.data.frame(X0[,-which(colnames(X0)%in%c("factor(xp1)0"))])
  colnames(Xobs1) <- c("x_1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
                       "x9", "x10", "x11", "x12", "x13", "x14")
  Ap.i <- inla.spde.make.A(mesh=meshfit, loc=coords)
  lp.i = rep(1,length(xp1))
  stk.point <- inla.stack(tag='point',
                          data=list(y=Numvacc,n=weights),
                          A=list(Ap.i,1,1,1),
                          effects=list(s=1:spde$n.spde, rr=1:length(weights), intercept=rep(1, length(xp1)), Xobs1))  
  
  # Validation data - covariates and stack object
  X0 <- model.matrix(~ -1 + factor(xv1) + xv2 + xv3 + xv4 + xv5 + xv6 + xv7 + xv8 +
                       xv9 + xv10 + xv11 + xv12 + xv13 + xv14)
  Xobs <-  as.data.frame(X0[,-which(colnames(X0)%in%c("factor(xv1)0"))])
  colnames(Xobs) <- c("x_1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
                      "x9", "x10", "x11", "x12", "x13", "x14")
  Aval <- inla.spde.make.A(mesh=meshfit, loc=coords.nc)
  lval = rep(1,length(xv1))
  stk.val <- inla.stack(tag='val',
                        data=list(y=yp.nc,n=np.nc),
                        A=list(Aval,1,1,1),
                        effects=list(s=1:spde$n.spde, rr=(length(weights)+1):(length(weights)+length(xv1)), intercept=rep(1, length(xv1)), Xobs)) 
  
  #Final model stack - can include prediction data
  stk.full <- inla.stack(stk.point)  
  
  # Fit model
  # Priors
  hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
  control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  
  
  #For smooth functions 
  ux6 <- sort(unique(vaxcov.all$n12_Mean_gp)); ux7 <- sort(unique(vaxcov.all$l_n51_Mean_gp))
  ux8 <- sort(unique(vaxcov.all$l_n28_Mean_gp)); ux9 <- sort(unique(vaxcov.all$n69_Mean_gp))
  ux10 <- sort(unique(vaxcov.all$logit_wealth.prop_gp)); ux11 <- sort(unique(vaxcov.all$logit_educ.prop_gp))
  ux12 <- sort(unique(vaxcov.all$health_card_doc.prop_gp)); ux13 <- sort(unique(vaxcov.all$logit_sba.prop_gp))
  ux14 <- sort(unique(vaxcov.all$l_n2_Mean_gp))
  
  #Model formula
  #Scaling the RW2 model is highly recommended
  #Using "values" in RW functions ensures that these are computed for all relevant values of the covariate
  #This is especially important when making predictions using posterior samples.
  formula  <- y ~ -1 + intercept + x_1 + x2 + x3 + x4 + x5 +
    f(x6, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux6) +
    f(x7, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux7) +
    f(x8, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux8) +
    f(x9, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux9) +
    f(x10, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux10) +
    f(x11, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux11) +
    f(x12, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux12) +
    f(x13, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux13) +
    f(x14, model = "rw2", constr=TRUE, scale.model=TRUE, values=ux14) +
    f(s, model=spde) + f(rr, model="iid", hyper = hyper.prec) 
  
  res <- inla(formula, data=inla.stack.data(stk.full), family="binomial", 
              Ntrials = stk.full$data$data$n,
              control.predictor=list(compute=TRUE, A=inla.stack.A(stk.full), link=1),
              control.compute=list(dic=TRUE, config = TRUE, waic=TRUE),
              control.fixed=control.fixed)
  
  spde.result <- inla.spde2.result(inla=res,name="s",spde=spde)
  
  
  ##################POSTERIOR SAMPLING
  nsamp <- 1000
  
  #Posterior sampling
  ps <- inla.posterior.sample(nsamp, res) 
  contents <- res$misc$configs$contents
  
  #ID for spatial random effect
  idSpace <- contents$start[which(contents$tag=="s")]-1 +
    (1:contents$length[which(contents$tag=="s")])
  
  #ID for iid effects
  idR <- contents$start[which(contents$tag=="rr")]-1 +
    (1:contents$length[which(contents$tag=="rr")])
  
  #ID for fixed effects
  idX <- contents$start[which(contents$tag=="intercept")]-1 + (1:6) # fixed effects, 6 = no of (fixed) regression coefficients
  
  #ID for x3
  #idx3 <- contents$start[which(contents$tag=="x3")]-1 +
   # (1:contents$length[which(contents$tag=="x3")])
  
  #ID for x4
  #idx4 <- contents$start[which(contents$tag=="x4")]-1 +
   # (1:contents$length[which(contents$tag=="x4")])
  
  #ID for x5
  #idx5 <- contents$start[which(contents$tag=="x5")]-1 +
   # (1:contents$length[which(contents$tag=="x5")])
  
  #ID for x6
  idx6 <- contents$start[which(contents$tag=="x6")]-1 +
    (1:contents$length[which(contents$tag=="x6")])
  
  #ID for x7
  idx7 <- contents$start[which(contents$tag=="x7")]-1 +
    (1:contents$length[which(contents$tag=="x7")])
  
  #ID for x8
  idx8 <- contents$start[which(contents$tag=="x8")]-1 +
    (1:contents$length[which(contents$tag=="x8")])
  
  
  #ID for x9
  idx9 <- contents$start[which(contents$tag=="x9")]-1 +
    (1:contents$length[which(contents$tag=="x9")])
  
  #ID for x10
  idx10 <- contents$start[which(contents$tag=="x10")]-1 +
    (1:contents$length[which(contents$tag=="x10")])
  
  #ID for x11
  idx11 <- contents$start[which(contents$tag=="x11")]-1 +
    (1:contents$length[which(contents$tag=="x11")])
  
  #ID for x12
  idx12 <- contents$start[which(contents$tag=="x12")]-1 +
    (1:contents$length[which(contents$tag=="x12")])
  
  #ID for x13
  idx13 <- contents$start[which(contents$tag=="x13")]-1 +
    (1:contents$length[which(contents$tag=="x13")])
  
  #ID for x14
  idx14 <- contents$start[which(contents$tag=="x14")]-1 +
    (1:contents$length[which(contents$tag=="x14")])
  
  # extract samples 
  xLatent <- matrix(0, nrow=length(ps[[1]]$latent), ncol=nsamp) 
  xHyper <- matrix(0, nrow=length(ps[[1]]$hyperpar), ncol=nsamp) 
  for(i in 1:nsamp){
    xLatent[,i] <- ps[[i]]$latent
    xHyper[,i] <- ps[[i]]$hyperpar
  }
  xSpace <- xLatent[idSpace,]
  XR <- xLatent[idR,]
  xX <- xLatent[idX,]
  #Xx3 <- xLatent[idx3,]
  #Xx4 <- xLatent[idx4,]
  #Xx5 <- xLatent[idx5,]
  Xx6 <- xLatent[idx6,]
  Xx7 <- xLatent[idx7,]
  Xx8 <- xLatent[idx8,]
  Xx9 <- xLatent[idx9,]
  Xx10 <- xLatent[idx10,]
  Xx11 <- xLatent[idx11,]
  Xx12 <- xLatent[idx12,]
  Xx13 <- xLatent[idx13,]
  Xx14 <- xLatent[idx14,]
  

  
  #Further processing for smooth functions
  #----------------------------------------------#
  ID6 <- ux6; ID7 <- ux7
  ID8 <- ux8; ID9 <- ux9; ID10 <- ux10
  ID11 <- ux11; ID12 <- ux12; ID13 <- ux13; ID14 <- ux14
  
  idx3 = idx4 = 0
  for (i in 1:nrow(vaxcov.nc)){
    #idx3[i] <- which(ID3==xv3[i]); idx4[i] <- which(ID4==xv4[i])
    #idx5[i] <- which(ID5==xv5[i]); 
    idx6[i] <- which(ID6==xv6[i])
    idx7[i] <- which(ID7==xv7[i]); idx8[i] <- which(ID8==xv8[i])
    idx9[i] <- which(ID9==xv9[i]); idx10[i] <- which(ID10==xv10[i])
    idx11[i] <- which(ID11==xv11[i]); idx12[i] <- which(ID12==xv12[i])
    idx13[i] <- which(ID13==xv13[i]); idx14[i] <- which(ID14==xv14[i])
  }
  fx6hat <- Xx6[idx6,]
  fx7hat <- Xx7[idx7,]; fx8hat <- Xx8[idx8,]
  fx9hat <- Xx9[idx9,]; fx10hat <- Xx10[idx10,]
  fx11hat <- Xx11[idx11,]; fx12hat <- Xx12[idx12,]
  fx13hat <- Xx13[idx13,]; fx14hat <- Xx14[idx14,]
  #----------------------------------------------#
  
  #Draw samples for IID term
  sample.IIDval <- matrix(0, length(xv1),  nsamp)
  for (i in 1:nsamp){
    ID.precision <- xHyper[3 + 9,i]   #9 is the number of smooth functions; 3 is (range for s, sd for s then prec for rr)
    ID.sigma <- ID.precision^-0.5
    sample.IIDval[, i] <- rnorm(length(xv1), sd=ID.sigma)
  }
  
  linpred <- as.matrix(Aval %*% xSpace + as.matrix(cbind(1, Xobs$x_1, Xobs$x2, Xobs$x3, Xobs$x4, Xobs$x5)) %*% xX + 
                         fx6hat + fx7hat + fx8hat +
                         fx9hat + fx10hat + fx11hat + fx12hat + fx13hat + fx14hat +
                         sample.IIDval) #Compute linear predictor
  inv.linpred <- inv.logit(linpred)    
  pred.obs <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ c(mean(x), sd(x), quantile(x, probs=c(0.025,0.1,0.5,0.9,0.975)))}))) 
  colnames(pred.obs) <- c("mean", "sd", "0.025", "0.1", "median", "0.9", "0.975")
  fitted.mean.val <- as.vector(data.matrix(as.vector(pred.obs[,"mean"])))
  fitted.sd.val   <- as.vector(data.matrix(as.vector(pred.obs[,"sd"])))
  fitted.low.val  <- as.vector(data.matrix(as.vector(pred.obs[,"0.025"])))
  fitted.up.val   <- as.vector(data.matrix(as.vector(pred.obs[,"0.975"])))
  fitted.low.val80  <- as.vector(data.matrix(as.vector(pred.obs[,"0.1"])))
  fitted.up.val80   <- as.vector(data.matrix(as.vector(pred.obs[,"0.9"])))

  #Calculate validation metrics
  prob.val <- Numvacc.nc/weights.nc
  corr <- cor(fitted.mean.val,prob.val)
  rsq.val  <- (cor(fitted.mean.val,prob.val))^2
  RMSE.val <- sum((fitted.mean.val-prob.val)^2)/length(prob.val)
  MAE <- sum(abs(fitted.mean.val-prob.val))/length(prob.val)
  perc_bias <- (sum(fitted.mean.val-prob.val)/sum(prob.val))*100
  avg_bias <- sum(fitted.mean.val-prob.val)/length(prob.val)
  #CRPS
  CRPS <- crps(inv.linpred, prob.val) 
  
  val.out[kk, ] <- c(corr, rsq.val, RMSE.val, MAE, perc_bias, avg_bias, CRPS)
  print(val.out[kk, ])

  qp <- data.frame(prob.val, fitted.mean.val, fitted.sd.val, fitted.low.val, fitted.up.val)
  obspred.out <- rbind(obspred.out, qp)
}

colnames(val.out) <- c("correl", "rsq.val", "RMSE.val", "MAE","perc_bias", "avg_bias", "CRPS")

#Write output files
write.csv(val.out, paste0(filePathData2, "val_out_dtp1_smooth_rand.csv"))
write.csv(obspred.out, paste0(filePathData2, "val_obsvpred_dtp1_smooth_rand.csv"))





