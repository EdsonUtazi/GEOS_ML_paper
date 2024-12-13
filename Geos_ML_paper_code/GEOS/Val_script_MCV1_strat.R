#Implemented on a HPC but can be run locally

#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(rgdal)
library(ggplot2)

#File paths
filePathData <- "/.../"
filePathData1 <- "/.../Shapefiles/"
filePathData2 <- "/.../MCV1/"

set.seed(500) #set seed - doesn't work for INLA functions

#Source CRPS code
source(paste0(filePathData, "CRPS_func.R"))

#Load the processed input data
vaxdata <- read.csv(paste0(filePathData,"NG_Vax_measles.csv"), header=TRUE) #Vaccination data
vaxcov  <- read.csv(paste0(filePathData,"Covariates_merged_selected_transformed.csv"), header = TRUE) #Covariate data

#Merge both data sets
data.merge <- merge(vaxcov, vaxdata, by = "DHSCLUST")
data.merge <- na.omit(data.merge)

#Split the merged data 
vaxcov.all  <- data.merge[,c(3,6:18,4,5)]   #includes coordinates as Longitude & Latitude
vaxdata.all <- data.merge[,20:31]

#Convert Urban-rural to 0's and 1's 
Urban_rural <- as.factor(vaxcov.all$URBAN_RURA)
Urban_rural <- as.numeric(Urban_rural)
Urban_rural[Urban_rural==1] <- 0
Urban_rural[Urban_rural==2] <- 1
vaxcov.all$URBAN_RURA <- Urban_rural

#Delete clusters where TotChild is <=1 
zero.clust <- which(is.na(vaxdata.all$Age9_35_Tot)|vaxdata.all$Age9_35_Tot<=1)
if (length(zero.clust)>0){
  vaxdata.all <- vaxdata.all[-zero.clust,]
  vaxcov.all  <- vaxcov.all[-zero.clust,]
}

#Numbers vaccinated and numbers sampled
Numvacc.all    <- vaxdata.all$Age9_35_Vax
weights.all    <- vaxdata.all$Age9_35_Tot

#Observation coordinates
coords.all    <- cbind(vaxcov.all$LONGNUM,vaxcov.all$LATNUM)


#######################################################
#Start k-fold cross-validation loop
cv <- 10   #i.e. 10-fold cross validation
lim <- floor(nrow(coords.all)/cv)

#Spatial stratification - cluster arrangement is spatial
strat <- 1:nrow(coords.all)

val.out <- matrix(0, cv, 7) #Output matrix for val metrics
obspred.out <- data.frame() #To store CV predictions

for (kk in 1:cv){
  if (kk < cv) {qq <- (((kk-1)*lim)+1):(lim*kk); samp.c <- strat[qq]}
  if (kk == cv) {qq <- (((kk-1)*lim)+1):nrow(coords.all); samp.c <- strat[qq]}
  
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
  xp2     <- vaxcov$n26_Mean
  xp3	    <- vaxcov$l_n2_Mean
  xp4	    <- vaxcov$l_n28_Mean
  xp5	    <- vaxcov$l_n51_Mean
  xp6     <- vaxcov$l_n58_Mean
  xp7     <- vaxcov$n69_Mean
  xp8     <- vaxcov$n12_Mean
  xp9     <- vaxcov$logit_wealth.prop
  xp10     <- vaxcov$logit_educ.prop
  xp11     <- vaxcov$health_card_doc.prop
  xp12     <- vaxcov$logit_sba.prop
  xp13     <- vaxcov$logit_phone_internet.prop
  xp14     <- vaxcov$media.prop
  
  #Covariates for validation
  xv1     <- vaxcov.nc$URBAN_RURA
  xv2     <- vaxcov.nc$n26_Mean
  xv3	    <- vaxcov.nc$l_n2_Mean
  xv4	    <- vaxcov.nc$l_n28_Mean
  xv5	    <- vaxcov.nc$l_n51_Mean
  xv6     <- vaxcov.nc$l_n58_Mean
  xv7     <- vaxcov.nc$n69_Mean
  xv8     <- vaxcov.nc$n12_Mean
  xv9     <- vaxcov.nc$logit_wealth.prop
  xv10     <- vaxcov.nc$logit_educ.prop
  xv11     <- vaxcov.nc$health_card_doc.prop
  xv12     <- vaxcov.nc$logit_sba.prop
  xv13     <- vaxcov.nc$logit_phone_internet.prop
  xv14     <- vaxcov.nc$media.prop
  
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
  
  # Observation data - covariates and stack object
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
                          effects=list(s=1:spde$n.spde, rr=1:length(weights), intercept=rep(1, length(xp1)), Xobs1))  #NOTE
  
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
  
  #Model formula
  formula  <- y ~ -1 + intercept + x_1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 +
    x9 + x10 + x11 + x12 + x13 + x14 + f(s, model=spde) + f(rr, model="iid", hyper = hyper.prec) 
  
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
  idX <- contents$start[which(contents$tag=="intercept")]-1 + (1:15) # fixed effects, 15 = no of (fixed) regression coefficients
  
  
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
  
  
  ####Prediction for validation locations
  #Draw samples for IID term
  sample.IIDval <- matrix(0, length(xv1),  nsamp)
  for (i in 1:nsamp){
    ID.precision <- xHyper[3,i]   #12 is the number of smooth functions; 3 is (range for s, sd for s then prec for rr)
    ID.sigma <- ID.precision^-0.5
    sample.IIDval[, i] <- rnorm(length(xv1), sd=ID.sigma)
  }
  
  linpred <- as.matrix(Aval %*% xSpace + as.matrix(cbind(1, Xobs$x_1, Xobs$x2, Xobs$x3, Xobs$x4, 
                                                         Xobs$x5, Xobs$x6, Xobs$x7, Xobs$x8, Xobs$x9, 
                                                         Xobs$x10, Xobs$x11, Xobs$x12, Xobs$x13, 
                                                         Xobs$x14)) %*% xX + sample.IIDval) #Compute linear predictor
  inv.linpred <- inv.logit(linpred)      #CHECK FUNCTION ON IRIDIS
  pred.obs <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975)))}))) 
  colnames(pred.obs) <- c("mean", "sd", "0.025", "median", "0.975")
  fitted.mean.val <- as.vector(data.matrix(as.vector(pred.obs[,"mean"])))
  fitted.sd.val   <- as.vector(data.matrix(as.vector(pred.obs[,"sd"])))
  fitted.low.val  <- as.vector(data.matrix(as.vector(pred.obs[,"0.025"])))
  fitted.up.val   <- as.vector(data.matrix(as.vector(pred.obs[,"0.975"])))
  
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

colnames(val.out) <- c("correl", "rsq.val", "RMSE.val", "MAE", "perc_bias", "avg_bias", "CRPS")

#Write output files
write.csv(val.out, paste0(filePathData2, "val_out_mcv1_strat.csv"))
write.csv(obspred.out, paste0(filePathData2, "val_obsvpred_mcv1_strat.csv"))
