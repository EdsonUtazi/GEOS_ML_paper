
#Load libraries
library(INLA)
#INLA:::inla.dynload.workaround()
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(ggplot2)
library(tictoc)

#Start timer
tic()
set.seed(500) #set seed - doesn't work for INLA functions

#Specify file paths
drive_path  <- "//patha/"
data_path   <- paste0(drive_path, "path_b/")
output_path <- paste0(drive_path, "path_c/")

#Load the processed input data
vaxdata <- read.csv(paste0(filePathData,"NG_Vax_measles.csv"), header=TRUE) #Vaccination data
vaxcov  <- read.csv(paste0(data_path,"Covariates_merged_selected_transformed.csv"), header = TRUE) #Covariate data

#Merge both data sets
data.merge <- merge(vaxcov, vaxdata, by = "DHSCLUST")
data.merge <- na.omit(data.merge) #Remove rows with missing data

#Split the merged data
vaxcov  <- data.merge[,c(3,6:18,4,5)] 
vaxdata <- data.merge[,20:31]

#Convert Urban-rural to 0's and 1's to conform to raster layer
Urban_rural <- as.factor(vaxcov$URBAN_RURA)
Urban_rural <- as.numeric(Urban_rural)
Urban_rural[Urban_rural==1] <- 0
Urban_rural[Urban_rural==2] <- 1
vaxcov$URBAN_RURA <- Urban_rural

#Delete clusters where TotChild is <=1 
zero.clust <- which(is.na(vaxdata$Age9_35_Tot)|vaxdata$Age9_35_Tot<=1)
if (length(zero.clust)>0){
  vaxdata <- vaxdata[-zero.clust,]
  vaxcov  <- vaxcov[-zero.clust,]
}

#Numbers vaccinated and numbers sampled
Numvacc    <- vaxdata$Age9_35_Vax
weights    <- vaxdata$Age9_35_Tot

#Observation coordinates
coords    <- cbind(vaxcov$LONGNUM,vaxcov$LATNUM)

#In-sample covariates
xp1     <- vaxcov$URBAN_RURA
xp2     <- vaxcov$l_n2_Mean
xp3	    <- vaxcov$n26_Mean
xp4a	    <- vaxcov$l_n28_Mean  
xp4b <- exp(xp4a) - 0.05
xp4  <- log(xp4b + 2.05)
xp5     <- vaxcov$media.prop
xp6	    <- vaxcov$n12_Mean       
xp7	    <- vaxcov$l_n51_Mean
xp8     <- vaxcov$l_n58_Mean
xp9     <- vaxcov$n69_Mean
xp10     <- vaxcov$logit_wealth.prop
xp11     <- vaxcov$logit_educ.prop
xp12     <- vaxcov$health_card_doc.prop
xp13     <- vaxcov$logit_sba.prop
xp14     <- vaxcov$logit_phone_internet.prop


#Prediction covariates
URBAN_RURA  	<- raster(paste0(data_path,"NGA_urban_rural_1km.tif"))
l_n2_Mean  	<- raster(paste0(data_path,"n2_cattle_FAO_2010.tif")) 
n26_Mean 	<- raster(paste0(data_path,"n26_guf_ghsl_dst_2014.tif"))
l_n28_Mean  	<- raster(paste0(data_path,"n28_viirs_2016_nightlights.tif"))
media.prop  	<- raster(paste0(data_path,"mediakrig.tif"))
n12_Mean	<- raster(paste0(data_path,"n12_ccilc_dst011_2015.tif"))
l_n51_Mean	<- raster(paste0(data_path,"n51_MODIS_lst.tif"))
l_n58_Mean	<- raster(paste0(data_path,"n58_TravelTime_to_HF_ver2.tif"))
n69_Mean	<- raster(paste0(data_path,"n69_WorldPop10povcons200.tif"))
logit_wealth.prop  	<- raster(paste0(data_path,"wealthkrig.tif"))
logit_educ.prop	<- raster(paste0(data_path,"educkrig.tif"))
health_card_doc	<- raster(paste0(data_path,"health_card_dockrig.tif"))
logit_sba.prop	<- raster(paste0(data_path,"sbakrig.tif"))
logit_phone_internet.prop	<- raster(paste0(data_path,"phone_internetkrig.tif"))


#Transform prediction covariates
x1gp  	<- getValues(URBAN_RURA)
x2gp 	<- log(getValues(l_n2_Mean) + 0.05)
x3gp 	<- getValues(n26_Mean)   			
x4gp	<- log(getValues(l_n28_Mean) + 2.05)
x5gp	<- getValues(media.prop)
x6gp	<- getValues(n12_Mean) 
x7gp	<- log(getValues(l_n51_Mean) + 0.05)
x8gp	<- log(getValues(l_n58_Mean) + 0.05) 
x9gp	<- getValues(n69_Mean) 
cc <- getValues(logit_wealth.prop); cc[cc==0] <- 0.01; cc[cc==1] <- 0.99
x10gp	<- logit(cc) 
cc <- getValues(logit_educ.prop); cc[cc==0] <- 0.01; cc[cc==1] <- 0.99
x11gp	<- logit(cc) 
x12gp	<- getValues(health_card_doc)
cc <- getValues(logit_sba.prop); cc[cc==0] <- 0.01; cc[cc==1] <- 0.99
x13gp	<- logit(cc)
cc <- getValues(logit_phone_internet.prop); cc[cc==0] <- 0.01; cc[cc==1] <- 0.99
x14gp	<- logit(cc)


#Extract prediction locations from a sample raster layer
n25.p      <- raster(paste0(data_path,"n69_WorldPop10povcons200.tif"))
Pred_grid2 <- coordinates(n25.p)

#Population data for population-weighted aggregation
popc     <- raster(paste0(data_path,"n105_WorldPop_under5s_2018.tif"))
popc	 <- getValues(popc) 

#Combine prediction grid, covariates and population data
pred.dat <- cbind(Pred_grid2, x1gp, x2gp, x3gp, x4gp, x5gp, x6gp, x7gp, x8gp, 
                  x9gp, x10gp, x11gp, x12gp, x13gp, x14gp, popc) 

#Determine grid cells with data and those with NAs  
ind <- apply(pred.dat, 1, function(x) any(is.na(x)))
miss    <- which(ind==TRUE)
nonmiss <- which(ind==FALSE)

#Prepare data to predict only where data are available
pred.dat.1 <- pred.dat[nonmiss, ]
pop <- pred.dat.1[,17]              
coord.p <- pred.dat.1[,1:2]
ypred=npred=rep(NA, nrow(pred.dat.1))

#Adm0 Shapefile
shp_ng  <- sf::st_read(paste0(data_path,"/Shapefiles/","gadm36_NGA_0.shp"))

#meshfit: fine triangulated mesh
shp_df <- sf::st_coordinates(shp_ng)
shp.bnd <- cbind(shp_df[,1], shp_df[,2])
c.bnd <- as.matrix(shp.bnd)	
meshfit <- inla.mesh.2d(loc=coords,loc.domain=c.bnd, max.edge=c(0.05, 0.6),cutoff=0.1)			
#plot(meshfit); plot(shp_ng, add=TRUE); points(coords)

#For priors
nu <- 1 #Matern smoothness parameter, redundant here as it implies alpha=2
alpha <- 2

r0 <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction (i.e. 0.05*(ymax-ymin))

#Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0, 0.01), prior.sigma=c(3, 0.01)) 


#For LGA-level estimates
coord.p <- pred.dat.1[,1:2]
spol1   <- readShapePoly(paste0(data_path,"Detailed_Boundary_ADM2.shp"))
spol    = as(spol1, "SpatialPolygons")   #NOTE - no data in spol
sp.1    <- rep(NA, nrow(coord.p))
for(i in 1:length(spol)){
  sp.1[as.vector(which(!is.na(over(SpatialPoints(coord.p), spol[i]))))] <- i
}

#For State-level estimates
spol2   <- readShapePoly(paste0(data_path,"Detailed_Boundary_ADM1.shp"))
spol    = as(spol2, "SpatialPolygons")   #NOTE - no data in spol
sp.2    <- rep(NA, nrow(coord.p))
for(i in 1:length(spol)){
  sp.2[as.vector(which(!is.na(over(SpatialPoints(coord.p), spol[i]))))] <- i
}

#For Regional-level estimates
spol3   <- readShapePoly(paste0(data_path,"sdr_subnational_boundaries.shp"))
spol    = as(spol3, "SpatialPolygons")   #NOTE - no data in spol
sp.3    <- rep(NA, nrow(coord.p))
for(i in 1:length(spol)){
  sp.3[as.vector(which(!is.na(over(SpatialPoints(coord.p), spol[i]))))] <- i
}


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

# Prediction covariate data frame
xpred1 <- pred.dat.1[,3]; xpred2 <- pred.dat.1[,4]; xpred3 <- pred.dat.1[,5]; xpred4 <- pred.dat.1[,6]; xpred5 <- pred.dat.1[,7]; 
xpred6 <- pred.dat.1[,8]; xpred7 <- pred.dat.1[,9]; xpred8 <- pred.dat.1[,10]
xpred9 <- pred.dat.1[,11]; xpred10 <- pred.dat.1[,12]; xpred11 <- pred.dat.1[,13]
xpred12 <- pred.dat.1[,14]; xpred13 <- pred.dat.1[,15]; xpred14 <- pred.dat.1[,16]

X0 <- model.matrix(~ -1 + factor(xpred1) + xpred2 + xpred3 + xpred4 + xpred5 + xpred6 + xpred7 + xpred8 +
                     xpred9 + xpred10 + xpred11 + xpred12 + xpred13 + xpred14)

Xpred <-  as.data.frame(X0[,-which(colnames(X0)%in%c("factor(xpred1)0"))])
colnames(Xpred) <- c("x_1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
                     "x9", "x10", "x11", "x12", "x13", "x14") 


# Final mdoel stack - can include prediction data
stk.full <- inla.stack(stk.point)  


# Fit model
#Priors
hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000) 

#Model formula
formula  <- y ~ -1 + intercept + x_1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + f(s, model=spde) + f(rr, model="iid", hyper = hyper.prec) #

res <- inla(formula, data=inla.stack.data(stk.full), family="binomial", 
            Ntrials = stk.full$data$data$n,
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk.full), link=1),
            control.compute=list(dic=TRUE, config = TRUE, waic=TRUE))

#save model
#save(res, file=paste0(output_path,"inla_model_dtp1_12_23.rda"))
spde.result <- inla.spde2.result(inla=res,name="s",spde=spde)

#Parameter estimates
#Regression coefficients
coeff.reg <- summary(res)$fixed[,1:5]

#Range for spatial RE
range.mean = inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]]); #range.mean
range.ex.sq = inla.emarginal(function(x) x^2, spde.result$marginals.range.nominal[[1]])
range.sd = sqrt(range.ex.sq-(range.mean^2)); #range.sd
range.quant = inla.qmarginal(c(0.025,0.5,0.975), spde.result$marginals.range.nominal[[1]]);# range.quant 
range <- c(range.mean, range.sd, range.quant)

#Variance for spatial RE
variance.mean = inla.emarginal(function(x) x, spde.result$marginals.variance.nominal[[1]]); #variance.mean
variance.ex.sq = inla.emarginal(function(x) x^2, spde.result$marginals.variance.nominal[[1]])
variance.sd = sqrt(variance.ex.sq-(variance.mean^2)); #variance.sd
variance.quant = inla.qmarginal(c(0.025,0.5,0.975), spde.result$marginals.variance.nominal[[1]]); #variance.quant 
variance <- c(variance.mean, variance.sd, variance.quant)

#Variance of IID random effect
var.ind      <- inla.tmarginal(function(x) 1/x, res$marginals.hyperpar[[3]])
var.iid      <- inla.zmarginal(var.ind,silent=TRUE)
variance.iid <- c(var.iid$mean, var.iid$sd, var.iid$quant0.025, var.iid$quant0.5, var.iid$quant0.975)

param.all <- rbind(coeff.reg,range,variance, variance.iid)
write.csv(param.all, paste0(output_path, "parameter_output_mcv1_9_35.csv"))


#In-sample prediction
index.pred.obs 	<- inla.stack.index(stk.full, tag = "point")$data
fitted.pred.all.obs = round(res$summary.fitted.values[index.pred.obs,1:5], 4)
fitted.pred.mean.obs1 = as.vector(data.matrix(as.vector(fitted.pred.all.obs[,"mean"])))
fitted.pred.sd.obs1 = as.vector(data.matrix(as.vector(fitted.pred.all.obs[,"sd"])))
fitted.pred.low.obs1 = as.vector(data.matrix(as.vector(fitted.pred.all.obs[,"0.025quant"])))
fitted.pred.up.obs1 = as.vector(data.matrix(as.vector(fitted.pred.all.obs[,"0.975quant"])))

prob.obs <- Numvacc/weights
ds <- data.frame(pred.prob=fitted.pred.mean.obs1, pred.obs=prob.obs, Vax=Numvacc, Tot=weights, 
                 pred.sd = fitted.pred.sd.obs1, pred.low=fitted.pred.low.obs1, pred.up=fitted.pred.up.obs1)
write.csv(ds, paste0(output_path, "obsvpred_mcv1_9_35.csv"))


##################POSTERIOR SAMPLING for out-of-sample prediction
nsamp <- 1000
ps <- inla.posterior.sample(nsamp, res) 
contents <- res$misc$configs$contents

ttp <- nrow(pred.dat.1)
rm(stk.point, stk.full)
rm(res)
rm(URBAN_RURA, l_n2_Mean, n26_Mean, l_n28_Mean, media.prop, n12_Mean, l_n51_Mean,	
   l_n58_Mean, n69_Mean, logit_wealth.prop, logit_educ.prop, health_card_doc, logit_sba.prop,	
   logit_phone_internet.prop)
rm(pred.dat)


#ID for spatial random effect
idSpace <- contents$start[which(contents$tag=="s")]-1 +
  (1:contents$length[which(contents$tag=="s")])

#ID for iid effects
idR <- contents$start[which(contents$tag=="rr")]-1 +
  (1:contents$length[which(contents$tag=="rr")])

#ID for fixed effects
idX <- contents$start[which(contents$tag=="intercept")]-1 + (1:15) # fixed effects, 15 = no of regression coefficients

#Extract samples 
xLatent <- matrix(0, nrow=length(ps[[1]]$latent), ncol=nsamp) 
xHyper <- matrix(0, nrow=length(ps[[1]]$hyperpar), ncol=nsamp) 
for(i in 1:nsamp){
  xLatent[,i] <- ps[[i]]$latent
  xHyper[,i] <- ps[[i]]$hyperpar
}
xSpace <- xLatent[idSpace,]
XR <- xLatent[idR,]
xX <- xLatent[idX,]


#Draw samples for IID term
sample.IIDpred <- matrix(0, nrow(pred.dat.1),  nsamp)
for (i in 1:nsamp){
  ID.precision <- xHyper[3,i]                         #the 3rd row contains precision for rr; same as ps[[i]]$hyperpar[3]
  ID.sigma <- ID.precision^-0.5
  sample.IIDpred[, i] <- rnorm(nrow(pred.dat.1), sd=ID.sigma)
}


#Prediction
pred.grid <- numeric()
inv.linpred.out <- numeric()
nsplits <- 10 #Split prediction data into subsets
splitsize <- floor(nrow(pred.dat.1)/nsplits) 
sind <- 1:nrow(pred.dat.1)

for (kk in 1:nsplits){
  print(kk)
  if (kk < nsplits) {qq <- (((kk-1)*splitsize)+1):(splitsize*kk)}
  if (kk == nsplits) {qq <- (((kk-1)*splitsize)+1):nrow(pred.dat.1)}
  
  A_pred = inla.spde.make.A(mesh = meshfit, loc = coord.p[qq,])
  
  linpred     <- as.matrix(A_pred %*% xSpace + as.matrix(cbind(1, Xpred$x_1[qq], Xpred$x2[qq], Xpred$x3[qq], Xpred$x4[qq], 
                                                               Xpred$x5[qq], Xpred$x6[qq], Xpred$x7[qq], Xpred$x8[qq], Xpred$x9[qq],
                                                               Xpred$x10[qq], Xpred$x11[qq], Xpred$x12[qq], Xpred$x13[qq], Xpred$x14[qq])) %*% xX + sample.IIDpred[qq,])  
  inv.linpred <- inv.logit(linpred) 
  inv.linpred.out <- rbind(inv.linpred.out, inv.linpred)
  
  pp.grid    <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975)))}))) 
  pred.grid <- rbind(pred.grid, pp.grid)
}#close of loop for splits

colnames(pred.grid) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
fitted.pred.mean   <- as.vector(data.matrix(as.vector(pred.grid[,"mean"])))
fitted.pred.sd     <- as.vector(data.matrix(as.vector(pred.grid[,"sd"])))
fitted.pred.median <- as.vector(data.matrix(as.vector(pred.grid[,"0.5quant"])))
fitted.pred.low    <- as.vector(data.matrix(as.vector(pred.grid[,"0.025quant"])))
fitted.pred.up     <- as.vector(data.matrix(as.vector(pred.grid[,"0.975quant"])))

#Sample raster layer used to form output raster layers
n25.p <- raster(paste0(data_path,"n69_WorldPop10povcons200.tif"))

#Mean
ll=1:length(ind); ll[nonmiss] = fitted.pred.mean; ll[miss] = NA
rr.mean = raster(n25.p); values(rr.mean) = ll

#sd
ll=1:length(ind); ll[nonmiss] = fitted.pred.sd; ll[miss] = NA
rr.sd = raster(n25.p); values(rr.sd) = ll

#low
ll=1:length(ind); ll[nonmiss] = fitted.pred.low; ll[miss] = NA
rr.low = raster(n25.p); values(rr.low) = ll

#up
ll=1:length(ind); ll[nonmiss] = fitted.pred.up; ll[miss] = NA
rr.up = raster(n25.p); values(rr.up) = ll

#median
ll=1:length(ind); ll[nonmiss] = fitted.pred.median; ll[miss] = NA
rr.med = raster(n25.p); values(rr.med) = ll

#Write output raster layers
writeRaster(rr.mean, paste0(output_path, "inla_vax_mcv1_9_35_mean.tif"), overwrite=TRUE)
writeRaster(rr.sd,   paste0(output_path, "inla_vax_mcv1_9_35_sd.tif"), overwrite=TRUE)
writeRaster(rr.low,  paste0(output_path, "inla_vax_mcv1_9_35_low.tif"), overwrite=TRUE)
writeRaster(rr.up,   paste0(output_path, "inla_vax_mcv1_9_35_up.tif"), overwrite=TRUE)
writeRaster(rr.med,  paste0(output_path, "inla_vax_mcv1_9_35_median.tif"), overwrite=TRUE)


#Calculate weighted population LGA, State and Regional estimates
inv.linpred <- inv.linpred.out
rm(inv.linpred.out)


#LGA estimates and uncertainty (sd) 
dd    <- 1:nrow(spol1)
dd.un <- unique(sp.1)
dmiss <- which(!dd%in%dd.un)

if (length(dmiss)>0) dd_num <- dd[-dmiss]
if (length(dmiss)==0) dd_num <- dd

dist_out <- matrix(0, length(dd_num), 5)
for (i in 1:length(dd_num)){
  if (length(which(sp.1==dd_num[i]))==1){ 
    pop.ext <- pop[which(sp.1==dd_num[i])] 
    ext <- as.vector(sapply(inv.linpred[which(sp.1==dd_num[i]),], FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE))) 
  }
  if (length(which(sp.1==dd_num[i]))>1){  
    pop.ext <- pop[which(sp.1==dd_num[i])]
    ext <- as.vector(apply(inv.linpred[which(sp.1==dd_num[i]),], 2, FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE)))
  }
  
  dist_out[i,] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs=c(0.025,0.5,0.975))))						
}

dist_out <- cbind(dd_num, dist_out)
colnames(dist_out) <- c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")

#The district-level estimates will have the same ordering as in the shapefile if they have the same no of areas
write.csv(dist_out, paste0(output_path, "District_estimates_mcv1_9_35.csv"))


#State estimates and uncertainty
dd    <- 1:nrow(spol2)
dd.un <- unique(sp.2)
dmiss <- which(!dd%in%dd.un)

if (length(dmiss)>0) dd_num <- dd[-dmiss]
if (length(dmiss)==0) dd_num <- dd

dist_out <- matrix(0, length(dd_num), 5)
for (i in 1:length(dd_num)){
  if (length(which(sp.2==dd_num[i]))==1){ 
    pop.ext <- pop[which(sp.2==dd_num[i])] 
    ext <- as.vector(sapply(inv.linpred[which(sp.2==dd_num[i]),], FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE))) 
  }
  if (length(which(sp.2==dd_num[i]))>1){  
    pop.ext <- pop[which(sp.2==dd_num[i])]
    ext <- as.vector(apply(inv.linpred[which(sp.2==dd_num[i]),], 2, FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE)))
  }
  
  dist_out[i,] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs=c(0.025,0.5,0.975))))						
}

dist_out <- cbind(dd_num, dist_out)
colnames(dist_out) <- c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")

write.csv(dist_out, paste0(output_path, "State_estimates_mcv1_9_35.csv"))


#Regional estimates and uncertainty
dd    <- 1:nrow(spol3)
dd.un <- unique(sp.3)
dmiss <- which(!dd%in%dd.un)

if (length(dmiss)>0) dd_num <- dd[-dmiss]
if (length(dmiss)==0) dd_num <- dd

dist_out <- matrix(0, length(dd_num), 5)
for (i in 1:length(dd_num)){
  if (length(which(sp.3==dd_num[i]))==1){ 
    pop.ext <- pop[which(sp.3==dd_num[i])] 
    ext <- as.vector(sapply(inv.linpred[which(sp.3==dd_num[i]),], FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE))) 
  }
  if (length(which(sp.3==dd_num[i]))>1){  
    pop.ext <- pop[which(sp.3==dd_num[i])]
    ext <- as.vector(apply(inv.linpred[which(sp.3==dd_num[i]),], 2, FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE)))
  }
  
  dist_out[i,] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs=c(0.025,0.5,0.975))))						
}

dist_out <- cbind(dd_num, dist_out)
colnames(dist_out) <- c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")

write.csv(dist_out, paste0(output_path, "Region_estimates_mcv1_9_35.csv"))


#Threshold probability calculations - 95%
ff1=function(x) length(which(x>=0.95))/nsamp
y.95 <- apply(inv.linpred, 1, ff1)   #Check me
ll=1:length(ind); ll[nonmiss] = y.95; ll[miss] = NA
rr.95 = raster(n25.p); values(rr.95) = ll
writeRaster(rr.95, paste0(output_path, "inla_vax_mcv1_9_35_95perc_thresh.tif"), overwrite=TRUE)


toc() #End timer



