
# A Stacked Generalization Modelling using Bayesian Geostatistical Model (STG)

library(INLA)
library(sf)
library(sp)
library(tidyverse)
library(terra)
library(tictoc)
library(inlabru)
library(groupdata2)

#Specify Drive Path
drive_path <- "//.../Working/"
input_path <- paste0(drive_path, "Methods_comparison_paper/.../MCV1/")
data_path <- paste0(drive_path, "Methods_comparison_paper/2018_NGA_DHS/")
output_path <- paste0(drive_path, "Methods_comparison_paper/.../MCV1/Predictions/")
shapefile_path <- paste0(drive_path, "Methods_comparison_paper/2018_NGA_DHS/Shapefiles/")

#Load data
gam_data <- read.csv(paste0(input_path, "Gam_Cubic_Predictions.csv"))
gbm_data <- read.csv(paste0(input_path, "GBM_Predictions.csv"))
lasso_data <- read.csv(paste0(input_path, "Lasso_Predictions.csv"))
mcv1 <- read.csv(paste0(input_path, "gam_data.csv"))
shp_ng  <- st_read(paste0(shapefile_path, "gadm36_NGA_0.shp"))


#Select variables and join
gam_data <- gam_data %>% 
  select(DHSCLUST, gam_predictions)

gbm_data <- gbm_data %>% 
  select(DHSCLUST, gbm_predictions)

lasso_data <- lasso_data %>% 
  select(DHSCLUST, lasso_predictions)

mcv1 <- mcv1 %>% 
  select(X, DHSCLUST, prop, DHSCLUST, LONGNUM, LATNUM, Age9_35_Vax, Age9_35_Tot)


#gam_gbm_join
gam_gbm_join <- full_join(gam_data, gbm_data, by = "DHSCLUST")

#join to lasso
lasso_join <- full_join(lasso_data, gam_gbm_join , by = "DHSCLUST")

#Join datasets
mcv1 <- full_join(mcv1, lasso_join, by = "DHSCLUST")


#Create a grouping variable sequentially to be used for spatial k-fold cross val
mcv1 <- mcv1 %>% 
  group(n = 133, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()


#Renames X as id
mcv1 <- mcv1 %>% 
  rename(id = X)

#Calculate mcv1 proportion
mcv1 <- mcv1 %>% 
  mutate(prop = Age9_35_Vax/Age9_35_Tot)

# Extract coordinates from your data frame
xy <- mcv1[, c("LONGNUM", "LATNUM")]

# Create the SpatialPointsDataFrame
dtp_sf <- SpatialPointsDataFrame(coords = xy, data = mcv1,
                                 proj4string = CRS("+proj=longlat +datum=WGS84")) 

#Inverse logit function for back transformation
inverse_logit <- function(aa) {
  
  out = exp(aa)/(1+exp(aa))
  return(out)
  
}


##############################################################################
# Geostatistical Model - INLABRU SPDE ----------------------------------------

#-Define the coordinates of centroids
coords <- cbind(dtp_sf$LONGNUM, dtp_sf$LATNUM) 

#measure distance between coordinates
#summary(dist(coords)) #summarizes the Euclidean distance between points in the spatial domain

#build mesh
mesh1 <- fm_mesh_2d_inla(loc=coords,boundary=shp_ng, max.edge=c(0.05, 0.6),cutoff=0.1)


plot(mesh1)
plot(mesh1, add=T)
points(dtp_sf, col="red", pch="*")

#Count of mesh nodes
mesh1$n

#spde parameters
alpha <- 2
r0 <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction (i.e. 0.05*(ymax-ymin))

#Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(mesh=mesh1, alpha=alpha, prior.range=c(r0, 0.01), prior.sigma=c(3, 0.01))


# Priors
hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  


# Define model components and constraint
cmp <- Age9_35_Vax ~ Intercept(1) + x(cbind(gam_predictions, 
                                            gbm_predictions, lasso_predictions), 
                                      model = "iid",
                                      mapper = bru_mapper_matrix(labels = c("gam_predictions", "gbm_predictions", "lasso_predictions")),
                                      hyper = list(prec = list(initial = -log(30), fixed = TRUE)), #can be initial = -log(30) or log(1/1000) #-2*log(30)
                                      extraconstr = list(A = cbind(1, 1, 1), e = 1)) +
  field(coordinates, model = spde) + f.iid(id, model="iid", hyper = hyper.prec)


#Fit model
res <- bru(components = cmp,  family = "binomial", 
           Ntrials = dtp_sf$Age9_35_Tot, data = dtp_sf,
           options = list(control.fixed=control.fixed))

summary(res)

#Get In-Sample Predictions
in_sample_predictions <- generate(res, newdata = dtp_sf, 
                                  formula = ~ Intercept + x + f.iid + 
                                    field_eval(cbind(LONGNUM, LATNUM)), n.samples = 1000)


#Get predictions as dataframe
in_sample_predictions <- data.frame(in_sample_predictions)

#back transform predictions
in_sample_predictions <- inverse_logit(in_sample_predictions) 

#summarize predictions
mean <- rowMeans(in_sample_predictions, na.rm = T)
std <- apply(in_sample_predictions, 1, sd)
lower_quantile <- apply(in_sample_predictions, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
upper_quantile <- apply(in_sample_predictions, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))

#cbind them as tibble
in_sample_predictions <- cbind(mean, std, lower_quantile, upper_quantile)


#Get prediction metrics
in_sample_predictions <- in_sample_predictions %>% 
  as_tibble() %>% 
  mutate(observed = dtp_sf$prop) %>% 
  rename(predicted = mean)


write.csv(in_sample_predictions, paste0(output_path, "STG_insample_pred_dtp1.csv"))

metrics <- in_sample_predictions %>% 
  mutate(residual = observed - predicted)%>% 
  
  summarise(
    Bias= mean(residual),
    Imprecision = sd(residual),
    mae = mean(abs(residual)),
    mse = mean((residual)^2),
    rmse = sqrt(mse),
    Corr = cor(observed, predicted))
metrics

####################################################################################################
################### Calculate CRPS ################################################################
#Generate samples from posteriors
smp <- generate(res, newdata = dtp_sf, 
                formula = ~ Intercept + x + f.iid + 
                  field_eval(cbind(LONGNUM, LATNUM)), n.samples = 1000)

#Cbind observed data to samples
smp_data <- cbind(dtp_sf$prop, smp) %>% 
  as_tibble()%>%
  drop_na()

#Number of iterations excluding the first column(which is the observed data)
its <- ncol(smp_data) - 1

#Function to calculate the crps
crps_fnc <- function(p, n) {
  out <- map_dbl(1:n, ~ abs(p[.x + 1] - p[1]))
  out1 <- map_dbl(1:n, ~ sum(abs(p[.x + 1] - p[-1])))
  
  out <- sum(out) / n
  out1 <- sum(out1) / (2 * n * n)
  out - out1
}

crps <- smp_data %>%
  rowwise() %>%
  mutate(crps_value = crps_fnc(c_across(everything()), its)) %>%
  ungroup() %>%
  summarise(mean_crps = mean(crps_value)) %>%
  pull(mean_crps)

crps

###########Additional CRPS Code ######################

tmp <- cbind(dtp_sf$prop, smp)
tmp<-na.omit(tmp)
its<-length(tmp[1,])-1

fnc<-function(p,n){
  out<-NULL
  out1<-NULL
  for(i in 1:n){
    out[i]<-abs(p[i+1]-p[1])
    #out1[i]<-abs(p[i+1]-sample(p[-c(1,i)],1))*abs(p[i+1]-sample(p[-c(1,i)],1))
    out1[i]<-sum(abs(p[i+1]-p[-1]))
  }
  out<-sum(out)/n
  out1<-sum(out1)/(2*n*n)
  out<-out-out1
  out
}
tmp<-apply(tmp,1,fnc,its)
tmp<-mean(tmp)
tmp

#Both CRPS functions give the same result
crps
tmp

########################################################################################################
########################### Kfold Cross Validation ####################################################

# K-Fold Cross Validation -------------------------------------------------

# function to calculate k-fold
kfold_cv <- function(data, k) {
  n <- nrow(data)
  fold_size <- n %/% k
  folds <- sample(rep(1:k, each = fold_size, length.out = n))
  
  #Train metrics placeholder
  
  # train_rmse_values <- numeric(k) # Placeholder for RMSE
  # train_pearson_values <- numeric(k) # Placeholder for corr
  # train_mae_values <- numeric(k)  # Placeholder for MAE
  # train_bias_values <- numeric(k)  # Placeholder for bias
  
  #Test metrics placeholder
  test_rmse_values <- numeric(k) # Placeholder for RMSE
  test_pearson_values <- numeric(k) # Placeholder for corr
  test_mae_values <- numeric(k)  # Placeholder for MAE
  test_bias_values <- numeric(k)  # Placeholder for bias
  test_APV_values <- numeric(k)  # Placeholder for Average Prediction Variance
  test_crps_values <- numeric(k)  # Placeholder for CRPS
  test_pred_covs50 <- numeric(k)  # Placeholder for 50% prediction covs
  test_pred_covs80 <- numeric(k)  # Placeholder for 80% prediction covs
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), observed = numeric(), inla_IHME_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    print(paste("Processing fold", i, "out of", k))
    
    #build non-convex hull mesh
    #non_convex_bdry <- inla.nonconvex.hull(train_data, -0.03, -0.05, resolution = c(100, 100))
    #mesh2 <- fm_mesh_2d_inla(boundary = non_convex_bdry, max.edge=c(0.1, 1), 
    # offset = c(0.05, 1),
    # cutoff = 0.003)
    
    #-Define the coordinates of centroids
    coords1 <- cbind(train_data$LONGNUM, train_data$LATNUM) 
    
    #build mesh
    mesh2 <- fm_mesh_2d_inla(loc=coords1,boundary=shp_ng, max.edge=c(0.05, 0.6),cutoff=0.1)
    
    
    plot(mesh2)
    #plot(mesh2, add=T)
    points(train_data, col="red", pch="*")
    
    #Count of mesh nodes
    mesh2$n
    
    #spde parameters
    alpha <- 2
    r0 <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction (i.e. 0.05*(ymax-ymin))
    
    #Matern SPDE model object using inla.pcmatern
    spde2 <- inla.spde2.pcmatern(mesh=mesh2, alpha=alpha, prior.range=c(r0, 0.01), prior.sigma=c(3, 0.01))
    
    
    # Priors
    hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
    control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  
    
    
    # Define model components and constraint
    cmp2 <- Age9_35_Vax ~ Intercept(1) + x(cbind(gam_predictions, 
                                                 gbm_predictions, lasso_predictions), 
                                           model = "iid",
                                           mapper = bru_mapper_matrix(labels = c("gam_predictions", "gbm_predictions", "lasso_predictions")),
                                           hyper = list(prec = list(initial = -log(30), fixed = TRUE)), #can be initial = -log(30) or log(1/1000) #-2*log(30)
                                           extraconstr = list(A = cbind(1, 1, 1), e = 1)) +
      field(coordinates, model = spde) + f.iid(id, model="iid", hyper = hyper.prec)
    
    #Fit model
    res2 <- bru(components = cmp2,  family = "binomial", 
                Ntrials = train_data$Age9_35_Tot, data = train_data,
                options = list(control.fixed=control.fixed))
    
    #Train data Predictions
    # train_predictions <- predict(res2, newdata = train_data, 
    #  formula = ~ Intercept + x + f.iid + 
    #   field_eval(cbind(LONGNUM, LATNUM)))
    
    #back transform predictions
    # train_predictions$predicted <- inverse_logit(train_predictions$mean) 
    
    #convert to tibble and get only observed and predicted values
    # train_predictions <- train_predictions %>% 
    # as_tibble() %>% 
    #  mutate(observed = prop,
    #        predicted = predicted)
    
    
    #Train data metrics
    # train_rmse_values[i] <- sqrt(mean((train_predictions$observed - train_predictions$predicted)^2))
    # train_pearson_values[i] <- cor(train_predictions$observed, train_predictions$predicted)
    # train_mae_values[i] <- mean(abs(train_predictions$observed - train_predictions$predicted))
    # train_bias_values[i] <- mean(train_predictions$observed - train_predictions$predicted)
    
    
    #Make predictions for test data
    test_samples <- generate(res2, newdata = test_data, 
                             formula = ~ Intercept + x + f.iid + 
                               field_eval(cbind(LONGNUM, LATNUM)), n.samples = 1000)
    
    
    #back transform predictions
    test_samples <- inverse_logit(test_samples) 
    
    
    #Get predictions as dataframe
    test_predictions <- data.frame(test_samples)
    
    
    #summarize predictions
    mean <- rowMeans(test_predictions, na.rm = T)
    std <- apply(test_predictions, 1, sd)
    lower_quantile <- apply(test_predictions, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
    upper_quantile <- apply(test_predictions, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))
    
    #cbind them as tibble
    test_predictions <- cbind(mean, std, lower_quantile, upper_quantile)
    
    
    #Get prediction metrics
    test_predictions <- test_predictions %>% 
      as_tibble() %>% 
      mutate(observed = test_data$prop, 
             DHSCLUST = test_data$DHSCLUST) %>% 
      rename(predicted = mean)
    
    
    test_rmse_values[i] <- sqrt(mean((test_predictions$observed - test_predictions$predicted)^2))
    test_pearson_values[i] <- cor(test_predictions$observed, test_predictions$predicted)
    test_mae_values[i] <- mean(abs(test_predictions$observed - test_predictions$predicted))
    test_bias_values[i] <- mean(test_predictions$observed - test_predictions$predicted)
    test_APV_values[i] <- mean(test_predictions$std)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_predictions$DHSCLUST, folds = i, observed = test_predictions$observed, inla_IHME_predictions = test_predictions$predicted)
    results_df <- rbind(results_df, fold_results)
    
    # CRPS calculation
    
    # Cbind observed data to samples
    smp_data <- cbind(test_data$prop, test_samples) %>% 
      as_tibble() %>%
      drop_na()
    
    # Number of iterations excluding the first column (which is the observed data)
    its <- ncol(smp_data) - 1
    
    # Function to calculate the crps
    crps_fnc <- function(p, n) {
      out <- map_dbl(1:n, ~ abs(p[.x + 1] - p[1]))
      out1 <- map_dbl(1:n, ~ sum(abs(p[.x + 1] - p[-1])))
      
      out <- sum(out) / n
      out1 <- sum(out1) / (2 * n * n)
      out - out1
    }
    
    crps <- smp_data %>%
      rowwise() %>%
      mutate(crps_value = crps_fnc(c_across(everything()), its)) %>%
      ungroup() %>%
      summarise(mean_crps = mean(crps_value)) %>%
      pull(mean_crps)
    
    test_crps_values[i] <- crps
    
  }
  
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "INLA_IHME_Kfold_Predictions.csv"))
  
  list(#train_rmse = mean(train_rmse_values),
    #train_corr = mean(train_pearson_values),
    #train_mae = mean(train_mae_values),
    #train_bias = mean(train_bias_values),
    test_rmse = mean(test_rmse_values),
    test_corr = mean(test_pearson_values),
    test_mae = mean(test_mae_values),
    test_bias = mean(test_bias_values),
    test_APV = mean(test_APV_values),
    test_crps = mean(test_crps_values),
    test_covs50 = mean(test_pred_covs50),
    test_covs80 = mean(test_pred_covs80))
  
}

tic()

# Apply function
result1 <- kfold_cv(data = dtp_sf, k = 10)
result1

toc()

##########################################################################################
###########################################################################################

# Spatial K-Fold Cross Validation -----------------------------------------

# function to calculate spatial k-fold

spatial_kfold_cv <- function(data) {
  fold <- data$Group_ID
  k <- length(unique(fold))
  
  #Train metrics placeholder
  
  #train_rmse_values <- numeric(k) # Placeholder for RMSE
  #train_pearson_values <- numeric(k) # Placeholder for corr
  # train_mae_values <- numeric(k)  # Placeholder for MAE
  # train_bias_values <- numeric(k)  # Placeholder for bias
  
  #Test metrics placeholder
  test_rmse_values <- numeric(k) # Placeholder for RMSE
  test_pearson_values <- numeric(k) # Placeholder for corr
  test_mae_values <- numeric(k)  # Placeholder for MAE
  test_bias_values <- numeric(k)  # Placeholder for bias
  test_APV_values <- numeric(k)  # Placeholder for Average Prediction Variance
  test_crps_values <- numeric(k)  # Placeholder for CRPS
  test_pred_covs50 <- numeric(k)  # Placeholder for 50% prediction covs
  test_pred_covs80 <- numeric(k)  # Placeholder for 80% prediction covs
  
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), observed = numeric(), inla_IHME_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(fold == i)
    train_indices <- which(fold != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    print(paste("Processing fold", i, "out of", k))
    
    #-Define the coordinates of centroids
    coords1 <- cbind(train_data$LONGNUM, train_data$LATNUM) 
    
    #build mesh
    mesh2 <- fm_mesh_2d_inla(loc=coords1,boundary=shp_ng, max.edge=c(0.05, 0.6),cutoff=0.1)
    
    plot(mesh2)
    #plot(mesh2, add=T)
    points(train_data, col="red", pch="*")
    
    #Count of mesh nodes
    mesh2$n
    
    #spde parameters
    alpha <- 2
    r0 <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction (i.e. 0.05*(ymax-ymin))
    
    #Matern SPDE model object using inla.pcmatern
    spde2 <- inla.spde2.pcmatern(mesh=mesh2, alpha=alpha, prior.range=c(r0, 0.01), prior.sigma=c(3, 0.01))
    
    # Priors
    hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
    control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  
    
    # Define model components and constraint
    cmp2 <- Age9_35_Vax ~ Intercept(1) + x(cbind(gam_predictions, 
                                                 gbm_predictions, lasso_predictions), 
                                           model = "iid",
                                           mapper = bru_mapper_matrix(labels = c("gam_predictions", "gbm_predictions", "lasso_predictions")),
                                           hyper = list(prec = list(initial = -log(30), fixed = TRUE)), #can be initial = -log(30) or log(1/1000) #-2*log(30)
                                           extraconstr = list(A = cbind(1, 1, 1), e = 1)) +
      field(coordinates, model = spde) + f.iid(id, model="iid", hyper = hyper.prec)
    
    #Fit model
    res2 <- bru(components = cmp2,  family = "binomial", 
                Ntrials = train_data$Age9_35_Tot, data = train_data,
                options = list(control.fixed=control.fixed))
    
    
    #Train data Predictions
    # train_predictions <- predict(res2, newdata = train_data, 
    # formula = ~ Intercept + x + f.iid + 
    #  field_eval(cbind(LONGNUM, LATNUM)))
    
    #back transform predictions
    # train_predictions$predicted <- inverse_logit(train_predictions$mean) 
    
    #convert to tibble and get only observed and predicted values
    # train_predictions <- train_predictions %>% 
    # as_tibble() %>% 
    # mutate(observed = prop,
    #       predicted = predicted)
    
    
    #Train data metrics
    # train_rmse_values[i] <- sqrt(mean((train_predictions$observed - train_predictions$predicted)^2))
    # train_pearson_values[i] <- cor(train_predictions$observed, train_predictions$predicted)
    # train_mae_values[i] <- mean(abs(train_predictions$observed - train_predictions$predicted))
    # train_bias_values[i] <- mean(train_predictions$observed - train_predictions$predicted)
    
    
    #Make predictions for test data
    test_samples <- generate(res2, newdata = test_data, 
                             formula = ~ Intercept + x + f.iid + 
                               field_eval(cbind(LONGNUM, LATNUM)), n.samples = 1000)
    
    #back transform predictions
    test_samples <- inverse_logit(test_samples) 
    
    #Get predictions as dataframe
    test_predictions <- data.frame(test_samples)
    
    #summarize predictions
    mean <- rowMeans(test_predictions, na.rm = T)
    std <- apply(test_predictions, 1, sd)
    lower_quantile <- apply(test_predictions, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
    upper_quantile <- apply(test_predictions, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))
    
    #cbind them as tibble
    test_predictions <- cbind(mean, std, lower_quantile, upper_quantile)
    
    
    #Get prediction metrics
    test_predictions <- test_predictions %>% 
      as_tibble() %>% 
      mutate(observed = test_data$prop, 
             DHSCLUST = test_data$DHSCLUST) %>% 
      rename(predicted = mean)
    
    test_rmse_values[i] <- sqrt(mean((test_predictions$observed - test_predictions$predicted)^2))
    test_pearson_values[i] <- cor(test_predictions$observed, test_predictions$predicted)
    test_mae_values[i] <- mean(abs(test_predictions$observed - test_predictions$predicted))
    test_bias_values[i] <- mean(test_predictions$observed - test_predictions$predicted)
    test_APV_values[i] <- mean(test_predictions$std)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_predictions$DHSCLUST, folds = i, observed = test_predictions$observed, inla_IHME_predictions = test_predictions$predicted)
    results_df <- rbind(results_df, fold_results)
    
    # CRPS calculation
    
    # Cbind observed data to samples
    smp_data <- cbind(test_data$prop, test_samples) %>% 
      as_tibble() %>%
      drop_na()
    
    # Number of iterations excluding the first column (which is the observed data)
    its <- ncol(smp_data) - 1
    
    # Function to calculate the crps
    crps_fnc <- function(p, n) {
      out <- map_dbl(1:n, ~ abs(p[.x + 1] - p[1]))
      out1 <- map_dbl(1:n, ~ sum(abs(p[.x + 1] - p[-1])))
      
      out <- sum(out) / n
      out1 <- sum(out1) / (2 * n * n)
      out - out1
    }
    
    crps <- smp_data %>%
      rowwise() %>%
      mutate(crps_value = crps_fnc(c_across(everything()), its)) %>%
      ungroup() %>%
      summarise(mean_crps = mean(crps_value)) %>%
      pull(mean_crps)
    
    test_crps_values[i] <- crps
    
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "INLA_IHME_Spatial_Kfold_Predictions.csv"))
  
  list(#train_rmse = mean(train_rmse_values),
    #train_corr = mean(train_pearson_values),
    #train_mae = mean(train_mae_values),
    #train_bias = mean(train_bias_values),
    test_rmse = mean(test_rmse_values),
    test_corr = mean(test_pearson_values),
    test_mae = mean(test_mae_values),
    test_bias = mean(test_bias_values),
    test_APV = mean(test_APV_values),
    test_crps = mean(test_crps_values),
    test_covs50 = mean(test_pred_covs50),
    test_covs80 = mean(test_pred_covs80))
  
}

# Apply function

tic()
result2 <- spatial_kfold_cv(data = dtp_sf)

result2

toc()

b <- cbind(result1, result2)
b

#########################################################################################################
################## Make Predictions #####################################################################


# Load Rasters for predictions --------------------------------------------

#import covaraites rasters

raster_list <-list.files(path=input_path, pattern= ".tif$", all.files=TRUE, full.names=FALSE)
raster_list

r1 <- rast(paste0(input_path, "GAM.tif"))
r2 <- rast(paste0(input_path, "GBM.tif"))
r3 <- rast(paste0(input_path, "lasso.tif"))


#Stack all rasters
stack_rasters<- c(r1, r2, r3)

#get raster values
raster_df <- terra::values(stack_rasters, dataframe = T)

# Get the xy coordinate of the centroid of each pixel as a dataframe
coord <- xyFromCell(r1, 1:ncell(r1))

#cbind coordinates to predictions
pred_covs <- cbind(raster_df, coord)

#Rename variables for consistency with model object
pred_covs <- pred_covs %>% 
  rename(LONGNUM = x, LATNUM = y,
         lasso_predictions = lasso) 

#create unique id for pixels to be used as a random effect component
pred_covs <- pred_covs %>% 
  tibble::rowid_to_column("id")

#Make predictions
mu <- generate(res, newdata = pred_covs,
               formula = ~ Intercept + x + f.iid+
                 field_eval(cbind(LONGNUM, LATNUM)), n.samples = 1000)


#Apply Inverse Logit function to posteriors
inv.linpred <- inverse_logit(mu) 


#summarize predictions
mean <- rowMeans(inv.linpred, na.rm = T)
std <- apply(inv.linpred, 1, sd)
lower_quantile <- apply(inv.linpred, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
upper_quantile <- apply(inv.linpred, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))

#cbind them as tibble
inv.linpred <- cbind(mean, std, lower_quantile, upper_quantile) %>% 
  as_tibble() %>% 
  mutate(LONGNUM = pred_covs$LONGNUM, LATNUM = pred_covs$LATNUM) # Add lat long to data


summary(inv.linpred$std)

# Rasterize Predictions ---------------------------------------------------

# convert pixel_predictions to an sf object
inv.linpred_sf <- st_as_sf(inv.linpred, coords = c("LONGNUM", "LATNUM"))

# set the CRS of the sf object
st_crs(inv.linpred_sf) <- 4326

# reproject to raster spatial reference
inv.linpred_sf <- st_transform(inv.linpred_sf, crs = st_crs(r1))


#Rasterize mean predictions
predictions_raster <- rasterize(inv.linpred_sf, r1, field = "mean")
plot(predictions_raster)


writeRaster(predictions_raster, paste0(input_path, "STG_mcv1_9_35_mean.tif"), 
            overwrite = T, names = "mcv1_Coverage")

#Rasterize sd
sd_raster <- rasterize(inv.linpred_sf, r1, field = "std")
plot(sd_raster)


writeRaster(sd_raster, paste0(input_path, "STG_mcv1_9_35_sd.tif"), 
            overwrite = T, names = "mcv1_sd")

#Rasterize upper bound
upper_raster <- rasterize(inv.linpred_sf, r1, field = "upper_quantile")
plot(upper_raster)


writeRaster(upper_raster, paste0(input_path, "STG_mcv1_9_35_upper.tif"), 
            overwrite = T, names = "mcv1_upper")


#Rasterize lower bound
lower_raster <- rasterize(inv.linpred_sf, r1, field = "lower_quantile")
plot(lower_raster)


writeRaster(lower_raster, paste0(input_path, "STG_mcv1_9_35_lower.tif"), 
            overwrite = T, names = "mcv1_lower")

###################### END #################################################################
############################################################################################
############################################################################################



