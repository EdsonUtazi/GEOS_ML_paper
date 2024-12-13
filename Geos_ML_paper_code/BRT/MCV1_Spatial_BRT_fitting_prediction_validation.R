## mcv1 Modelling using XGBoost Regression

#Load libraries
library(caret)
library(raster)
library(tidyverse)
library(gbm)
library(terra)
library(sf)
library(secr)
library(groupdata2)


options(scipen = 999) # turn off scientific notation for all variables

#specify drive path
drive_path <- "/.../Working/"
data_path <- paste0(drive_path, "Methods_comparison_paper/2018_NGA_DHS/")
output_path <- paste0(drive_path, "Methods_comparison_paper/.../MCV1/")

#read data
mcv1 <- read.csv(paste0(data_path, "NG_Vax_measles.csv"))
covs <- read.csv(paste0(data_path, "Covariates_merged_selected.csv"))

#select covs and add to response variable
mcv1 <- mcv1 %>% 
  select(Age9_35_Vax, Age9_35_Tot, DHSCLUST) %>% 
  inner_join(covs, by = "DHSCLUST") %>% 
  filter(Age9_35_Tot > 1) %>%    #Remove observations where total is not 0
  mutate(prop = (Age9_35_Vax/Age9_35_Tot), logit_prop = logit(abs(prop-0.000001)))

#Create a grouping variable sequentially to be used for spatial k-fold cross val
mcv1 <- mcv1 %>% 
  group(n = 133, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()


#Recode Rural & Urban as dummy variables
mcv1 <- mcv1 %>% 
  mutate(urban_rural = case_when(URBAN_RURA == "U" ~ 1,   
                                 URBAN_RURA == "R" ~ 0))

mod <-  gbm::gbm(
  formula = logit_prop ~ n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
    educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop + 
    n2_Mean + media.prop + n12_Mean + urban_rural+ LATNUM + LONGNUM,
  distribution = "gaussian",
  data = mcv1,
  n.trees = 10000,
  shrinkage = 0.01)

# Get the predictions
predictions <- predict(mod, n.trees = 10000)

# Calculate the standard deviation
std_dev <- sd(predictions)

# Calculate the variance
variance <- var(predictions)


# In-Sample Predictions ---------------------------------------------------
#Select covs
covs1 <- mcv1 %>% 
  select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
         educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
         n2_Mean,media.prop, n12_Mean, urban_rural, LATNUM, LONGNUM)


#Make predictions using the fitted model
in_sample_predictions <- gbm::predict.gbm(
  object = mod, 
  newdata = covs1, 
  n.trees = 10000, 
  type = "response")

in_sample_predictions <- as.numeric(in_sample_predictions) 

#Inverse logit function
inverse_logit <- function(aa) {
  
  out = exp(aa)/(1+exp(aa))
  return(out)
  
}

#Apply function
in_sample_predictions <- inverse_logit(in_sample_predictions) %>% 
  as_tibble() %>% 
  rename(gbm_predictions = value)

#cbind cluster id and response variables
gbm_data <- mcv1 %>% 
  select(DHSCLUST, LATNUM, LONGNUM, Age9_35_Vax, Age9_35_Tot, prop) %>% 
  cbind(in_sample_predictions)

write.csv(gbm_data, paste0(output_path, "gbm_data_spatial.csv"))

#In sample metrics
metrics <- gbm_data  %>% 
  rename(observed = prop,
         predicted = gbm_predictions) %>% 
  mutate(residual = observed - predicted)%>% 
  
  summarise(
    Bias= mean(residual),
    Imprecision = sd(residual),
    mae = mean(abs(residual)),
    mse = mean((residual)^2),
    rmse = sqrt(mse),
    Corr = cor(observed, predicted))
metrics


# K-Fold Cross Validation -------------------------------------------------

# function to calculate k-fold
kfold_cv <- function(data, formula, k = 10) {
  n <- nrow(data)
  fold_size <- n %/% k
  folds <- sample(rep(1:k, each = fold_size, length.out = n))
  
  rmse_values <- numeric(k) # Placeholder for RMSE
  pearson_values <- numeric(k) # Placeholder for corr
  mae_values <- numeric(k)  # Placeholder for MAE
  bias_values <- numeric(k)  # Placeholder for bias
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), gbm_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    #Model
    model <- gbm::gbm(formula, data = train_data, distribution = "gaussian", n.trees = 10000, shrinkage = 0.01) 
    
    #Predictions
    predictions <- gbm::predict.gbm(model, newdata = test_data, n.trees = 10000, type = "response")
    
    predictions <- as.numeric(predictions) 
    
    #Apply inverse logit function
    predictions <- inverse_logit(predictions) %>% 
      as_tibble() %>% 
      rename(gbm_predictions = value)
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$gbm_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$gbm_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$gbm_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$gbm_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, gbm_predictions = predictions$gbm_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Spatial_GBM_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}

# Apply function
result1 <- kfold_cv(data = mcv1, formula = logit_prop ~ n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                      educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop + 
                      n2_Mean + media.prop + n12_Mean + urban_rural+ LATNUM + LONGNUM)

result1

# Spatial K-Fold Cross Validation -----------------------------------------

# function to calculate spatial k-fold

spatial_kfold_cv <- function(data, formula) {
  fold <- data$Group_ID
  k <- length(unique(fold))
  
  rmse_values <- numeric(k) # Placeholder for RMSE
  pearson_values <- numeric(k) # Placeholder for corr
  mae_values <- numeric(k)  # Placeholder for MAE
  bias_values <- numeric(k)  # Placeholder for bias
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), gbm_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(fold == i)
    train_indices <- which(fold != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    #Model
    model <- gbm::gbm(formula, data = train_data, distribution = "gaussian", n.trees = 10000, shrinkage = 0.01) 
    
    #Predictions
    predictions <- gbm::predict.gbm(model, newdata = test_data, n.trees = 10000, type = "response")
    
    predictions <- as.numeric(predictions) 
    
    #Apply inverse logit function
    predictions <- inverse_logit(predictions) %>% 
      as_tibble() %>% 
      rename(gbm_predictions = value)
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$gbm_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$gbm_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$gbm_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$gbm_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, gbm_predictions = predictions$gbm_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Spatial_GBM_Spatial_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}


# Apply function
result2 <- spatial_kfold_cv(data = mcv1, formula = logit_prop ~ n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                              educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop + 
                              n2_Mean + media.prop + n12_Mean + urban_rural+ LATNUM + LONGNUM)
result2


# Load Rasters for predictions --------------------------------------------
#import covaraites rasters

raster_list <-list.files(path=data_path, pattern= ".tif$", all.files=TRUE, full.names=FALSE)
raster_list

raster_covs <- raster_list[-c(7, 12, 14, 17)]
raster_covs


#Stack all rasters
stack_rasters<- rast(paste0(data_path, c(raster_covs)))

#get raster values
raster_df <- terra::values(stack_rasters, dataframe = T)

#Read raster and get xy values
r1 <- rast(paste0(data_path, "n28_viirs_2016_nightlights.tif"))

# Get the xy coordinate of the centroid of each pixel as a dataframe
coord <- xyFromCell(r1, 1:ncell(r1))

#Rename coord lat-long
coord <- coord %>% 
  as_tibble() %>% 
  rename(LONGNUM = x, LATNUM = y)

#cbind coordinates to raster
raster_df <- cbind(raster_df, coord)

#get covs names in order for  model predictions
pred_covs <- raster_df %>% 
  rename(n26_Mean = n26_guf_ghsl_dst_2014,
         n28_Mean = n28_viirs_2016_nightlights,
         n51_Mean = n51_MODIS_lst,
         n58_Mean = n58_TravelTime_to_HF_ver2,
         n69_Mean = n69_WorldPop10povcons200,
         wealth.prop = wealthkrig,
         educ.prop = educkrig,
         health_card_doc.prop = health_card_dockrig,
         sba.prop = sbakrig, 
         phone_internet.prop = phone_internetkrig,
         n2_Mean = n2_cattle_FAO_2010,
         media.prop = mediakrig,
         n12_Mean = n12_ccilc_dst011_2015,
         urban_rural = NGA_urban_rural_1km)

#Make predictions using the fitted model
raster_predictions <- gbm::predict.gbm(
  object = mod, 
  newdata = pred_covs, 
  n.trees = 10000, 
  type = "response")

raster_predictions <- as.numeric(raster_predictions) 


#Apply inverse logit transformation function
raster_predictions <- inverse_logit(raster_predictions) %>% 
  as_tibble() %>% 
  rename(predictions = value)


# Rasterize predictions ---------------------------------------------------

#Read raster and get xy values
r1 <- rast(paste0(data_path, "n28_viirs_2016_nightlights.tif"))

# Get the xy coordinate of the centroid of each pixel as a dataframe
coord <- xyFromCell(r1, 1:ncell(r1))

#cbind coordinates to predictions
raster_predictions <- cbind(raster_predictions, coord)

# convert the dataframe to an sf object
raster_predictions1 <- st_as_sf(raster_predictions, coords = c("x", "y"))

# set the CRS of the sf object
st_crs(raster_predictions1) <- 4326

# set the CRS of the sf object
raster_predictions1 <- st_transform(raster_predictions1, crs = st_crs(r1))

#Rasterize
GBM_raster <- rasterize(raster_predictions1, r1, field = "predictions")

plot(GBM_raster)

writeRaster(GBM_raster , paste0(output_path, "GBM_spatial.tif"), 
            overwrite = T, names = "gbm_spatial_predictions")


###################################################################################################
###################################################################################################
output_path1 <- paste0(drive_path, "Methods_comparison_paper/.../MCV1/JackKnife/")

mcv_full <- mcv1  #NOTE - Run this once before starting the jackknife part

# JackKnife Predictions to Get Uncertainty Calculations -------------------

### Start of jackknife section: Analyzing with different replicates/block sizes

# Define number of replicates (iterations)
nrast <- 100  # Number of replicates to run (e.g., 10, 30, 50, 100, 200)

# Calculate block size based on replicates and data size
bsize <- floor(nrow(mcv_full) / nrast)  # Number of data points per replicate

# Set random seed for reproducibility
srand <- sample(1:nrow(mcv_full), nrow(mcv_full), replace = FALSE)

# Loop through each replicate
for (kk in 1:nrast) {
  
  # Define data points to exclude for this replicate
  if (kk < nrast) {
    qq <- (((kk - 1) * bsize) + 1):(bsize * kk)
  } else {
    qq <- (((kk - 1) * bsize) + 1):nrow(mcv_full)
  }
  samp.k <- srand[qq]  # Data points to be excluded
  
  # Create a subset of the full data excluding the current replicate's points
  mcv1 <- mcv_full[-samp.k, ]
  
  # Print the number of data points in the current subset
  print(nrow(mcv1))
  
  mod2 <-  gbm::gbm(
    formula = logit_prop ~ n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
      educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop + 
      n2_Mean + media.prop + n12_Mean + urban_rural+ LATNUM + LONGNUM,
    distribution = "gaussian",
    data = mcv1,
    n.trees = 10000,
    shrinkage = 0.01)
  
  #Make predictions using the fitted model
  raster_predictions <- gbm::predict.gbm(
    object = mod2, 
    newdata = pred_covs, 
    n.trees = 10000, 
    type = "response")
  
  raster_predictions <- as.numeric(raster_predictions) 
  
  # Transform predictions to probability
  raster_predictions <- inverse_logit(raster_predictions)
  
  
  r1 <- raster(paste0(data_path, "n28_viirs_2016_nightlights.tif"))
  values(r1) <- as.numeric(raster_predictions)
  
  writeRaster(r1 , paste0(output_path1, "spatial_jack_gbm_", kk, ".tif"), 
              overwrite = T)
}

#Calculate jackknife mean and std. dev.
raster_list <-list.files(path=output_path1,pattern="spatial_jack_gbm.*\\.tif$" , all.files=TRUE, full.names=FALSE)
raster_list

stack_rasters<- rast(paste0(output_path1, c(raster_list)))

#Jackknife mean
mean.rast <- app(stack_rasters, mean, na.rm=TRUE)
plot(mean.rast)

#Mask Rasters before exporting
r3 <- rast(paste0(output_path, "GAM.tif"))

mean.rast <- terra::mask(mean.rast, r3)
plot(mean.rast)

#Jackknife sd
fun=function(i) i*i
sum <- 0
for (i in 1:nrast){
  s <- subset(stack_rasters, i)
  difr <- s - mean.rast 
  sum <- sum + app(difr, fun)
}

sd.rast <- sqrt(((nrast-1)/nrast)*sum)
#plot(sd.rast, range=c(0,0.5))
plot(sd.rast)

sd.rast <- terra::mask(sd.rast, r3)
plot(sd.rast)

sd.rast <- raster(sd.rast)
summary(getValues(sd.rast))

writeRaster(mean.rast, paste0(output_path1, "Spatia_GBM_mean_raster_mcv1_", nrast,".tif"), 
            overwrite = T)
writeRaster(sd.rast, paste0(output_path1, "Spatial_GBM_stddev_raster_mcv1_", nrast,".tif"), 
            overwrite = T)

##################### END OF SCRIPT ####################################################
########################################################################################























