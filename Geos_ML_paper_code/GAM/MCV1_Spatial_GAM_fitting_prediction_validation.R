#Generalize Additive Model (GAM) including Lat Long as Variables
#Load packages

library(tidyverse)
library(mgcv)
library(terra)
library(sf)
library(groupdata2)

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
  mutate(prop = Age9_35_Vax/Age9_35_Tot)

#Create a grouping variable sequentially to be used for spatial k-fold cross val
mcv1 <- mcv1 %>% 
  group(n = 133, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()

#Recode Rural & Urban as dummy variables
mcv1 <- mcv1 %>% 
  mutate(urban_rural = case_when(URBAN_RURA == "U" ~ 1,   
                                 URBAN_RURA == "R" ~ 0))

#fit GAM model
fit_GAM <- mgcv::gam(
  formula = Age9_35_Vax/Age9_35_Tot ~ s(n26_Mean, bs = "cs")  + s(n28_Mean, bs = "cs") + s(n51_Mean, bs = "cs") +
    s(n58_Mean, bs = "cs") + s(n69_Mean, bs = "cs")  + s(wealth.prop, bs = "cs") + s(educ.prop, bs = "cs") +
    s(health_card_doc.prop, bs = "cs") + s(sba.prop, bs = "cs") + s(phone_internet.prop, bs = "cs") + 
    s(n2_Mean, bs = "cs") + s(media.prop, bs = "cs") +s(n12_Mean, bs = "cs")+
    urban_rural + s(LATNUM, LONGNUM, bs="sos",m=2,k=100),
  family = binomial("logit"),
  weights = Age9_35_Tot,
  data = mcv1)

summary(fit_GAM)

#Get the variance
gam.vcomp(fit_GAM)

# Get In Sample Predictions ---------------------------------------------------
#Select covs
covs1 <- mcv1 %>% 
  select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
         educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
         n2_Mean,media.prop, n12_Mean, urban_rural, LATNUM, LONGNUM)


#Make predictions using the fitted model
in_sample_predictions <- mgcv::predict.gam(object = fit_GAM,  newdata = covs1, type = "response")

in_sample_predictions <- as_tibble(in_sample_predictions) %>% 
  rename(gam_predictions = value)


#cbind cluster id and response variables
gam_data <- mcv1 %>% 
  select(DHSCLUST, LATNUM, LONGNUM, Age9_35_Vax, Age9_35_Tot, prop) %>% 
  cbind(in_sample_predictions)

write.csv(gam_data, paste0(output_path, "spatial_gam_data.csv"))


#In sample metrics
metrics <- gam_data %>% 
  rename(observed = prop,
         predicted = gam_predictions) %>% 
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
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), gam_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    #Model
    model <- mgcv::gam(formula, data = train_data, weights = Age9_35_Tot, family = binomial("logit"))
    
    #Predictions
    predictions <- predict(model, newdata = test_data)
    predictions <- mgcv::predict.gam(object = model,  newdata = test_data, type = "response")
    
    predictions <- as_tibble(predictions) %>% 
      rename(gam_predictions = value)
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$gam_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$gam_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$gam_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$gam_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, gam_predictions = predictions$gam_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Spatial_Gam_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}

# Apply function
result1 <- kfold_cv(data = mcv1,  formula = Age9_35_Vax/Age9_35_Tot ~ s(n26_Mean, bs = "cs")  + s(n28_Mean, bs = "cs") + s(n51_Mean, bs = "cs") +
                      s(n58_Mean, bs = "cs") + s(n69_Mean, bs = "cs")  + s(wealth.prop, bs = "cs") + s(educ.prop, bs = "cs") +
                      s(health_card_doc.prop, bs = "cs") + s(sba.prop, bs = "cs") + s(phone_internet.prop, bs = "cs") + 
                      s(n2_Mean, bs = "cs") + s(media.prop, bs = "cs") +s(n12_Mean, bs = "cs")+
                      urban_rural + s(LATNUM, LONGNUM, bs="sos",m=2,k=100))

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
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), gam_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(fold == i)
    train_indices <- which(fold != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    #Model
    model <- mgcv::gam(formula, data = train_data, weights = Age9_35_Tot,
                       family = binomial("logit"))
    
    #Predictions
    predictions <- predict(model, newdata = test_data)
    predictions <- mgcv::predict.gam(object = model,  newdata = test_data, type = "response")
    
    predictions <- as_tibble(predictions) %>% 
      rename(gam_predictions = value)
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$gam_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$gam_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$gam_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$gam_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, gam_predictions = predictions$gam_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Spatial_Gam_Spatial_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}

# Apply function
result2 <- spatial_kfold_cv(data = mcv1,  formula = Age9_35_Vax/Age9_35_Tot ~ s(n26_Mean, bs = "cs")  + s(n28_Mean, bs = "cs") + s(n51_Mean, bs = "cs") +
                              s(n58_Mean, bs = "cs") + s(n69_Mean, bs = "cs")  + s(wealth.prop, bs = "cs") + s(educ.prop, bs = "cs") +
                              s(health_card_doc.prop, bs = "cs") + s(sba.prop, bs = "cs") + s(phone_internet.prop, bs = "cs") + 
                              s(n2_Mean, bs = "cs") + s(media.prop, bs = "cs") +s(n12_Mean, bs = "cs")+
                              urban_rural + s(LATNUM, LONGNUM, bs="sos",m=2,k=100))
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
#Remember predictions are on the log scale
predictions <- mgcv::predict.gam(object = fit_GAM,  newdata = pred_covs, type = "response", se.fit = TRUE) 


# Calculate the lower and upper bounds on the log-odds scale
lower_bound_log_odds <- predictions$fit - 1.96 * predictions$se.fit
upper_bound_log_odds <- predictions$fit + 1.96 * predictions$se.fit

# Back-transform the bounds to the probability scale
lower_bound_prob <- plogis(lower_bound_log_odds)
upper_bound_prob <- plogis(upper_bound_log_odds)

#Back transform predicted values
predicted_values <- plogis(predictions$fit)

# Transform the standard errors using the delta method
se_probabilities <- predictions$se.fit * predicted_values * (1 - predicted_values)


# Combine the results into a data frame for easier viewing
raster_predictions <- data.frame(
  #predictions = predicted_values,
  standard_erors = predictions$se.fit
  #lower_bound = lower_bound_prob,
  #upper_bound = upper_bound_prob
)

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
GAM_raster <- rasterize(raster_predictions1, r1, field = "predictions")

plot(GAM_raster)

writeRaster(GAM_raster , paste0(output_path, "Spatial_GAM.tif"), 
            overwrite = T, names = "gam_predictions")

#Rasterize lower bound
lower_raster <- rasterize(raster_predictions1, r1, field = "lower_bound")

plot(lower_raster)

writeRaster(lower_raster , paste0(output_path, "Spatial_GAM_lower.tif"), 
            overwrite = T, names = "lower_bound")

#Rasterize upper bound
upper_raster <- rasterize(raster_predictions1, r1, field = "upper_bound")

plot(upper_raster)

writeRaster(upper_raster , paste0(output_path, "Spatial_GAM_upper.tif"), 
            overwrite = T, names = "upper_bound")


#Rasterize standard errors
se_raster <- rasterize(raster_predictions1, r1, field = "standard_erors")

plot(se_raster)

writeRaster(se_raster , paste0(output_path, "Spatial_GAM_SE.tif"), 
            overwrite = T, names = "Standard Errors")
