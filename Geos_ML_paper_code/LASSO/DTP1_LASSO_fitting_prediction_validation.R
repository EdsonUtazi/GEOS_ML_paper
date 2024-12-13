## dtp1 Modelling using Lasso Regression

#Load libraries
library(caret)
library(raster)
library(tidyverse)
library(glmnet)
library(terra)
library(sf)
library(groupdata2)

#specify drive path
drive_path <- "/.../Working/"
data_path <- paste0(drive_path, "Methods_comparison_paper/2018_NGA_DHS/")
output_path <- paste0(drive_path, "Methods_comparison_paper/.../DTP1/")

#read data
dtp1 <- read.csv(paste0(data_path, "NG_Vax_DTP1.csv"))
covs <- read.csv(paste0(data_path, "Covariates_merged_selected.csv"))


#select covs and add to response variable
dtp1 <- dtp1 %>% 
  select(Age12_23_Vax, Age12_23_Tot, DHSCLUST) %>% 
  inner_join(covs, by = "DHSCLUST") %>% 
  filter(Age12_23_Tot > 1) %>%  
  mutate(prop = Age12_23_Vax/Age12_23_Tot)

#Create a grouping variable sequentially to be used for spatial k-fold cross val
dtp1 <- dtp1 %>% 
  group(n = 116, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()

#Recode Rural & Urban as dummy variables
dtp1 <- dtp1 %>% 
  mutate(urban_rural = case_when(URBAN_RURA == "U" ~ 1,   
                                 URBAN_RURA == "R" ~ 0))

#Fit a model using the LASSO with the glmnet package

#Select covariates
#x <- dtp1 %>% 
# dplyr::select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
#   educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop,
# n2_Mean, media.prop, n12_Mean, urban_rural) %>%  
# as.matrix()



#Select response
y <- dtp1 %>% 
  dplyr::select(success = Age12_23_Vax, Age12_23_Tot) %>%   #count = number of success and total = total cases
  mutate(failure = Age12_23_Tot - success) %>% 
  select(failure, success) %>% 
  as.matrix()

##################NEW - create a model matrix for x
xx1 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                      educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                      n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = dtp1)

xx <-  as.data.frame(xx1[,-which(colnames(xx1)%in%c("factor(urban_rural)0"))])
xx <- as.matrix(xx)

# Perform cross-validation to find optimal lambda
cvfit <- cv.glmnet(x = xx, y = y, alpha = 1, family = "binomial")

#  optimal lambda
best_lambda <- cvfit$lambda.min

#Fit lasso model with optimal lambda

fit_lasso <- glmnet::glmnet(x=xx, y=y, family = "binomial", alpha = 1, lambda = best_lambda, intercept = TRUE)

#Find coefficient estimates
coef(fit_lasso)


# In-Sample Predictions ---------------------------------------------------
#Select covs
#covs1 <- dtp1 %>% 
#  select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
#         educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
#         n2_Mean,media.prop, n12_Mean, urban_rural) %>% 
#  as.matrix()

covs1 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                        educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                        n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = dtp1)

covs <-  as.data.frame(covs1[,-which(colnames(covs1)%in%c("factor(urban_rural)0"))])
covs <- as.matrix(covs)


#Make predictions using the fitted model
in_sample_predictions <-  glmnet::predict.glmnet(fit_lasso, newx = covs)

in_sample_predictions <- as.numeric(in_sample_predictions) 

#Inverse logit function
inverse_logit <- function(aa) {
  
  out = exp(aa)/(1+exp(aa))
  return(out)
  
}
#Apply function
in_sample_predictions <- inverse_logit(in_sample_predictions) %>% 
  as_tibble() %>% 
  rename(lasso_predictions = value)

#cbind cluster id and response variables
lasso_data <- dtp1 %>% 
  select(DHSCLUST, LATNUM, LONGNUM, Age12_23_Vax, Age12_23_Tot, prop) %>% 
  cbind(in_sample_predictions)

write.csv(lasso_data, paste0(output_path, "lasso_data.csv"))

#In sample metrics
metrics <- lasso_data  %>% 
  rename(observed = prop,
         predicted = lasso_predictions) %>% 
  mutate(residual = observed - predicted)%>% 
  
  summarise(
    Bias= mean(residual),
    Imprecision = sd(residual),
    mae = mean(abs(residual)),
    mse = mean((residual)^2),
    rmse = sqrt(mse),
    Corr = cor(observed, predicted))  #Changed something here
metrics

# K-Fold Cross Validation -------------------------------------------------

# function to calculate k-fold
kfold_cv <- function(data, k = 10) {
  n <- nrow(data)
  fold_size <- n %/% k
  folds <- sample(rep(1:k, each = fold_size, length.out = n))
  
  rmse_values <- numeric(k) # Placeholder for RMSE
  pearson_values <- numeric(k) # Placeholder for corr
  mae_values <- numeric(k)  # Placeholder for MAE
  bias_values <- numeric(k)  # Placeholder for bias
  
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), lasso_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    y <- train_data %>% 
      dplyr::select(success = Age12_23_Vax, Age12_23_Tot) %>%   
      mutate(failure = Age12_23_Tot - success) %>% 
      select(failure, success) %>% 
      as.matrix()
    
    #x <- train_data %>% 
    # dplyr::select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
    #               educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
    #               n2_Mean,media.prop,n12_Mean, urban_rural) %>% 
    #as.matrix()
    
    xx1 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                          educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                          n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = train_data)
    
    xx <-  as.data.frame(xx1[,-which(colnames(xx1)%in%c("factor(urban_rural)0"))])
    xx <- as.matrix(xx)
    
    
    
    #Model
    model <-  glmnet::glmnet(y=y, x=xx, family = "binomial", alpha = 1, lambda = best_lambda, intercept = T)
    
    #test_data1 <- test_data %>% 
    #  dplyr::select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
    #               educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
    #               n2_Mean,media.prop,n12_Mean, urban_rural) %>% 
    # as.matrix()
    
    xx2 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                          educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                          n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = test_data)
    
    xxa <-  as.data.frame(xx2[,-which(colnames(xx2)%in%c("factor(urban_rural)0"))])
    test_data1 <- as.matrix(xxa)
    
    #Predictions
    predictions <- glmnet::predict.glmnet(model, newx = test_data1)
    
    predictions <- as.numeric(predictions) 
    
    predictions <- inverse_logit(predictions) %>% 
      as_tibble() %>% 
      rename(lasso_predictions = value) 
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$lasso_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$lasso_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$lasso_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$lasso_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, lasso_predictions = predictions$lasso_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Lasso_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}

# Apply function
result1 <- kfold_cv(data = dtp1)

result1 


# Spatial K-Fold Cross Validation -----------------------------------------

# function to calculate spatial k-fold

spatial_kfold_cv <- function(data) {
  fold <- data$Group_ID
  k <- length(unique(fold))
  
  rmse_values <- numeric(k) # Placeholder for RMSE
  pearson_values <- numeric(k) # Placeholder for corr
  mae_values <- numeric(k)  # Placeholder for MAE
  bias_values <- numeric(k)  # Placeholder for bias
  
  # Initialize data frame to store prop and predictions outside the loop
  results_df <- data.frame(DHSCLUST = numeric(), folds = numeric(), prop = numeric(), lasso_predictions = numeric())
  
  for (i in 1:k) {
    test_indices <- which(fold == i)
    train_indices <- which(fold != i)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    y <- train_data %>% 
      dplyr::select(success = Age12_23_Vax, Age12_23_Tot) %>%   
      mutate(failure = Age12_23_Tot - success) %>% 
      select(failure, success) %>% 
      as.matrix()
    
    #x <- train_data %>% 
    # dplyr::select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
    #               educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
    #               n2_Mean,media.prop,n12_Mean, urban_rural) %>% 
    # as.matrix()
    
    xx1 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                          educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                          n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = train_data)
    
    xx <-  as.data.frame(xx1[,-which(colnames(xx1)%in%c("factor(urban_rural)0"))])
    xx <- as.matrix(xx)
    
    #Model
    model <-  glmnet::glmnet(y=y, x=xx, family = "binomial", alpha = 1, lambda = best_lambda, intercept = T)
    
    #test_data1 <- test_data %>% 
    #  dplyr::select(n26_Mean, n28_Mean, n51_Mean, n58_Mean, n69_Mean, wealth.prop,
    #                educ.prop, health_card_doc.prop, sba.prop, phone_internet.prop, 
    #                n2_Mean,media.prop,n12_Mean, urban_rural) %>% 
    # as.matrix()
    
    xx2 <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                          educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                          n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = test_data)
    
    xxa <-  as.data.frame(xx2[,-which(colnames(xx2)%in%c("factor(urban_rural)0"))])
    test_data1 <- as.matrix(xxa)
    
    #Predictions
    predictions <- glmnet::predict.glmnet(model, newx = test_data1)
    
    predictions <- as.numeric(predictions) 
    
    predictions <- inverse_logit(predictions) %>% 
      as_tibble() %>% 
      rename(lasso_predictions = value) 
    
    rmse_values[i] <- sqrt(mean((test_data$prop - predictions$lasso_predictions)^2))
    pearson_values[i] <- cor(test_data$prop, predictions$lasso_predictions)
    mae_values[i] <- mean(abs(test_data$prop - predictions$lasso_predictions))
    bias_values[i] <- mean(test_data$prop - predictions$lasso_predictions)
    
    # Append prop and predictions to results_df with fold information
    fold_results <- data.frame(DHSCLUST = test_data$DHSCLUST, folds = i, prop = test_data$prop, lasso_predictions = predictions$lasso_predictions)
    results_df <- rbind(results_df, fold_results)
  }
  
  # Write results to CSV
  write.csv(results_df, paste0(output_path, "Lasso_Spatial_Kfold_Predictions.csv"))
  
  list(rmse = mean(rmse_values),
       corr = mean(pearson_values),
       mae = mean(mae_values),
       bias = mean(bias_values))
}

# Apply function
result2 <- spatial_kfold_cv(data = dtp1)
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
         urban_rural = NGA_urban_rural_1km) #%>% 
#as.matrix()


xxp <- model.matrix(~ -1 + n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean + wealth.prop +
                      educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                      n2_Mean + media.prop + n12_Mean + factor(urban_rural), data = pred_covs)

xxp <- xxp[match(rownames(pred_covs),rownames(xxp)),] #To match the no of rows to pred_covs
rownames(xxp) <- rownames(pred_covs)

xxb <-  as.data.frame(xxp[,-which(colnames(xxp)%in%c("factor(urban_rural)0", "factor(urban_rural)NaN"))])
pred_covs <- as.matrix(xxb)


#Make predictions using the fitted model
raster_predictions <- glmnet::predict.glmnet(fit_lasso, newx = pred_covs)
raster_predictions <- as.numeric(raster_predictions) 

#Apply function
#raster_predictions <- inverse_logit(raster_predictions) %>% 
#  as_tibble() %>% 
#  rename(predictions = value)

raster_predictions <- inverse_logit(raster_predictions)

# Rasterize predictions ---------------------------------------------------

#Read raster and get xy values
#r1 <- rast(paste0(data_path, "n28_viirs_2016_nightlights.tif"))

#library(raster)
r1 <- raster(paste0(data_path, "n28_viirs_2016_nightlights.tif"))
values(r1) <- as.numeric(raster_predictions)
plot(r1)

writeRaster(r1 , paste0(output_path, "lasso.tif"), 
            overwrite = T, names = "lasso_predictions")


# Get the xy coordinate of the centroid of each pixel as a dataframe
#coord <- xyFromCell(r1, 1:ncell(r1))

#cbind coordinates to predictions
#raster_predictions <- cbind(raster_predictions, coord)

# convert the dataframe to an sf object
#raster_predictions1 <- st_as_sf(raster_predictions, coords = c("x", "y"))

# set the CRS of the sf object
#st_crs(raster_predictions1) <- 4326

# set the CRS of the sf object
#raster_predictions1 <- st_transform(raster_predictions1, crs = st_crs(r1))

#Rasterize
#lasso_raster <- rasterize(raster_predictions1, r1, field = "predictions")

#plot(lasso_raster)

#writeRaster(lasso_raster , paste0(output_path, "lasso.tif"), 
#            overwrite = T, names = "lasso_predictions")

#####################################################################################
###################################################################################
output_path <- paste0(drive_path, "Methods_comparison_paper/.../DTP1/JackKnife/")


dtp1_full <- dtp1  #NOTE - Run this once before starting the jackknife part

# JackKnife Predictions to Get Uncertainty Calculations -------------------

### Start of jackknife section: Analyzing with different replicates/block sizes

# Define number of replicates (iterations)
nrast <- 100  # Number of replicates to run (e.g., 10, 30, 50, 100, 200)

# Calculate block size based on replicates and data size
bsize <- floor(nrow(dtp1_full) / nrast)  # Number of data points per replicate

# Set random seed for reproducibility
srand <- sample(1:nrow(dtp1_full), nrow(dtp1_full), replace = FALSE)

# Loop through each replicate
for (kk in 1:nrast) {
  
  # Define data points to exclude for this replicate
  if (kk < nrast) {
    qq <- (((kk - 1) * bsize) + 1):(bsize * kk)
  } else {
    qq <- (((kk - 1) * bsize) + 1):nrow(dtp1_full)
  }
  samp.k <- srand[qq]  # Data points to be excluded
  
  # Create a subset of the full data excluding the current replicate's points
  dtp1 <- dtp1_full[-samp.k, ]
  
  # Print the number of data points in the current subset
  print(nrow(dtp1))
  
  # Select response variables
  y <- dtp1 %>%
    dplyr::select(success = Age12_23_Vax, Age12_23_Tot) %>%  # Success and total counts
    mutate(failure = Age12_23_Tot - success) %>%           # Calculate failures
    dplyr::select(failure, success) %>%                   # Select only needed columns
    as.matrix()                                            # Convert to matrix
  
  ################## Create model matrix for explanatory variables
  xx1 <- model.matrix(~ -1 +  # Intercept
                        n26_Mean + n28_Mean + n51_Mean + n58_Mean + n69_Mean +  # Add your explanatory variables here
                        wealth.prop + educ.prop + health_card_doc.prop + sba.prop + phone_internet.prop +
                        n2_Mean + media.prop + n12_Mean + factor(urban_rural),  # Add your explanatory variables here
                      data = dtp1)
  
  # Exclude the intercept column from the model matrix
  xx <- as.data.frame(xx1[,-which(colnames(xx1) %in% c("factor(urban_rural)0"))])
  xx <- as.matrix(xx)
  
  # Perform cross-validation to find the optimal regularization parameter (lambda)
  cvfit <- cv.glmnet(x = xx, y = y, alpha = 1, family = "binomial")
  best_lambda <- cvfit$lambda.min  # Get the optimal lambda
  
  # Fit the Lasso model with the optimal lambda
  fit_lasso <- glmnet::glmnet(x = xx, y = y, family = "binomial", alpha = 1, lambda = best_lambda, intercept = T)
  
  # Find coefficient estimates (optional, commented out)
  # coef(fit_lasso)
  
  
  # Make predictions using the fitted model
  raster_predictions <- glmnet::predict.glmnet(fit_lasso, newx = pred_covs)  # Replace pred_covs with your prediction data
  raster_predictions <- as.numeric(raster_predictions)
  
  # Transform predictions to probability
  raster_predictions <- inverse_logit(raster_predictions)
  
  # Rasterize predictions (assumed code for reading and writing raster data exists)
  # Read raster and get xy values
  # r1 <- rast(paste0(data_path, "n28_viirs_2016_nightlights.tif"))
  
  library(raster)
  r1 <- raster(paste0(data_path, "n28_viirs_2016_nightlights.tif"))
  values(r1) <- as.numeric(raster_predictions)

  writeRaster(r1 , paste0(output_path, "jack_lasso_", kk, ".tif"), 
              overwrite = T)
}

#Calculate jackknife mean and std. dev.
raster_list <-list.files(path=output_path,pattern="^jack_lasso_.*\\.tif$" , all.files=TRUE, full.names=FALSE)
raster_list

stack_rasters<- rast(paste0(output_path, c(raster_list)))

#Jackknife mean
mean.rast <- app(stack_rasters, mean, na.rm=TRUE)
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

sd.rast <- raster(sd.rast)
summary(getValues(sd.rast))

writeRaster(mean.rast, paste0(output_path, "LASSO_mean_raster_dtp1_", nrast,".tif"), 
            overwrite = T)
writeRaster(sd.rast, paste0(output_path, "LASSO_stddev_raster_dtp1_", nrast,".tif"), 
            overwrite = T)

##################### END OF SCRIPT ####################################################
########################################################################################



