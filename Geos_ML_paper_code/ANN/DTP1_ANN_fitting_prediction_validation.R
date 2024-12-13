#Load libraries
library(caret)
library(ggplot2)
library(gbm)
library(gridExtra)
library(groupdata2) # To create Group_ID
library(secr)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(viridis)

library(h2o)
library(data.table)
library(terra)
library(DMwR)

filePathData <- "Data/"
filePathData1 <- "Data/Raster Files/"
filePathData2 <- "Output/"

# Import data
dtp1 <- read.csv(paste0(filePathData,"NGA_Vax_DTP1.csv"), header=TRUE)
covs <- read.csv(paste0(filePathData,"Covariates_merged_selected.csv"),header=TRUE)
grp <- read.csv(paste0(filePathData, "NGA_Groups.csv"))


#select covs and add to response variable
dtp1 <- dtp1 %>% 
  select(Age12_23_Vax, Age12_23_Tot, DHSCLUST) %>% 
  inner_join(covs, by = "DHSCLUST") %>% 
  inner_join(grp, by = "DHSCLUST") %>%
  # filter(Age12_23_Tot != 0) %>%    # Keep values where total is not Zero 
  filter(Age12_23_Tot > 1) %>%    # Keep values where total is greater than 1
  mutate(prop = (Age12_23_Vax/Age12_23_Tot), logit_prop = logit(abs(prop-0.000001)))

#Create a grouping variable sequentially to be used for spatial k-fold cross val
dtp1 <- dtp1 %>% 
  group(n = 116, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()

data <- dtp1
names(data)

# Initialize the H2O cluster
h2o.init()

# Convert the data to an H2O frame
df <- data[,c(3, 6:7, 1:2, 8:20, 5, 23:24, 25)]
names(df)
data_h2o <- as.h2o(df)

# Define the response variable and predictor variables
response <- "logit_prop"
predictors <- c(
  "n2_Mean", "n12_Mean",
  "n26_Mean", "n28_Mean", 
  "n51_Mean", "n58_Mean", 
  "n69_Mean", "wealth.prop", 
  "educ.prop", "health_card_doc.prop", 
  "sba.prop",  "media.prop",
  "phone_internet.prop", "urban1_rural0",
  "LATNUM", "LONGNUM"
)

# Split the data into training and testing sets
splits <- h2o.splitFrame(data_h2o, ratios = 0.99, seed = 1234)
train <- splits[[1]]
test <- splits[[2]]

# Define the neural network model
nn_model <- h2o.deeplearning(
  x = predictors,
  y = response,
  training_frame = train,
  validation_frame = test,
  hidden = c(100, 100),       # Two hidden layers with 100 neurons each
  epochs = 100,               # Number of epochs
  distribution = "gaussian",
  activation = "Rectifier",        # Activation function
  stopping_metric = "RMSE",   # Metric to evaluate the model performance
  stopping_tolerance = 0.001, # Stopping tolerance
  stopping_rounds = 5,        # Stopping rounds
  seed = 1234                 # Seed for reproducibility
)

# Print the model summary
summary(nn_model)

# Predict on the test set
predictions <- h2o.predict(nn_model, test) 

predictions <- as.numeric(predictions) 

#Inverse logit function
inverse_logit <- function(aa) {
  out = exp(aa)/(1+exp(aa))
  return(out)
  
}

predictions <- inverse_logit(predictions) %>%
  as_tibble()

# Evaluate the model performance
performance <- h2o.performance(nn_model, test)
print(performance)

# Convert predictions and actual values to data frame for plotting
pred_df <- as.data.frame(predictions)
actual_df <- as.data.frame(test[20])

# Combine predictions and actual values
result_df <- data.frame(Actual = actual_df$prop, Predicted = pred_df$exp.predict.)
summary(result_df)


# Plot predictions vs. actual values
ggplot(result_df, aes(x = Actual, y = Predicted)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  ggtitle("Predictions vs. Actual Values") +
  xlab("Actual Values") +
  ylab("Predicted Values") +
  xlim(0, max(result_df$Actual, result_df$Predicted)) +
  ylim(0, max(result_df$Actual, result_df$Predicted)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()

# Calculate evaluation metrics
h2o_metrics <- list(
  RMSE = sqrt(mean((result_df$Predicted - result_df$Actual)^2)),
  Correlation = cor(result_df$Predicted, result_df$Actual),
  MAE = mean(abs(result_df$Predicted - result_df$Actual)),
  Bias = mean(result_df$Predicted - result_df$Actual)
)
h2o_metrics


#===============================================================================
# In-Sample Prediction
#===============================================================================
in_sample_pred <- h2o.predict(nn_model, data_h2o)

in_sample_predictions <- inverse_logit(in_sample_pred) %>%
  as_tibble() %>% 
  rename(nn_predictions = exp.predict.)

# Cbind cluster id and response variables
nn_insample_df <- as.data.frame(data_h2o) %>% 
  select(DHSCLUST, LATNUM, LONGNUM, Age12_23_Vax, Age12_23_Tot, prop) %>% 
  mutate(nn_predictions = in_sample_predictions$nn_predictions)

summary(nn_insample_df[,6:7])

write.csv(nn_insample_df, paste0(filePathData2, "dtp1_nn_insample_result.csv"))

# Plot predictions vs. actual values
ggplot(nn_insample_df, aes(x = prop, y = nn_predictions)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  ggtitle("Predictions vs. Actual Values") +
  xlab("Actual Values") +
  ylab("Predicted Values") +
  xlim(0, max(nn_insample_df$prop, nn_insample_df$nn_predictions)) +
  ylim(0, max(nn_insample_df$prop, nn_insample_df$nn_predictions)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()

# Calculate evaluation metrics
in_sample_metrics <- list(
  RMSE = sqrt(mean((nn_insample_df$nn_predictions - nn_insample_df$prop)^2)),
  Correlation = cor(nn_insample_df$nn_predictions, nn_insample_df$prop),
  MAE = mean(abs(nn_insample_df$nn_predictions - nn_insample_df$prop)),
  Bias = mean(nn_insample_df$nn_predictions - nn_insample_df$prop)
)

in_sample_metrics

#===============================================================================
# Random K-Fold Cross Validation
#===============================================================================

# Function to perform k-fold cross-validation with H2O and save predictions

kfold_rand_h2o <- function(data = data_h2o, predictors, response, k = 10) {
  
  n <- nrow(data_h2o)
  fold_size <- n %/% k
  folds <- rep(1:k, length.out = n)
  
  # Initialize vectors to store metrics
  rmse_values <- numeric(k)
  pearson_values <- numeric(k)
  mae_values <- numeric(k)
  bias_values <- numeric(k)
  
  # Data frame to store predictions and actual values
  prediction_data <- data.frame(Fold = integer(), DHSCLUST = integer(), Actual = numeric(), Predicted = numeric(), Prop = numeric())
  
  for (i in 1:k) {
    
    test_indices <- which(folds == i)
    train_indices <- which(folds != i)
    
    # Filter data for train and test sets based on fold indices
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Check if validation frame has > 0 rows
    if (nrow(test_data) > 0) {
      # Define the neural network model
      nn_model <- h2o.deeplearning(
        x = predictors,
        y = response,
        training_frame = train_data,
        validation_frame = test_data,
        hidden = c(100, 100),         # Two hidden layers with 100 neurons each
        epochs = 100,                 # Number of epochs
        distribution = "gaussian",
        activation = "Rectifier",     # Activation function
        stopping_metric = "RMSE",     # Metric to evaluate the model performance
        stopping_tolerance = 0.001,   # Stopping tolerance
        seed = 1234                   # Seed for reproducibility
      )
      print(i)
      
      # Predict on the test set
      predictions <- h2o.predict(nn_model, test_data)
      predictions <- as.numeric(predictions)
      print(predictions)
      
      # Apply inverse logit function
      kfold_predictions <- inverse_logit(predictions) %>% 
        as_tibble() %>% 
        rename(nn_predictions = exp.predict.)
      
      # Retrieve actual values and prop values from the test set
      DHSCLUST <- as.vector(test_data[[1]])
      actual <- as.vector(test_data[[20]])
      prop_values <- as.vector(test_data$prop)
      predicted <- as.vector(kfold_predictions$nn_predictions)
      
      # Store predictions, actual values, and prop values for the fold
      fold_data <- data.frame(Fold = rep(i, length(actual)), DHSCLUST = DHSCLUST, Actual = actual, Predicted = predicted, Prop = prop_values)
      prediction_data <- rbind(prediction_data, fold_data)
      
      # Compute evaluation metrics for the fold
      rmse_values[i] <- sqrt(mean((actual - predicted)^2))
      pearson_values[i] <- cor(actual, predicted)
      mae_values[i] <- mean(abs(actual - predicted))
      bias_values[i] <- mean(actual - predicted)
    } else {
      # If test set is empty, skip this fold
      cat("Skipping fold", i, "due to empty test set.\n")
    }
  }
  
  # Compute mean metrics across all non-empty folds
  valid_folds <- which(rmse_values > 0)  # Exclude folds with empty test sets
  metrics <- list(
    RMSE = mean(rmse_values[valid_folds]),
    Correlation = mean(pearson_values[valid_folds]),
    MAE = mean(mae_values[valid_folds]),
    Bias = mean(bias_values[valid_folds])
  )
  
  # Return metrics and prediction data
  return(list(metrics = metrics, prediction_data = prediction_data))
}

# Perform k-fold cross-validation
result_rand <- kfold_rand_h2o(data_h2o, predictors, response)

# Extract metrics and prediction data
# rand_metrics <- result_rand$metrics
rand_prediction_data <- result_rand$prediction_data


# Save results in a list
Correlation <- cor(rand_prediction_data$Prop, rand_prediction_data$Predicted)
RMSE <- sqrt(mean((rand_prediction_data$Prop - rand_prediction_data$Predicted)^2))
MAE <- mean(abs(rand_prediction_data$Prop - rand_prediction_data$Predicted))
Bias <- mean(rand_prediction_data$Prop - rand_prediction_data$Predicted)

rand_metrics <- data.frame(RMSE, Correlation, MAE, Bias)

# Print cross-validation metrics
print("Cross-Validation Metrics:")
print(rand_metrics)


kfold_rand_df <- as.data.frame(data_h2o) %>%
  select(DHSCLUST, LATNUM, LONGNUM, Age12_23_Vax, Age12_23_Tot, prop) %>%
  mutate(Fold = rand_prediction_data$Fold,
         Actual = rand_prediction_data$Prop,
         Predicted = rand_prediction_data$Predicted)

# summary(kfold_rand_df[,2:3])
write.csv(rand_prediction_data, paste0(filePathData2, "dtp1_nn_rand_kfold_results.csv"))



# Create ggplot comparing observed and predicted values for n folds
ggplot(kfold_rand_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  facet_wrap(~ Fold, scales = "free") +
  ggtitle("Observed vs Predicted Values for Each Fold") +
  xlab("Observed Values") +
  ylab("Predicted Values") +
  xlim(0, max(kfold_rand_df$Actual, kfold_rand_df$Predicted)) +
  ylim(0, max(kfold_rand_df$Actual, kfold_rand_df$Predicted)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()


# Create ggplot comparing observed and predicted values for all folds combined
ggplot(kfold_rand_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  ggtitle("Observed vs Predicted Values (All Folds Combined)") +
  xlab("Observed Values") +
  ylab("Predicted Values") +
  xlim(0, max(kfold_rand_df$Actual, kfold_rand_df$Predicted)) +
  ylim(0, max(kfold_rand_df$Actual, kfold_rand_df$Predicted)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()


# Summary of all metrics
summary_metrics <- data.frame(
  RMSE = c(h2o_metrics$RMSE, in_sample_metrics$RMSE, rand_metrics$RMSE),
  Correlation = c(h2o_metrics$Correlation, in_sample_metrics$Correlation, rand_metrics$Correlation),
  MAE = c(h2o_metrics$MAE, in_sample_metrics$MAE, rand_metrics$MAE),
  Bias = c(h2o_metrics$Bias, in_sample_metrics$Bias, rand_metrics$Bias)
)

# Add row names
row.names(summary_metrics) <- c("Test", "In_Sample", "Rand_K-fold")

# Print the summary metrics
print(summary_metrics)



#===============================================================================
# Spatial K-Fold Cross Validation
#===============================================================================

# Function to perform spatial k-fold cross-validation with H2O and save predictions

kfold_spatial_h2o <- function(data = data_h2o, predictors, response) {
  
  # Extract fold assignments from Group_ID
  fold <- as.numeric(as.vector(data$Group_ID))
  k <- length(unique(fold))  # Number of unique folds
  
  # Initialize vectors to store metrics
  rmse_values <- numeric(k)
  pearson_values <- numeric(k)
  mae_values <- numeric(k)
  bias_values <- numeric(k)
  
  # Data frame to store predictions and actual values
  prediction_data <- data.frame(Fold = integer(), DHSCLUST = integer(), Actual = numeric(), Predicted = numeric(), Prop = numeric())
  
  # Iterate over each unique fold
  for (i in seq_len(k)) {
    
    # Identify test and train indices based on the fold assignment
    test_indices <- which(fold == i)
    train_indices <- which(fold != i)
    
    # Filter data for train and test sets based on fold indices
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Check if validation frame has > 0 rows
    if (nrow(test_data) > 0) {
      # Define the neural network model
      nn_model <- h2o.deeplearning(
        x = predictors,
        y = response,
        training_frame = train_data,
        validation_frame = test_data,
        hidden = c(100, 100),         # Two hidden layers with 100 neurons each
        epochs = 100,                 # Number of epochs
        distribution = "gaussian",
        activation = "Rectifier",     # Activation function
        stopping_metric = "RMSE",     # Metric to evaluate the model performance
        stopping_tolerance = 0.001,   # Stopping tolerance
        seed = 1234                   # Seed for reproducibility
      )
      print(paste("Processing Fold:", i))
      
      # Predict on the test set
      predictions <- h2o.predict(nn_model, test_data)
      predictions <- as.numeric(predictions)
      print(predictions)
      
      # Apply inverse logit function
      kfold_predictions <- inverse_logit(predictions) %>% 
        as_tibble() %>% 
        rename(nn_predictions = exp.predict.)
      
      # Retrieve actual values and prop values from the test set
      DHSCLUST <- as.vector(test_data$DHSCLUST)
      actual <- as.vector(test_data$logit_prop)
      prop_values <- as.vector(test_data$prop)
      predicted <- as.vector(kfold_predictions$nn_predictions)
      
      # Store predictions, actual values, and prop values for the fold
      fold_data <- data.frame(Fold = rep(i, length(actual)), DHSCLUST = DHSCLUST, Actual = actual, Predicted = predicted, Prop = prop_values)
      prediction_data <- rbind(prediction_data, fold_data)
      
      # Compute evaluation metrics for the fold
      rmse_values[i] <- sqrt(mean((actual - predicted)^2))
      pearson_values[i] <- cor(actual, predicted)
      mae_values[i] <- mean(abs(actual - predicted))
      bias_values[i] <- mean(actual - predicted)
    } else {
      # If test set is empty, skip this fold
      cat("Skipping fold", i, "due to empty test set.\n")
    }
  }
  
  # Compute mean metrics across all non-empty folds
  valid_folds <- which(rmse_values > 0)  # Exclude folds with empty test sets
  metrics <- list(
    RMSE = mean(rmse_values[valid_folds]),
    Correlation = mean(pearson_values[valid_folds]),
    MAE = mean(mae_values[valid_folds]),
    Bias = mean(bias_values[valid_folds])
  )
  
  # Return metrics and prediction data
  return(list(metrics = metrics, prediction_data = prediction_data))
}

# Perform spatial k-fold cross-validation using Group_ID
result_spatial <- kfold_spatial_h2o(data_h2o, predictors, response)


# Extract metrics and prediction data
# spatial_metrics <- result_spatial$metrics
spatial_prediction_data <- result_spatial$prediction_data


# Save results in a list
Correlation <- cor(spatial_prediction_data$Prop, spatial_prediction_data$Predicted)
RMSE <- sqrt(mean((spatial_prediction_data$Prop - spatial_prediction_data$Predicted)^2))
MAE <- mean(abs(spatial_prediction_data$Prop - spatial_prediction_data$Predicted))
Bias <- mean(spatial_prediction_data$Prop - spatial_prediction_data$Predicted)

spatial_metrics <- data.frame(RMSE, Correlation, MAE, Bias)



# Print cross-validation metrics
print("Cross-Validation Metrics:")
print(spatial_metrics)


     
summary(spatial_prediction_data[,4:5])
write.csv(spatial_prediction_data, paste0(filePathData2, "dtp1_nn_spatial_kfold_results.csv"))


kfold_spatial_df <- as.data.frame(data_h2o) %>%
  select(DHSCLUST, LATNUM, LONGNUM, Age12_23_Vax, Age12_23_Tot) %>%
  mutate(Fold = spatial_prediction_data$Fold,
         Actual = spatial_prediction_data$Prop,
         Predicted = spatial_prediction_data$Predicted)




# Create ggplot comparing observed and predicted values for n folds
ggplot(kfold_spatial_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  facet_wrap(~ Fold, scales = "free") +
  ggtitle("Observed vs Predicted Values for Each Fold") +
  xlab("Observed Values") +
  ylab("Predicted Values") +
  xlim(0, max(kfold_spatial_df$Actual, kfold_spatial_df$Predicted)) +
  ylim(0, max(kfold_spatial_df$Actual, kfold_spatial_df$Predicted)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()


# Create ggplot comparing observed and predicted values for all folds combined
ggplot(kfold_spatial_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  ggtitle("Observed vs Predicted Values (All Folds Combined)") +
  xlab("Observed Values") +
  ylab("Predicted Values") +
  xlim(0, max(kfold_spatial_df$Actual, kfold_spatial_df$Predicted)) +
  ylim(0, max(kfold_spatial_df$Actual, kfold_spatial_df$Predicted)) +
  # xlim(0, 1) +
  # ylim(0, 1) +
  theme_minimal()


# Summary of all metrics
summary_metrics <- data.frame(
  RMSE = c(h2o_metrics$RMSE, in_sample_metrics$RMSE, rand_metrics$RMSE, spatial_metrics$RMSE),
  Correlation = c(h2o_metrics$Correlation, in_sample_metrics$Correlation, rand_metrics$Correlation, spatial_metrics$Correlation),
  MAE = c(h2o_metrics$MAE, in_sample_metrics$MAE, rand_metrics$MAE, spatial_metrics$MAE),
  Bias = c(h2o_metrics$Bias, in_sample_metrics$Bias, rand_metrics$Bias, spatial_metrics$Bias)
)

# Add row names
row.names(summary_metrics) <- c("Test", "In_Sample", "Rand_K-fold", "Spatial_K-fold")

# Print the summary metrics
print(summary_metrics)
write.csv(summary_metrics, paste0(filePathData2, "dtp1_summary_metrics.csv"))

#===============================================================================

# Define the response variable and predictor variables
response <- "logit_prop"
predictors <- c(
  "n2_Mean", "n12_Mean",
  "n26_Mean", "n28_Mean", 
  "n51_Mean", "n58_Mean", 
  "n69_Mean", "wealth.prop", 
  "educ.prop", "health_card_doc.prop", 
  "sba.prop",  "media.prop",
  "phone_internet.prop", "urban1_rural0",
  "LATNUM", "LONGNUM"
)

# Split the data into training and testing sets
splits <- h2o.splitFrame(data_h2o, ratios = 0.99, seed = 1234)
train <- splits[[1]]
test <- splits[[2]]

# Define the neural network model
nn_model_rast <- h2o.deeplearning(
  x = predictors,
  y = response,
  training_frame = train,
  validation_frame = test,
  hidden = c(100, 100),       # Two hidden layers with 100 neurons each
  epochs = 100,               # Number of epochs
  distribution = "gaussian",
  activation = "Rectifier",        # Activation function
  stopping_metric = "RMSE",   # Metric to evaluate the model performance
  stopping_tolerance = 0.001, # Stopping tolerance
  stopping_rounds = 5,        # Stopping rounds
  seed = 1234                 # Seed for reproducibility
)





# Load rasters for predictions
raster_list <- list.files(path = filePathData1, pattern = ".tif$", all.files = TRUE, full.names = FALSE)
raster_covs <- raster_list[-c(7, 12, 14, 17)]
stack_rasters <- rast(paste0(filePathData1, raster_covs))

# Extract raster values and create a dataframe
raster_df <- terra::values(stack_rasters, dataframe = TRUE)

# Rename columns to match covariate names
pred_covs <- raster_df %>%
  rename(
    n2_Mean = n2_cattle_FAO_2010,
    n12_Mean = n12_ccilc_dst011_2015,
    n26_Mean = n26_guf_ghsl_dst_2014,
    n28_Mean = n28_viirs_2016_nightlights,
    n51_Mean = n51_MODIS_lst,
    n58_Mean = n58_TravelTime_to_HF_ver2,
    n69_Mean = n69_WorldPop10povcons200,
    wealth.prop = wealthkrig,
    educ.prop = educkrig,
    health_card_doc.prop = health_card_dockrig,
    sba.prop = sbakrig, 
    media.prop = mediakrig,
    phone_internet.prop = phone_internetkrig,
    urban1_rural0 = NGA_urban_rural_1km
  )

# Convert covariate data to h2o frame
pred_covs_h2o <- as.h2o(pred_covs)

# Predict using the h2o Deep Learning model
raster_predictions <- h2o.predict(nn_model_rast, pred_covs_h2o) %>% as.data.frame()
raster_predictions <- raster_predictions$predict


# Apply inverse logit transformation
raster_predictions <- inverse_logit(raster_predictions)

# Load reference raster and extract coordinates
r1 <- rast(paste0(filePathData1, "n28_viirs_2016_nightlights.tif"))
coord <- xyFromCell(r1, 1:ncell(r1))


# Combine predictions with coordinates
raster_predictions <- cbind(raster_predictions, coord)

# Convert to sf object and set CRS
raster_predictions1 <- st_as_sf(as.data.frame(raster_predictions), coords = c("x", "y"))

# Update the values greater than 1 to 1
# raster_predictions1 <- raster_predictions1 %>%
#   mutate(raster_predictions = ifelse(raster_predictions > 1, 1, raster_predictions))

st_crs(raster_predictions1) <- 4326
raster_predictions1 <- st_transform(raster_predictions1, crs = st_crs(r1))

# Rasterize predictions
DL_raster <- rasterize(raster_predictions1, r1, field = "raster_predictions")

# Plot raster
plot(DL_raster)

# Crop raster
crop_raster <- terra::mask(DL_raster, r1)
plot(crop_raster)

# Define breaks for visualization
breaks <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1)

# Visualize using ggplot2
p2 <- ggplot() +
  geom_spatraster(data = crop_raster) +
  scale_fill_viridis_c(option = "D", breaks = breaks, 
                       guide = guide_legend(reverse = TRUE)) +
  labs(fill = "Proportion") +
  ggtitle("ANN (h2o)") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 16))

p2

# Interactive mapview
mapview::mapview(crop_raster)

# Save raster output
writeRaster(DL_raster, paste0(filePathData2, "dtp1_spatial.tif"), overwrite = TRUE, names = "dl_predictions")

# Shut down h2o
h2o.shutdown(prompt = FALSE)
