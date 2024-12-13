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

# Create factor for Urban-Rural Binary values
dtp1$urban1 <- as.factor(dtp1$urban1_rural0) ###### Additional Change

#Create a grouping variable sequentially to be used for spatial k-fold cross val
dtp1 <- dtp1 %>% 
  group(n = 134, method = "greedy", col_name = "Group_ID") %>% 
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
  "phone_internet.prop",
  "urban1" ###### Additional Change
  # "urban1_rural0"
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

# #cbind cluster id and response variables
# nn_data <- test %>% 
#   select(DHSCLUST, LATNUM, LONGNUM, Age12_23_Vax, Age12_23_Tot, prop) %>% 
#   cbind(predictions)
# 
# write.csv(gbm_data, paste0(output_path, "gbm_data.csv"))

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
# Raster Prediction
#===============================================================================
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
    urban1 = NGA_urban_rural_1km ###### Additional Change
    # urban1_rural0 = NGA_urban_rural_1km
  )

# Convert covariate data to h2o frame
pred_covs_h2o <- as.h2o(pred_covs)

# Predict using the h2o Deep Learning model
raster_predictions <- h2o.predict(nn_model, pred_covs_h2o) %>% as.data.frame()
raster_predictions <- raster_predictions$predict


# Apply inverse logit transformation
raster_predictions <- inverse_logit(raster_predictions)

# Load reference raster and extract coordinates
r1 <- rast(paste0(filePathData1, "n28_viirs_2016_nightlights.tif"))
coord <- xyFromCell(r1, 1:ncell(r1))



#===============================================================================
# Jackknife Delete-a-Spatial-Block Analysis for h2o Model
#===============================================================================

# Function for inverse logit transformation
inverse_logit <- function(aa) {
  exp(aa) / (1 + exp(aa))
}

dtp_full <- dtp1 #Run this line only once before starting the jackknife analysis

# Define block size and create random order of rows
nrast <- 10 # 10 blocks
bsize <- floor(nrow(dtp_full) / nrast) # Block size
srand <- sample(1:nrow(dtp_full), nrow(dtp_full), replace = FALSE) # Random order of rows

# Initialize an empty list to store predictions for each block
jackknife_preds <- list()

# Loop through blocks
for (kk in 1:nrast) {
  
  # Define which rows to remove in this iteration
  if (kk < nrast) {
    qq <- (((kk - 1) * bsize) + 1):(bsize * kk)
  } else {
    qq <- (((kk - 1) * bsize) + 1):nrow(dtp_full)
  }
  
  samp.k <- srand[qq]
  
  # Remove the selected block of observations from the dataset
  dtp1 <- dtp_full[-samp.k, ]
  print(paste("Block", kk, "size:", nrow(dtp1)))
  
  # Convert the block data to h2o frame
  dtp1_h2o_block <- as.h2o(dtp1[, c(3, 6:7, 1:2, 8:20, 5, 23:24, 25)])
  
  # Train the model on the reduced data (excluding the selected block)
  nn_model_block <- h2o.deeplearning(
    x = predictors,
    y = response,
    training_frame = dtp1_h2o_block,
    hidden = c(100, 100),
    epochs = 100,
    #distribution = "gaussian",
    activation = "Rectifier",
    stopping_metric = "RMSE",
    stopping_tolerance = 0.001,
    #stopping_rounds = 5,
    seed = 1234
  )
  
  # Predict using the raster covariates for this block
  raster_preds_block <- h2o.predict(nn_model_block, pred_covs_h2o) %>% as.data.frame()
  raster_preds_block <- inverse_logit(raster_preds_block$predict)
  
  # Create the raster and assign predictions
  block_raster <- raster(paste0(filePathData1, "n28_viirs_2016_nightlights.tif"))
  values(block_raster) <- as.numeric(raster_preds_block)
  
  # Save block-specific raster
  writeRaster(block_raster, paste0(filePathData2, "dtp1_jackknife_block_", kk, ".tif"), overwrite = TRUE)
}

#===============================================================================
# Calculate Jackknife Mean
#===============================================================================

# Load all saved block rasters
block_raster_list <- list.files(path = filePathData2, pattern = "dtp1_jackknife_block_.*.tif", full.names = TRUE)
stack_block_rasters <- rast(block_raster_list)

# Jackknife mean across blocks
jackknife_mean_raster <- app(stack_block_rasters, mean, na.rm = TRUE)
plot(jackknife_mean_raster)

# Crop raster
crop_raster_mean <- terra::mask(jackknife_mean_raster, r1)
plot(crop_raster_mean)

# Save the jackknife mean raster
writeRaster(crop_raster_mean, paste0(filePathData2, "jackknife_mean_dtp1_.tif"), overwrite = TRUE)

#===============================================================================
# Calculate Jackknife Standard Deviation
#===============================================================================

# # Jackknife standard deviation across blocks
# fun_sq <- function(i) i * i
# sum_diff_sq <- 0
# 
# for (i in 1:nrast) {
#   block_raster <- subset(stack_block_rasters, i)
#   diff_raster <- block_raster - jackknife_mean_raster
#   sum_diff_sq <- sum_diff_sq + app(diff_raster, fun_sq)
# }
# 
# jackknife_sd_raster <- sqrt(((nrast - 1) / nrast) * sum_diff_sq)
# plot(jackknife_sd_raster)
# 
# 
# # Check values
# sd.rast <- raster(jackknife_sd_raster)
# summary(raster::getValues(sd.rast))
# 
# plot(sd.rast)

#===============================================================================

# Step 1: Initialize a raster for storing the sum of squared deviations
squared_deviation_sum <- 0  

# Step 2: Loop over each jackknife raster and calculate the squared deviation from the overall mean
for (i in 1:length(nrast)) { #### nrast replace by stack_block_rasters
  if (i %in% c(2, 4, 8)){next} #### Block this also with the previous line update
  
  # Calculate deviation (difference from overall mean raster)
  deviation <- stack_block_rasters[[i]] - jackknife_mean_raster
  
  # Square the deviation
  squared_deviation <- deviation^2
  
  # Sum the squared deviations
  squared_deviation_sum <- squared_deviation_sum + squared_deviation
}

# Step 3: Calculate the mean of the squared deviations (divide by the number of jackknife rasters, nrast)
mean_squared_deviation <- squared_deviation_sum / 7 #### 7 replace by stack_block_rasters

# Step 4: Calculate the standard deviation (square root of the mean squared deviation)
jackknife_sd_raster <- sqrt(mean_squared_deviation)

#===============================================================================

# Step 5: Save or plot the SD raster
# Save the SD raster to a file
# writeRaster(jackknife_sd_raster, "jackknife_sd_raster_loop.tif", format = "GTiff", overwrite = TRUE)

# Plot the SD raster
plot(jackknife_sd_raster, zlim = c(0, 0.5), col = rev(terrain.colors(100)))

# Crop raster
crop_raster_sd <- terra::mask(jackknife_sd_raster, r1)
plot(crop_raster_sd, zlim = c(0, 0.5), col = rev(terrain.colors(100)))


sd.rast1 <- raster(squared_deviation_sum)
summary(raster::getValues(sd.rast1))


plot(crop_raster_sd)


# Save the jackknife standard deviation raster
writeRaster(crop_raster_sd, paste0(filePathData2, "jackknife_sd_dtp1_.tif"), overwrite = TRUE)

#===============================================================================
# Visualize Jackknife Results
#===============================================================================

# Define breaks for visualization
breaks <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1)

# Plot mean raster
p_mean <- ggplot() +
  geom_spatraster(data = crop_raster_mean) +
  scale_fill_viridis_c(option = "D", breaks = breaks, guide = guide_legend(reverse = TRUE)) +
  labs(fill = "Mean Proportion") +
  ggtitle("Jackknife Mean Predictions") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 16))

p_mean

# Plot sd raster
p_sd <- ggplot() +
  geom_spatraster(data = crop_raster_sd) +
  scale_fill_viridis_c(option = "D", breaks = breaks, guide = guide_legend(reverse = TRUE)) +
  labs(fill = "Standard Deviation") +
  ggtitle("Jackknife Standard Deviation") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 16))

p_sd


# save.image("jackknife_DTP1.RData")
# 
# load("jackknife_DTP1.RData")
