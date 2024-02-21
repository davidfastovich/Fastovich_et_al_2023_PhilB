################################################################################
# THIS SCRIPT DEMONSTREATES THE VARIABILITY BETWEEN THE ALTERNATE PALEOCLIMATE
# ESTIMATES USED IN FASTOVICH ET AL., (2023) WITH CORRELOGRAMS.
#
#   LAND MASK
#     * MODERN_DATA/POSTPROCESSED/GLDAS_50KM_ALBERS.NC
# 
#   PALEOPRECIPITATION
#     * MODEL_PRECIPITATION_ANOMALY/CCSM-MARUM_PR_ANOM_50KM_ENA.NC
#     * MODEL_PRECIPITATION_ANOMALY/CM2MC_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/COSMOS-W_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/IPSL_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/MIROC-W_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/TRACE_MWF_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/CCSM-NCAR_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/COSMOS-S_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/HADCM3_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/MIROC-S_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/TRACE_LORENZ_PR_ANOM.NC
#     * MODEL_PRECIPITATION_ANOMALY/TRACE_PR_ANOM.NC
# 
#   PALEOTEMPERATURE
#     * MODEL_TEMPERATURE_ANOMALY/CCSM-MARUM_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/CM2MC_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/COSMOS-W_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/IPSL_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/MIROC-W_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/TRACE_MWF_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/CCSM-NCAR_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/COSMOS-S_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/HADCM3_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/MIROC-S_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/TRACE_LORENZ_TAS_ANOM.NC
#     * MODEL_TEMPERATURE_ANOMALY/TRACE_TAS_ANOM.NC
# 
# OUTPUT:
#
#   FIGURE SX
#     * FIGURES/CLIMATE_MODEL_CORRELOGRAM.PDF
#     * FIGURES/CLIMATE_MODEL_CORRELOGRAM.PNG
#     
# ANACONDA R ENVIRONMENT WITH:
#
# R-BASE: 4.1.0
# R-TIDYVERSE: 1.3.1
#
# ALL PACKAGES WITH AN `R-` PREFIX WERE INSTALLED FROM THE CONDA-FORGE CHANNEL
#
# AUTHOR: DAVID FASTOVICH
# CONCTACT: FASTOVICH@WISC.EDU
################################################################################

#####################
# PREPARE ENVIRONMENT
#####################

library(raster)
library(tidyverse)
library(parallel)

# Set working directory to project
setwd("~/Fastovich_et_al_2023_PhilB")

########################
# LOAD IN NECESSARY DATA
########################

# Land sea mask
mask <- raster("modern_data/postprocessed/GLDAS_50km_albers.nc")

# Modern temperature correlate
modern_tave <- raster("modern_data/postprocessed/wc_tave_annual_50km_ena.tif")

# Modern precipitation correlate
modern_pr <- raster("modern_data/postprocessed/wc_pr_annual_50km_ena.tif")

# Loop through precipitation files and perform climate velocity calculation
tas_files <- list.files("model_temperature_anomaly_50km_ena/")
pr_files <- list.files("model_precipitation_anomaly_50km_ena/")

model_tave <- sapply(tas_files, function(x)
  raster(
    paste0("model_temperature_anomaly_50km_ena/", x),
    varname = "tas"
  )
)

model_pr <- sapply(pr_files, function(x)
  raster(
    paste0("model_precipitation_anomaly_50km_ena/", x),
    varname = "pr"
  )
)

# Apply land mask to each of the raster
model_tave <- lapply(
  model_tave, function(x)
    raster::mask(
      x,
      mask = mask,
      maskvalue = 0
    )
)

model_pr <- lapply(
  model_pr, function(x)
    raster::mask(
      x,
      mask = mask,
      maskvalue = 0
    )
)

# Crop to the study region
model_tave <- lapply(
  model_tave, function(x)
    raster::crop(
      x,
      extent(320000, 2850000, -1550000, 1150000)
    )
)

model_pr <- lapply(
  model_pr, function(x)
    raster::crop(
      x,
      extent(320000, 2850000, -1550000, 1150000)
    )
)

#######################
# CALCULATE CORRELOGRAM
#######################

model_tave_correl <- mclapply(
  model_tave,
  function (x) {
    df = as.data.frame(rasterToPoints(x))
    # return(pgirmess(x = df[,1], y = df[,2], z = df[,3], increment = 100*1000, resamp = 100))
    return(pgirmess::correlog(coords = df[,1:2], z = df[,3]))
  },
  mc.cores = 13
)

model_pr_correl <- mclapply(
  model_pr,
  function (x) {
    df = as.data.frame(rasterToPoints(x))
    # return(pgirmess(x = df[,1], y = df[,2], z = df[,3], increment = 100*1000, resamp = 100))
    return(pgirmess::correlog(coords = df[,1:2], z = df[,3]))
  },
  mc.cores = 13
)

##############################################
# EXTRACT RESULTS INTO DATA FRAME FOR PLOTTING
##############################################

# TAVE -------------------------------------------------------------------------
# Create data frame of relevant correlogram results
model_tave_corel_df <- lapply(
  model_tave_correl,
  as.data.frame
)

# Tidy up the data frame
model_tave_corel_df <- bind_rows(model_tave_corel_df, .id = 'Model')
model_tave_corel_df$Model <- gsub(
  "_tas_anom_50km_ena.nc",
  "",
  model_tave_corel_df$Model
)

# Rename models to make them appear cleaner for the ones that have informal
# names
model_tave_corel_df$Model[model_tave_corel_df$Model == "PROXY_KRIGING"] <- "Proxy Kriging"
model_tave_corel_df$Model[model_tave_corel_df$Model == "TRACE_LORENZ"] <- "TraCE-21ka\n(Statistcally Downscaled)"
model_tave_corel_df$Model[model_tave_corel_df$Model == "TRACE_MWF"] <- "TraCE-MWF (Single Forcing)"
model_tave_corel_df$Model[model_tave_corel_df$Model == "TRACE"] <- "TraCE-21ka"

# Add ID
model_tave_corel_df$clim_var <- "Temperature Anomaly"

# PR ---------------------------------------------------------------------------
# Create data frame of relevant correlogram results
model_pr_corel_df <- lapply(
  model_pr_correl,
  as.data.frame
)

# Tidy up the data frame
model_pr_corel_df <- bind_rows(model_pr_corel_df, .id = 'Model')
model_pr_corel_df$Model <- gsub(
  "_pr_anom_50km_ena.nc",
  "",
  model_pr_corel_df$Model
)

# Rename models to make them appear cleaner for the ones that have informal
# names
model_pr_corel_df$Model[model_pr_corel_df$Model == "PROXY_KRIGING"] <- "Proxy Kriging"
model_pr_corel_df$Model[model_pr_corel_df$Model == "TRACE_LORENZ"] <- "TraCE-21ka\n(Statistcally Downscaled)"
model_pr_corel_df$Model[model_pr_corel_df$Model == "TRACE_MWF"] <- "TraCE-MWF (Single Forcing)"
model_pr_corel_df$Model[model_pr_corel_df$Model == "TRACE"] <- "TraCE-21ka"

# Add ID
model_pr_corel_df$clim_var <- "Precipitation Anomaly"

# Combine PR and TAVE data frames ----------------------------------------------
model_corel <- rbind(model_tave_corel_df, model_pr_corel_df)

# Make significance a factor
model_corel$Significance <- ifelse(model_corel$p.value > 0.05, "Not Significant", "Significant")

######
# PLOT
######

# Setting theme
theme <- 
  theme_bw() + 
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.key = element_blank(),
        panel.spacing.y = unit(-.1, "cm"))
theme_set(theme)

corel_plot <- ggplot(data = model_corel, aes(x = dist.class/1000, y = coef, shape = Significance, color = Model, group = Model)) + 
  geom_point() + 
  geom_line(size = 0.25) + 
  scale_shape_manual(values = c(1, 19)) +
  xlab("Distance (km)") +
  ylab("Moran's I") + 
  facet_wrap(~ clim_var, ncol = 2, scales = "free_y")

ggsave(
  filename = "figures/climate_model_correlogram.pdf",
  corel_plot,
  dpi = 300,
  height = 4.5,
  width = 8.5
)

ggsave(
  filename = "figures/climate_model_correlogram.png",
  corel_plot,
  dpi = 300,
  height = 4.5,
  width = 8.5
)
