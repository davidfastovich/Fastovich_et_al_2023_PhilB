################################################################################
########### THESE ARE SPATIAL GRAIN SENSITIVITY TESTS ON A 100KM GRID ##########
################################################################################
# THIS SCRIPT PERFORMS THE STATISTICAL ANALYSES FOR FASTOVICH ET AL., (2023,
# PHIL B) AND GENERATES RELEVANT DIAGNOSTICS. WE USE SPATIAL ERROR MODELS TO
# BUILD RELATIONSHIPS BETWEEN SPECIES RICHNESS FOR AMPHIBIANS, BIRDS, MAMMALS,
# REPTILES, AND TREES AND MODERN CLIMATE AND PALEOCLIMATE CORRELATES. WE ASSESS
# VARIOUS NEIGHBORHOOD ADJACENCY MATRICES AND FIND THAT THE ROOKS CASE MOST
# EFFECTIVELY HANDLES RESIDUAL SPATIAL AUTOCORRELATION BUT CODE FOR ALL REMAINS.
# THE PRODUCTS OF THIS SCRIPT ARE SEVERAL OF THE FIGURES IN THE FINAL
# PUBLICATION AND DATA IN TABLE 1. ALL SPATIAL ERROR MODELS ARE OF THE FORM:
# 
# SPECIES RICHNESS ~ MODERN MEAN ANNUAL TEMPERATURE + MODERN MEAN ANNUAL + 
# PRECIPITATION + PALEOTEMPERATURE + PALEOPRECIPITATION + SPATIAL ERROR
# 
# INPUT:
#   RICHNESS ESTIMATES
#     * RESOLUTION_SENSITIVITY/100KM/RICHNESS/AMPH_RICHNESS_EPSG_102008.TIF
#     * RESOLUTION_SENSITIVITY/100KM/RICHNESS/BIRDS_PASSERIFORMES_RICHNESS_EPSG_102008.TIF
#     * RESOLUTION_SENSITIVITY/100KM/RICHNESS/MAMMS_RICHNESS_EPSG_102008.TIF
#     * RESOLUTION_SENSITIVITY/100KM/RICHNESS/REPS_RICHNESS_EPSG_102008.TIF
#     * RESOLUTION_SENSITIVITY/100KM/RICHNESS/TREE_RICHNESS_EPSG_102008.TIF
# 
#   CONTEMPORARY MEAN ANNUAL TEMPERATURE
#     * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/WC_TAVE_ANNUAL_50KM_ENA.TIF
#
#   CONTEMPORARY MEAN ANNUAL PRECIPITATION
#     * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/WC_PR_ANNUAL_50KM_ENA.TIF
#
#   LAND MASK
#     * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/GLDAS_50KM_ALBERS.NC
# 
#   PALEOPRECIPITATION
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CCSM-MARUM_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CM2MC_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/COSMOS-W_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/IPSL_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/MIROC-W_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_MWF_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CCSM-NCAR_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/COSMOS-S_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/HADCM3_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/MIROC-S_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_LORENZ_PR_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_PR_100KM_ENA.NC
# 
#   PALEOTEMPERATURE
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CCSM-MARUM_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CM2MC_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/COSMOS-W_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/IPSL_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/MIROC-W_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_MWF_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CCSM-NCAR_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/COSMOS-S_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/HADCM3_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/MIROC-S_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_LORENZ_TAS_100KM_ENA.NC
#     * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_TAS_100KM_ENA.NC
# 
# OUTPUT:
#   FIGURE S3
#     * RESOLUTION_SENSITIVITY/100KM/FIGURES/SAR_BIODIVERSITY_TEMPERATURE_MODERN_ENA_MODERN_CORRELATES.PDF
#     
#   FIGURE S4
#     * RESOLUTION_SENSITIVITY/100KM/FIGURES/SAR_BIODIVERSITY_PRECIPITATION_MODERN_ENA_MODERN_CORRELATES.PDF
#     
# ANACONDA R ENVIRONMENT WITH:
#
# R-BASE: 4.1.0
# R-SPDEP: 1.2_2
# R-SF: 1.0_5
# SPATIALREG: 1.2-3
# NCF: 1.3-2
# R-TIDYVERSE: 1.3.1
# PARALLEL: 4.1.0
# USDM: 1.1-18
# PATCHWORK: 1.1.2
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
library(spdep)
library(sf)
library(spatialreg)
library(tidyverse)

# Set working directory to project
setwd("~/Fastovich_et_al_2023_PhilB")

#############
# TEMPERATURE
#############

# Loop through organisms biodiversity data
biodiversity_files <- list(
  "amph_richness_epsg_102008.tif",
  "birds_Passeriformes_richness_epsg_102008.tif",
  "mamms_richness_epsg_102008.tif",
  "reps_richness_epsg_102008.tif",
  "tree_richness_epsg_102008.tif"
)

# Land sea mask
mask <- raster("resolution_sensitivity/100km/modern_data/postprocessed/GLDAS_100km_albers.nc")

# Modern temperature correlate
modern_tave <- raster("resolution_sensitivity/100km/modern_data/postprocessed/wc_tave_annual_100km_ena.tif")

# Modern precipitation correlate
modern_pr <- raster("resolution_sensitivity/100km/modern_data/postprocessed/wc_pr_annual_100km_ena.tif")

# Preparing some data holding lists
tave_fitted <- list()
pr_fitted <- list()

for(i in seq(biodiversity_files)) {
  # Loop through precipitation files and perform climate velocity calculation
  tas_files <- list.files("resolution_sensitivity/100km/model_temperature_anomaly_100km_ena/")
  pr_files <- list.files("resolution_sensitivity/100km/model_precipitation_anomaly_100km_ena/")
  
  for (j in seq(tas_files)) {
    
    model_tave <- raster(
      paste0("resolution_sensitivity/100km/model_temperature_anomaly_100km_ena/", tas_files[j]),
      varname = "tas"
    )
    
    model_pr <- raster(
      paste0("resolution_sensitivity/100km/model_precipitation_anomaly_100km_ena/", pr_files[j]),
      varname = "pr"
    )
    
    richness <- raster(
      paste0("resolution_sensitivity/100km/richness/", biodiversity_files[i])
    )
    
    # Aptly land mask to each of the raster
    # The mask layer is a bit larger than the richness by design and
    # needs to be cropped first
    model_tave <- raster::mask(
      model_tave,
      mask = mask,
      maskvalue = 0
    )
    model_pr <- raster::mask(
      model_pr,
      mask = mask,
      maskvalue = 0
    )
    
    richness <- raster::mask(
      raster::crop(
        richness,
        extent(-4025000, 3975000, -2975000, 2975000)
      ),
      mask = raster::crop(
        mask,
        extent(-4025000, 3975000, -2975000, 2975000)
      ),
      maskvalue=0
    )
    
    # Build data frame with corresponding richness and model output
    model_input <- as.data.frame(rasterToPoints(richness))[,c("x","y")]
    model_input$richness <- raster::extract(richness, model_input[,c("x","y")], simple = TRUE)
    model_input$model_tave <- raster::extract(model_tave, model_input[,c("x","y")], simple = TRUE)
    model_input$model_pr <- raster::extract(model_pr, model_input[,c("x","y")], simple = TRUE)
    model_input$modern_tave <- raster::extract(modern_tave, model_input[,c("x","y")], simple = TRUE)
    model_input$modern_pr <- raster::extract(modern_pr, model_input[,c("x","y")], simple = TRUE)
    model_input$model_pr2 <- model_input$model_pr^2
    
    # Clipping raster to eastern North America since thats where we have proxy
    # defined here as north of 24°24'2.458", south of 49°40'35.204", east of 
    # -91°18'54.247" and west of -56°39'55.054" in EPSG:4326 WGS 84.
    # data to verify the climate models then removing NAs which has the 
    # This also has the effect of removing model temperatures above the ocean
    model_input <- model_input %>% 
      filter(x >= 320000 & x <= 2850000 & y >= -1550000 & y <= 1150000) %>% 
      na.omit()
    
    # Convert to SpatialPointsDataFrame
    model_input_sf <- st_as_sf(model_input, coords = c("x","y"))
    model_input_sf <- st_set_crs(model_input_sf, CRS("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) 
    
    # Create null model
    fit_null <- lm(richness ~  modern_tave + modern_pr + model_tave + model_pr, data = model_input)
    
    # Create adjacency matrix
    neighbors_100 <- dnearneigh(model_input_sf, 0, 100*1000)
    W_100 <- nb2listw(neighbors_100)
    fit_100 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE)

    # Predict values to plot GLM lines with error
    # Source: https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
    # Modern + paleoclimate parameters parameters tave effect
    fit <- fit_100
    W <- W_100
    
    fit_modern_paleo_modern_tave_effect <- with(
      model_input, 
      tibble(
        model_tave = mean(model_tave),
        model_pr = mean(model_pr),
        model_pr2 = mean(model_pr2),
        modern_tave = seq(min(modern_tave), max(modern_tave), length = 200),
        modern_pr = mean(modern_pr)
      )
    ) %>% 
      add_column(fit = predict(fit, listw = W, newdata = ., type = "TC")[,1]) %>% 
      add_column(model = gsub("_tas_anom_100km_ena.nc", "", tas_files[j])) %>% 
      add_column(species  = case_when(i == 1 ~ "Amphibians",
                                      i == 2 ~ "Birds",
                                      i == 3 ~ "Mammals",
                                      i == 4 ~ "Reptiles",
                                      i == 5 ~ "Trees"))
    
    # Modern + paleoclimate parameters parameters pr effect
    fit_modern_paleo_modern_pr_effect <- with(model_input, 
                                             tibble(
                                               model_pr = mean(model_pr),
                                               model_pr2 = mean(model_pr2),
                                               model_tave = mean(model_tave),
                                               modern_tave = mean(modern_tave),
                                               modern_pr = seq(min(modern_pr), max(modern_pr), length = 200)
                                             )
    ) %>% 
      add_column(fit = predict(fit, listw = W, newdata = ., type = "TC")[,1]) %>% 
      add_column(model = gsub("_tas_anom_100km_ena.nc", "", tas_files[j])) %>% 
      add_column(species  = case_when(i == 1 ~ "Amphibians",
                                      i == 2 ~ "Birds",
                                      i == 3 ~ "Mammals",
                                      i == 4 ~ "Reptiles",
                                      i == 5 ~ "Trees"))
    
    # Conditional to add significance identifier
    model_summary_fit <- summary(fit_100)
    
    if(model_summary_fit$Coef["modern_tave","Pr(>|z|)"] <= 0.05) {
      fit_modern_paleo_modern_tave_effect$sig <- "Significant"
    } else {
      fit_modern_paleo_modern_tave_effect$sig <- "Not Significant"
    }
    
    if(model_summary_fit$Coef["modern_pr","Pr(>|z|)"] <= 0.05) {
      fit_modern_paleo_modern_pr_effect$sig <- "Significant"
    } else {
      fit_modern_paleo_modern_pr_effect$sig <- "Not Significant"
    }
    
    # Retain results in premade lists
    tave_fitted <- c(tave_fitted, list(fit_modern_paleo_modern_tave_effect))
    pr_fitted <- c(pr_fitted, list(fit_modern_paleo_modern_pr_effect))
  }
}

######
# PLOT
######

# Getting plotting settings in order ------------------------------------------

# Reordering the facet wrap plots from alphabetical to model skill order
tave_fitted_tibble <- tave_fitted %>% 
  bind_rows() %>% 
  mutate(
    model = case_when(
      model == "PROXY_KRIGING" ~ "Proxy Kriging",
      model == "TRACE_LORENZ" ~ "TraCE-21ka\n(Statistcally Downscaled)",
      model == "TRACE_MWF" ~ "TraCE-MWF (Single Forcing)",
      model == "TRACE" ~ "TraCE-21ka",
      TRUE ~ model
    )
  ) %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "Proxy Kriging",
        "TraCE-21ka\n(Statistcally Downscaled)",
        "TraCE-21ka",
        "MIROC-S",
        "CM2Mc",
        "CCSM-NCAR",
        "TraCE-MWF (Single Forcing)",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  ) %>% 
  select(-c(model_tave, model_pr))

pr_fitted_tibble <- pr_fitted %>% 
  bind_rows() %>% 
  mutate(
    model = case_when(
      model == "PROXY_KRIGING" ~ "Proxy Kriging",
      model == "TRACE_LORENZ" ~ "TraCE-21ka\n(Statistcally Downscaled)",
      model == "TRACE_MWF" ~ "TraCE-MWF (Single Forcing)",
      model == "TRACE" ~ "TraCE-21ka",
      TRUE ~ model
    )
  ) %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "Proxy Kriging",
        "TraCE-21ka\n(Statistcally Downscaled)",
        "TraCE-21ka",
        "MIROC-S",
        "CM2Mc",
        "CCSM-NCAR",
        "TraCE-MWF (Single Forcing)",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  ) %>% 
  select(-c(model_tave, model_pr))

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

# We're plotting the raw data against the model fits which means creating a
# large data frame of every lat, lon, richness estimate, temperature, and
# precipitation pairing. The code below accomplishes that.

# Richness --------------------------------------------------------------------
# Same order as above
richness_raster_list <- lapply(
  biodiversity_files,
  function(x) raster(
    paste0("resolution_sensitivity/100km/richness/", x)
  )
) %>% 
  purrr::map(function(x) raster::mask(
    raster::crop(
      x,
      extent(-4025000, 3975000, -2975000, 2975000)
    ),
    mask = raster::crop(
      mask,
      extent(-4025000, 3975000, -2975000, 2975000)
    ),
    maskvalue=0
  )) %>% 
  purrr::map(function(x) raster::extract(x, model_input[,c("x","y")], simple = TRUE)) %>% 
  purrr::map(as.data.frame) %>% 
  purrr:::map(function(x) add_column(x, x = model_input$x, y = model_input$y)) %>% 
  purrr:::map(function(x) rename(x, richness = `.x[[i]]`)) %>% 
  purrr::map(function(x) right_join(x, model_input[,c("x", "y", "modern_tave", "modern_pr")], by = c("x", "y"))) # Can do this join because the modern climate variables do not change across different paleoclimate estimates

# Add species identifier to each tibble with the list
richness_raster_list[[1]] <- add_column(richness_raster_list[[1]], species = "Amphibians")
richness_raster_list[[2]] <- add_column(richness_raster_list[[2]], species = "Birds")
richness_raster_list[[3]] <- add_column(richness_raster_list[[3]], species = "Mammals")
richness_raster_list[[4]] <- add_column(richness_raster_list[[4]], species = "Reptiles")
richness_raster_list[[5]] <- add_column(richness_raster_list[[5]], species = "Trees")

# Convert to a single data frame
# https://stackoverflow.com/a/68181517
richness_df <- richness_raster_list %>% 
  bind_rows() %>% 
  add_column(sig = NA) %>% 
  crossing(
    model = c(
      "Proxy Kriging",
      "TraCE-21ka\n(Statistcally Downscaled)",
      "TraCE-21ka",
      "MIROC-S",
      "CM2Mc",
      "CCSM-NCAR",
      "TraCE-MWF (Single Forcing)",
      "CCSM-MARUM", 
      "COSMOS-S",
      "IPSL",
      "COSMOS-W",
      "MIROC-W",
      "HadCM3"
    )
  )

# PLOT 4 - TEMPERATURE FOR ALL MODELS PALEOCLIMATE + MODERN ONLY ---------------

supp_modern_biodiversity_temperature_models <- tave_fitted_tibble %>% 
  ggplot(aes(x = modern_tave, y = fit, color = species, linetype = sig)) + 
  geom_point(data = richness_df, aes(x = modern_tave, y = richness, color = species), alpha = 0.025) +
  geom_path() + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species") + 
  scale_linetype_manual(values = c("dotted", "solid"), name = "Model Significance") + 
  facet_wrap(~model, scales = "free", ncol = 3) +
  ylab("Species Richness") + 
  xlab("Modern Mean Annual Temperature (K)")

ggsave(supp_modern_biodiversity_temperature_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_temperature_modern_ena_modern_correlates_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(supp_modern_biodiversity_temperature_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_temperature_modern_ena_modern_correlates_quadratic.png", height = 10, width = 8, dpi = 300)

# PLOT 6 - PRECIPITATION FOR ALL MODELS PALEOCLIMATE + MODERN ONLY -------------

supp_modern_biodiversity_precipitation_models <- pr_fitted_tibble %>% 
  ggplot(aes(x = modern_pr, y = fit, color = species, linetype = sig)) + 
  geom_point(data = richness_df, aes(x = modern_pr, y = richness, color = species), alpha = 0.025) +
  geom_path() + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species") + 
  scale_linetype_manual(values = c("dotted", "solid"), name = "Model Significance") + 
  facet_wrap(~model, scales = "free", ncol = 3) +
  ylab("Species Richness") + 
  xlab("Modern Mean Annual Precipitation (mm/day)")

ggsave(supp_modern_biodiversity_precipitation_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_precipitation_modern_ena_modern_correlates_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(supp_modern_biodiversity_precipitation_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_precipitation_modern_ena_modern_correlates_quadratic.png", height = 10, width = 8, dpi = 300)

