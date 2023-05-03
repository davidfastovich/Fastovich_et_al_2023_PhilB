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
# BECAUSE THIS IS ONLY A SENSITITY TEST, NO DIAGNOSTICS ARE KEPT AS IN THE MAIN
# SCRIPT.
# 
# ALL STATISTICAL ANALYSES ARE PERFORMED FOR ALL CLASSES OF ORGANISMS IN A LOOP
# TO PRODUCE FIGURES WHERE ALL MODEL FITS ARE COMPARED.
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
#   FIGURE S5
#     * RESOLUTION_SENSITIVITY/100KM/FIGURES/SAR_BIODIVERSITY_TEMPERATURE_MODEL_ENA_MODERN_CORRELATES.PDF
#     
#   FIGURE S6
#     * RESOLUTION_SENSITIVITY/100KM/FIGURES/SAR_BIODIVERSITY_PRECIPITATION_MODEL_ENA_MODERN_CORRELATES.PDF
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
library(ncf)
library(tidyverse)
library(parallel)
library(usdm)

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
diag_list <- list()
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
    
    # Create adjacency matrix
    fit_null <- lm(richness ~  modern_tave + modern_pr + model_tave + model_pr, data = model_input)

    neighbors_100 <- dnearneigh(model_input_sf, 0, 100*1000)
    W_100 <- nb2listw(neighbors_100)
    fit_100 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE)

    # Need to force lamdba value - the interval command tells the model optimizer
    # what values to consider. We're tricking the optimizer a bit here by giving
    # it a 0.0000002 window of values to consider, thereby forcing the optimizer
    # to converge on the lamdda value to that from fit_50. This makes the models
    # identical and comparable.
    fit_100_null <- errorsarlm(richness ~  1, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + .0000001))
    fit_100_no_modern_tave <- errorsarlm(richness ~  modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + 0.0000001))
    fit_100_no_modern_pr <- errorsarlm(richness ~  modern_tave + model_tave + model_pr + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + 0.0000001))
    fit_100_no_paleo_tave <- errorsarlm(richness ~  modern_tave + modern_pr + model_pr + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + 0.0000001))
    fit_100_no_paleo_pr <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr2, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + 0.0000001))
    fit_100_no_paleo_pr2 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr, data = model_input, listw = W_100, zero.policy=TRUE, interval = c(fit_100$lambda - 0.0000001, fit_100$lambda + 0.0000001))

    neighbors_200 <- dnearneigh(model_input_sf, 0, 200*1000)
    W_200 <- nb2listw(neighbors_200)
    fit_200 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_200, zero.policy=TRUE)

    neighbors_400 <- dnearneigh(model_input_sf, 0, 400*1000)
    W_400 <- nb2listw(neighbors_400)
    fit_400 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_400, zero.policy=TRUE)

    neighbors_600 <- dnearneigh(model_input_sf, 0, 600*1000)
    W_600 <- nb2listw(neighbors_600)
    fit_600 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_600, zero.policy=TRUE)

    neighbors_800 <- dnearneigh(model_input_sf, 0, 800*1000)
    W_800 <- nb2listw(neighbors_800)
    fit_800 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_800, zero.policy=TRUE)

    neighbors_1000 <- dnearneigh(model_input_sf, 0, 1000*1000)
    W_1000 <- nb2listw(neighbors_1000)
    fit_1000 <- errorsarlm(richness ~  modern_tave + modern_pr + model_tave + model_pr + model_pr2, data = model_input, listw = W_1000, zero.policy=TRUE)

    # Calculate global Morans' I using `raster`
    res_raster_fit_null <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_null)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_100 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_100)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_200 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_200)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_400 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_400)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_600 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_600)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_800 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_800)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    res_raster_fit_1000 <- rasterFromXYZ(
      data.frame(
        x = model_input$x,
        y = model_input$y,
        z = residuals(fit_1000)
      ),
      res = c(50000, 50000),
      crs = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    )

    # Legacy portion of the code that I have not bothered to clean. Prior to
    # using Passeriformes data for birds the rooks case was not always the best
    # fit model and this enabled me to adjust for that
    fit <- fit_100
    fit_null <- fit_100_null
    fit_no_paleo_tave <- fit_100_no_paleo_tave
    fit_no_paleo_pr <- fit_100_no_paleo_pr
    fit_no_paleo_pr2 <- fit_100_no_paleo_pr2
    fit_no_modern_tave <- fit_100_no_modern_tave
    fit_no_modern_pr <- fit_100_no_modern_pr
    W <- W_100

    # Data frame of diagnostics - saving as list then csv
    # Psuedo-R2 explanation: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
    diag_list <- c(
      diag_list,
      list(
        tibble(
          paleo_tave_deviance = ((deviance(fit_no_paleo_tave) - deviance(fit))/deviance(fit_null)) * 100,
          paleo_pr_deviance = ((deviance(fit_no_paleo_pr) - deviance(fit))/deviance(fit_null)) * 100,
          paleo_pr2_deviance = ((deviance(fit_no_paleo_pr2) - deviance(fit))/deviance(fit_null)) * 100,
          modern_tave_deviance = ((deviance(fit_no_modern_tave) - deviance(fit))/deviance(fit_null)) * 100,
          modern_pr_deviance = ((deviance(fit_no_modern_pr) - deviance(fit))/deviance(fit_null)) * 100,
          AIC_fit_null = AIC(fit_null),
          AIC_fit_100 = AIC(fit_100),
          AIC_fit_200 = AIC(fit_200),
          AIC_fit_400 = AIC(fit_400),
          AIC_fit_600 = AIC(fit_600),
          AIC_fit_800 = AIC(fit_800),
          AIC_fit_1000 = AIC(fit_1000),
          pseudo_r2_fit_null = summary(fit_null, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_100 = summary(fit_100, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_200 = summary(fit_200, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_400 = summary(fit_400, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_600 = summary(fit_600, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_800 = summary(fit_800, Nagelkerke = TRUE)$NK,
          pseudo_r2_fit_1000 = summary(fit_1000, Nagelkerke = TRUE)$NK,
          moran_fit_null = raster::Moran(res_raster_fit_null),
          moran_fit_100 = raster::Moran(res_raster_fit_100),
          moran_fit_200 = raster::Moran(res_raster_fit_200),
          moran_fit_400 = raster::Moran(res_raster_fit_400),
          moran_fit_600 = raster::Moran(res_raster_fit_600),
          moran_fit_800 = raster::Moran(res_raster_fit_800),
          moran_fit_1000 = raster::Moran(res_raster_fit_1000),
          climate_simulation = gsub("_tas_anom_100km_ena.nc", "", tas_files[j]),
          species = str_sub(biodiversity_files[i], 0, 4)
        )
      )
    )
    
    # Predict values to plot GLM lines with error
    # Source: https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
    # Modern + paleoclimate parameters parameters tave effect
    fit_modern_paleo_model_tave_effect <- with(
      model_input, 
      tibble(
        model_tave = seq(min(model_tave), max(model_tave), length = 200),
        model_pr = mean(model_pr),
        model_pr2 = mean(model_pr)^2,
        modern_tave = mean(modern_tave),
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
    fit_modern_paleo_model_pr_effect <- with(model_input, 
                                             tibble(
                                               model_pr = seq(min(model_pr), max(model_pr), length = 200),
                                               model_pr2 = seq(min(model_pr), max(model_pr), length = 200)^2,
                                               model_tave = mean(model_tave),
                                               modern_tave = mean(modern_tave),
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
    
    # Conditional to add significance identifier
    model_summary_fit <- summary(fit)
    
    if(model_summary_fit$Coef["model_tave","Pr(>|z|)"] <= 0.05) {
      fit_modern_paleo_model_tave_effect$sig <- "Significant"
    } else {
      fit_modern_paleo_model_tave_effect$sig <- "Not Significant"
    }
    
    if(model_summary_fit$Coef["model_pr","Pr(>|z|)"] <= 0.05 & model_summary_fit$Coef["model_pr2","Pr(>|z|)"] <= 0.05) {
      fit_modern_paleo_model_pr_effect$sig <- "Significant"
    } else if(model_summary_fit$Coef["model_pr2","Pr(>|z|)"] <= 0.05  & model_summary_fit$Coef["model_pr","Pr(>|z|)"] > 0.05) {
      fit_modern_paleo_model_pr_effect$sig <- "Only Quadratic Term Significant"
    } else if(model_summary_fit$Coef["model_pr","Pr(>|z|)"] <= 0.05 & model_summary_fit$Coef["model_pr2","Pr(>|z|)"] > 0.05) {
      fit_modern_paleo_model_pr_effect$sig <- "Only Linear Term Significant"
    } else {
      fit_modern_paleo_model_pr_effect$sig <- "Not Significant"
    }
    
    # Retain results in premade lists
    tave_fitted <- c(tave_fitted, list(fit_modern_paleo_model_tave_effect))
    pr_fitted <- c(pr_fitted, list(fit_modern_paleo_model_pr_effect))
  }
}

#########################
# SAVE DIAGNOSTICS AS CSV
#########################

diag_list %>%
  bind_rows() %>% 
  write_csv("resolution_sensitivity/100km/skill_scores/sar_diag.csv")

######
# PLOT
######

# Getting plotting settings in order ------------------------------------------

# Reordering the facet wrap plots from alphabetical to model skill order
tave_fitted_tibble <- tave_fitted %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "PROXY_KRIGING", # proxy
        "TRACE_LORENZ",
        "TRACE", # sig+agree
        "MIROC-S", # sig+oppo
        "CM2Mc",
        "CCSM-NCAR", # insig
        "TRACE_MWF",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  ) %>% 
  select(-c(modern_tave, modern_pr))

pr_fitted_tibble <- pr_fitted %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "PROXY_KRIGING", # proxy
        "TRACE_LORENZ",
        "TRACE", # sig+agree
        "MIROC-S", # sig+oppo
        "CM2Mc",
        "CCSM-NCAR", # insig
        "TRACE_MWF",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  ) %>% 
  select(-c(modern_tave, modern_pr))

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

# Temperature ------------------------------------------------------------------

# We're plotting the raw data against the model fits which means creating a
# large data frame of every lat, lon, richness estimate, paleoprecipitation, and
# paleotemperature pairing. The code below accomplishes that.

# Read in all temperature files
model_tave_raster_list <- lapply(
  tas_files,
  function(x) raster(
    paste0("resolution_sensitivity/100km/model_temperature_anomaly_100km_ena/", x),
    varname = "tas"
  )
) %>% 
  # Apply mask
  purrr::map(function(x) raster::mask(
    x,
    mask = mask,
    maskvalue = 0
  )) %>% 
  # Extract output at model_input coordinates
  purrr::map(function(x) raster::extract(x, model_input[,c("x","y")], simple = TRUE)) %>%
  # Convert to data frame
  purrr::map(as.data.frame) %>% 
  # Add x/y which will be used as a key later
  purrr:::map(function(x) add_column(x, x = model_input$x, y = model_input$y)) %>% 
  # Rename column
  purrr:::map(function(x) rename(x, model_tave = `.x[[i]]`))

# Add model identifier to each tibble with the list
model_tave_raster_list[[1]] <- add_column(model_tave_raster_list[[1]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[1]))
model_tave_raster_list[[2]] <- add_column(model_tave_raster_list[[2]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[2]))
model_tave_raster_list[[3]] <- add_column(model_tave_raster_list[[3]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[3]))
model_tave_raster_list[[4]] <- add_column(model_tave_raster_list[[4]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[4]))
model_tave_raster_list[[5]] <- add_column(model_tave_raster_list[[5]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[5]))
model_tave_raster_list[[6]] <- add_column(model_tave_raster_list[[6]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[6]))
model_tave_raster_list[[7]] <- add_column(model_tave_raster_list[[7]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[7]))
model_tave_raster_list[[8]] <- add_column(model_tave_raster_list[[8]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[8]))
model_tave_raster_list[[9]] <- add_column(model_tave_raster_list[[9]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[9]))
model_tave_raster_list[[10]] <- add_column(model_tave_raster_list[[10]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[10]))
model_tave_raster_list[[11]] <- add_column(model_tave_raster_list[[11]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[11]))
model_tave_raster_list[[12]] <- add_column(model_tave_raster_list[[12]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[12]))
model_tave_raster_list[[13]] <- add_column(model_tave_raster_list[[13]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[13]))

# Convert to data frame
model_tave_df <- do.call(rbind, model_tave_raster_list)

# Precipitation ----------------------------------------------------------------
# Same order as above
# Read in all temperature files
model_pr_raster_list <- lapply(
  pr_files,
  function(x) raster(
    paste0("resolution_sensitivity/100km/model_precipitation_anomaly_100km_ena/", x),
    varname = "pr"
  )
) %>% 
  # Apply mask
  purrr::map(function(x) raster::mask(
    x,
    mask = mask,
    maskvalue = 0
  )) %>% 
  # Extract output at model_input coordinates
  purrr::map(function(x) raster::extract(x, model_input[,c("x","y")], simple = TRUE)) %>%
  # Convert to data frame
  purrr::map(as.data.frame) %>% 
  # Add x/y which will be used as a key later
  purrr:::map(function(x) add_column(x, x = model_input$x, y = model_input$y)) %>% 
  # Rename column
  purrr:::map(function(x) rename(x, model_pr = `.x[[i]]`))

# Add model identifier to each tibble with the list
model_pr_raster_list[[1]] <- add_column(model_pr_raster_list[[1]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[1]))
model_pr_raster_list[[2]] <- add_column(model_pr_raster_list[[2]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[2]))
model_pr_raster_list[[3]] <- add_column(model_pr_raster_list[[3]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[3]))
model_pr_raster_list[[4]] <- add_column(model_pr_raster_list[[4]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[4]))
model_pr_raster_list[[5]] <- add_column(model_pr_raster_list[[5]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[5]))
model_pr_raster_list[[6]] <- add_column(model_pr_raster_list[[6]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[6]))
model_pr_raster_list[[7]] <- add_column(model_pr_raster_list[[7]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[7]))
model_pr_raster_list[[8]] <- add_column(model_pr_raster_list[[8]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[8]))
model_pr_raster_list[[9]] <- add_column(model_pr_raster_list[[9]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[9]))
model_pr_raster_list[[10]] <- add_column(model_pr_raster_list[[10]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[10]))
model_pr_raster_list[[11]] <- add_column(model_pr_raster_list[[11]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[11]))
model_pr_raster_list[[12]] <- add_column(model_pr_raster_list[[12]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[12]))
model_pr_raster_list[[13]] <- add_column(model_pr_raster_list[[13]], model = gsub("_tas_anom_100km_ena.nc", "", tas_files[13]))

# Convert to data frame
model_pr_df <- do.call(rbind, model_pr_raster_list)

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
  purrr:::map(function(x) rename(x, richness = `.x[[i]]`))

# Add species identifier to each tibble with the list
richness_raster_list[[1]] <- add_column(richness_raster_list[[1]], species = "Amphibians")
richness_raster_list[[2]] <- add_column(richness_raster_list[[2]], species = "Birds")
richness_raster_list[[3]] <- add_column(richness_raster_list[[3]], species = "Mammals")
richness_raster_list[[4]] <- add_column(richness_raster_list[[4]], species = "Reptiles")
richness_raster_list[[5]] <- add_column(richness_raster_list[[5]], species = "Trees")

# Convert to a single data frame
richness_df <- richness_raster_list %>% 
  bind_rows() %>% 
  add_column(sig = NA)

# Combine model tave with species richness to get the desired output data frame
model_tave_richness_df <- model_tave_raster_list %>% 
  purrr::map(function(x) right_join(x, richness_df, by = c("x", "y"))) %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "PROXY_KRIGING", # proxy
        "TRACE_LORENZ",
        "TRACE", # sig+agree
        "MIROC-S", # sig+oppo
        "CM2Mc",
        "CCSM-NCAR", # insig
        "TRACE_MWF",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  )

# Combine model pr with species richness to get the desired output data frame
model_pr_richness_df <- model_pr_raster_list %>% 
  purrr::map(function(x) right_join(x, richness_df, by = c("x", "y"))) %>% 
  bind_rows() %>% 
  mutate(
    model = factor(
      model, 
      levels = c(
        "PROXY_KRIGING", # proxy
        "TRACE_LORENZ",
        "TRACE", # sig+agree
        "MIROC-S", # sig+oppo
        "CM2Mc",
        "CCSM-NCAR", # insig
        "TRACE_MWF",
        "CCSM-MARUM", 
        "COSMOS-S",
        "IPSL",
        "COSMOS-W",
        "MIROC-W",
        "HadCM3"
      )
    )
  )

# PLOT 1 - TEMPERATURE FOR ALL MODELS PALEOCLIMATE + MODERN ONLY ---------------

supp_modern_biodiversity_temperature_models <- tave_fitted_tibble %>% 
  ggplot(aes(x = model_tave, y = fit, color = species, linetype = sig)) + 
  geom_point(data = model_tave_richness_df, aes(x = model_tave, y = richness, color = species), alpha = 0.025) +
  geom_path() + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species") + 
  scale_linetype_manual(values = c("dotted", "dashed", "twodash", "solid"), name = "Model Significance", breaks=c("Not Significant", "Only Linear Term Significant", "Only Quadratic Term Significant", "Significant")) + 
  facet_wrap(~model, scales = "free", ncol = 3) +
  ylab("Fitted Species Richness") + 
  xlab("Temperature Anomaly (K)")

ggsave(supp_modern_biodiversity_temperature_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_temperature_model_ena_modern_correlates_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(supp_modern_biodiversity_temperature_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_temperature_model_ena_modern_correlates_quadratic.png", height = 10, width = 8, dpi = 300)

# PLOT 2 - PRECIPITATION FOR ALL MODELS PALEOCLIMATE + MODERN ONLY -------------

supp_modern_biodiversity_precipitation_models <- pr_fitted_tibble %>% 
  ggplot(aes(x = model_pr, y = fit, color = species, linetype = sig)) + 
  geom_point(data = model_pr_richness_df, aes(x = model_pr, y = richness, color = species), alpha = 0.025) +
  geom_path() + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species") + 
  scale_linetype_manual(values = c("dotted", "dashed", "twodash", "solid"), name = "Model Significance", breaks=c("Not Significant", "Only Linear Term Significant", "Only Quadratic Term Significant", "Significant")) + 
  facet_wrap(~model, scales = "free", ncol = 3) +
  ylab("Fitted Species Richness") + 
  xlab("Precipitation Anomaly (mm/day)")

ggsave(supp_modern_biodiversity_precipitation_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_precipitation_model_ena_modern_correlates_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(supp_modern_biodiversity_precipitation_models, filename = "resolution_sensitivity/100km/figures/sar_biodiversity_precipitation_model_ena_modern_correlates_quadratic.png", height = 10, width = 8, dpi = 300)
