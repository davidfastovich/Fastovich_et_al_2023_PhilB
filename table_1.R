################################################################################
# THIS SCRIPT PRODUCES PIVOT TABLES THAT AVERAGE THE RESULTS OF SPATIAL ERROR
# MODEL FITS ON THE BASIS OF CLIMATE MODELS AND TAXANOMIC CLASS. IN THE
# MANUSCRIPT, WE ONLY USE THE AVERAGES BY TAXAONOMMIC CLASS BUT RETAIN BOTH
# FOR COMPLETENESSS.
# 
# INPUT:
# * SKILL_SCORES/SAR_DIAG.CSV
# 
# OUTPUT:
#   TABLE 1:
#     * FIGURES/TABLE_1.CSV
#   
#   EXTRA
#     * FIGURES/TABLE_1.CSV
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

library(tidyverse)

# Set working directory to project
setwd("~/Fastovich_et_al_2023_PhilB")

##############
# READ IN DATA
##############

# AIC scores calculated by me and are the result of spatial error models
aic <- read_csv("skill_scores/sar_diag.csv") %>%
  select(
    AIC_fit_50,
    AIC_fit_100,
    paleo_tave_deviance,
    paleo_pr_deviance,
    paleo_pr2_deviance,
    modern_tave_deviance,
    modern_pr_deviance,
    climate_simulation,
    pseudo_r2_fit_50,
    species
  ) %>% 
  rename(
    `Climate Model` = climate_simulation,
    Species = species,
  )

# Retain 50 km neighborhood AIC for all taxonomic groups except birds
aic <- aic %>% 
  mutate(AIC =
           case_when(
             Species != "bird" ~ AIC_fit_50,
             Species == "bird" ~ AIC_fit_100,
           )
  )


#######################
# RENAME SPECIES COLUMN
#######################

# The species identifiers are shorthand and need to be converted into proper
# names. This helpfer function takes care of that.

rename_species <- function(x) {
  if(x == "mamm") {
    return("Mammals") 
  } else if (x == "amph") {
    return("Amphibians")
  } else if (x == "bird") {
    return("Birds")
  } else if (x == "reps") {
    return("Reptiles")
  } else {
    return("Trees")
  }
}

# Apply function to the AIC data frames
aic$Species <- sapply(aic$Species, rename_species, USE.NAMES = FALSE)

####################
# SUMMARIZE THE DATA
####################

species_pivot <- aic %>% 
  group_by(Species) %>% 
  summarise(`Nagelkerke's Psuedo R2` = mean(pseudo_r2_fit_50),
            `Paleotemperature Deviance Explained (%)` = mean(paleo_tave_deviance),
            `Paleoprecipitation Deviance Explained (%)` = mean(paleo_pr_deviance),
            `Paleoprecipitation^2 Deviance Explained (%)` = mean(paleo_pr2_deviance),
            `Modern Temperature Deviance Explained (%)` = mean(modern_tave_deviance),
            `Modern Precipitation Deviance Explained (%)` = mean(modern_pr_deviance)) %>% 
  mutate(
    `Paleoclimate Deviance Explained (%)` = `Paleotemperature Deviance Explained (%)` + `Paleoprecipitation Deviance Explained (%)` + `Paleoprecipitation^2 Deviance Explained (%)`,
    `Modern Climate Deviance Explained (%)` = `Modern Temperature Deviance Explained (%)` + `Modern Precipitation Deviance Explained (%)`
  ) %>% 
  write_csv("figures/table_1.csv")

model_pivot <- aic %>% 
  group_by(`Climate Model`) %>% 
  summarise(`Nagelkerke's Psuedo R2` = mean(pseudo_r2_fit_50),
            `Paleotemperature Deviance Explained (%)` = mean(paleo_tave_deviance),
            `Paleoprecipitation Deviance Explained (%)` = mean(paleo_pr_deviance),
            `Paleoprecipitation^2 Deviance Explained (%)` = mean(paleo_pr2_deviance),
            `Modern Temperature Deviance Explained (%)` = mean(modern_tave_deviance),
            `Modern Precipitation Deviance Explained (%)` = mean(modern_pr_deviance)) %>% 
  mutate(
    `Paleoclimate Deviance Explained (%)` = `Paleotemperature Deviance Explained (%)` + `Paleoprecipitation Deviance Explained (%)`  + `Paleoprecipitation^2 Deviance Explained (%)`,
    `Modern Climate Deviance Explained (%)` = `Modern Temperature Deviance Explained (%)` + `Modern Precipitation Deviance Explained (%)`
  ) %>% 
  write_csv("figures/table_1_model.csv")
