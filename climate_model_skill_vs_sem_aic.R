################################################################################
# THIS SCRIPT PRODUCES FIGURE S10 AND S11 IN FASTOVICH ET AL., (2023, PHIL B) 
# THAT COMPARES SPATIAL ERROR MODEL INFORMATION LOSS AGAINST CLIMATE MODEL SKILL
# SCORE
# 
# INPUT:
# * SKILL_SCORES/TAS_SKILL.CSV
# * SKILL_SCORES/PR_SKILL.CSV
# * SKILL_SCORES/SAR_DIAG.CSV
# 
# OUTPUT:
#   FIGURE S10
#     * FIGURES/SYNTHESIS_TEMPERATURE_PALEO_MODERN_SAR.PDF
#     
#   FIGURE S11:
#     * FIGURES/SYNTHESIS_PRECIPITATION_PALEO_MODERN_SAR.PDF
#     
# ANACONDA R ENVIRONMENT WITH:
#
# R-BASE: 4.1.0
# R-TIDYVERSE: 1.3.1
#
# ALL PACKAGES WITH AN `R-` PREFIX WERE INSTALLED FROM THE CONDA-FORGE CHANNEL
#
# AUTHOR: DAVID FASTOVICH
# CONTACT: FASTOVICH@WISC.EDU
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

# Calculated by me from tas_skill.py
tas_skill <- read.csv("skill_scores/tas_skill.csv")

# Calculated by me from pr_skill.py 
pr_skill <- read.csv("skill_scores/pr_skill.csv")

# AIC from the SEM models calculated here
aic <- read_csv("skill_scores/sar_diag_nonlinear_again.csv") %>%
  select(AIC_fit_50, AIC_fit_100, climate_simulation, species) %>% 
  rename(model = climate_simulation)

aic <- aic %>% 
  mutate(
    AIC = case_when(
      species != "bird" ~ AIC_fit_50,
      species == "bird" ~ AIC_fit_100,
    )
  )

# Removing proxy AIC because there is no skill score for proxy-proxy comparisons
# its strictly a model skill metric. However, we still want this value because
# its a useful metric for compariing AIC between models and the proxy record.
aic_proxy <- aic[aic$model == "PROXY_KRIGING",]
aic <- aic[aic$model != "PROXY_KRIGING",]

####################
# CLEAN MODEL COLUMN
####################

# The model column has extra information about the files that we are removing,
# to retain the model only

tas_skill$model <- gsub("_tas_anom_50km_ena.nc", "", tas_skill$model)
tas_skill$model <- gsub("model_temperature_anomaly_50km_ena/", "", tas_skill$model)

pr_skill$model <- gsub("_pr_anom_50km_ena.nc", "", pr_skill$model)
pr_skill$model <- gsub("model_precipitation_anomaly_50km_ena/", "", pr_skill$model)

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
aic$species <- sapply(aic$species, rename_species, USE.NAMES = FALSE)
aic_proxy$species <- sapply(aic_proxy$species, rename_species, USE.NAMES = FALSE)

#########################
# MERGE ALL DATA TOGETHER
#########################

tas_skill_aic <- full_join(aic, tas_skill, by = "model")
pr_skill_aic <- left_join(aic, pr_skill, by = "model")


######
# PLOT
######

# Change the default ordering of the three legend boxes
# https://stackoverflow.com/a/10037881
guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

# AIC vs Temperature skill -----------------------------------------------------
tas_plot = ggplot(data = tas_skill_aic, aes(x = median_skill, y = AIC, shape = model, color = species, fill = species)) + 
  # geom_ribbon(data = regression_line_tas, aes(x = median_skill, y = fit, ymin = lwr, ymax = upr), color = "transparent", alpha = 0.2) +
  # geom_line(data = regression_line_tas, aes(x = median_skill, y = fit, linetype = sig)) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "Model Significance") + 
  geom_point() + 
  scale_shape_manual(values = 1:12, name = "Model") + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species", guide = "none") +
  # geom_hline(data = aic_proxy, aes(yintercept = AIC), color = "grey50") + 
  facet_wrap(~species+null, scales = "free_y", nrow = 5, labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))
  }) + 
  xlab("Temperature Climate Model Skill Score") + 
  scale_y_reverse() +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(vjust = -1),
        panel.spacing.y = unit(-.1, "cm"),
        legend.key = element_blank())

ggsave(tas_plot, filename = "figures/synthesis_temperature_paleo_modern_sar_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(tas_plot, filename = "figures/synthesis_temperature_paleo_modern_sar_quadratic.png", height = 10, width = 8, dpi = 300)

# AIC vs Precipitation skill ---------------------------------------------------
pr_plot = ggplot(data = pr_skill_aic, aes(x = median_skill, y = AIC, shape = model, color = species, fill = species)) + 
  # geom_ribbon(data = regression_line_pr, aes(x = median_skill, y = fit, ymin = lwr, ymax = upr), color = "transparent", alpha = 0.2) +
  # geom_line(data = regression_line_pr, aes(x = median_skill, y = fit, linetype = sig)) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "Model Significance") + 
  geom_point() + 
  scale_shape_manual(values = 1:12, name = "Model") + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species", guide = "none") +
  # geom_hline(data = aic_proxy, aes(yintercept = AIC), color = "grey50") + 
  facet_wrap(~species+null, scales = "free_y", nrow = 5, labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))
  }) + 
  xlab("Precipitation Climate Model Skill Score") + 
  scale_y_reverse() +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(vjust = -1),
        panel.spacing.y = unit(-.1, "cm"),
        legend.key = element_blank())

ggsave(pr_plot, filename = "figures/synthesis_precipitation_paleo_modern_sar_quadratic.pdf", height = 10, width = 8, dpi = 300)
ggsave(pr_plot, filename = "figures/synthesis_precipitation_paleo_modern_sar_quadratic.png", height = 10, width = 8, dpi = 300)
