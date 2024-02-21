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

# loglik from the SEM models calculated here
loglik <- read_csv("skill_scores/sar_diag.csv") %>%
  select(logLik_fit_50, logLik_fit_100, climate_simulation, species) %>% 
  rename(model = climate_simulation)

loglik <- loglik %>% 
  mutate(
    loglik = case_when(
      species != "bird" ~ logLik_fit_50,
      species == "bird" ~ logLik_fit_100,
    )
  )

# Removing proxy loglik because there is no skill score for proxy-proxy comparisons
# its strictly a model skill metric. However, we still want this value because
# its a useful metric for compariing loglik between models and the proxy record.
loglik_proxy <- loglik[loglik$model == "PROXY_KRIGING",]
loglik <- loglik[loglik$model != "PROXY_KRIGING",]

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

# Apply function to the loglik data frames
loglik$species <- sapply(loglik$species, rename_species, USE.NAMES = FALSE)

#########################
# MERGE ALL DATA TOGETHER
#########################

tas_skill_loglik <- full_join(loglik, tas_skill, by = "model")
pr_skill_loglik <- left_join(loglik, pr_skill, by = "model")

# Renaming PROXY_KRIGING, TRACE_LORENZ, TRACE, and TRACE-MWF
tas_skill_loglik$model[tas_skill_loglik$model == "PROXY_KRIGING"] <- "Proxy Kriging"
tas_skill_loglik$model[tas_skill_loglik$model == "TRACE_LORENZ"] <- "TraCE-21ka\n(Statistcally Downscaled)"
tas_skill_loglik$model[tas_skill_loglik$model == "TRACE_MWF"] <- "TraCE-MWF (Single Forcing)"
tas_skill_loglik$model[tas_skill_loglik$model == "TRACE"] <- "TraCE-21ka"

pr_skill_loglik$model[pr_skill_loglik$model == "PROXY_KRIGING"] <- "Proxy Kriging"
pr_skill_loglik$model[pr_skill_loglik$model == "TRACE_LORENZ"] <- "TraCE-21ka\n(Statistcally Downscaled)"
pr_skill_loglik$model[pr_skill_loglik$model == "TRACE_MWF"] <- "TraCE-MWF (Single Forcing)"
pr_skill_loglik$model[pr_skill_loglik$model == "TRACE"] <- "TraCE-21ka"

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

# loglik vs Temperature skill -----------------------------------------------------
tas_plot = ggplot(data = tas_skill_loglik, aes(x = median_skill, y = loglik, shape = model, color = species, fill = species)) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "Model Significance") + 
  geom_point() + 
  scale_shape_manual(values = 1:12, name = "Model") + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species", guide = "none") +
  facet_wrap(~species+null, scales = "free_y", nrow = 5, labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))
  }) + 
  xlab("Temperature Climate Model Skill Score") + 
  ylab("Log-Likelihood") + 
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

# loglik vs Precipitation skill ---------------------------------------------------
pr_plot = ggplot(data = pr_skill_loglik, aes(x = median_skill, y = loglik, shape = model, color = species, fill = species)) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "Model Significance") + 
  geom_point() + 
  scale_shape_manual(values = 1:12, name = "Model") + 
  scale_color_brewer(palette = "Dark2", name = "Species") + 
  scale_fill_brewer(palette = "Dark2", name = "Species", guide = "none") +
  facet_wrap(~species+null, scales = "free_y", nrow = 5, labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))
  }) + 
  xlab("Precipitation Climate Model Skill Score") + 
  ylab("Log-Likelihood") + 
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
