################################################################################
# THIS SCRIPT PRODUCES FIGURE 3 IN THE FASTOVICH ET AL., (2023, PHIL B) THAT
# COMPARES SPATIAL ERROR MODEL INFORMATION LOSS AGAINST CHANGE IN THE STRENGTH
# OF THE AMOC. ESTIMATES FROM CHANGE IN THE STRENGTH OF THE AMOC ARE TAKEN
# DIRECTLY FROM THE LITERATURE THAT AGGREGATED THESE MODELS (FASTOVICH ET AL.,
# 2022, QSR AND KAGEYAMA ET AL., 2013, CLIMATE OF THE PAST). ESTIMATES OF AMOC
# REDUCTION IN THE PROXY RECORD ARE BASED ON ESTIMATES FROM RITZ ET AL., WHERE
# AMOC STRENGTH WAS CALCULATED FROM PROXY BASED SEA SURFACE TEMPERATURES AND
# CLIMATE SIMULATIONS. RITZ ET AL. (2013) ASSUMED THAT PROXY SEA SURFACE
# TEMPERATURES ARE A LINEAR COMBINATION OF AMOC STRENGTH AND GLOBAL MEAN
# ATMOSPHERIC TEMPERATURE.
# 
# INPUT:
# * SKILL_SCORES/AMOC.CSV
# 
# OUTPUT:
#   FIGURE 3:
#     * FIGURES/SYNTHESIS_AMOC_SAR.PDF
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

# Read in AMOC data taken directly from Kageyama et al., or calculated by me
# for TraCE-21ka, TraCE-MWF, HadCM3 , and CM2Mc
amoc <- read.csv("skill_scores/amoc.csv")

# AIC from the SEM models calculated here
aic <- read_csv("skill_scores/sar_diag.csv") %>%
  select(AIC_fit_50, AIC_fit_100, climate_simulation, species) %>% 
  rename(model = climate_simulation)

# Retain 50 km neighborhood AIC for all taxonomic groups except birds
aic <- aic %>% 
  mutate(AIC =
           case_when(
             species != "bird" ~ AIC_fit_50,
             species == "bird" ~ AIC_fit_100,
           )
  )


#######################
# RENAME SPECIES COLUMN
#######################

# Small function to make the species names nicer for plotting purposes
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

# Apply function
aic$species <- sapply(aic$species, rename_species, USE.NAMES = FALSE)

###########################
# COMBINE ALL DATA TOGETHER
###########################

# Join AMOC and AIC df's together
amoc_aic <- left_join(aic, amoc, by = "model")
amoc_aic$amoc_change <- amoc_aic$amoc_ref - amoc_aic$amoc_pert

# Setting value for AMOC change manually with data from Ritz et al., 2013, NGS
# Using manual approach because data is only given as anomalies from mean so
# while I dont konw ref state or perturb state I do know the change between the
# two.
amoc_aic$amoc_change[amoc_aic$model == "PROXY_KRIGING"] <- 11.013636

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

# Generate plot
amoc_plot <- ggplot(data = amoc_aic, aes(x = amoc_change, y = AIC, shape = model, color = species, fill = species)) + 
  # geom_ribbon(data = regression_line, aes(x = amoc_change, y = fit, ymin = lwr, ymax = upr), color = "transparent", alpha = 0.2) +
  # geom_line(data = regression_line, aes(x = amoc_change, y = fit, linetype = sig)) + 
  # scale_linetype_manual(values = c("dashed", "solid"), name = "Model Significance") + 
  geom_point() + 
  scale_shape_manual(values = 1:13, name = "Model") + 
  scale_color_brewer(palette = "Dark2", name = "Species") +
  scale_fill_brewer(palette = "Dark2", name = "Species", guide = "none") +
  facet_wrap(.~species, scales = "free_y", nrow = 5) +
  xlab("AMOC Reduction (Sv)") + 
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  scale_y_reverse() +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(vjust = -1),
        panel.spacing.y = unit(-.1, "cm"),
        legend.key = element_blank())

# Save plot as PDF and PNG
ggsave(amoc_plot, filename = "figures/synthesis_amoc_sar_quadratic.pdf", height = 8, width = 4.5, dpi = 300)
ggsave(amoc_plot, filename = "figures/synthesis_amoc_sar_quadratic.png", height = 8, width = 4.5, dpi = 300)
