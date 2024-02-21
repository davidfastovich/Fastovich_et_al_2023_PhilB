#####################
# PREPARE ENVIRONMENT
#####################

# Install micromamba
curl micro.mamba.pm/install.sh | zsh

# Create environment from provided environment.yml
micromamba env create --file=~/Fastovich_et_al_2023_PhilB/environment/environment_combined.yml

# Install older versions of two packages not stored in conda-forge
Rscript environment/require.R

##########
# ANALYSES
##########

echo "Figure 1"

python model_0.5degree_bilinear_interpolation.py
python figure_1_fingerprint_fulcrum_on_bottom.py

echo "Figure 2, Figure S2, Figure S3"

# This code chunk will fail because I am not allowed to redistribute range
# maps.
python amphibian_richness.py
python bird_richness.py
python mammalian_richness.py
python reptilian_richness.py
python tree_richness.py
python model_50km_bilinear_interpolation.py
python proxy_kriging.py
python modern_data_reproject.py
Rscript sem_paleoclimate_effect.R

echo "Figure 3"
Rscript amoc_vs_sem_loglik.R

echo "Figure S1"
python richness_maps.py

echo "Figure S4, Figure S5"

python tas_model_skill.py
python pr_model_skill.py
Rscript climate_model_skill_vs_sem_loglik.R

echo "Figure S6, Figure S7"

Rscript sem_modern_climate_effect.R

echo "Figure S8"

# This code chunk will fail because I am not allowed to redistribute range
# maps.
python resolution_sensitivity/100km/amphibian_richness.py
python resolution_sensitivity/100km/bird_richness.py
python resolution_sensitivity/100km/mammalian_richness.py
python resolution_sensitivity/100km/reptilian_richness.py
python resolution_sensitivity/100km/tree_richness.py
python resolution_sensitivity/100km/richness_maps.py #DONE

echo "Figure S9, Figure S10"

python resolution_sensitivity/100km/proxy_kriging.py
python resolution_sensitivity/100km/model_100km_bilinear_interpolation.py
python resolution_sensitivity/100km/modern_data_reproject.py
Rscript resolution_sensitivity/100km/sem_paleoclimate_effect.R

echo "Figure S11, Figure S12"

Rscript resolution_sensitivity/100km/sem_modern_climate_effect.R

echo "Figure S13"

Rscript resolution_sensitivity/100km/amoc_vs_sem_loglik.R

echo "Figure S14, Figure S15"

python plot_paleoclimate.py

echo "Figure S16"

Rscript model_correlogram.R
