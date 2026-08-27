## Diversity is positively correlated with macroalgal abundance and stability on coral reefs across decadal time scales

Data analysis and visualization in support of for Enright et al. -  Diversity is positively correlated with macroalgal abundance and stability on coral reefs across decadal time scales (in review at Ecology)

### Description 
This project analyzes data from a 17-year time series of benthic macroalgae on coral reefs in Moorea, French Polynesia. Scripts are numbered in the order in which they should be run.

Code is achieved at  [will place link to zenodo here]

### Contents 

- `00_functions_and_aes.R` - defines functions and figure aesthetics used throughout the analyses
- `01_data_diversitycalculations.R` -  calculates diversity metrics
- `02_time_series_plot.R` - generates Figure 2 and Figure S2
- `03a_diversity_cover_models.R` - analyzes relationships between diversity and abundance at the plot level
- `03b_diversity_cover_figures.R` - generates Figure 3
- `03c_diversity_cover_site_models.R` - analyzes relationships between diversity and abundance at the site level
- `03d_diversity_cover_site_figures.R` - generates Figure S3
- `04a_calculate_synchrony.R` - calculates synchrony metrics 
- `04b_dss_models.R` - analyzes relationships among diversity, species synchrony, and stability at the plot level
- `04c_dss_figures.R` - generates Figure 4 and Figure S4
- `04d_site-level_dss_models.R` - analyzes relationships among diversity, species synchrony, and stability at the site level
- `04e_site-level_dss_figures.R` - generates Figure S5
- `05a_calculate_spatial_synchrony_and_bray.R` - calculates spatial synchrony and Bray-Curtis dissimilarity metrics 
- `05b_spatial_synchrony_and_bray_models.R` - analyzes relationships among spatial synchrony, site-level stability, and beta diversity
- `05c_spatial_synchrony_and_bray_figures.R` -  generates Figure 5 and Figure S6
- `06a_supplements.R` - generates Figure S1

### Data
Data can be found at [will place link to EDI here]

Column Names

- taxa: taxonomy of the observed algae, or name of the observed substrate or functional group.
- year: year of observation.
- location: name of the unique location combines the site, habitat, transect and quadrat.
- site: name of the LTER location/LTER site.
- habitat: name of the habitat.
- site_habitat: name of the site habitat combines the LTER location/LTER site and habitat, for 24 unique sites
- transect: transect number.
- quadrat: quadrat number identifies the 0.25 square meter quadrat.
- prop_cover: percent of area covered in a quadrat, expressed as a proportion. For example, 1 = 100 % cover and 0.4 = 40% cover.
- functional_group: functional grouping of the taxa, based on existing literature or AlgaeBase.

Additional metadata is available by following the EDI link. 

