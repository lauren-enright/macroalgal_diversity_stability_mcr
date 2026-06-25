rm(list = ls())
library(glmmTMB)
library(emmeans)
library(car)
library(effects)
library(DHARMa)
library(codyn)
library(multcomp)
library(here)

source("00_functions_and_aes.R")
#source("04a_calculate_synchrony.R")

#read in data produced in 04a_calculate_synchrony script
diversity_stability_synchrony <- read.csv(here::here("data", "full_plot_level_dss_09262025.csv")) 
diversity_stability_synchrony_site <- read.csv(here::here("data", "diversity_stability_synchrony_site_09262025.csv"))

diversity_stability_synchrony$site_habitat <- paste(diversity_stability_synchrony$site,
                                                    diversity_stability_synchrony$habitat,
                                                    sep = "_")

#want habitats to be in order.
diversity_stability_synchrony$habitat <- factor(diversity_stability_synchrony$habitat,
                                             levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

# EXTRACT RANGES
dss_ranges <- extract_ranges(diversity_stability_synchrony, "habitat",
                             c("richness_mean", "functional_richness_mean", "cover_stability", "synchrony"))


#### STABILITY ~ RICHNESS ####
#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
rich_stab_mod <- glmmTMB(cover_stability~ richness_mean*habitat + (1|site_habitat), family = Gamma(link = "log"), data = diversity_stability_synchrony)
performance::check_singularity(rich_stab_mod)
# (1|site/site_habitat) --> singular
# (1|site_habitat) --> ok 
summary(rich_stab_mod)
car::Anova(rich_stab_mod)
#Response: cover_stability
#Chisq Df Pr(>Chisq)    
#richness_mean         1427.967  1  < 2.2e-16 ***
#habitat                 31.588  3   6.39e-07 ***
#richness_mean:habitat   39.739  3   1.21e-08 ***
performance::r2(rich_stab_mod) 
# Conditional R2: 0.840
#Marginal R2: 0.793
performance::check_model(rich_stab_mod)

em_rich_stab_mod <- emtrends(rich_stab_mod, pairwise ~ habitat, var = "richness_mean")
cld_rich_stab_mod <- multcomp::cld(em_rich_stab_mod, Letters = letters, sort = FALSE)

#habitat      richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                   0.906 0.0423 Inf     0.823     0.989  a    
#Backreef                   0.666 0.0309 Inf     0.606     0.727   bc  
#Forereef 10m               0.552 0.0401 Inf     0.474     0.631   b   
#Forereef 17m               0.750 0.0403 Inf     0.671     0.829    c  

#### STABILITY ~ FUNCTIONAL RICHNESS ####
  
#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
rich_stab_mod_fg <- glmmTMB(cover_stability ~ functional_richness_mean*habitat + (1|site_habitat), family = Gamma(link = "log"), data = diversity_stability_synchrony)
car::Anova(rich_stab_mod_fg)
#Response: cover_stability
#Chisq Df Pr(>Chisq)    
#functional_richness_mean         1639.776  1  < 2.2e-16 ***
#habitat                            12.746  3    0.00522 ** 
#functional_richness_mean:habitat   42.057  3  3.902e-09 ***



performance::check_singularity(rich_stab_mod_fg) #this is not singular! 
performance::check_model(rich_stab_mod_fg)
performance::r2(rich_stab_mod_fg) 
#  Conditional R2: 0.855
#Marginal R2: 0.807

em_rich_stab_mod_fg <- emtrends(rich_stab_mod_fg, pairwise ~ habitat, var = "functional_richness_mean") 
#$contrasts
#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef           0.2621 0.0598 Inf   4.383 <0.0001
#Fringing - Forereef 10m       0.4090 0.0669 Inf   6.116 <0.0001
#Fringing - Forereef 17m       0.3277 0.0627 Inf   5.230 <0.0001
#Backreef - Forereef 10m       0.1469 0.0585 Inf   2.511  0.0583
#Backreef - Forereef 17m       0.0656 0.0536 Inf   1.224  0.6115
#Forereef 10m - Forereef 17m  -0.0813 0.0606 Inf  -1.343  0.5355

cld_rich_stab_mod_fg <- multcomp::cld(em_rich_stab_mod_fg, Letters = letters, sort = FALSE)


#habitat      functional_richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                              1.099 0.0478 Inf     1.005     1.192  a    
#Backreef                              0.836 0.0359 Inf     0.766     0.907   b   
#Forereef 10m                          0.690 0.0462 Inf     0.599     0.780   b   
#Forereef 17m                          0.771 0.0398 Inf     0.693     0.849   b   


#### STABILITY ~ SYNCHRONY ####
#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
synch_stab_mod <- glmmTMB(cover_stability ~ synchrony*habitat + (1|site_habitat), family = Gamma("log"), data = diversity_stability_synchrony)
performance::check_singularity(synch_stab_mod) #false

performance::check_model(synch_stab_mod)
summary(synch_stab_mod)
car::Anova(synch_stab_mod)
#                      Chisq Df Pr(>Chisq)    
#synchrony         1080.0358  1    < 2e-16 ***
#habitat              7.2496  3    0.06435 .  
#synchrony:habitat  163.9107  3    < 2e-16 ***

performance::r2(synch_stab_mod)
#Conditional R2: 0.815
#Marginal R2: 0.575

em_synch_stab_mod <- emtrends(synch_stab_mod, pairwise ~ habitat, var = "synchrony")


#Fringing - Backreef            0.415 0.0857 Inf   4.842 <0.0001
#Fringing - Forereef 10m        1.124 0.1450 Inf   7.764 <0.0001
#Fringing - Forereef 17m        1.605 0.1440 Inf  11.185 <0.0001
#Backreef - Forereef 10m        0.709 0.1490 Inf   4.760 <0.0001
#Backreef - Forereef 17m        1.190 0.1480 Inf   8.058 <0.0001
#Forereef 10m - Forereef 17m    0.481 0.1880 Inf   2.556  0.0518

cld_synch_stab_mod <- multcomp::cld(em_synch_stab_mod, Letters = letters, sort = FALSE)

#habitat      synchrony.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing              -0.895 0.0554 Inf     -1.00    -0.786  a    
#Backreef              -1.310 0.0655 Inf     -1.44    -1.182   b   
#Forereef 10m          -2.019 0.1340 Inf     -2.28    -1.757    c  
#Forereef 17m          -2.500 0.1320 Inf     -2.76    -2.241    c  


#pull out data from models to make supplemental figures:

# get estimates for taxonomic richness
t.s4.a <- as.tibble(cld_rich_stab_mod) %>%
  mutate(Predictor = "Taxonomic richness") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

# get estimates for fucntional richness
t.s4.b <- as.tibble(cld_rich_stab_mod_fg) %>%
  mutate(Predictor = "Functional richness") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

# get estimates for synchrony
t.s4.c <- as.tibble(cld_synch_stab_mod) %>% 
  mutate(Predictor = "Synchrony") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

rbind(t.s4.a, t.s4.b, t.s4.c) %>% 
  mutate(`Spatial scale` = "Plot-level", 
         Table = "S4") %>% 
  dplyr::select(-SE, -df) -> table_s4_plotlevel

supplemental_tables_tableS4_plot <- table_s4_plotlevel %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
                Lower_CI = asymp.LCL,
                Upper_CI = asymp.UCL)


