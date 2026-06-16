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
rich_stab_mod <- glmmTMB(cover_stability~ richness_mean*habitat + (1|site), family = Gamma(link = "log"), data = diversity_stability_synchrony)
summary(rich_stab_mod)
car::Anova(rich_stab_mod)

#great, matches
#                        Chisq Df Pr(>Chisq)    
#richness_mean         2721.440  1  < 2.2e-16 ***
#habitat                325.876  3  < 2.2e-16 ***
#richness_mean:habitat   33.926  3  2.054e-07 ***

performance::check_model(rich_stab_mod)

hist(residuals(rich_stab_mod)) # looks good
plot(residuals(rich_stab_mod) ~ fitted(rich_stab_mod)) # heteroskedastic - bring to group
#r.squaredGLMM(rich_stab_mod) # 0.803, 0.808
performance::r2(rich_stab_mod) # marginal: 0.812, conditional: 0.816 - yes, matches
em_rich_stab_mod <- emtrends(rich_stab_mod, pairwise ~ habitat, var = "richness_mean") # fringing not diff from outer17, backreef not different from outer10

#contrasts
#contrast                     estimate     SE  df z.ratio p.value
#Fringing - Backreef          0.186671 0.0414 Inf   4.512  <.0001
#Fringing - Forereef 10m      0.187196 0.0455 Inf   4.114  0.0002
#Fringing - Forereef 17m      0.027111 0.0463 Inf   0.586  0.9364
#Backreef - Forereef 10m      0.000525 0.0405 Inf   0.013  1.0000
#Backreef - Forereef 17m     -0.159559 0.0411 Inf  -3.885  0.0006
#Forereef 10m - Forereef 17m -0.160084 0.0439 Inf  -3.647  0.0015

cld_rich_stab_mod <- multcomp::cld(em_rich_stab_mod, Letters = letters, sort = FALSE)
# habitat      richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                   0.877 0.0326 Inf     0.813     0.941  a    
#Backreef                   0.690 0.0252 Inf     0.641     0.740   b   
#Forereef 10m               0.690 0.0315 Inf     0.628     0.752   b   
#Forereef 17m               0.850 0.0317 Inf     0.788     0.912  a    


#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
rich_stab_mod_NEW <- glmmTMB(cover_stability~ richness_mean*habitat + (1|site_habitat), family = Gamma(link = "log"), data = diversity_stability_synchrony)
performance::check_singularity(rich_stab_mod_NEW)
# (1|site/site_habitat) --> singular
# (1|site_habitat) --> ok 
summary(rich_stab_mod_NEW)
car::Anova(rich_stab_mod_NEW)
#Response: cover_stability
#Chisq Df Pr(>Chisq)    
#richness_mean         1427.967  1  < 2.2e-16 ***
#habitat                 31.588  3   6.39e-07 ***
#richness_mean:habitat   39.739  3   1.21e-08 ***
performance::r2(rich_stab_mod_NEW) 
# Conditional R2: 0.840
#Marginal R2: 0.793
performance::check_model(rich_stab_mod_NEW)

em_rich_stab_mod_NEW <- emtrends(rich_stab_mod_NEW, pairwise ~ habitat, var = "richness_mean")
cld_rich_stab_mod_NEW <- multcomp::cld(em_rich_stab_mod_NEW , Letters = letters, sort = FALSE)

#habitat      richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                   0.906 0.0423 Inf     0.823     0.989  a    
#Backreef                   0.666 0.0309 Inf     0.606     0.727   bc  
#Forereef 10m               0.552 0.0401 Inf     0.474     0.631   b   
#Forereef 17m               0.750 0.0403 Inf     0.671     0.829    c  

#### STABILITY ~ FUNCTIONAL RICHNESS ####
rich_stab_mod_fg <- glmmTMB(cover_stability ~ functional_richness_mean*habitat + (1|site), family = Gamma(link = "log"), data = diversity_stability_synchrony)

summary(rich_stab_mod_fg)
car::Anova(rich_stab_mod_fg)

#new, with 8 functional groups

#Chisq Df Pr(>Chisq)    
#functional_richness_mean         3104.718  1  < 2.2e-16 ***
#habitat                           119.466  3  < 2.2e-16 ***
# functional_richness_mean:habitat   24.858  3  1.653e-05 ***

performance::check_model(rich_stab_mod_fg)

hist(residuals(rich_stab_mod_fg)) # looks fine
plot(residuals(rich_stab_mod_fg) ~ fitted(rich_stab_mod_fg)) # same as usual 
performance::r2(rich_stab_mod_fg) 
#new for 8 functional groups: #Conditional R2: 0.833, Marginal R2: 0.824
em_rich_stab_mod_fg <- emtrends(rich_stab_mod_fg, pairwise ~ habitat, var = "functional_richness_mean") 
#  fringing different from all others

#new, with 8 functional groups
#$contrasts
#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef           0.1791 0.0467 Inf   3.836  0.0007
#Fringing - Forereef 10m       0.2337 0.0517 Inf   4.521  <.0001
#Fringing - Forereef 17m       0.1979 0.0488 Inf   4.056  0.0003
#Backreef - Forereef 10m       0.0546 0.0471 Inf   1.159  0.6525
#Backreef - Forereef 17m       0.0189 0.0439 Inf   0.429  0.9735
#Forereef 10m - Forereef 17m  -0.0357 0.0475 Inf  -0.753  0.8755

cld_rich_stab_mod_fg <- multcomp::cld(em_rich_stab_mod_fg, Letters = letters, sort = FALSE)

#habitat      functional_richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                              1.064 0.0359 Inf     0.994     1.134  a    
#Backreef                              0.885 0.0296 Inf     0.827     0.943   b   
#Forereef 10m                          0.830 0.0370 Inf     0.758     0.903   b   
#Forereef 17m                          0.866 0.0314 Inf     0.804     0.928   b   

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
rich_stab_mod_fg_new <- glmmTMB(cover_stability ~ functional_richness_mean*habitat + (1|site_habitat), family = Gamma(link = "log"), data = diversity_stability_synchrony)
car::Anova(rich_stab_mod_fg_new)
#Response: cover_stability
#Chisq Df Pr(>Chisq)    
#functional_richness_mean         1639.776  1  < 2.2e-16 ***
#habitat                            12.746  3    0.00522 ** 
#functional_richness_mean:habitat   42.057  3  3.902e-09 ***



performance::check_singularity(rich_stab_mod_fg_new) #this is not singular! 
performance::check_model(rich_stab_mod_fg_new)
performance::r2(rich_stab_mod_fg_new) 
#  Conditional R2: 0.855
#Marginal R2: 0.807

em_rich_stab_mod_fg_new <- emtrends(rich_stab_mod_fg_new, pairwise ~ habitat, var = "functional_richness_mean") 
cld_rich_stab_mod_fg_new <- multcomp::cld(em_rich_stab_mod_fg_new, Letters = letters, sort = FALSE)


#habitat      functional_richness_mean.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                              1.099 0.0478 Inf     1.005     1.192  a    
#Backreef                              0.836 0.0359 Inf     0.766     0.907   b   
#Forereef 10m                          0.690 0.0462 Inf     0.599     0.780   b   
#Forereef 17m                          0.771 0.0398 Inf     0.693     0.849   b   

#check again for NOAM
rich_stab_mod_fg_new2 <- glmmTMB(cover_stability ~ functional_richness_mean*habitat + (1|site/habitat), family = Gamma(link = "log"), data = diversity_stability_synchrony)
performance::check_singularity(rich_stab_mod_fg_new2)
#(1|site/habitat) is singular

#Looking at these 3 model structures -- AIC VS BIC
AIC(rich_stab_mod_fg, rich_stab_mod_fg_new, rich_stab_mod_fg_new2)
#df       AIC
#rich_stab_mod_fg      10 -491.7148
#rich_stab_mod_fg_new  10 -678.4380 #rich_stab_mod_fg_new IS BEST
#rich_stab_mod_fg_new2 11 -676.4380

BIC(rich_stab_mod_fg, rich_stab_mod_fg_new, rich_stab_mod_fg_new2)
#df       BIC
#rich_stab_mod_fg      10 -440.9482
#rich_stab_mod_fg_new  10 -627.6715 #rich_stab_mod_fg_new IS BEST
#rich_stab_mod_fg_new2 11 -620.5948


#### STABILITY ~ SYNCHRONY ####
synch_stab_mod <- glmmTMB(cover_stability ~ synchrony*habitat + (1|site), family = Gamma("log"), data = diversity_stability_synchrony)
summary(synch_stab_mod)
car::Anova(synch_stab_mod)

#great, matches.  
#synchrony         1308.33  1  < 2.2e-16 ***
#  habitat            103.13  3  < 2.2e-16 ***
#  synchrony:habitat  327.97  3  < 2.2e-16 ***

hist(residuals(synch_stab_mod)) # looks good
plot(residuals(synch_stab_mod) ~ fitted(synch_stab_mod)) # same issue as always

#r.squaredGLMM(synch_stab_mod) # 0.63, 0.71
performance::r2(synch_stab_mod) # marginal: 0.66, conditional: 0.74 - matches
em_synch_stab_mod <- emtrends(synch_stab_mod, pairwise ~ habitat, var = "synchrony") # fringing and backreef not different
# contrast                    estimate    SE  df z.ratio p.value
# Fringing - Backreef            0.181 0.107 Inf   1.695  0.3262
# Fringing - Forereef 10m        1.492 0.142 Inf  10.487  <.0001
# Fringing - Forereef 17m        2.436 0.149 Inf  16.336  <.0001
# Backreef - Forereef 10m        1.311 0.157 Inf   8.346  <.0001
# Backreef - Forereef 17m        2.255 0.161 Inf  14.008  <.0001
# Forereef 10m - Forereef 17m    0.944 0.181 Inf   5.222  <.0001
#matches

cld_synch_stab_mod <- multcomp::cld(em_synch_stab_mod, Letters = letters, sort = FALSE)
# habitat      synchrony.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing               -1.22 0.0638 Inf     -1.35     -1.10  a    
#Backreef               -1.41 0.0849 Inf     -1.57     -1.24  a    
#Forereef 10m           -2.72 0.1300 Inf     -2.97     -2.46   b   
#Forereef 17m           -3.66 0.1380 Inf     -3.93     -3.39    c  

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
synch_stab_mod_NEW <- glmmTMB(cover_stability ~ synchrony*habitat + (1|site_habitat), family = Gamma("log"), data = diversity_stability_synchrony)
performance::check_singularity(synch_stab_mod_NEW) #false

performance::check_model(synch_stab_mod_NEW)
summary(synch_stab_mod_NEW)
car::Anova(synch_stab_mod_NEW)
#                      Chisq Df Pr(>Chisq)    
#synchrony         1080.0358  1    < 2e-16 ***
#habitat              7.2496  3    0.06435 .  
#synchrony:habitat  163.9107  3    < 2e-16 ***

performance::r2(synch_stab_mod_NEW)
#Conditional R2: 0.815
#Marginal R2: 0.575

em_synch_stab_mod_NEW <- emtrends(synch_stab_mod_NEW, pairwise ~ habitat, var = "synchrony")


#Fringing - Backreef            0.415 0.0857 Inf   4.842 <0.0001
#Fringing - Forereef 10m        1.124 0.1450 Inf   7.764 <0.0001
#Fringing - Forereef 17m        1.605 0.1440 Inf  11.185 <0.0001
#Backreef - Forereef 10m        0.709 0.1490 Inf   4.760 <0.0001
#Backreef - Forereef 17m        1.190 0.1480 Inf   8.058 <0.0001
#Forereef 10m - Forereef 17m    0.481 0.1880 Inf   2.556  0.0518

cld_synch_stab_mod_new <- multcomp::cld(em_synch_stab_mod_NEW, Letters = letters, sort = FALSE)

#habitat      synchrony.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing              -0.895 0.0554 Inf     -1.00    -0.786  a    
#Backreef              -1.310 0.0655 Inf     -1.44    -1.182   b   
#Forereef 10m          -2.019 0.1340 Inf     -2.28    -1.757    c  
#Forereef 17m          -2.500 0.1320 Inf     -2.76    -2.241    c  

#check for Noam
synch_stab_mod_NEW2 <- glmmTMB(cover_stability ~ synchrony*habitat + (1|site/habitat), family = Gamma("log"), data = diversity_stability_synchrony)
performance::check_singularity(synch_stab_mod_NEW2)
#this is singular

AIC(synch_stab_mod, synch_stab_mod_NEW, synch_stab_mod_NEW2)
#                    df       AIC
#synch_stab_mod      10  137.7269
#synch_stab_mod_NEW  10 -457.2483 #synch_stab_mod_NEW is BEST
#synch_stab_mod_NEW2 11 -455.2483
BIC(synch_stab_mod, synch_stab_mod_NEW, synch_stab_mod_NEW2)
#                    df       BIC
#synch_stab_mod      10  188.4935
#synch_stab_mod_NEW  10 -406.4818 #synch_stab_mod_NEW is BEST
#synch_stab_mod_NEW2 11 -399.4051


#pull out data from models to make supplemental figures:

# get estimates for taxonomic richness
t.s4.a <- as.tibble(cld_rich_stab_mod_NEW) %>%
  mutate(Predictor = "Taxonomic richness") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

# get estimates for fucntional richness
t.s4.b <- as.tibble(cld_rich_stab_mod_fg_new) %>%
  mutate(Predictor = "Functional richness") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

# get estimates for synchrony
t.s4.c <- as.tibble(cld_synch_stab_mod_new) %>% 
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


