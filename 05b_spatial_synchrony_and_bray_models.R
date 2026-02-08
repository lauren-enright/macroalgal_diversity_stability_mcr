rm(list = ls())
library(glmmTMB)
library(emmeans)
library(car)
library(effects)
library(DHARMa)
library(codyn)

source("05a_calculate_spatial_synchrony_and_bray.R")


#remind r of the correct order 
#want habitats to be in order.
dss_spatial_2$habitat <- factor(dss_spatial_2$habitat,
                                                levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

#### REGIONAL STABILITY ~ LOCAL STABILITY ####
alpha_gamma_stab_mod <- glmmTMB(cover_stability ~ stability_mean*habitat, family = Gamma("log"), data = dss_spatial_2)
summary(alpha_gamma_stab_mod)
car::Anova(alpha_gamma_stab_mod)
#                        Chisq Df Pr(>Chisq)    
#stability_mean         61.2725  1   4.97e-15 ***
#habitat                12.9323  3   0.004785 ** 
#stability_mean:habitat  9.4832  3   0.023511 *
performance::r2(alpha_gamma_stab_mod) # 0.838 #matches

  
em_alpha_gamma_stab_mod <- emtrends(alpha_gamma_stab_mod, pairwise ~ habitat, var = "stability_mean") # backreef different from fringe and forereef 10

# $contrasts
#  contrast                    estimate    SE  df z.ratio p.value
#   Fringing - Backreef            1.078 0.381 Inf   2.828  0.0242
#   Fringing - Forereef 10m        0.175 0.345 Inf   0.509  0.9571
#   Fringing - Forereef 17m        0.380 0.305 Inf   1.248  0.5965
#   Backreef - Forereef 10m       -0.903 0.349 Inf  -2.589  0.0475
#   Backreef - Forereef 17m       -0.698 0.309 Inf  -2.256  0.1085
#   Forereef 10m - Forereef 17m    0.205 0.263 Inf   0.779  0.8638



cld_em_alpha_gamma_stab_mod <- multcomp::cld(em_alpha_gamma_stab_mod, Letters = letters, sort = FALSE)

#### REGIONAL STABILITY ~ SPATIAL SYNCHRONY ####
spatial_synchrony_mod <- glmmTMB(cover_stability ~ spatial_synchrony*habitat, family = Gamma("log"), data = dss_spatial_2)
hist(residuals(spatial_synchrony_mod)) # a little skew
plot(residuals(spatial_synchrony_mod) ~ predict(spatial_synchrony_mod))
summary(spatial_synchrony_mod)
car::Anova(spatial_synchrony_mod)
  
#spatial_synchrony         16.2914  1  5.431e-05 ***
#habitat                   19.5303  3  0.0002124 ***
#spatial_synchrony:habitat  6.7239  3  0.0812368 .  
performance::r2(spatial_synchrony_mod) # 0.64 #this matches

## interaction is not signficant, re run 

spatial_synchrony_mod2 <- glmmTMB(cover_stability ~ spatial_synchrony + habitat, family = Gamma("log"), data = dss_spatial_2)
hist(residuals(spatial_synchrony_mod2)) # a little skew
plot(residuals(spatial_synchrony_mod2) ~ predict(spatial_synchrony_mod2))
summary(spatial_synchrony_mod2)
performance::r2(spatial_synchrony_mod2) # 0.54 #this matches
car::Anova(spatial_synchrony_mod2)

#Response: cover_stability
#Chisq Df Pr(>Chisq)    
#spatial_synchrony 11.483  1  0.0007024 ***
#habitat           16.655  3  0.0008321 ***


#BUT SHOULDN'T IT BE SINCE NO SIGN INTERACTION? 
# Yes, use this 
em_spatial_synchrony_mod2 <- emmeans(spatial_synchrony_mod2, pairwise ~ habitat, type = "response")

#$contrasts
#contrast                    ratio     SE  df null z.ratio p.value
#Fringing / Backreef         0.827 0.1294 Inf    1  -1.215  0.6172
#Fringing / Forereef 10m     0.655 0.1197 Inf    1  -2.314  0.0949
#Fringing / Forereef 17m     0.540 0.0849 Inf    1  -3.921  0.0005
#Backreef / Forereef 10m     0.792 0.1441 Inf    1  -1.279  0.5762
#Backreef / Forereef 17m     0.653 0.1026 Inf    1  -2.716  0.0335
#Forereef 10m / Forereef 17m 0.823 0.1441 Inf    1  -1.110  0.6833

cld_spatial_synchrony_mod2 <- multcomp::cld(em_spatial_synchrony_mod2, Letters = letters, sort = FALSE)

#### REGIONAL STABILITY : LOCAL STABILITY ~ SPATIAL SYNCHRONY ####

#run model, without logging response 
ratio_mod2 <- glmmTMB(ratio ~ spatial_synchrony, family = Gamma("log"), data = dss_spatial_2)

hist(residuals(ratio_mod2)) # THIS LOOKS WEIRD?
plot(residuals(ratio_mod2) ~ predict(ratio_mod2)) # looks fine I think?
summary(ratio_mod2)
car::Anova(ratio_mod2)

#Response: ratio
#Chisq Df Pr(>Chisq)    
#spatial_synchrony   126  1  < 2.2e-16 ***

#log to improve residuals
#this is logged - GOING TO USE THIS MODEL
ratio_mod3 <- glmmTMB(ratio ~ log(spatial_synchrony), family = Gamma("log"), data = dss_spatial_2)

hist(residuals(ratio_mod3)) # this looks better
plot(residuals(ratio_mod3) ~ predict(ratio_mod3)) # looks fine I think?
summary(ratio_mod3)
car::Anova(ratio_mod3)

#Response: ratio
#                       Chisq Df Pr(>Chisq)    
#log(spatial_synchrony) 774.21  1  < 2.2e-16 ***

performance::r2(ratio_mod3)  #0.97


#### SPATIAL SYNCHRONY ~ BRAY ####

bray_mod2 <- glmmTMB(spatial_synchrony ~ mean_bray, family = beta_family(), data = beta_spatialsync)
hist(residuals(bray_mod2)) # looks fine
plot(residuals(bray_mod2) ~ predict(bray_mod2)) # looks good
summary(bray_mod2)
car::Anova(bray_mod2)
#Response: spatial_synchrony
#Chisq Df Pr(>Chisq)    
#mean_bray 13.832  1  0.0001999 ***

performance::r2(bray_mod2) # 0.40 --> this got way worse.. 


#pull out data from models to make supplemental figures:

# get estimates for mean plot level stability
t.s6.a <- as_tibble(cld_em_alpha_gamma_stab_mod) %>%
  mutate(Predictor = "Mean plot-level stability") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

# get estimates for spatial synchrony
t.s6.b <- as_tibble(cld_spatial_synchrony_mod2) %>%
  mutate(Predictor = "Spatial synchrony") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean") %>%
  dplyr::rename(Mean = response) #note these are different things, but just want for plotting.. 


rbind(t.s6.a, t.s6.b) %>% 
  mutate(Table = "S6") %>% 
  dplyr::select(-SE, -df) -> table_s6

supplemental_tableS6 <- table_s6 %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
                Lower_CI = asymp.LCL,
                Upper_CI = asymp.UCL)
