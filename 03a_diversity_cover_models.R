rm(list = ls())
library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(effects)
library(DHARMa)
library(multcomp)
library(multcompView)

#Scripts starting with 03 focus on diversity/cover relationships. Need to have run script 01 to get necessary files. 

source("00_functions_and_aes.R")

alpha_diversity_quad_macro <- read.csv(here::here("data", "alpha_diversity_quad_macro_09262025.csv"))
#alpha_diversity_site_macro <- read.csv(here::here("data", "alpha_diversity_site_macro.csv"))

alpha_diversity_quad_macro$habitat <- factor(alpha_diversity_quad_macro$habitat,
                     levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

cor(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$functional_richness, use = "complete.obs", method = "pearson") #0.9584012
cor.test(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$functional_richness, method = "pearson")   # or "spearman"


#### COVER ~ RICHNESS ####

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
cover_mod_NEW <- glmmTMB(cover_trans ~ richness*habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro)
summary(cover_mod_NEW)
car::Anova(cover_mod_NEW)

#Chisq Df Pr(>Chisq)    
#richness         28840.4013  1    < 2e-16 ***
#habitat              8.0026  3    0.04596 *  
#richness:habitat   327.8860  3    < 2e-16 ***

performance::r2(cover_mod_NEW) 
#Conditional R2: 0.648
#Marginal R2: 0.626

#need to check singularity 
performance::check_singularity(cover_mod_NEW) 
#(1|site/site_habitat/location) = singular
#(1|site_habitat/location) = not singular, so confirming this is the random effects structure that we want 

em_cover_mod_NEW <- emtrends(cover_mod_NEW, pairwise ~ habitat, var = "richness") 
em_cover_mod_NEW
cld_cover_mod_new <- multcomp::cld(em_cover_mod_NEW$emtrends, Letters = letters, sort = FALSE)
#Fringing - Backreef          0.35748 0.0247 Inf  14.475 <0.0001
#Fringing - Forereef 10m      0.36451 0.0249 Inf  14.630 <0.0001
#Fringing - Forereef 17m      0.44403 0.0256 Inf  17.313 <0.0001
#Backreef - Forereef 10m      0.00704 0.0206 Inf   0.342  0.9863
#Backreef - Forereef 17m      0.08656 0.0214 Inf   4.051  0.0003
#Forereef 10m - Forereef 17m  0.07952 0.0213 Inf   3.741  0.0011

# habitat      richness.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing               1.73 0.0204 Inf      1.69      1.76  a    
#Backreef               1.37 0.0145 Inf      1.34      1.40   b   
#Forereef 10m           1.36 0.0148 Inf      1.33      1.39   b   
#Forereef 17m           1.28 0.0158 Inf      1.25      1.31    c  


#### COVER ~ FUNCTIONAL RICHNESS ####
#check number of functional richness categories
unique(alpha_diversity_quad_macro$functional_richness)

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
cover_mod_fg_new <- glmmTMB(cover_trans ~ functional_richness*habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro)

summary(cover_mod_fg_new)
car::Anova(cover_mod_fg_new)
#Response: cover_trans
#Chisq Df Pr(>Chisq)    
#functional_richness         30210.2372  1     <2e-16 ***
#habitat                         3.0836  3     0.3789    
#functional_richness:habitat   502.7920  3     <2e-16 ***

#want to check these against each other
hist(residuals(cover_mod_fg)) # good
hist(residuals(cover_mod_fg_new)) # good

#need to check singularity 
performance::check_singularity(cover_mod_fg_new)
##(1|site/site_habitat/location) = singular
#(1|site_habitat/location) = not singular 

performance::r2(cover_mod_fg_new)
#Conditional R2: 0.642
#Marginal R2: 0.622

em_cover_mod_fg_new <- emtrends(cover_mod_fg_new, pairwise ~ habitat, var = "functional_richness") 
#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef            0.236 0.0269 Inf   8.748 <0.0001
#Fringing - Forereef 10m        0.358 0.0268 Inf  13.372 <0.0001
#Fringing - Forereef 17m        0.575 0.0268 Inf  21.494 <0.0001
#Backreef - Forereef 10m        0.122 0.0231 Inf   5.285 <0.0001
#Backreef - Forereef 17m        0.340 0.0231 Inf  14.698 <0.0001
#Forereef 10m - Forereef 17m    0.217 0.0225 Inf   9.667 <0.0001

cld_cover_mod_fg_new <- multcomp::cld(em_cover_mod_fg_new$emtrends, Letters = letters, sort = FALSE)

#habitat      functional_richness.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                          1.88 0.0216 Inf      1.84      1.93  a    
#Backreef                          1.65 0.0167 Inf      1.62      1.68   b   
#Forereef 10m                      1.53 0.0163 Inf      1.49      1.56    c  
#Forereef 17m                      1.31 0.0162 Inf      1.28      1.34     d 


#pull out data from models to make supplemental figures:

#richness by cover
t.s2.a <- as_tibble(cld_cover_mod_new) %>%
  mutate(Predictor = "Taxonomic richness",
         `Spatial scale` = "Plot-level") %>% 
  # rename this so it matches for joining purposes
  rename_if(str_detect(names(.), ".trend"), ~"Mean")
#functional richness by cover

t.s2.b <-  as_tibble(cld_cover_mod_fg_new) %>%
  mutate(Predictor = "Functional richness",
         `Spatial scale` = "Plot-level") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

#joining together

rbind(t.s2.a, t.s2.b) %>% 
  mutate(Table = "S2_plotlevel") %>% 
  dplyr::select(-SE, -df) -> table_s2_plot

supplemental_tables_tableS2_plot <- table_s2_plot %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
         Lower_CI = asymp.LCL,
         Upper_CI = asymp.UCL)



