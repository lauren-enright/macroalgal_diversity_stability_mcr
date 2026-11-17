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
library(MuMIn)

#Scripts starting with 03 focus on diversity/cover relationships. Need to have run script 01 to get necessary files. 

source("00_functions_and_aes.R")

#### Read in Data #### 

alpha_diversity_quad_macro <- read.csv(here::here("data", "alpha_diversity_quad_macro_09262025.csv"))
#alpha_diversity_site_macro <- read.csv(here::here("data", "alpha_diversity_site_macro.csv"))

unique(alpha_diversity_quad_macro$habitat)
unique(alpha_diversity_quad_macro$site_habitat) #24, this is our site-level random effect

#### Remove plots with 0s and log data ####

alpha_diversity_quad_macro_present <- alpha_diversity_quad_macro %>%
  dplyr::filter(macroalgal_prop_cover > 0) %>%
  mutate(cover_trans_no0s = sv_trans(macroalgal_prop_cover), #slightly different numbers than if transforming, because lower sample size 
         log_rich = log(richness),
         log_fun  = log(functional_richness))

alpha_diversity_quad_macro_present$habitat <- factor(alpha_diversity_quad_macro_present$habitat,
                                             levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

#### correlations for methods/discussion ####

#all plots
cor.test(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$functional_richness, method = "pearson")   # or "spearman"
#0.9584012 
#zeros removed
cor.test(alpha_diversity_quad_macro_present$richness, alpha_diversity_quad_macro_present$functional_richness, method = "pearson")   # or "spearman"
#0.8869669 

unique(alpha_diversity_quad_macro_present$cover_trans_no0s) #make sure no exact 1s!


#### Taxonomic richness and cover ####


#####  M1: interaction, no log transformation of richness, no dispformula = habitat ##### 

#first trying linear model
taxo_richness_cover_m1 <- glmmTMB(
  cover_trans_no0s ~ richness * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)


##### M2: interaction, no log transformation of richness, yes dispformula = habitat #####
#adding disformula = cover be more variable in some habitats than others 

taxo_richness_cover_m2 <- glmmTMB(
  cover_trans_no0s ~ richness * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(taxo_richness_cover_m2)

res_taxo_lin <- simulateResiduals(taxo_richness_cover_m2, n = 1000)
plotResiduals(taxo_richness_cover_m2, form = alpha_diversity_quad_macro_present$richness, rank = FALSE)
#looks hump shaped? overiftting at low values and underfiting and high values

##### M3: interaction, log transformation of richness, no dispformula = habitat #####  

#trying logged richness without habitat disp 
taxo_richness_cover_m3 <- glmmTMB(
  cover_trans_no0s ~ log_rich * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(taxo_richness_cover_m3)
testDispersion(taxo_richness_cover_m3)
res_taxo_log_NODISP <- simulateResiduals(taxo_richness_cover_m3, n = 1000)
plot(res_taxo_log_NODISP)

plotResiduals(res_taxo_log_NODISP, form = alpha_diversity_quad_macro_present$habitat)

#####   M4: interaction, log transformation of richness, yes dispformula = habitat #####  

#trying logged richness with habitat disp 
taxo_richness_cover_m4 <- glmmTMB(
  cover_trans_no0s ~ log_rich * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

unique(alpha_diversity_quad_macro_present$site_habitat)

##### Let's see which model is best #####
AIC(taxo_richness_cover_m1, taxo_richness_cover_m2, taxo_richness_cover_m3, taxo_richness_cover_m4)
#                       df       AIC
#taxo_richness_cover_m1 12 -24299.40
#taxo_richness_cover_m2 15 -25421.75
#taxo_richness_cover_m3 12 -24426.14
#taxo_richness_cover_m4 15 -25575.60


#taxo_richness_cover_m4 (logged with habitat dispersion is best!)

#R2 does not work with dispersion sub model (disp habitat) 
MuMIn::r.squaredGLMM(model_occ_log_taxo)
summary(model_occ_log_taxo)
performance::r2_nakagawa(model_occ_log_taxo)
performance::model_performance(model_occ_log_taxo)

##### checks for best model (taxo_richness_cover_m4):##### 
taxo_richness_cover_m4$fit$convergence           # should be 0 - yes
taxo_richness_cover_m4$sdr$pdHess #TRUE
diagnose(taxo_richness_cover_m4) #Unusually large Z-statistics, but ok 
VarCorr(taxo_richness_cover_m4) # all signficantly above 0 
#Groups                Name        Std.Dev.
#location:site_habitat (Intercept) 0.15807 
#site_habitat          (Intercept) 0.23880 
#year                  (Intercept) 0.13890 
testUniformity(taxo_richness_cover_m4) #looks ok 
res_taxo_log <- simulateResiduals(taxo_richness_cover_m4, n = 1000)
plot(res_taxo_log) #looks ok 
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$richness, rank = FALSE) #not hump shaped
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$richness,
              rank = FALSE, asFactor = FALSE) #ok
#throws warning but just bc we only have 6 richness 
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$habitat) #ok
testDispersion(res_taxo_log)  #slightly high but ok?  #make sure to run this on the residuals

performance::check_model(taxo_richness_cover_m4)

hist(residuals(taxo_richness_cover_m4)) #ok
hist(residuals(res_taxo_log), breaks = 25) #blocky, but expected

testOutliers(res_taxo_log, type = "bootstrap")

testDispersion(recalculateResiduals(res_taxo_log, group = interaction(alpha_diversity_quad_macro_present$site_habitat, alpha_diversity_quad_macro_present$location)))

##### Anova and PostHoc for Best Model #####
car::Anova(taxo_richness_cover_m4, type = "II")
#Chisq Df Pr(>Chisq)    
#log_rich         4236.3045  1    < 2e-16 ***
#habitat             5.8679  3    0.11822    
#log_rich:habitat   11.0501  3    0.01146 *  
car::Anova(taxo_richness_cover_m3, type = "II") #just to check, is also significant interaction w/o dispersion


#Run likelihood ratio test of whether the richness × habitat interaction improves the taxonomic richness model. 
anova(update(taxo_richness_cover_m4, . ~ . - log_rich:habitat), taxo_richness_cover_m4)
#also suggests keeping interaction, but not overwhelming different

pairs(emtrends(taxo_richness_cover_m4, ~ habitat, var = "log_rich"))

#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef          -0.1780 0.0564 Inf  -3.157  0.0087
#Fringing - Forereef 10m      -0.1212 0.0535 Inf  -2.267  0.1059
#Fringing - Forereef 17m      -0.1454 0.0499 Inf  -2.917  0.0186
#Backreef - Forereef 10m       0.0568 0.0460 Inf   1.236  0.6041
#Backreef - Forereef 17m       0.0326 0.0418 Inf   0.781  0.8631
#Forereef 10m - Forereef 17m  -0.0242 0.0371 Inf  -0.652  0.9150

em_cover_mod <- emtrends(taxo_richness_cover_m4, pairwise ~ habitat, var = "log_rich") 

em_cover_mod #yep, same as above

# will need this for supplemental figure 
cld_cover_mod_richness <- multcomp::cld(em_cover_mod$emtrends, Letters = letters, sort = FALSE)

#### Functional Richness and Cover ####

###### M1: interaction, no log transformation of richness, no dispformula = habitat ######
#first trying linear with no dispersion
functional_richness_cover_m1 <- glmmTMB(
  cover_trans_no0s ~ functional_richness * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(functional_richness_cover_m1)

###### M2: interaction, no log transformation of richness, yes dispformula = habitat ######
#then trying linear with WITH dispersion
functional_richness_cover_m2 <- glmmTMB(
  cover_trans_no0s ~ functional_richness * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(functional_richness_cover_m2)
res_taxo_lin_fun <- simulateResiduals(functional_richness_cover_m2, n = 1000)
plotResiduals(res_taxo_lin_fun, form = alpha_diversity_quad_macro_present$functional_richness, rank = FALSE)
#humped shaped.... 

##### M3: interaction, log transformation of richness, no dispformula = habitat #####  

#now trying logging functiona richness WITHOUT DISPERSION
functional_richness_cover_m3 <- glmmTMB(
  cover_trans_no0s ~ log_fun * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

######   M4: interaction, log transformation of richness, yes dispformula = habitat #####  

#now trying logging functional richness WITH DISPERSION
functional_richness_cover_m4 <- glmmTMB(
  cover_trans_no0s ~ log_fun * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

##### want to test which model is best #####
AIC(functional_richness_cover_m1, functional_richness_cover_m2, functional_richness_cover_m3, functional_richness_cover_m4)
#better with dispformula and logged functional richness 

#df       AIC
#functional_richness_cover_m1 12 -23376.21
#functional_richness_cover_m2 15 -24693.52
#functional_richness_cover_m3 12 -23443.74
#functional_richness_cover_m4 15 -24786.20

##### checks for best model: functional_richness_cover_m4 ##### 
performance::check_singularity(functional_richness_cover_m4) 
#FALSE

performance::check_model(functional_richness_cover_m4)
### this takes awhile

functional_richness_cover_m4$fit$convergence           # should be 0 - yes
functional_richness_cover_m4$sdr$pdHess #TRUE

diagnose(functional_richness_cover_m4) #Unusually large Z-statistics, but ok 
VarCorr(functional_richness_cover_m4) # all signficantly above 0 
#Conditional model:
# Groups                Name        Std.Dev.
#location:site_habitat (Intercept) 0.16251 
#site_habitat          (Intercept) 0.25236 
#year                  (Intercept) 0.13403 

testUniformity(functional_richness_cover_m4) #looks ok-ish
res_taxo_log_fun <- simulateResiduals(functional_richness_cover_m4, n = 1000)
plotResiduals(res_taxo_log_fun, form = alpha_diversity_quad_macro_present$functional_richness, rank = FALSE)
plotResiduals(res_taxo_log_fun, form = alpha_diversity_quad_macro_present$functional_richness,
              rank = FALSE, asFactor = FALSE)
#throws warning but just bc we only have 4 functional richness 
plotResiduals(res_taxo_log_fun, form = alpha_diversity_quad_macro_present$habitat) #ok?

testDispersion(res_taxo_log_fun)  #slightly high but ok?  #make sure to run this on the residuals


hist(residuals(functional_richness_cover_m4))
hist(residuals(functional_richness_cover_m4), breaks = 25)

testOutliers(model_occ_log_fun, type = "bootstrap")


##### Anova and PostHoc for Best Model #####
car::Anova(functional_richness_cover_m4, type = "II")
#Chisq Df Pr(>Chisq)    
#log_fun         3231.6259  1  < 2.2e-16 ***
#habitat            8.5129  3    0.03652 *  
#log_fun:habitat   25.9497  3  9.771e-06 ***


#make sure interaction between log functional richness and habitat is important in our model
anova(update(functional_richness_cover_m4, . ~ . - log_fun:habitat), functional_richness_cover_m4)
#yes!

pairs(emtrends(functional_richness_cover_m4, ~ habitat, var = "log_fun"))
# contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef          -0.2730 0.0660 Inf  -4.138  0.0002
#Fringing - Forereef 10m      -0.1624 0.0607 Inf  -2.673  0.0377
#Fringing - Forereef 17m      -0.2584 0.0555 Inf  -4.656 <0.0001
#Backreef - Forereef 10m       0.1107 0.0543 Inf   2.038  0.1742
#Backreef - Forereef 17m       0.0146 0.0484 Inf   0.302  0.9904
#Forereef 10m - Forereef 17m  -0.0960 0.0408 Inf  -2.354  0.0862


em_cover_mod_fg <- emtrends(functional_richness_cover_m4, pairwise ~ habitat, var = "log_fun") 

pairs(em_cover_mod_fg)
#yep, same as other way

#need this for plotting Figure S2
cld_cover_mod_fg <- multcomp::cld(em_cover_mod_fg$emtrends, Letters = letters, sort = FALSE)


#### Data wrangling from models to make supplemental figures: ####

#richness by cover

t.s2.a <- as_tibble(cld_cover_mod_richness) %>%
  mutate(Predictor = "Taxonomic richness",
         `Spatial scale` = "Plot-level") %>% 
  # rename this so it matches for joining purposes
  rename_if(str_detect(names(.), ".trend"), ~"Mean")
#functional richness by cover

t.s2.b <-  as_tibble(cld_cover_mod_fg) %>%
  mutate(Predictor = "Functional richness",
         `Spatial scale` = "Plot-level") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

summary(cld_cover_mod_fg, infer = TRUE)

#joining together

rbind(t.s2.a, t.s2.b) %>% 
  mutate(Table = "S2_plotlevel") %>% 
  dplyr::select(-SE, -df) -> table_s2_plot

supplemental_tables_tableS2_plot <- table_s2_plot %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
         Lower_CI = asymp.LCL,
         Upper_CI = asymp.UCL)



