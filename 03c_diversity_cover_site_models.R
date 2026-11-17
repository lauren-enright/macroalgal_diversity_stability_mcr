#### Diversity Cover Site Models ####

source("00_functions_and_aes.R")

# read in data

alpha_diversity_site_macro <- read.csv(here::here("data", "alpha_diversity_site_macro_09262025.csv"))

alpha_diversity_site_macro_above0 <- alpha_diversity_site_macro %>%
  filter(richness > 0) %>%
  mutate(cover_trans_no0s = sv_trans(prop_cover), #slightly different numbers than if transforming, because lower sample size 
  log_rich_site = log(richness),
  log_fun_site  = log(functional_richness))

alpha_diversity_site_macro_above0$habitat <- factor(alpha_diversity_site_macro_above0$habitat,
                                                     levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

unique(alpha_diversity_site_macro_above0$site_habitat)  #24, this is our site-level random effect

cor(alpha_diversity_site_macro_above0$richness, alpha_diversity_site_macro_above0$functional_richness) #0.73 


min(alpha_diversity_site_macro_above0$richness) #1
max(alpha_diversity_site_macro_above0$richness) #14
str(alpha_diversity_site_macro_above0)

alpha_diversity_site_macro$habitat <- as.factor(alpha_diversity_site_macro$habitat)
alpha_diversity_site_macro$site <- as.factor(alpha_diversity_site_macro$site)
alpha_diversity_site_macro$site_habitat <- as.factor(alpha_diversity_site_macro$site_habitat)

#### Taxonomic Diversity ####

#####   M1: interaction, no log transformation of richness, no dispformula = habitat #####  

site_cover_mod <- glmmTMB(cover_trans_no0s ~ richness*habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)
summary(site_cover_mod)

res_taxo_lin_site <- simulateResiduals(site_cover_mod, n = 1000)
plotResiduals(res_taxo_lin_site, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
#humped shaped 
#this model sometimes throws a warning
#In (function (start, objective, gradient = NULL, hessian = NULL,  : NA/NaN function evaluation
#if choosing this model, would need to carefuly check that final fit is ok 
#If pdHess is TRUE, the standard errors look reasonable and diagnose() reports no problems, can generally treat the fit as fine.

#####  M2: interaction, no log transformation of richness, yes dispformula = habitat #####  

site_cover_mod2 <- glmmTMB(cover_trans_no0s ~ richness*habitat + (1|site_habitat) + (1|year), 
                           dispformula = ~ habitat,
                           family = beta_family(), data = alpha_diversity_site_macro_above0)

res_taxo_lin_site_m2 <- simulateResiduals(site_cover_mod2, n = 1000)
plotResiduals(res_taxo_lin_site_m2, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
#also humped shaped 

##### M3: interaction, log transformation of richness, no dispformula = habitat #####  

site_cover_mod3 <- glmmTMB(
  cover_trans ~ log_rich_site* habitat + (1 | site_habitat) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

res_taxo_log_site_m3 <- simulateResiduals(site_cover_mod3, n = 1000)
plotResiduals(res_taxo_log_site_m3, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)

#better - less humped shaped! 


######   M4: interaction, log transformation of richness, yes dispformula = habitat #####  

site_cover_mod4 <- glmmTMB(
  cover_trans ~ log_rich_site * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)
#no warning

res_taxo_log_site_m4 <- simulateResiduals(site_cover_mod4, n = 1000)
plotResiduals(res_taxo_log_site_m4, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)

##### want to test which model is best #####
AIC(site_cover_mod, site_cover_mod2, site_cover_mod3, site_cover_mod4)
#                df       AIC
#site_cover_mod  11 -1151.550
#site_cover_mod2 14 -1182.796
#site_cover_mod3 11 -1157.059
#site_cover_mod4 14 -1189.225

# model 4 is the best!

##### checks for best model: functional_richness_cover_m4 ##### 
performance::check_singularity(site_cover_mod4) 
#FALSE

performance::check_model(site_cover_mod4)

site_cover_mod4$fit$convergence           # should be 0 - yes
site_cover_mod4$sdr$pdHess #TRUE

diagnose(site_cover_mod4) #Unusually large Z-statistics, but ok 
VarCorr(site_cover_mod4) # all signficantly above 0 
#Conditional model:
#Groups       Name        Std.Dev.
#site_habitat (Intercept) 0.52502 
#year         (Intercept) 0.14667 

testUniformity(res_taxo_log_site_m4) #feel good about this -- run on residuals 

#glance at residuals again
plotResiduals(res_taxo_log_site_m4, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
plotResiduals(res_taxo_log_site_m4, form = alpha_diversity_site_macro_above0$habitat) #ok 


testDispersion(res_taxo_log_site_m4)  #looks good!

hist(residuals(site_cover_mod4))
hist(residuals(site_cover_mod4), breaks = 25)

testOutliers(site_cover_mod4, type = "bootstrap") #no outlier problems!
testOutliers(site_cover_mod4) ##no outlier problems!

##### Anova and PostHoc for Best Model #####

car::Anova(site_cover_mod4, type = "II")
#                        Chisq Df Pr(>Chisq)    
#log(richness)         58.2655  1   2.29e-14 ***
#habitat                5.3609  3   0.147197    
#log(richness):habitat 12.0441  3   0.007233 ** 

em_cover_mod_disp <- emtrends(site_cover_mod4, pairwise ~ habitat, var = "log_rich_site") 
cld_cover_mod_disp <- multcomp::cld(em_cover_mod_disp$emtrends, Letters = letters, sort = FALSE)

#### Functional Richness ####

#####   M1: interaction, no log transformation of richness, no dispformula = habitat #####  

site_cover_mod_fg <- glmmTMB(cover_trans_no0s ~ functional_richness*habitat + (1|site_habitat) + (1|year), family = beta_family(), 
                             data = alpha_diversity_site_macro_above0)
summary(site_cover_mod_fg)

res_fun_lin_site <- simulateResiduals(site_cover_mod_fg, n = 1000)
plotResiduals(res_fun_lin_site, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)
#hump shaped  

#####  M2: interaction, no log transformation of richness, yes dispformula = habitat #####  

site_cover_mod_fg2 <- glmmTMB(cover_trans_no0s ~ functional_richness*habitat + (1|site_habitat) + (1|year), 
                           dispformula = ~ habitat,
                           family = beta_family(), data = alpha_diversity_site_macro_above0)

res_fun_lin_site_m2 <- simulateResiduals(site_cover_mod_fg2, n = 1000)
plotResiduals(res_fun_lin_site_m2, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)
##hump shaped  

##### M3: interaction, log transformation of richness, no dispformula = habitat #####  

site_cover_mod_fg3 <- glmmTMB(
  cover_trans ~ log_fun_site* habitat + (1 | site_habitat) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

res_fun_log_site_m3 <- simulateResiduals(site_cover_mod_fg3, n = 1000)
plotResiduals(res_fun_log_site_m3, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)

#looks slightly better than none-log transformed functional richness! 


#####   M4: interaction, log transformation of richness, yes dispformula = habitat #####  

site_cover_mod_fg4 <- glmmTMB(
  cover_trans ~ log_fun_site * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

res_fun_log_site_m4 <- simulateResiduals(site_cover_mod_fg4, n = 1000)
plotResiduals(res_fun_log_site_m4, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)
plotResiduals(res_fun_log_site_m4, form = alpha_diversity_site_macro_above0$functional_richness,
              rank = FALSE, asFactor = FALSE)

car::Anova(site_cover_mod_fg4, type = "II")
#interaction is not signficant....

##### want to test which model is best #####
AIC(site_cover_mod_fg, site_cover_mod_fg2, site_cover_mod_fg3, site_cover_mod_fg4)

#df       AIC
#site_cover_mod_fg  11 -1121.788
#site_cover_mod_fg2 14 -1154.881
#site_cover_mod_fg3 11 -1123.081
#site_cover_mod_fg4 14 -1155.479

# model 4 is the best! but we want to remove the interaction since it is not signficant

######   M5: NO interaction, log transformation of richness, yes dispformula = habitat #####  

site_cover_mod_fg5 <- glmmTMB(
  cover_trans ~ log_fun_site + habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

res_fun_log_site_m5 <- simulateResiduals(site_cover_mod_fg5, n = 1000)
plotResiduals(res_fun_log_site_m5, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)
plotResiduals(res_fun_log_site_m5, form = alpha_diversity_site_macro_above0$functional_richness,
              rank = FALSE, asFactor = FALSE)

#looks ok, such little data at 6 and 7, not suprising it is predicting poorly... 


##### checks for best model: functional_richness_cover_m5 ##### 
performance::check_singularity(site_cover_mod_fg5) 
#FALSE

performance::check_model(site_cover_mod_fg5)

site_cover_mod_fg5$fit$convergence           # should be 0 - yes
site_cover_mod_fg5$sdr$pdHess #TRUE

diagnose(site_cover_mod_fg5) #Unusually large Z-statistics, but ok
VarCorr(site_cover_mod_fg5) # all signficantly above 0 
#Conditional model:
# Groups       Name        Std.Dev.
#site_habitat (Intercept) 0.53813 
#year         (Intercept) 0.14347 

testUniformity(res_fun_log_site_m5) #feel good about this -- run on residuals 

#glance at residuals again
plotResiduals(res_fun_log_site_m5, form = alpha_diversity_site_macro_above0$functional_richness, rank = FALSE)
plotResiduals(res_fun_log_site_m5, form = alpha_diversity_site_macro_above0$habitat) #ok 


testDispersion(res_fun_log_site_m5)  #looks good!

hist(residuals(site_cover_mod_fg5))
hist(residuals(site_cover_mod_fg5), breaks = 25)

testOutliers(site_cover_mod_fg5, type = "bootstrap") #no outlier problems!
testOutliers(site_cover_mod_fg5) ##no outlier problems!


#### Anova and PostHoc for Best Model ####

car::Anova(site_cover_mod_fg5, type = "II")
#Chisq Df Pr(>Chisq)    
#log_fun_site 26.8108  1  2.244e-07 ***
#habitat       3.9536  3     0.2665  

em_cover_mod_disp_fg <- emtrends(site_cover_mod_fg5, pairwise ~ habitat, var = "log_fun_site") 
cld_cover_mod_disp_fg <- multcomp::cld(em_cover_mod_disp_fg$emtrends, Letters = letters, sort = FALSE)
#yep -- all the same letter (makes sense since habitat is not significant)

