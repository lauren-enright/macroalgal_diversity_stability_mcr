#### Diversity Cover Site Models ####

source("00_functions_and_aes.R")

# read in data

alpha_diversity_site_macro <- read.csv(here::here("data", "alpha_diversity_site_macro_09262025.csv"))

cor.test(alpha_diversity_site_macro$richness, alpha_diversity_site_macro$functional_richness)
plot(alpha_diversity_site_macro$richness, alpha_diversity_site_macro$functional_richness)

table(alpha_diversity_site_macro$richness, alpha_diversity_site_macro$cover_trans)

min(alpha_diversity_site_macro$richness)

str(alpha_diversity_site_macro)

alpha_diversity_site_macro$habitat <- as.factor(alpha_diversity_site_macro$habitat)
alpha_diversity_site_macro$site <- as.factor(alpha_diversity_site_macro$site)
alpha_diversity_site_macro$site_habitat <- as.factor(alpha_diversity_site_macro$site_habitat)

#### Taxonomic Diversity ####

site_cover_mod <- glmmTMB(cover_trans ~ richness*habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)
summary(site_cover_mod)
car::Anova(site_cover_mod)

#taxo richness
#                   Chisq Df Pr(>Chisq)    
#richness         58.9195  1  1.643e-14 ***
#  habitat           6.6247  3    0.08487 .  
#richness:habitat  3.2967  3    0.34811  


performance::check_singularity(site_cover_mod)
#site/site_habitat = singular
#site_habitat = NOT singular

#this matches
hist(residuals(site_cover_mod)) # some positive skew, not the worst thing ever
plot(residuals(site_cover_mod) ~ fitted(site_cover_mod)) # a little heteroscedastic but can't be less than 0
#r.squaredGLMM(site_cover_mod) # 0.2212534, 0.5990825
performance::r2(site_cover_mod) # marginal: 0.22, conditional: 0.59; matches

# do not need to do a posthoc test here because no interaction and habitat was also not signficant. 

#running again without the interaction, this is what we will report 
site_cover_mod2 <- glmmTMB(cover_trans ~ richness + habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)
summary(site_cover_mod2)
car::Anova(site_cover_mod2)
#           Chisq Df Pr(>Chisq)     #same
#richness 58.0746  1  2.524e-14 ***
# habitat   6.7036  3    0.08197 .  

site_cover_mod2$fit$convergence           # should be 0 - yes
site_cover_mod2$sdr$pdHess    #TRUE

diagnose(site_cover_mod2) #Unusually large Z-statistics, but ok bc SD is above 0
VarCorr(site_cover_mod2) # all signficantly above 0 

hist(residuals(site_cover_mod2)) # some positive skew, not the worst thing ever
plot(residuals(site_cover_mod2) ~ fitted(site_cover_mod2)) # a little heteroscedastic but can't be less than 0
performance::r2(site_cover_mod2) # marginal: 0.21, conditional: 0.58

#trying log 

alpha_diversity_site_macro_above0 <- alpha_diversity_site_macro %>%
  filter(richness > 0)

#lose 2 sites
# year   site  habitat    site_habitat shannon richness prop_cover  cover_trans functional_richness
#1 2018 lter_5 Fringing lter_5_fringing       0        0          0 1.302083e-08                   0
#2 2023 lter_5 Fringing lter_5_fringing       0        0          0 1.302083e-08                   0

str(alpha_diversity_site_macro_above0)


site_cover_mod3 <- glmmTMB(
  cover_trans ~ log(richness) * habitat + (1 | site_habitat) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

car::Anova(site_cover_mod3)

site_cover_mod4 <- glmmTMB(
  cover_trans ~ log(richness) * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

car::Anova(site_cover_mod4)

car:

em_cover_mod_disp <- emtrends(site_cover_mod4, pairwise ~ habitat, var = "richness") 
cld_cover_mod_disp <- multcomp::cld(em_cover_mod_disp$emtrends, Letters = letters, sort = FALSE)

pred_test <- ggpredict(site_cover_mod4, terms = c("richness [all]", "habitat"))
plot(pred_test)

emmip(site_cover_mod4, habitat ~ richness,
      at = list(richness = 1:10), type = "response", CIs = TRUE)

alpha_diversity_site_macro_above0$log_richness <- log(alpha_diversity_site_macro_above0$richness)
alpha_diversity_site_macro_above0$log_fun_richness <- log(alpha_diversity_site_macro_above0$functional_richness)

site_cover_mod5 <- glmmTMB(
  cover_trans ~ log_richness * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)

site_cover_mod6 <- glmmTMB(
  cover_trans ~ log_richness * habitat + (1 | site_habitat) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_site_macro_above0
)




car::Anova(site_cover_mod5, type = "II")

AIC(site_cover_mod3, site_cover_mod4, site_cover_mod5)

performance::r2(site_cover_mod4)

testDispersion(res_taxo_log_site)

site_cover_mod5_lin <- update(site_cover_mod5, . ~ . - log(richness) + richness)
AIC(site_cover_mod5, site_cover_mod5_lin)

plotResiduals(simulateResiduals(site_cover_mod5),     form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
plotResiduals(simulateResiduals(site_cover_mod3), form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
plotResiduals(simulateResiduals(site_cover_mod5_lin), form = alpha_diversity_site_macro_above0$richness, rank = FALSE)

site_cover_log <- glmmTMB(
  cover_trans ~ log(richness) * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family = beta_family(link = "logit"),
  data = alpha_diversity_site_macro_above0
)

car::Anova(site_cover_log)



site_cover_lin <- glmmTMB(
  cover_trans ~ richness * habitat + (1 | site_habitat) + (1 | year),
  dispformula = ~ habitat,
  family = beta_family(link = "logit"),
  data = alpha_diversity_site_macro_above0
)

AIC(site_cover_log, site_cover_lin)

anova(update(site_cover_log, . ~ . - log(richness):habitat), site_cover_log)

diagnose(site_cover_mod3) #Unusually large Z-statistics, but ok bc SD is above 0
VarCorr(site_cover_mod3) # all signficantly above 0 
# Groups       Name        Std.Dev.
#site_habitat (Intercept) 0.50202 
#year         (Intercept) 0.16699 
testUniformity(site_cover_mod2)

res_taxo_log_site <- simulateResiduals(site_cover_mod3, n = 1000)
plotResiduals(res_taxo_log_site, form = alpha_diversity_site_macro_above0$richness, rank = FALSE)
plotResiduals(res_taxo_log_site, form = alpha_diversity_site_macro_above0$richness,
              rank = FALSE, asFactor = FALSE)
#throws warning but just bc we only have 6 richness 
testDispersion(res_taxo_log_site)  #slightly high but ok?  #make sure to run this on the residuals
performance::check_model(res_taxo_log_site)

hist(residuals(model_occ_log))
hist(residuals(res_taxo_log), breaks = 25)


car::Anova(site_cover_mod3)

AIC(site_cover_mod3, site_cover_mod2, site_cover_mod)

anova(site_cover_mod2, site_cover_mod3)

pairs(emtrends(site_cover_mod5, ~ habitat, var = "log_richness"))

### running checks on oct 21 -- v good flow###

d <- alpha_diversity_site_macro_above0

models <- list(
  log_disp   = site_cover_log,                             # log, habitat dispersion
  lin_disp   = site_cover_lin,                             # linear, habitat dispersion
  log_nodisp = update(site_cover_log, dispformula = ~ 1)   # log, constant dispersion
)

# 1. Model fit (lower is better; same data and response required)
sapply(models, nobs)       # must all match
sapply(models, AIC)

# 2. Residuals for each model
res <- lapply(models, simulateResiduals, n = 1000)

par(mfrow = c(1, 2))
for (nm in names(res)) { plot(res[[nm]], main = nm) }

plotResiduals(res$log_disp,   form = d$habitat)
plotResiduals(res$log_nodisp, form = d$habitat)

plotResiduals(res$log_disp, form = d$richness, rank = FALSE)
plotResiduals(res$lin_disp, form = d$richness, rank = FALSE)

sapply(res, function(r) c(
  KS_p         = testUniformity(r, plot = FALSE)$p.value,
  dispersion   = unname(testDispersion(r, plot = FALSE)$statistic),
  dispersion_p = testDispersion(r, plot = FALSE)$p.value
))



#### Functional Diversity ####

site_cover_mod_fg <- glmmTMB(cover_trans ~ functional_richness*habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)

summary(site_cover_mod_fg)
car::Anova(site_cover_mod_fg)

#with 8 functional groups
#Response: cover_trans
#Chisq Df Pr(>Chisq)    
#functional_richness         24.1673  1  8.832e-07 ***
#  habitat                      4.1544  3    0.24527    
#functional_richness:habitat  6.5034  3    0.08953 .  

performance::check_singularity(site_cover_mod_fg)
#site/site_habitat = singular
#site_habitat = OK

hist(residuals(site_cover_mod_fg)) # positive skew
plot(residuals(site_cover_mod_fg) ~ fitted(site_cover_mod_fg)) # a little heteroscedastic but can't be less than 0
performance::r2(site_cover_mod_fg) 
#new with 8 funtional groups:Conditional R2: 0.566, Marginal R2: 0.169

#don't need a post hoc test here because no interaction and habitat was also not signficant. 

#running again without the interaction, this is what we will report 
site_cover_mod_fg2 <- glmmTMB(cover_trans ~ functional_richness + habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)

summary(site_cover_mod_fg2)
car::Anova(site_cover_mod_fg2)

#new -- 8 functional groups # okay cool 
#Chisq Df Pr(>Chisq)    
#functional_richness 23.1372  1  1.508e-06 ***
# habitat              3.8533  3     0.2777  

hist(residuals(site_cover_mod_fg2)) # positive skew
plot(residuals(site_cover_mod_fg2) ~ fitted(site_cover_mod_fg2)) # a little heteroscedastic but can't be less than 0
performance::r2(site_cover_mod_fg2) 
#new - 8 fg =   Conditional R2: 0.563; Marginal R2: 0.145

site_cover_mod_fg2 <- glmmTMB(cover_trans ~ log(functional_richness)*habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)
car::Anova(site_cover_mod_fg2)


AIC(site_cover_mod_fg2,site_cover_mod_fg)

site_cover_mod_fg3 <- glmmTMB(cover_trans ~ log(functional_richness)*habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)

site_cover_mod_fg4 <- glmmTMB(cover_trans ~ log(functional_richness)+habitat + (1|site_habitat) + (1|year), family = beta_family(), data = alpha_diversity_site_macro_above0)

AIC(site_cover_mod_fg4,site_cover_mod_fg2)

anova(site_cover_mod_fg2, site_cover_mod_fg4)

site_cover_mod_fg5 <- glmmTMB(cover_trans ~ log(functional_richness)*habitat + (1|site_habitat) + (1|year), 
                              dispformula = ~ habitat,
                              family = beta_family(), data = alpha_diversity_site_macro_above0)

car::Anova(site_cover_mod_fg5)

site_cover_mod_fg6 <- glmmTMB(cover_trans ~ log_fun_richness + habitat + (1|site_habitat) + (1|year), 
                              dispformula = ~ habitat,
                              family = beta_family(), data = alpha_diversity_site_macro_above0)

site_cover_mod_fg7 <- glmmTMB(cover_trans ~ log_fun_richness + habitat + (1|site_habitat) + (1|year), 
                              family = beta_family(), data = alpha_diversity_site_macro_above0)

site_cover_mod_fg6 <- glmmTMB(cover_trans ~ log() + habitat + (1|site_habitat) + (1|year), 
                              dispformula = ~ habitat,
                              family = beta_family(), data = alpha_diversity_site_macro_above0)
car::Anova(site_cover_mod_fg6)

site_cover_mod_fg7 <- glmmTMB(cover_trans ~ functional_richness + habitat + (1|site_habitat) + (1|year), 
                            #  dispformula = ~ habitat,
                              family = beta_family(), data = alpha_diversity_site_macro_above0)

car::Anova(site_cover_mod_fg7)

plotResiduals(simulateResiduals(site_cover_mod_fg2), form = alpha_diversity_site_macro_above0$habitat)   # without
plotResiduals(simulateResiduals(site_cover_mod_fg3), form = alpha_diversity_site_macro_above0$habitat)
plotResiduals(simulateResiduals(site_cover_mod_fg4), form = alpha_diversity_site_macro_above0$habitat)   # without
plotResiduals(simulateResiduals(site_cover_mod_fg5), form = alpha_diversity_site_macro_above0$habitat)   # with

perform
AIC(site_cover_mod_fg4,site_cover_mod_fg2, site_cover_mod_fg5, site_cover_mod_fg6)

car::Anova(site_cover_mod_fg5, type = "II")

emmip(site_cover_mod_fg5, habitat ~ functional_richness,
      at = list(functional_richness = 1:10), type = "response", CIs = TRUE)

cover_trans ~ log(richness) * habitat + (1 | site_habitat) + (1 | year),
dispformula = ~ habitat,
family      = beta_family(link = "logit"),
data        = alpha_diversity_site_macro_above0

### running checks on oct 21 -- v good flow###

d <- alpha_diversity_site_macro_above0

models_func <- list(
  log_disp   = site_cover_mod_fg6,                             # log, habitat dispersion
  lin_disp   = site_cover_mod_fg7,                             # linear, habitat dispersion
  log_nodisp = update(site_cover_mod_fg6, dispformula = ~ 1)   # log, constant dispersion
)

# 1. Model fit (lower is better; same data and response required)
sapply(models_func, nobs)       # must all match
sapply(models_func, AIC)

# 2. Residuals for each model
res <- lapply(models_func, simulateResiduals, n = 1000)

par(mfrow = c(1, 2))
for (nm in names(res)) { plot(res[[nm]], main = nm) }

plotResiduals(res$log_disp,   form = d$habitat)
plotResiduals(res$log_nodisp, form = d$habitat)

plotResiduals(res$log_disp, form = d$functional_richness, rank = FALSE)
plotResiduals(res$lin_disp, form = d$functional_richness, rank = FALSE)

sapply(res, function(r) c(
  KS_p         = testUniformity(r, plot = FALSE)$p.value,
  dispersion   = unname(testDispersion(r, plot = FALSE)$statistic),
  dispersion_p = testDispersion(r, plot = FALSE)$p.value
))


em_cover_mod_Site <- emtrends(site_cover_mod6, pairwise ~ habitat, var = "log_richness") 

em_cover_mod_Site #yep, same as above

cld_cover_mod_richness_site <- multcomp::cld(em_cover_mod_Site$emtrends, Letters = letters, sort = FALSE)


t.s2.a_site <- as_tibble(cld_cover_mod_richness_site) %>%
  mutate(Predictor = "Taxonomic richness",
         `Spatial scale` = "Site-level") %>% 
  # rename this so it matches for joining purposes
  rename_if(str_detect(names(.), ".trend"), ~"Mean")
