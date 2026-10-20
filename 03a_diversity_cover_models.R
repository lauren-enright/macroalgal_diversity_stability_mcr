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

unique(alpha_diversity_quad_macro$habitat)

alpha_diversity_quad_macro %>%
#filter(habitat == "Fringing")
  filter(location == "lter_4_backreef_algae_transect_5_quad_10")
 # filter(location == "lter_1_fringing_reef_algae_transect_3_quad_10")
  #filter(location == "lter_2_outer_17_m_algae_transect_5_quad_10")
  filter(location == "LTER 2 Outer 17 m Algae Transect 5 Quad 10")
  
alpha_diversity_quad_macro %>%
  filter(site_habitat == "lter_5_fringing") %>%
  filter(year == 2018) %>%
  filter(macroalgal_prop_cover > 0)


unique(alpha_diversity_quad_macro$cover_trans)

alpha_diversity_quad_macro$habitat <- factor(alpha_diversity_quad_macro$habitat,
                     levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

cor(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$functional_richness, use = "complete.obs", method = "pearson") #0.9584012
cor.test(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$functional_richness, method = "pearson")   # or "spearman"
#0.9584012 
cor.test(alpha_diversity_quad_macro_present$richness, alpha_diversity_quad_macro_present$functional_richness, method = "pearson")   # or "spearman"
#0.8869669 

alpha_diversity_quad_macro_present


## Need to separate out non-0 data first

alpha_diversity_quad_macro_present <- alpha_diversity_quad_macro %>%
  dplyr::filter(macroalgal_prop_cover > 0) %>%
  mutate(cover_trans_no0s = sv_trans(macroalgal_prop_cover), #slightly different numbers than if transforming, because lower sample size 
        log_rich = log(richness),
        log_fun  = log(functional_richness))

#also want to turn algae in present absence for step 1 of hurdle mdoel 
alpha_diversity_quad_macro$macro_present <- ifelse(alpha_diversity_quad_macro$macroalgal_prop_cover > 0, 1, 0)

#just making sure the code is doing what it should 
table(alpha_diversity_quad_macro$habitat, alpha_diversity_quad_macro$macro_present)

#0    1
#Fringing     2575 2214
#Backreef     2307 2486
#Forereef 10m 1186 3613
#Forereef 17m  917 3875

#perfect this matches exactly, as expected

#1 Fringing. 2575
#2 Backreef 2307
#3 Forereef 10m 1186
#4 Forereef 17m 917


#### Hurdle model ####

# Part 1: Does the probability of macroalgae occurring differ by habitat?
#this uses ALL the data

#random effects of site_habitat aka site ID within location aka plot and then also year.. 

presence_model <- glmmTMB(macro_present ~ habitat + (1|site_habitat/location)  + (1 | year),
                  family = binomial, data = alpha_diversity_quad_macro)

performance::check_model(presence_model)
#looks ok overall
summary(presence_model)

#running specific model checks
diagnose(presence_model)       # flags huge estimates or near-zero variances -- looks ok
VarCorr(presence_model)        # all random-effect SDs clearly above 0 -- looks ok 
#Conditional model:
#Groups                Name        Std.Dev.
#location:site_habitat (Intercept) 0.68913 
#site_habitat          (Intercept) 1.07002 
#year                  (Intercept) 0.38117 
table(alpha_diversity_quad_macro$habitat, alpha_diversity_quad_macro$macro_present)
#no habitat all 0s or all 1s

res_presence <- simulateResiduals(presence_model, n = 1000)
plot(res_presence)
#this looks kinda amazing? 

plotResiduals(res_presence, form = alpha_diversity_quad_macro$habitat)
#also looks amazing

testDispersion(res_presence)
#data:  simulationOutput
#dispersion = 1.009, p-value = 0.866
#alternative hypothesis: two.sided

#residuals by groups
res_loc <- recalculateResiduals(res_presence, group = interaction(alpha_diversity_quad_macro$site_habitat,
                                                           alpha_diversity_quad_macro$location))
testDispersion(res_loc)
#data:  simulationOutput
#dispersion = 1.0249, p-value = 0.8
#alternative hypothesis: two.sided
#feel good about this as well 

#The quadrat totals vary just as much as the model predicts, so the location random effect captures the clustering of presence within quadrats
#aka good 2 go

##also trying LRTs 
#make null model with no habitat
m_pres_null <- update(presence_model, . ~ . - habitat)
#now compare my model with null model
anova(m_pres_null, presence_model)
#m_pres_null     4 19660 19692 -9826.2    19652                            
#presence_model  7 19654 19709 -9819.8    19640 12.697      3   0.005341 **
#this tell us that habitat is important :) 

## getting stats to report ##

car::Anova(presence_model, type = "II")
#Response: macro_present
#Chisq Df Pr(>Chisq)    
#habitat 16.653  3   0.000833 ***

performance::r2(presence_model) 
#  Conditional R2: 0.439
#Marginal R2: 0.138

#for presence model, now calculate post-hoc
pres_emm <- emmeans(presence_model, ~ habitat, type = "response")
pres_emm                                   # probability of presence per habitat
pairs(pres_emm)                            # pairwise odds ratios, Tukey-adjusted
# contrast                    odds.ratio     SE  df null z.ratio p.value
#Fringing / Backreef              0.752 0.4680 Inf    1  -0.459  0.9680
#Fringing / Forereef 10m          0.203 0.1260 Inf    1  -2.560  0.0513
#Fringing / Forereef 17m          0.114 0.0710 Inf    1  -3.484  0.0028
#Backreef / Forereef 10m          0.270 0.1680 Inf    1  -2.102  0.1525
#Backreef / Forereef 17m          0.151 0.0944 Inf    1  -3.027  0.0132
#Forereef 10m / Forereef 17m      0.560 0.3500 Inf    1  -0.929  0.7893

pairs(emmeans(presence_model, ~ habitat))  # same comparisons on the log-odds scale -- these look exactly the same

# Compact letter display (habitats sharing a letter don't differ)
library(multcomp)
cld(emmeans(presence_model, ~ habitat), Letters = letters)

em_cover_pres_emm <- emmeans(presence_model, pairwise ~ habitat, type = "response") 

em_cover_pres_emm
#yep, same as other way

cld_presence_model <- multcomp::cld(em_cover_pres_emm$emmeans, Letters = letters, sort = FALSE)


#### Taxonomic richness and cover ####

#okay now let's do the "occupancy" model for taxonomic richness 

#first trying linear
model_occ_lin_taxo <- glmmTMB(
  cover_trans_no0s ~ richness * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(model_occ_lin_taxo)

res_taxo_lin <- simulateResiduals(model_occ_lin_taxo, n = 1000)
plotResiduals(res_taxo_lin, form = alpha_diversity_quad_macro_present$richness, rank = FALSE)


#looks log? overiftting at low values and underfiting and high values

model_occ_log_taxo <- glmmTMB(
  cover_trans_no0s ~ log_rich * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(model_occ_log_taxo)

# checks for this model:
diagnose(model_occ_log_taxo) #Unusually large Z-statistics, but ok bc SD is above 0
VarCorr(model_occ_log_taxo) # all signficantly above 0 
#Groups                Name        Std.Dev.
#location:site_habitat (Intercept) 0.15807 
#site_habitat          (Intercept) 0.23880 
#year                  (Intercept) 0.13890 
testUniformity(model_occ_log_taxo)
res_taxo_log <- simulateResiduals(model_occ_log, n = 1000)
plot(res_taxo_log)
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$richness, rank = FALSE)
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$richness,
              rank = FALSE, asFactor = FALSE)
#throws warning but just bc we only have 6 richness 
testDispersion(res_taxo_log)  #slightly high but ok?  #make sure to run this on the residuals
performance::check_model(res_taxo_log)

hist(residuals(model_occ_log))
hist(residuals(res_taxo_log), breaks = 25)

#want to test which model is best
AIC(model_occ_lin, model_occ_log)
#              df       AIC
#model_occ_lin 15 -25421.75
#model_occ_log 15 -25575.60

testOutliers(res_taxo_log, type = "bootstrap")

testDispersion(recalculateResiduals(res_taxo_log, group = interaction(d$site_habitat, d$location)))

#making sure dispformula = ~ habitat improve the mdoel

model_occ_log_taxo_NODISP <- glmmTMB(
  cover_trans_no0s ~ log_rich * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

AIC(model_occ_lin, model_occ_log, model_occ_log_taxo_NODISP)


summary(model_occ_log_taxo_NODISP)
summary(model_occ_log_taxo)
testDispersion(res_taxo_log)
testDispersion(res_taxo_log_NODISP)

res_taxo_log_NODISP <- simulateResiduals(model_occ_log_taxo_NODISP, n = 1000)
plot(res_taxo_log_NODISP)

plotResiduals(res_taxo_log_NODISP, form = alpha_diversity_quad_macro_present$habitat)
plotResiduals(res_taxo_log, form = alpha_diversity_quad_macro_present$habitat)
plotResiduals(res, form = alpha_diversity_quad_macro_present$habitat)

### now do the Anova
car::Anova(model_occ_log_taxo, type = "II")
car::Anova(model_occ_log_taxo_NODISP, type = "II")
#Chisq Df Pr(>Chisq)    
#Chisq Df Pr(>Chisq)    
#log_rich         4236.3045  1    < 2e-16 ***
#habitat             5.8679  3    0.11822    
#log_rich:habitat   11.0501  3    0.01146 *  

#Run likelihood ratio test of whether the richness × habitat interaction improves the taxonomic richness model. 
anova(update(model_occ_log_taxo, . ~ . - log_rich:habitat), model_occ_log_taxo)
#suggests keeping interaction, but not overwhelming different

pairs(emtrends(model_occ_log_taxo, ~ habitat, var = "log_rich"))
pairs(emtrends(model_occ_log_taxo_NODISP, ~ habitat, var = "log_rich"))

#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef          -0.1780 0.0564 Inf  -3.157  0.0087
#Fringing - Forereef 10m      -0.1212 0.0535 Inf  -2.267  0.1059
#Fringing - Forereef 17m      -0.1454 0.0499 Inf  -2.917  0.0186
#Backreef - Forereef 10m       0.0568 0.0460 Inf   1.236  0.6041
#Backreef - Forereef 17m       0.0326 0.0418 Inf   0.781  0.8631
#Forereef 10m - Forereef 17m  -0.0242 0.0371 Inf  -0.652  0.9150

em_cover_mod <- emtrends(model_occ_log, pairwise ~ habitat, var = "log_rich") 

em_cover_mod #yep, same as above

cld_cover_mod_richness <- multcomp::cld(em_cover_mod$emtrends, Letters = letters, sort = FALSE)


#:trying a shape-free check. If the log model were the wrong shape, this more flexible model would fit noticeably better.
m_fac_int <- update(model_occ_log, . ~ . - log_rich * habitat + factor(richness) * habitat)
AIC(model_occ_log, m_fac_int)
#the log model was better by about 14 -- So the log curve captures the relationship adequately.
emmeans(m_fac_int, ~ richness | habitat, type = "response")

#### Functional Richness ####

#first trying linear
model_occ_lin_functional <- glmmTMB(
  cover_trans_no0s ~ functional_richness * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

summary(model_occ_lin_functional)

#now trying logistic  
model_occ_log_fun <- glmmTMB(
  cover_trans_no0s ~ log_fun * habitat + (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)

model_occ_log_NODISP <- glmmTMB(
  cover_trans_no0s ~ log_fun * habitat + (1 | site_habitat/location) + (1 | year),
  family      = beta_family(link = "logit"),
  data        = alpha_diversity_quad_macro_present
)
testDispersion(model_occ_log_NODISP)

AIC(model_occ_log_fun,model_occ_log_NODISP)
#better with dispformula 

# checks for this model:
diagnose(model_occ_log_fun) #Unusually large Z-statistics, but ok bc SD is above 0
VarCorr(model_occ_log_fun) # all signficantly above 0 
#Conditional model:
# Groups                Name        Std.Dev.
#location:site_habitat (Intercept) 0.16251 
#site_habitat          (Intercept) 0.25236 
#year                  (Intercept) 0.13403 

testUniformity(model_occ_log_fun) #looks ok-ish
res_taxo_log_fun <- simulateResiduals(model_occ_log_fun, n = 1000)
plotResiduals(res_taxo_log_fun, form = alpha_diversity_quad_macro_present$functional_richness, rank = FALSE)
plotResiduals(res_taxo_log_fun, form = alpha_diversity_quad_macro_present$functional_richness,
              rank = FALSE, asFactor = FALSE)
#throws warning but just bc we only have 4 functional richness 

testDispersion(res_taxo_log_fun)  #slightly high but ok?  #make sure to run this on the residuals


hist(residuals(model_occ_log_fun))
hist(residuals(res_taxo_log_fun), breaks = 25)

#want to test which model is best
AIC(model_occ_lin_functional, model_occ_log_fun)

#df       AIC
#model_occ_lin_functional 15 -24693.52
#model_occ_log_fun        15 -24786.20

testOutliers(model_occ_log_fun, type = "bootstrap")

testDispersion(recalculateResiduals(res_taxo_log_fun, group = interaction(d$site_habitat, d$location)))
#high disperison

res_taxo_lin_fun <- simulateResiduals(model_occ_lin_functional, n = 1000)
plotResiduals(res_taxo_lin_fun, form = alpha_diversity_quad_macro_present$functional_richness, rank = FALSE)

performance::check_singularity(model_occ_log_fun) 
#FALSE

performance::check_model(model_occ_log_fun)
###


diagnose(model_occ_log_fun2)   # large Z for theta_ (log SDs) is fine if SDs are clearly > 0
VarCorr(model_occ_log_fun) 
plotResiduals(simulateResiduals(model_occ_log_fun2), form = alpha_diversity_quad_macro_present$functional_richness, rank = FALSE)

AIC(model_occ_log_fun2, model_occ_log_fun)

### now do that Anova
car::Anova(model_occ_log_fun, type = "II")
#Chisq Df Pr(>Chisq)    
#log_fun         3231.6259  1  < 2.2e-16 ***
#habitat            8.5129  3    0.03652 *  
#log_fun:habitat   25.9497  3  9.771e-06 ***

performance::r2(model_occ_log_fun2) 

#make sure interaction between log functional richness and habitat is important in our model
anova(update(model_occ_log_fun, . ~ . - log_fun:habitat), model_occ_log_fun)
#yeppp

pairs(emtrends(model_occ_log_fun, ~ habitat, var = "log_fun"))
# contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef          -0.2730 0.0660 Inf  -4.138  0.0002
#Fringing - Forereef 10m      -0.1624 0.0607 Inf  -2.673  0.0377
#Fringing - Forereef 17m      -0.2584 0.0555 Inf  -4.656 <0.0001
#Backreef - Forereef 10m       0.1107 0.0543 Inf   2.038  0.1742
#Backreef - Forereef 17m       0.0146 0.0484 Inf   0.302  0.9904
#Forereef 10m - Forereef 17m  -0.0960 0.0408 Inf  -2.354  0.0862


em_cover_mod_fg <- emtrends(model_occ_log_fun, pairwise ~ habitat, var = "log_fun") 

pairs(em_cover_mod_fg)
#yep, same as other way

cld_cover_mod_fg <- multcomp::cld(em_cover_mod_fg$emtrends, Letters = letters, sort = FALSE)


#### trying the betabinomial model ###

d <- alpha_diversity_quad_macro_present %>%
mutate(
  habitat  = factor(habitat, levels = c("Fringing", "Backreef",
                                        "Forereef 10m", "Forereef 17m")),
  year     = factor(year),
  hits     = round(macroalgal_prop_cover * 25),
  extra    = hits - 1,          # hits beyond the first
  log_rich = log(richness),
  log_fun  = log(functional_richness)
)

# Check: cover must be whole multiples of 1/25 (should be TRUE)
all.equal(d$hits, d$macroalgal_prop_cover * 25)
table(d$hits)
range(d$hits)                     # should be 1 to 25

# ---------------------------------------------------------------
# 1. Taxonomic richness
# ---------------------------------------------------------------
bb_taxo <- glmmTMB(
  cbind(extra, 24 - extra) ~ log_rich * habitat +
    (1 | site_habitat/location) + (1 | year),
  dispformula = ~ habitat,
  family      = betabinomial(link = "logit"),
  data        = d
)

# site-specific year effects (compare by AIC, keep the better one)
bb_taxo_sy <- update(bb_taxo, . ~ . + (1 | site_habitat:year))
AIC(bb_taxo, bb_taxo_sy)

# Does habitat-varying dispersion help? (compare by AIC)
AIC(update(bb_taxo, dispformula = ~ 1), bb_taxo)

# --- Convergence ---
diagnose(bb_taxo)
VarCorr(bb_taxo)

# --- Residual checks (should now be much closer to uniform) ---
res_bb_taxo <- simulateResiduals(bb_taxo, n = 1000)
plot(res_bb_taxo)
hist(residuals(res_bb_taxo), breaks = 40)        # want roughly flat
plotResiduals(res_bb_taxo, form = d$richness, rank = FALSE)
plotResiduals(res_bb_taxo, form = d$habitat)
plotResiduals(res_bb_taxo, form = d$year)
testDispersion(res_bb_taxo)
testOutliers(res_bb_taxo, type = "bootstrap")
testDispersion(recalculateResiduals(res_bb_taxo,
                                    group = interaction(d$site_habitat, d$location)))
hist(residuals(model_occ_log_fun))
hist(residuals(res_bb_taxo), breaks = 25)

# --- Tests ---
car::Anova(bb_taxo, type = "II")   
#Chisq Df Pr(>Chisq)    
#log_rich         4055.4618  1     <2e-16 ***
#habitat             2.3163  3     0.5094    
#log_rich:habitat   14.5926  3     0.0022 ** 

# Wald chi-square
anova(update(bb_taxo, . ~ . - log_rich:habitat), bb_taxo)        # LRT for interaction

# --- Slopes by habitat (logit scale, per unit log richness) ---
slopes_bb_taxo <- emtrends(bb_taxo, ~ habitat, var = "log_rich")
slopes_bb_taxo
pairs(slopes_bb_taxo)

#very similar to other model

#Fringing - Backreef          -0.2929 0.0774 Inf  -3.786  0.0009
#Fringing - Forereef 10m      -0.1913 0.0749 Inf  -2.556  0.0518
#Fringing - Forereef 17m      -0.2117 0.0714 Inf  -2.967  0.0159
#Backreef - Forereef 10m       0.1015 0.0623 Inf   1.629  0.3624
#Backreef - Forereef 17m       0.0811 0.0583 Inf   1.393  0.5038
#Forereef 10m - Forereef 17m  -0.0204 0.0538 Inf  -0.379  0.9814

bb_taxo_25 <- update(bb_taxo, cbind(hits, 25 - hits) ~ .)

AIC(bb_taxo_25, bb_taxo)

res_bb_taxo_25 <- simulateResiduals(bb_taxo_25)
testOutliers(res_bb_taxo_25, type = "bootstrap")
plotResiduals(res_bb_taxo_25, form = d$richness, rank = FALSE)
plotResiduals(res_bb_taxo_25, form = d$habitat)
plotResiduals(res_bb_taxo_25, form = d$year)
testDispersion(res_bb_taxo_25)
plot()

plot(simulateResiduals(bb_taxo_25))                     # look for a low-end deficit
plot()
pairs(emtrends(bb_taxo_25, ~ habitat, var = "log_rich"))
pairs(emtrends(bb_taxo_sy, ~ habitat, var = "log_rich"))

# ---------------------------------------------------------------
# 2. Functional richness (same structure)
# ---------------------------------------------------------------
bb_fun <- update(bb_taxo, . ~ . - log_rich * habitat + log_fun * habitat)

diagnose(bb_fun)
res_bb_fun <- simulateResiduals(bb_fun, n = 1000)
plot(res_bb_fun)
hist(residuals(res_bb_fun), breaks = 40)
plotResiduals(res_bb_fun, form = d$functional_richness, rank = FALSE)
plotResiduals(res_bb_fun, form = d$habitat)
testDispersion(res_bb_fun)


car::Anova(bb_fun, type = "II")
#                    Chisq Df Pr(>Chisq)    
#log_fun         3073.5611  1  < 2.2e-16 ***
#habitat            3.4286  3     0.3301    
#log_fun:habitat   28.1933  3  3.308e-06 ***
anova(update(bb_fun, . ~ . - log_fun:habitat), bb_fun)

slopes_bb_fun <- emtrends(bb_fun, ~ habitat, var = "log_fun")
slopes_bb_fun
pairs(slopes_bb_fun)
# contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef          -0.4246 0.0886 Inf  -4.793 <0.0001
#Fringing - Forereef 10m      -0.2373 0.0832 Inf  -2.852  0.0225
#Fringing - Forereef 17m      -0.3554 0.0777 Inf  -4.571 <0.0001
#Backreef - Forereef 10m       0.1873 0.0721 Inf   2.597  0.0464
#Backreef - Forereef 17m       0.0692 0.0659 Inf   1.050  0.7198
#Forereef 10m - Forereef 17m  -0.1181 0.0581 Inf  -2.034  0.1756

# ---------------------------------------------------------------
# 3. Compare with the beta models (robustness check)
# ---------------------------------------------------------------
# Same direction and ranking of slopes, and same significance of the
# interaction = conclusions are robust to the distribution choice.
cmp <- bind_rows(
  as.data.frame(emtrends(model_occ_log, ~ habitat, var = "log_rich")) %>%
    transmute(habitat, model = "beta", measure = "taxonomic", slope = log_rich.trend),
  as.data.frame(slopes_bb_taxo) %>%
    transmute(habitat, model = "beta-binomial", measure = "taxonomic", slope = log_rich.trend),
  as.data.frame(emtrends(model_occ_log_fun, ~ habitat, var = "log_fun")) %>%
    transmute(habitat, model = "beta", measure = "functional", slope = log_fun.trend),
  as.data.frame(slopes_bb_fun) %>%
    transmute(habitat, model = "beta-binomial", measure = "functional", slope = log_fun.trend)
) %>%
  pivot_wider(names_from = model, values_from = slope)
cmp
# Note: slopes are not numerically identical between models (the
# beta-binomial models hits beyond the first), so compare ranking and
# significance, not exact values.



#### Below here is the code from old verison of the model ####

#### COVER ~ RICHNESS ####

#6985 - richness = 0 OR macroalgal_prop_cover =0 
#6985 - also when richness = 0 AND macroalgal_prop_cover =0 
#1 Fringing. 2575
#2 Backreef 2307
#3 Forereef 10m 1186
#4 Forereef 17m 917

zeros_richness <- alpha_diversity_quad_macro %>%
  filter(macroalgal_prop_cover == 0) %>%
  filter(richness ==  0) %>%
  group_by(habitat, site) %>%
  tally()

zeros_richness %>%
ggplot(aes(x = site, y = n)) +
  geom_col() +
  facet_wrap(~ habitat)

sum(zeros_richness$n)


alpha_diversity_quad_macro_no0s <- alpha_diversity_quad_macro %>%
  filter(macroalgal_prop_cover != 0) %>%
  filter(richness !=  0) 

alpha_diversity_quad_macro_no0s$cover_trans_no0s <- sv_trans(alpha_diversity_quad_macro_no0s$macroalgal_prop_cover)

#okay 12188 = checks out 

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
cover_mod <- glmmTMB(cover_trans ~ richness*habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro)
summary(cover_mod)
car::Anova(cover_mod)

#Chisq Df Pr(>Chisq)    
#richness         28840.4013  1    < 2e-16 ***
#habitat              8.0026  3    0.04596 *  
#richness:habitat   327.8860  3    < 2e-16 ***

testDispersion(cover_mod)
res_cover_mod <- simulateResiduals(cover_mod, n = 1000)
hist(residuals(res_cover_mod))

performance::r2(cover_mod) 
#Conditional R2: 0.648
#Marginal R2: 0.626

#need to check singularity 
performance::check_singularity(cover_mod) 
#(1|site/site_habitat/location) = singular
#(1|site_habitat/location) = not singular, so confirming this is the random effects structure that we want 

em_cover_mod <- emtrends(cover_mod, pairwise ~ habitat, var = "richness") 
em_cover_mod
cld_cover_mod <- multcomp::cld(em_cover_mod$emtrends, Letters = letters, sort = FALSE)
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


## trying with NO 0s 
#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
cover_mod_nozeros <- glmmTMB(cover_trans_no0s ~ richness*habitat + (1|site/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro_no0s)
summary(cover_mod_nozeros)
car::Anova(cover_mod_nozeros)
performance::r2(cover_mod_nozeros) 
#  Conditional R2: 0.469
# Marginal R2: 0.260

cover_mod_nozeros <- glmmTMB(cover_trans_no0s ~ richness + habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro_no0s)
summary(cover_mod_nozeros)
car::Anova(cover_mod_nozeros)
performance::r2(cover_mod_nozeros) 

cover_mod_nozeros <- glmmTMB(macroalgal_prop_cover ~ richness + habitat + (1|site_habitat/location) + (1|year), family = ordbeta(), data = alpha_diversity_quad_macro_no0s)
summary(cover_mod_nozeros)
car::Anova(cover_mod_nozeros)
performance::r2(cover_mod_nozeros) 
#conditional 

#### COVER ~ FUNCTIONAL RICHNESS ####
#check number of functional richness categories
unique(alpha_diversity_quad_macro$functional_richness)

#random effects as site_habitat (aka site ID) instead of site (coarser LTER locations)
cover_mod_fg <- glmmTMB(cover_trans ~ functional_richness*habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro)

summary(cover_mod_fg)
car::Anova(cover_mod_fg)
#Response: cover_trans
#Chisq Df Pr(>Chisq)    
#functional_richness         30210.2372  1     <2e-16 ***
#habitat                         3.0836  3     0.3789    
#functional_richness:habitat   502.7920  3     <2e-16 ***

testDispersion(cover_mod_fg)   
hist(residuals(cover_mod_fg))

#want to check these against each other
# good
hist(residuals(cover_mod_fg)) # good

#need to check singularity 
performance::check_singularity(cover_mod_fg)
##(1|site/site_habitat/location) = singular
#(1|site_habitat/location) = not singular 

performance::r2(cover_mod_fg)
#Conditional R2: 0.642
#Marginal R2: 0.622

em_cover_mod_fg <- emtrends(cover_mod_fg, pairwise ~ habitat, var = "functional_richness") 
#contrast                    estimate     SE  df z.ratio p.value
#Fringing - Backreef            0.236 0.0269 Inf   8.748 <0.0001
#Fringing - Forereef 10m        0.358 0.0268 Inf  13.372 <0.0001
#Fringing - Forereef 17m        0.575 0.0268 Inf  21.494 <0.0001
#Backreef - Forereef 10m        0.122 0.0231 Inf   5.285 <0.0001
#Backreef - Forereef 17m        0.340 0.0231 Inf  14.698 <0.0001
#Forereef 10m - Forereef 17m    0.217 0.0225 Inf   9.667 <0.0001

cld_cover_mod_fg <- multcomp::cld(em_cover_mod_fg$emtrends, Letters = letters, sort = FALSE)

#habitat      functional_richness.trend     SE  df asymp.LCL asymp.UCL .group
#Fringing                          1.88 0.0216 Inf      1.84      1.93  a    
#Backreef                          1.65 0.0167 Inf      1.62      1.68   b   
#Forereef 10m                      1.53 0.0163 Inf      1.49      1.56    c  
#Forereef 17m                      1.31 0.0162 Inf      1.28      1.34     d 

cover_mod_fg_nozeros_beta<- glmmTMB(macroalgal_prop_cover ~ functional_richness*habitat + (1|site_habitat/location) + (1|year), family = ordbeta(), data = alpha_diversity_quad_macro_no0s)
cover_mod_fg_nozeros_ordbeta <- glmmTMB(cover_trans_no0s ~ functional_richness*habitat + (1|site_habitat/location) + (1|year), family = beta_family(), data = alpha_diversity_quad_macro_no0s)

summary(cover_mod_fg_nozeros)
car::Anova(cover_mod_fg_nozeros_beta)

hist(residuals(cover_mod_fg_nozeros))

table(alpha_diversity_quad_macro$richness, alpha_diversity_quad_macro$macroalgal_prop_cover)
table(alpha_diversity_quad_macro$site_habitat, alpha_diversity_quad_macro$richness)
table(alpha_diversity_quad_macro_no0s$habitat, alpha_diversity_quad_macro_no0s$macroalgal_prop_cover)
table(alpha_diversity_quad_macro_no0s$habitat, alpha_diversity_quad_macro_no0s$cover_trans)

table(alpha_diversity_quad_macro_no0s$richness, alpha_diversity_quad_macro_no0s$macroalgal_prop_cover)

table(alpha_diversity_quad_macro$habitat, alpha_diversity_quad_macro$functional_richness)

performance::r2(cover_mod_fg_nozeros)
#Conditional R2: 0.642
#Marginal R2: 0.622

plot(fitted(cover_mod_fg_nozeros), residuals(cover_mod_fg_nozeros))   # look for a fan or trend
plot(fitted(cover_mod_fg_nozeros), cover_mod_fg_nozeros$y)            # or your response column; observed vs fitted
abline(0, 1)

summary(alpha_diversity_quad_macro_no0s$cover_trans_no0s)          # raw column
summary(cover_mod_fg_nozeros$cover_trans_no0s)   #

em_cover_mod_fg_nozeros <- emtrends(cover_mod_fg_nozeros_beta, pairwise ~ habitat, var = "functional_richness") 

cor.test(alpha_diversity_quad_macro_no0s$richness, alpha_diversity_quad_macro_no0s$functional_richness)

#### pull out data from models to make supplemental figures: ####

#don't need to pull out anything from presence model, since Figure 3a plots the model output

t.s2.a <-  as_tibble(cld_presence_model) %>%
  mutate(Predictor = "presence",
         `Spatial scale` = "Plot-level") %>% 
  rename_if(str_detect(names(.), "prob"), ~"Mean")

#richness by cover

t.s2.b <- as_tibble(cld_cover_mod_richness) %>%
  mutate(Predictor = "Taxonomic richness",
         `Spatial scale` = "Plot-level") %>% 
  # rename this so it matches for joining purposes
  rename_if(str_detect(names(.), ".trend"), ~"Mean")
#functional richness by cover

t.s2.c <-  as_tibble(cld_cover_mod_fg) %>%
  mutate(Predictor = "Functional richness",
         `Spatial scale` = "Plot-level") %>% 
  rename_if(str_detect(names(.), ".trend"), ~"Mean")

#joining together

rbind(t.s2.b, t.s2.c) %>% 
  mutate(Table = "S2_plotlevel") %>% 
  dplyr::select(-SE, -df) -> table_s2_plot

supplemental_tables_tableS2_plot <- table_s2_plot %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
         Lower_CI = asymp.LCL,
         Upper_CI = asymp.UCL)



