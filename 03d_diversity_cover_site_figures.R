
rm(list = ls())
source("03c_diversity_cover_site_models.R")
library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(effects)
library(DHARMa)
library(ggplot2)
library(ggpubr)

#### General Setup ####

diversity_ranges_site <- extract_ranges(alpha_diversity_site_macro_above0, "habitat", c("richness", "cover_trans_no0s", "functional_richness"))

#### Taxonomic Richness at Site - Figure S3a ####

#using integers (predict at whole number)
rich_grid_site <- 1:max(alpha_diversity_site_macro_above0$richness)

pred_taxo_richness_site <- emmeans(site_cover_mod4, ~ log_rich_site| habitat,
                              #at = list(log_rich = log(1:6)),
                              at = list(log_rich_site = log(rich_grid_site)),
                              component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(taxo_richness_site = round(exp(log_rich_site), 6)) %>%
  left_join(diversity_ranges_site, by = "habitat") %>%
  filter(taxo_richness_site >= min_richness, taxo_richness_site <= max_richness)

table(pred_taxo_richness_site$habitat, pred_taxo_richness_site$taxo_richness_site)
#make sure that max richness matches data

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_taxo_richness_site)[names(pred_taxo_richness_site) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_taxo_richness_site)[names(pred_taxo_richness_site) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

#maintain habitat order
pred_taxo_richness_site$habitat <- factor(pred_taxo_richness_site$habitat,
                                                    levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

figure_S3a <- ggplot() +
  # geom_jitter(data = occ, #don't want raw data
  #            aes(x = richness, y = cover_p * 100, colour = habitat),
  #           width = 0.15, height = 0, alpha = 0.15, size = 1) +
  geom_ribbon(data = pred_taxo_richness_site,
              aes(x = taxo_richness_site, ymin = lcl, ymax = ucl,
                  fill = habitat),
              alpha = 0.25) +
  geom_line(data = pred_taxo_richness_site,
            aes(x = taxo_richness_site, y = response, colour = habitat),
            linewidth = 1) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0.1, 0.4, by = 0.1)) +
  #ylim(0,0.35) +
  #ylim(0.05,0.55) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_site_macro_above0$richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic(base_size = 12) +
  labs(y = "Cover", x = "Taxonomic richness", title = "") +
  model_themes +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)

figure_S3a


#### Functional Richness at Site - Figure S3b ####

#using integers (predict at whole number)
functional_rich_grid_site <- 1:max(alpha_diversity_site_macro_above0$functional_richness)

pred_functional_richness_site <- emmeans(site_cover_mod_fg5, ~ log_fun_site| habitat,
                                   #at = list(log_rich = log(1:6)),
                                   at = list(log_fun_site = log(functional_rich_grid_site)),
                                   component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(functional_richness_site = round(exp(log_fun_site), 6)) %>%
  left_join(diversity_ranges_site, by = "habitat") %>%
  filter(functional_richness_site >= min_functional_richness, functional_richness_site <= max_functional_richness)

#maintain habitat order
pred_functional_richness_site$habitat <- factor(pred_functional_richness_site$habitat,
                                          levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))

table(pred_functional_richness_site$habitat, pred_functional_richness_site$functional_richness_site)
#make sure that max richness matches data

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_functional_richness_site)[names(pred_functional_richness_site) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_functional_richness_site)[names(pred_functional_richness_site) %in% c("asymp.UCL", "upper.CL")] <- "ucl"


figure_S3b <- ggplot() +
  geom_ribbon(data = pred_functional_richness_site,
              aes(x = functional_richness_site, ymin = lcl, ymax = ucl,
                  fill = habitat),
              alpha = 0.25) +
  geom_line(data = pred_functional_richness_site,
            aes(x = functional_richness_site, y = response, colour = habitat),
            linewidth = 1) +
 # ylim(0,0.35) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_site_macro_above0$functional_richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic(base_size = 12) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0.1, 0.4, by = 0.1)) +
  labs(y = "Cover", x = "Functional richness", title = "") +
  model_themes +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)

figure_S3b

#### Combine Panels Together to Make Figure S3 ####

(figure_S3 <- ggarrange(figure_S3a, figure_S3b,  
                       ncol = 2, nrow = 1, 
                       common.legend = TRUE,
                       labels = c("a.", "b."),
                       legend = "bottom", 
                       font.label = list(size = 26, color = "black", face = "plain")))

#save

#ggsave(filename = "figures/Supp_FigS3_10272026.jpg", figure_S3, height = 10, width = 18)

