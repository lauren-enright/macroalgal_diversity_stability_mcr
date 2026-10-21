
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

unique(alpha_diversity_site_macro$habitat)

ranges_site <- alpha_diversity_site_macro_above0 %>%
  group_by(habitat) %>%
  summarise(rmin = min(richness), rmax = max(richness), .groups = "drop")

rich_grid_site <- seq(1, max(alpha_diversity_site_macro_above0$richness), length.out = 100)
max(alpha_diversity_site_macro_above0$richness)

pred_taxo_richness_site <- emmeans(site_cover_mod5, ~ log_richness| habitat,
                              at = list(log_richness = log(rich_grid_site)),
                              component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(taxo_richness_site = exp(log_richness)) %>%
  left_join(ranges_site, by = "habitat") %>%
  filter(taxo_richness_site >= rmin, taxo_richness_site <= rmax)

max(pred_taxo_richness_site$rmax)
# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_taxo_richness_site)[names(pred_taxo_richness_site) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_taxo_richness_site)[names(pred_taxo_richness_site) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

cover_emm_site <- emmip(site_cover_mod5, habitat ~ log_richness,
                   at = list(log_richness = 1:14),
                   type = "response", plotit = FALSE, CIs = TRUE)


figure_S3a2 <- ggplot() +
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
  ylim(0,0.3) +
  #ylim(0.05,0.55) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_site_macro_above0$richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic(base_size = 12) +
  labs(y = "Cover", x = "Taxonomic richness", title = "") +
  model_themes +
  theme(axis.title.x   = element_text(size = 24,
                                      margin = margin(t = -7)),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)

figure_S3a2

#### Taxonomic Richness Setup ####

diversity_ranges_site <- extract_ranges(alpha_diversity_site_macro, "habitat", c("richness", "functional_richness", "cover_trans"))

site_cover_emm <- emmip(site_cover_mod5,
                        habitat ~ log_richness,
                        at = list(log_richness = seq(0,17,1)), plotit = F, CIs = T)

filtered_site_cover_effects <- filter_ranges(site_cover_emm, diversity_ranges_site, "habitat", "log_richness")

#since we are removing raw data, need to remind R what order we want the legend in

filtered_site_cover_effects$habitat <- factor(filtered_site_cover_effects$habitat,
                                             levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))


### Taxonomic Richness Figure - Figure S3a ####

(figure_S3a <- 
   ggplot() +
  # geom_point(data = alpha_diversity_site_macro,
      #        aes(x = richness, y = cover_trans, colour = habitat), size = 1, alpha = 0.3) + 
   geom_line(data = filtered_site_cover_effects, aes(x = richness, y = inv_logit(yvar), colour = habitat), size = 1.5) +
   geom_ribbon(data = filtered_site_cover_effects, aes(x = richness, ymin = inv_logit(LCL), ymax = inv_logit(UCL), fill = habitat),
               alpha = 0.3) +
   scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
   scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
   labs(y = "Cover\n", x = "Taxonomic richness", title = "") +
   scale_y_continuous(
     limits = c(0, 0.4),               # force axis range
     breaks = c(0, 0.2, 0.4)) + 
   trends_theme +
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
)
# ggsave(filename = "output/figure_3d.png", figure_3d, height = 10, width = 14)




#### Functional Richness Setup ####

d_site <- alpha_diversity_site_macro_above0

# Observed range of functional richness in each habitat
ranges_site_fun <- d_site %>%
  group_by(habitat) %>%
  summarise(rmin = min(functional_richness), rmax = max(functional_richness), .groups = "drop")

# Predictions across the observed range
grid_fun <- seq(1, max(alpha_diversity_site_macro_above0$functional_richness), length.out = 100)

pred_fun_site <- emmeans(site_cover_mod_fg6, ~ log_fun_richness | habitat,
                         at = list(log_fun_richness = log(grid_fun)),
                         type = "response") %>%
  summary() %>%
  data.frame() %>%
  mutate(functional_richness_site = exp(log_fun_richness)) %>%
  left_join(ranges_site_fun, by = "habitat") %>%
  filter(functional_richness_site >= rmin, functional_richness_site <= rmax)

summary(alpha_diversity_site_macro_above0$cover_trans)
summary(pred_fun_site$response)

range(pred_fun_site$functional_richness)
range(alpha_diversity_site_macro_above0$functional_richness)

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_fun_site)[names(pred_fun_site) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_fun_site)[names(pred_fun_site) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

fg_max_site <- max(alpha_diversity_site_macro_above0$functional_richness)

site_cover_fg_emm <- emmip(site_cover_mod_fg6,     # your final model
                           habitat ~ log_fun_richness,
                           at = list(log_fun_richness = 1:fg_max_site),
                           type = "response", plotit = FALSE, CIs = TRUE)

figure_S3b2 <- ggplot() +
  # geom_jitter(data = occ, #don't want raw data
  #            aes(x = richness, y = cover_p * 100, colour = habitat),
  #           width = 0.15, height = 0, alpha = 0.15, size = 1) +
  geom_ribbon(data = pred_fun_site,
              aes(x = functional_richness_site, ymin = lcl, ymax = ucl,
                  fill = habitat),
              alpha = 0.25) +
  geom_line(data = pred_fun_site,
            aes(x = functional_richness_site, y = response, colour = habitat),
            linewidth = 1) +
  ylim(0,0.3) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_site_macro_above0$functional_richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic(base_size = 12) +
  labs(y = "Cover", x = "Functional richness", title = "") +
  model_themes +
  theme(axis.title.x   = element_text(size = 24,
                                      margin = margin(t = -7)),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)

figure_S3b2

## old 
site_cover_fg_emm <- emmip(site_cover_mod_fg2,
                           habitat ~ functional_richness,
                           at = list(functional_richness = seq(0,17,1)), plotit = F, CIs = T)


filtered_site_fg_cover_effects <- filter_ranges(site_cover_fg_emm, diversity_ranges_site, "habitat", "log_fun_richness")

#since we are removing raw data, need to remind R what order we want the legend in

filtered_site_fg_cover_effects$habitat <- factor(filtered_site_fg_cover_effects$habitat,
                                              levels = c("Fringing", "Backreef", "Forereef 10m", "Forereef 17m"))


### Functional Richness Figure -  Figure S3b ####

(figure_S3b <- 
    ggplot() +
    #geom_point(data = alpha_diversity_site_macro,
           #    aes(x = functional_richness, y = cover_trans, colour = habitat), size = 1, alpha = 0.3) + 
    geom_line(data = filtered_site_fg_cover_effects, aes(x = functional_richness, y = inv_logit(yvar), colour = habitat),
              linewidth = 1.5) +
    geom_ribbon(data = filtered_site_fg_cover_effects, aes(x = functional_richness, ymin = inv_logit(LCL), ymax = inv_logit(UCL), fill = habitat),
                alpha = 0.3) +
    scale_colour_manual(values = habitat_colours) +
    scale_fill_manual(values = habitat_colours) +
    labs(y = "", x = "Functional richness", title = "") +
   scale_y_continuous(
     limits = c(0, 0.4),               # force axis range
     breaks = c(0, 0.2, 0.4)) + 
    trends_theme +
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
)
#ggsave(filename = "output/figure_3e_v2.png", figure_3e, height = 10, width = 14)

(figure_S3 <- ggarrange(figure_S3a2, figure_S3b2,  
                       ncol = 2, nrow = 1, 
                       common.legend = TRUE,
                       labels = c("a.", "b."),
                       legend = "bottom", 
                       font.label = list(size = 26, color = "black", face = "plain")))

#ggsave(filename = "figures/Supp_FigS3_10202026.jpg", figure_S3, height = 10, width = 18)


t.s2.a_site

supplemental_tables_tableS3_plot <- t.s2.a_site %>% 
  # rename to match previous version/code
  dplyr::rename(Habitat = habitat,
                Lower_CI = asymp.LCL,
                Upper_CI = asymp.UCL)

supplemental_tables_tableS3_plot %>% 
  #filter(Table == "S2_plotlevel") %>% 
  #filter(Predictor == "Taxonomic richness") %>% 
  # order the habitats by distance from shore
  mutate(Habitat = factor(Habitat, levels = c("Forereef 17m", "Forereef 10m",  "Backreef", "Fringing"))) %>% 
  ggplot(aes(x = Mean, y = Habitat, color = Habitat)) +
  geom_point(size = 5) +
  # add confidence intervals:
  geom_linerange(aes(xmin = Lower_CI, xmax = Upper_CI), linewidth = 3) +
  scale_colour_manual(values = habitat_colours, label = habitat_labels) +
  xlab("Estimated slope") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  ylab("") +
  #xlim(-0.05, 1.3) +
  ggtitle("a. Site-level Taxonomic richness") +
  model_themes + 
  # add letters denoting significance
  geom_text(aes(x = Upper_CI + 0.10 , y = Habitat, 
                label = .group), colour = "black", size = 10) +
  scale_y_discrete(labels = c(
    "Fringing"     = "Fringing reef",
    "Backreef"     = "Back reef",
    "Forereef 10m" = "Fore reef 10 m",
    "Forereef 17m" = "Forereef 17 m"
  )) +
  guides(color = guide_legend(reverse = TRUE)) -> s2.taxo_site

s2.taxo_site
