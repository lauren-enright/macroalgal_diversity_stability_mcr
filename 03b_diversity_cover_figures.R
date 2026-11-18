rm(list = ls())
source("03a_diversity_cover_models.R")
library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(effects)
library(DHARMa)
library(ggplot2)
library(patchwork)
library(ggpubr)


diversity_ranges_quad <- extract_ranges(alpha_diversity_quad_macro_present, "habitat", c("richness", "cover_trans_no0s", "functional_richness"))

#### Figure 3A: COVER ~ RICHNESS ####
# SETUP

#using integers (predict at whole number)
rich_grid <- 1:max(alpha_diversity_quad_macro_present$richness)
#smooth lines
#rich_grid <- seq(1, max(alpha_diversity_quad_macro_present$richness), length.out = 100)

pred_taxo_richness <- emmeans(taxo_richness_cover_m4, ~ log_rich| habitat,
                              #at = list(log_rich = log(1:6)),
                               at = list(log_rich = log(rich_grid)),
                                    component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(taxo_richness = exp(log_rich)) %>%
  left_join(diversity_ranges_quad, by = "habitat") %>%
  filter(taxo_richness >= min_richness, taxo_richness <= max_richness)

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_taxo_richness)[names(pred_taxo_richness) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_taxo_richness)[names(pred_taxo_richness) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

table(pred_taxo_richness$habitat, pred_taxo_richness$taxo_richness)


Figure3a <- ggplot() +
  # geom_jitter(data = occ, #don't want raw data
  #            aes(x = richness, y = cover_p * 100, colour = habitat),
  #           width = 0.15, height = 0, alpha = 0.15, size = 1) +
  geom_ribbon(data = pred_taxo_richness,
              aes(x = taxo_richness, ymin = lcl, ymax = ucl,
                  fill = habitat),
              alpha = 0.25) +
  geom_line(data = pred_taxo_richness,
            aes(x = taxo_richness, y = response, colour = habitat),
            linewidth = 1) +
  ylim(0.05,0.55) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_quad_macro_present$richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic(base_size = 12) +
  labs(y = "Cover", x = "Taxonomic richness", title = "") +
  model_themes +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)

Figure3a

#looks good! 

#### Figure 3B: COVER ~ FUNCTIONAL RICHNESS ####
# SETUP

functional_rich_grid <- 1:max(alpha_diversity_quad_macro_present$functional_richness)

pred_functional_richness <- emmeans(functional_richness_cover_m4, ~ log_fun| habitat,
                                   # at = list(log_fun = log(1:4)),
                                    at = list(log_fun = log(functional_rich_grid)),
                component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(functional_richness = exp(log_fun)) %>%
  left_join(diversity_ranges_quad, by = "habitat") %>%
  filter(functional_richness >= min_functional_richness, functional_richness <= max_functional_richness)

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_functional_richness)[names(pred_functional_richness) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_functional_richness)[names(pred_functional_richness) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

str(pred_functional_richness)
table(pred_functional_richness$habitat, pred_functional_richness$functional_richness)

Figure3b <- ggplot() +
  # geom_jitter(data = occ, #don't want raw data
  #            aes(x = richness, y = cover_p * 100, colour = habitat),
  #           width = 0.15, height = 0, alpha = 0.15, size = 1) +
  geom_ribbon(data = pred_functional_richness,
              aes(x = functional_richness, ymin = lcl, ymax = ucl,
                  fill = habitat),
              alpha = 0.25) +
  geom_line(data = pred_functional_richness,
            aes(x = functional_richness, y = response, colour = habitat),
            linewidth = 1) +
  ylim(0.05,0.55) +
  # facet_wrap(~ habitat, nrow = 1) +
  scale_x_continuous(breaks = 1:max(alpha_diversity_quad_macro_present$functional_richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic() +
  labs(y = "Cover", x = "Functional richness", title = "") +
  model_themes  +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        aspect.ratio = 1)


Figure3b

### COMBINING TOGETHER ####

#### PANEL GRAPH ####

pfigure_3a <- Figure3a + theme(plot.margin = margin(10, 10, 10, 10))
pfigure_3b <- Figure3b + theme(plot.margin = margin(10, 10, 10, 10))


(figure_3 <- ggarrange(pfigure_3a, pfigure_3b,
                       ncol = 2, nrow = 1, 
                       common.legend = TRUE,
                       align = "hv",
                       labels = c("a.", "b."),
                       label.x = 0.02,   # left padding (0 = hard left, 1 = hard right)
                       label.y = 0.98,   # top padding (0 = bottom, 1 = top)
                       hjust = 0, vjust = 1,
                       legend = "bottom", 
                       font.label = list(size = 26, color = "black", face = "plain")) )

#ggsave(filename = "figures/figure3_10272026.jpg", figure_3, height = 8, width = 16)


#### FIGURE S2   #### 
# Adding supplemental figures - Figure S2
#want taxonomic richness as first panel to match

supplemental_tables_tableS2_plot <- supplemental_tables_tableS2_plot %>%
  mutate(Predictor = factor(Predictor, 
                          levels = c("Taxonomic richness", "Functional richness")))

supplemental_tables_tableS2_plot %>% 
  filter(Table == "S2_plotlevel") %>% 
  filter(Predictor == "Taxonomic richness") %>% 
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
  xlim(-0.01, 1.3) +
  ggtitle("a. Taxonomic richness") +
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
  guides(color = guide_legend(reverse = TRUE)) -> s2.taxo

supplemental_tables_tableS2_plot %>% 
  filter(Table == "S2_plotlevel") %>% 
  filter(Predictor == "Functional richness") %>% 
  # order the habitats by distance from shore
  mutate(Habitat = factor(Habitat, levels = c("Forereef 17m", "Forereef 10m",  "Backreef", "Fringing"))) %>% 
  ggplot(aes(x = Mean, y = Habitat, color = Habitat)) +
  geom_point(size = 5) +
  # add confidence intervals:
  geom_linerange(aes(xmin = Lower_CI, xmax = Upper_CI), linewidth = 3) +
  scale_color_manual(values = habitat_colours, label = habitat_labels) +
  xlab("Estimated slope") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  ylab("") +
  xlim(-0.01, 1.3) +
  ggtitle("b. Functional richness") +
  scale_y_discrete(labels = c(
    "Fringing" = "Fringing reef",
    "Backreef" = "Back reef",
    "Forereef 10m" = "Fore reef 10 m",
    "Forereef 17m" = "Fore reef 17 m"
  )) +
  model_themes + 
 # theme(legend.position = "none") +
  # add letters denoting significance
  geom_text(aes(x = Upper_CI + 0.08, y = Habitat, 
                label = .group), colour = "black", size = 10) +
  guides(color = guide_legend(reverse = TRUE)) -> s2.functional

s2_figure <- s2.taxo + s2.functional +  plot_layout(guides = "collect") &          # collect legends into one
  theme(legend.position = "bottom")  

#ggsave(filename = "figures/Supp_FigS2_10272026_final.jpg", s2_figure, height = 10, width = 22)


