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

#### Figure 3A: Presence plot ####

# Model estimates: probability of an empty plot by habitat

pres_pred <- emmeans(presence_model, ~ habitat, type = "response") %>%
  as.data.frame()

# Standardize column names for plotting
names(pres_pred)[names(pres_pred) %in% c("prob", "response", "emmean")] <- "fit"
names(pres_pred)[names(pres_pred) %in% c("asymp.LCL", "lower.CL")]      <- "lcl"
names(pres_pred)[names(pres_pred) %in% c("asymp.UCL", "upper.CL")]      <- "ucl"

# Observed proportion of occupied plots at each site
pres_obs <- df %>%
  group_by(habitat, site_habitat) %>%
  summarise(prop_present = mean(present), n = n(), .groups = "drop")

habitat_labels_wrapped <- c("Fringing"     = "Fringing\nreef",
                    "Backreef"     = "Back\nreef",
                    "Forereef 10m" = "Fore reef\n10 m",
                    "Forereef 17m" = "Fore reef\n17 m")


figure3a <- ggplot() +
  #geom_jitter(data = pres_obs,
         #    aes(x = habitat, y = prop_present * 100, colour = habitat),
          #  width = 0.12, height = 0, alpha = 0.5, size = 2) +
  geom_errorbar(data = pres_pred,
                aes(x = habitat, ymin = lcl * 100, ymax = ucl * 100, color = habitat),
                width = 0.15, linewidth = 0.8) +
  geom_point(data = pres_pred,
             aes(x = habitat, y = fit * 100, color = habitat, fill = habitat),
             size = 3.5, shape = 21, stroke = 1.2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "Probability of macroalgal presence (%)") +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  scale_x_discrete(labels = habitat_labels_wrapped) +
  model_themes +
  theme(legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=22))
       
        
        # axis.text.x = element_text(angle = 15, hjust = 1))

#### Figure 3B: COVER ~ RICHNESS ####
# SETUP

rich_grid <- seq(1, max(alpha_diversity_quad_macro_present$richness), length.out = 100)

pred_taxo_richness <- emmeans(model_occ_log, ~ log_rich| habitat,
                                    at = list(log_rich = log(rich_grid)),
                                    component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(taxo_richness = exp(log_rich)) %>%
  left_join(ranges, by = "habitat") %>%
  filter(taxo_richness >= rmin, taxo_richness <= rmax)

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_taxo_richness)[names(pred_taxo_richness) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_taxo_richness)[names(pred_taxo_richness) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

cover_emm <- emmip(model_occ_log, habitat ~ log_rich,
                   at = list(richness = 1:6),
                   type = "response", plotit = FALSE, CIs = TRUE)

Figure3b <- ggplot() +
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
  scale_x_continuous(breaks = 1:max(occ$richness)) +
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

Figure3b

#### Figure 3C: COVER ~ FUNCTIONAL RICHNESS ####
# SETUP
pred_functional_richness <- emmeans(m_occ_log_functional, ~ log_function_rich| habitat,
                at = list(log_function_rich = log(rich_grid_fun)),
                component = "cond", type = "response") %>%
  as.data.frame() %>%
  mutate(functional_richness = exp(log_function_rich)) %>%
  left_join(ranges, by = "habitat") %>%
  filter(functional_richness >= rmin, functional_richness <= rmax)

# emmeans names the CI columns asymp.LCL/asymp.UCL or lower.CL/upper.CL
# depending on version; standardize them
names(pred_functional_richness)[names(pred_functional_richness) %in% c("asymp.LCL", "lower.CL")] <- "lcl"
names(pred_functional_richness)[names(pred_functional_richness) %in% c("asymp.UCL", "upper.CL")] <- "ucl"

str(pred_functional_richness)

Figure3c <- ggplot() +
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
  scale_x_continuous(breaks = 1:max(occ$richness)) +
  scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
  scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
  theme_classic() +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(y = "Cover", x = "Functional richness", title = "") +
  model_themes +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.title.x   = element_text(size = 24,
        margin = margin(t = -7)),
   # axis titles
  axis.text    = element_text(size = 8),    # tick labels
  legend.text  = element_text(size = 8),    # legend labels
  legend.title = element_text(size = 9),    # legend title
  plot.title   = element_text(size = 11)) # plot title)


Figure3c

### COMBINING TOGETHER ####

#### PANEL GRAPH ####

pfigure_3a <- figure3a + theme(plot.margin = margin(10, 10, 10, 10))
pfigure_3b <- Figure3b + theme(plot.margin = margin(10, 10, 10, 10))
pfigure_3c <- Figure3c + theme(plot.margin = margin(10, 10, 10, 10))

leg <- get_legend(pfigure_3b)

(figure_3 <- ggarrange(pfigure_3a, pfigure_3b, pfigure_3c,  
                       ncol = 3, nrow = 1, 
                       legend.grob = leg,
                       common.legend = TRUE,
                       align = "hv",
                       labels = c("a.", "b.", "c."),
                       label.x = 0.02,   # left padding (0 = hard left, 1 = hard right)
                       label.y = 0.98,   # top padding (0 = bottom, 1 = top)
                       hjust = 0, vjust = 1,
                       legend = "bottom", 
                       font.label = list(size = 26, color = "black", face = "plain")) )

ggsave(filename = "figures/figure3_10182026.jpg", figure_3, height = 8, width = 20)

### OLD BELOW ####

#### COVER ~ RICHNESS ####
# SETUP

diversity_ranges_quad <- extract_ranges(alpha_diversity_quad_macro_no0s, "habitat", c("richness", "cover_trans_no", "functional_richness"))

cover_emm <- emmip(cover_mod,
                   habitat ~ richness,
                   at = list(richness = seq(0,6,1)), plotit = F, CIs = T)

filtered_cover_effects <- filter_ranges(trend = cover_emm, range_obj = diversity_ranges_quad, group = "habitat", value = "richness")

range(cover_emm$yvar)            # if this is already within 0–1, don't inv_logit it
range(inv_logit(cover_emm$yvar))
range(inv_logit(cover_emm$UCL))

# PLOT
(figure_2a <- 
    ggplot() +
#    geom_point(data = alpha_diversity_quad_macro,
#               aes(x = richness, y = cover_trans, colour = habitat), size = 1, alpha = 0.75) + 
    geom_line(data = filtered_cover_effects, aes(x = richness, y = inv_logit(yvar), colour = habitat)) +
    geom_ribbon(data = filtered_cover_effects, aes(x = richness, ymin = inv_logit(LCL), ymax = inv_logit(UCL), fill = habitat),
                alpha = 0.5) +
    scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
    scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
    labs(y = "Cover\n", x = "Taxonomic richness", title = "") +
    model_themes 
)

#### COVER ~ FUNCTIONAL RICHNESS ####
# SETUP
cover_emm_fg <- emmip(cover_mod_fg_nozeros,
                      habitat ~ functional_richness,
                      at = list(functional_richness = seq(0,10,1)), plotit = F, CIs = T)

filtered_cover_effects_fg <- filter_ranges(cover_emm_fg, diversity_ranges_quad, "habitat", "functional_richness")

# PLOT
(figure_2b <- 
    ggplot() +
    #geom_jitter(data = alpha_diversity_quad_macro,
             #  aes(x = functional_richness, y = cover_trans, colour = habitat), width = 0.2, size = 1, alpha = 0.1) + 
    geom_line(data = filtered_cover_effects_fg, aes(x = functional_richness, y = inv_logit(yvar), colour = habitat)) +
    geom_ribbon(data = filtered_cover_effects_fg, aes(x = functional_richness, ymin = inv_logit(LCL), ymax = inv_logit(UCL), fill = habitat),
                alpha = 0.5) +
    scale_colour_manual(values = habitat_colours, labels = habitat_labels) +
    scale_fill_manual(values = habitat_colours, labels = habitat_labels) +
    labs(y = "", x = "Functional richness", title = "") +
    model_themes
)
#ggsave(filename = "output/figure_3b_quad_v2.png", figure_3b, height = 10, width = 14)

#### PANEL GRAPH ####

pfigure_2a <- figure_2a + theme(plot.margin = margin(10, 10, 10, 10))
pfigure_2b <- figure_2b + theme(plot.margin = margin(10, 10, 10, 10))

(figure_2 <- ggarrange(pfigure_2a, pfigure_2b,  
                       ncol = 2, nrow = 1, 
                       common.legend = TRUE,
                       labels = c("a.", "b."),
                       label.x = 0.02,   # left padding (0 = hard left, 1 = hard right)
                       label.y = 0.98,   # top padding (0 = bottom, 1 = top)
                       hjust = 0, vjust = 1,
                       legend = "bottom", 
                       font.label = list(size = 26, color = "black", face = "plain")) )

figure2_attempt2 <- pfigure_2a + pfigure_2b +  plot_layout(guides = "collect") &          # collect legends into one
  theme(legend.position = "bottom")  

#this ended up being figure 3 in the MS, so saving it as so... 
#ggsave(filename = "figures/figure3_06252026_final.jpg", figure_2, height = 8, width = 16)

#adding supplemental figures

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

#ggsave(filename = "figures/Supp_FigS2_10192026_final.jpg", s2_figure, height = 10, width = 22)


