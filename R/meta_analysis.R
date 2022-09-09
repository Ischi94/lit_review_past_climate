library(merTools)
library(tidyverse)
library(here)
library(distributional)
library(ggdist)
library(lme4)



# general setting ---------------------------------------------------------


# set ggplot output
theme_set(ggplot2::theme_minimal(base_size = 11))

theme_update(text = element_text(colour = "grey20", size = 11),
             legend.text = element_text(colour = "grey20", size = 11),
             legend.title = element_text(colour = "grey20", size = 11),
             axis.text = element_text(colour = "grey20", size = 11),
             axis.ticks = element_blank(),
             strip.text = element_text(colour = "grey20", size = 11),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             plot.title = element_text(colour = "grey20", size = 11))

# define output sizes
image_width <- 183
image_height <- 100
image_units <- "mm"




# read in data and get overview -------------------------------------------


# read in data
dat_meta <- read_csv(here("data",
                          "meta_analysis_legacies.csv"))

# how many different biotic units
dat_meta %>% 
  count(biotic_unit)

# how many different kingdoms
dat_meta %>% 
  count(kingdom)

# how many animals
str_count(dat_meta$kingdom, "animalia") %>% sum(na.rm = TRUE)

# plants
str_count(dat_meta$kingdom, "plantae") %>% sum(na.rm = TRUE)

# and fungi
str_count(dat_meta$kingdom, "fungi") %>% sum(na.rm = TRUE)


# how many studies include climate legacies
dat_meta %>% 
  drop_na(past_climate, quantitative) %>% 
  count(past_climate, quantitative)



# coverage of scales ------------------------------------------------------

# visualize
plot1 <- dat_meta %>% 
  drop_na(scale_in_years) %>% 
  ggplot(aes(scale_in_years, fill = past_climate, 
             colour = past_climate)) +
  geom_dots(
    aes(y = as.numeric(past_climate == "yes"), 
        side = ifelse(past_climate == "yes", "bottom", "top")),
    scale = 0.8) +
  labs(x = "Years", 
       y = NULL) +
  scale_x_continuous(trans = "log10", 
                    breaks = c(10e-8, 10e-4, 
                               10e0, 10e4), 
                    labels = c(expression("10e"^{"-8"}), 
                               expression("10e"^{"-4"}),
                               expression("10e"^{"0"}),
                               expression("10e"^{"4"}))) +
  annotate(geom = "curve", 
           x = 7e-3, xend = 8e-1, 
           y = 0.77, yend = 0.9, 
           curvature = 0.25, 
           arrow = arrow(length = unit(0.05, "inch"), 
                         ends = "last"), 
           colour = "#de970bff", lwd = 0.3) +
  annotate(geom = "text", 
           x = 10e-5, y = 0.78, 
           label = "Studies including\nclimate legacies", 
           colour = "#de970bff",
           size = 11/.pt) +
  scale_fill_manual(values = alpha(c("grey60", "#de970bff"), 0.8)) +
  scale_colour_manual(values = c("grey60", "#de970bff")) +
  theme(axis.text.y = element_blank(), 
        legend.position = "none")


# save scale plot
ggsave(plot1, filename = here("figures",
                              "fig2_meta_scales.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)





# logistic regression for time effect -------------------------------------


model_time <- glm(past_climate == "yes" ~ publish_year, data = dat_meta, family = binomial)


plot2 <- tibble(publish_year = seq(min(dat_meta$publish_year), 
                         max(dat_meta$publish_year),
                         by = 1)) %>% 
  bind_cols(predict(model_time, ., se.fit = TRUE)) %>% 
  mutate(
    # distribution describing uncertainty in log odds
    log_odds = dist_normal(fit, se.fit),
    # inverse-logit transform the log odds to get
    # distribution describing uncertainty in Pr(sex == "male")
    p_legacy = dist_transformed(log_odds, plogis, qlogis)
  ) %>%
  ggplot(aes(x = publish_year)) +
  geom_dots(
    aes(y = as.numeric(past_climate == "yes"),
        side = ifelse(past_climate == "yes", "bottom", "top"), 
        fill = past_climate, 
        colour = past_climate),
    scale = 0.7,
    data = dat_meta %>% drop_na(past_climate), 
    dotsize = 1.2
  ) +
  stat_lineribbon(
    aes(ydist = p_legacy), 
    alpha = 0.35, 
    fill = "#de970bff", 
    colour = "grey20"
  ) +
  labs(
    x = "Date published",
    y = "Pr (Climate Legacy)"
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_fill_manual(values = alpha(c(rep("grey10", 3), 
                               "grey60", "#de970bff"), 0.8)) +
  scale_colour_manual(values = c("grey60", "#de970bff")) +
  theme(legend.position = "none")
  
# save logistic plot
ggsave(plot2, filename = here("figures",
                              "fig3_meta_year.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)


# get model summaries
confint(model_time)
coef(model_time)

# the probability of climate legacies being included at a particular year
tibble(publish_year = 2023) %>% 
  bind_cols(predict(model_time, ., se.fit = TRUE, type = "response"))




# conduct meta-analysis ---------------------------------------------------


# read in effect size data
dat_cohens_d <- read_csv(here("data", 
                              "cohens_d_effect.csv"))

# model as a mixed effect model
model_meta <- lmer(cohens_d ~ 1 + (1 | study),
                   weights = 1 / cohens_d_var, 
                   data = dat_cohens_d) 

# extract model estimates
dat_cohen_res <- predictInterval(model_meta, 
                newdata = dat_cohens_d %>% 
                  distinct(study) %>% 
                  as.data.frame(),
                which = "full", 
                include.resid.var = FALSE) %>% 
  as_tibble() %>% 
  add_column(study = dat_cohens_d %>% 
               distinct(study) %>% 
               pull()) %>% 
  select(study, mean_est = fit, lwr, upr)

# calculate overall effect
dat_cohen_ov <- FEsim(model_meta, 1e4) %>% 
  as_tibble() %>% 
  mutate(lwr = mean - 1.96*sd, 
         upr = mean + 1.96*sd) %>% 
  add_column(study = " ") %>% 
  select(study, mean_est = mean, lwr, upr) %>% 
  bind_rows(dat_cohen_res) %>% 
  # add scale of the focal climate legacy from spreadsheet
  add_column(scale = c(4000.5, 1, 1, 
                       8000, 10000000, 
                       1, 6000000))

# visualize as a forest plot
plot3 <- dat_cohen_ov %>% 
  mutate(col_lev = if_else(study == " ", 
                           "yes", "no")) %>% 
  ggplot() +
  geom_linerange(aes(xmin = lwr, 
                     xmax = upr, 
                     y = scale, 
                     colour = col_lev), 
                 position = position_nudge(y = c(0, -0.1, 0, 0, 
                                                 0, 0.1, 0)), 
                 alpha = 0.5) +
  geom_point(aes(mean_est, scale, 
                 fill = col_lev), 
             position = position_nudge(y = c(0, -0.1, 0, 0, 
                                             0, 0.1, 0)), 
             shape = 21, 
             size = 3, 
             colour = "grey30") +
  geom_text(aes(mean_est + 0.85, scale, label = study), 
            position = position_nudge(y = c(0, -0.4, 0, 0, 
                                 0.1, 0.4, -0.1)), 
            colour = "grey30",
            size = 11/.pt) +
  annotate("text", x = 1.02, y = 1600, 
           label = "Overall", size = 11/.pt,
           colour = "#de970bff") +
  annotate(geom = "curve",
           x = c(0, 0.51, 0.81, 1.21), 
           xend = c(0.49, 0.79, 1.19, 2),
           y = 0.095, yend = 0.095,
           curvature = 0,
           arrow = arrow(length = unit(0.05, "inch"),
                         ends = "both"),
           colour = "grey30", lwd = 0.3) +
  annotate("label", x = c(0.25, 0.65,
                         1, 1.58), y = 0.1, 
           label = c("small", "medium", 
                     "large", "very large"),
           size = 10/.pt, label.size = 0,
           label.padding = unit(0.1, "lines"),
           colour = "grey30") +
  labs(y = "Temporal scale of the climate legacy\nin years", 
       x = "Effect size expressed as Cohen's d") +
  coord_cartesian(xlim = c(0, 2.2)) +
  scale_color_manual(values = c("grey60", "#de970bff")) +
  scale_fill_manual(values = c("grey70", "#de970bff")) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(10e0, 10e2, 
                                10e4, 10e6), 
                     labels = c(expression("10e"^{"0"}), 
                                expression("10e"^{"2"}),
                                expression("10e"^{"4"}),
                                expression("10e"^{"6"}))) +
  theme(legend.position = "none")

# save forest plot
ggsave(plot3, filename = here("figures",
                              "fig4_meta_effect.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)



# funnel plot -------------------------------------------------------------

# calculate average effect from meta-model
av_effect <- FEsim(model_meta, 1e4)

# create funnel plot
plot4 <- av_effect %>% 
  as_tibble() %>% 
  expand_grid(effect_se = seq(0, 0.28, by = 0.01)) %>% 
  mutate(lwr = mean - 1.96*effect_se, 
         upr = mean + 1.96*effect_se) %>% 
  ggplot(aes(y = effect_se)) +
  geom_line(aes(x = lwr), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_line(aes(x = upr), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_line(aes(x = mean, y = effect_se), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_point(aes(mean_est, effect_se), 
             shape = 21, 
             size = 4,
             fill = "coral2", 
             colour = "grey20",
             alpha = 0.8,
             data = dat_cohen_res %>% 
               mutate(effect_se = (upr - lwr) / 3.92)) +
  scale_y_reverse() +
  labs(y = "Standard error", 
       x = "Effect size expressed as Cohen's d")



# kendall´s rank correlation ------------------------------------------------

dat_cohen_res %>% 
  mutate(effect_se = (upr - lwr) / 3.92) %>% 
  { cor.test(.$mean_est, .$effect_se, method="kendall") }
