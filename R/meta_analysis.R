library(tidyverse)
library(here)
library(distributional)
library(ggdist)
library(metafor)
library(lme4)
library(patchwork)


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


# add colour indicator column
dat_meta_plot <- dat_meta %>%
  drop_na(scale_in_years) %>%
  mutate(
    colour_ind = case_when(
      past_climate == "no" ~ "a",
      past_climate == "yes" & quantitative == "no" ~ "b",
      past_climate == "yes" & quantitative == "yes" ~ "c"
    )
  )


# visualize
plot1 <- dat_meta_plot %>%
  ggplot(aes(scale_in_years, fill = past_climate, 
             colour = colour_ind)) +
  geom_dots(
    aes(y = as.numeric(past_climate == "yes"), 
        side = ifelse(past_climate == "yes", "bottom", "top")),
    scale = 3, 
    position = position_dodgejust(width = -0.2)) +
  labs(x = "Temporal scale in years", 
       y = NULL) +
  scale_x_continuous(trans = "log10", 
                    breaks = c(10e-8, 10e-4, 
                               10e0, 10e4), 
                    labels = c(expression("10e"^{"-8"}), 
                               expression("10e"^{"-4"}),
                               expression("10e"^{"0"}),
                               expression("10e"^{"4"}))) +
  annotate(geom = "curve", 
           x = 7e-3, xend = 9e0, 
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
  annotate(geom = "curve", 
           x = 3e5, xend = 25e3, 
           y = 0.75, yend = 0.95, 
           curvature = -0.2, 
           arrow = arrow(length = unit(0.05, "inch"), 
                         ends = "last"), 
           colour = "grey10", lwd = 0.3) +
  annotate(geom = "text",
           x = 30e4, y = 0.65,
           label = "Studies quantifying\nclimate legacies",
           colour = "grey10",
           size = 11/.pt) +
  scale_fill_manual(values = alpha(c("grey60", "#de970bff"), 0.8)) +
  scale_colour_manual(values = c("grey60", "#de970bff", "grey10")) +
  theme(axis.text.y = element_blank(), 
        legend.position = "none")


# # save scale plot
# ggsave(plot1, filename = here("figures",
#                               "fig2_meta_scales.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)





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
        colour = colour_ind),
    scale = 0.7,
    data = dat_meta_plot, 
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
  scale_colour_manual(values = c("grey60", "#de970bff", "grey10")) +
  theme(legend.position = "none")
  
# # save logistic plot
# ggsave(plot2, filename = here("figures",
#                               "fig3_meta_year.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)


# get model summaries
confint(model_time)
coef(model_time)

# the probability of climate legacies being included at a particular year
tibble(publish_year = 2023) %>% 
  bind_cols(predict(model_time, ., se.fit = TRUE, type = "response"))




# conduct meta-analysis ---------------------------------------------------


# read in effect size data
dat_cohens_d <- read_csv(here("data", 
                              "cohens_d_effect.csv")) %>% 
  # add id col
  rownames_to_column("id")

# model as a mixed effect model
model_metafor <- rma.mv(
  yi = cohens_d,
  V = cohens_d_var,
  data = dat_cohens_d,
  random = list( ~  1 | study, 
                 ~ 1 | id))

model_metafor %>% 
  ranef() %>% 
  pluck("study") %>% 
  as_tibble(rownames = "study") %>% 
  mutate(across(c(intrcpt, pi.lb, pi.ub), 
                ~ .x + coef(model_metafor)))


# extract estimates
dat_metafor_ov  <- model_metafor %>% 
  predict() %>% 
  as_tibble() %>% 
  add_column(study = " ") %>% 
  select(mean_est = pred, 
         ci.lb, ci.ub, study) %>% 
  # add study wise effect
  bind_rows(model_metafor %>% 
              ranef() %>% 
              pluck("study") %>% 
              as_tibble(rownames = "study") %>% 
              mutate(across(c(intrcpt, pi.lb, pi.ub), 
                            ~ .x + coef(model_metafor))) %>% 
              select(mean_est = intrcpt, 
                     ci.lb = pi.lb, 
                     ci.ub = pi.ub, 
                     study)) %>% 
      # add scale of the focal climate legacy from spreadsheet
      add_column(scale = c(1400, 37, 57, 
                           21000, 6000000, 
                           1, 6000000))



# visualize as a forest plot
plot3 <- dat_metafor_ov %>%
  # add reference system from proceedings B
  left_join(tibble(study = dat_metafor_ov$study, 
                   study_ref = c("", 
                                 "(29)", 
                                 "(31)", 
                                 "(32)", 
                                 "(34)", 
                                 "(30)", 
                                 "(33)"))) %>% 
  mutate(col_lev = if_else(study == " ",
                           "yes", "no")) %>%
  ggplot() +
  geom_linerange(aes(xmin = ci.lb, 
                      xmax = ci.ub, 
                      y = scale, 
                      colour = col_lev), 
                  position = position_nudge(y = c(0, -0.2, 0, 0, 
                                                  0, 0.2, 0.2)), 
                  alpha = 0.5) +
  geom_point(aes(mean_est, scale, 
                 fill = col_lev), 
             position = position_nudge(y = c(0, -0.2, 0, 0, 
                                             0, 0.2, 0.2)), 
             shape = 21, 
             size = 3, 
             colour = "grey30") +
  geom_text(aes(mean_est, scale, label = study_ref), 
            position = position_nudge(y = c(0.4, -0.5, 0.4, 0.4,
                                            -0.3, 0.6,0.5)), 
            colour = "grey30",
            size = 10/.pt) +
  annotate("text", x = 1.5, y = 4000, 
           label = "Overall", size = 11/.pt,
           colour = "#de970bff") +
  annotate(geom = "curve",
           x = c(0, 0.51, 0.81, 1.21, 2.01), 
           xend = c(0.49, 0.79, 1.19, 2, 3.3),
           y = 0.095, yend = 0.095,
           curvature = 0,
           arrow = arrow(length = unit(0.05, "inch"),
                         ends = c(rep("both", 4), "first")),
           colour = "grey30", lwd = 0.3) +
  annotate("label", x = c(0.25, 0.65,
                         1, 1.58, 2.6), y = 0.1, 
           label = c("small", "medium", 
                     "large", "very large", 
                     "huge"),
           size = 10/.pt, label.size = 0,
           label.padding = unit(0.1, "lines"),
           colour = "grey30") +
  labs(y = "Temporal scale of the climate legacy\nin years", 
       x = "Effect size expressed as Cohen's d") +
  coord_cartesian(xlim = c(0, 3)) +
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

# # saveplot
# ggsave(plot3, filename = here("figures",
#                               "fig4_meta_effect.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)



# funnel plot -------------------------------------------------------------

dat_cohens_d %>%
  mutate(cohens_d_sd = abs(sqrt(cohens_d_var))) %>% 
  ggplot(aes(cohens_d, 
             -cohens_d_sd)) +
  geom_point() +
  geom_vline(xintercept = mean(dat_cohens_d$cohens_d)) +
  geom_abline(intercept = -2, slope = )


funnel(model_metafor) %>% 
  as_tibble()

# create funnel plot
plot4 <- coef(model_metafor) %>%
  as_tibble() %>% 
  expand_grid(effect_se = seq(0, 2, by = 0.5)) %>% 
  mutate(lwr = value - 1.96*effect_se, 
         upr = value + 1.96*effect_se) %>% 
  ggplot(aes(y = effect_se)) +
  geom_line(aes(x = lwr), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_line(aes(x = upr), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_line(aes(x = value, y = effect_se), 
            linetype = "dotted", 
            colour = "grey20") +
  geom_point(aes(mean_est, effect_se, 
                 fill = study), 
             shape = 21, 
             size = 4,
             colour = "grey20",
             alpha = 0.8,
             data = dat_cohens_d %>% 
               mutate(effect_se = abs(sqrt(cohens_d_var))) %>% 
               rename(mean_est = cohens_d) ) +
  scale_y_reverse(limits = c(2, 0), 
                  breaks = c(0, 1, 2)) +
  scale_fill_brewer(type = "qual", 
                    palette = 8, 
                    name = "Study") +
  labs(y = "Standard error", 
       x = "Effect size expressed as Cohen's d")

# # save plot
# ggsave(plot4, filename = here("figures",
#                               "figs1_funnel_plot.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)



# kendall´s rank correlation 

ranktest(model_metafor)



# save plots --------------------------------------------------------------

fig2 <- plot1 /
  plot2 +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")") +
  plot_layout(heights = c(1.5, 1))


# save plot
ggsave(fig2, filename = here("figures",
                              "fig2.png"),
       width = image_width, height = image_height*1.5,
       units = image_units,
       bg = "white", device = ragg::agg_png)


# save plot
ggsave(plot3, filename = here("figures",
                             "fig3.png"),
       width = image_width, height = image_height,
       units = image_units,
       bg = "white", device = ragg::agg_png)


# save plot
ggsave(plot4, filename = here("figures",
                              "figs1_funnel_plot.png"),
       width = image_width, height = image_height,
       units = image_units,
       bg = "white", device = ragg::agg_png)
