library(tidyverse)
library(here)
library(distributional)
library(ggdist)


#set ggplot output
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11))

ggplot2::theme_update(text = element_text(colour = "grey20", size = 11), 
                      legend.text = element_text(colour = "grey20", size = 11),
                      legend.title = element_text(colour = "grey20", size = 11),
                      axis.text = element_text(colour = "grey20", size = 11),
                      axis.ticks = element_blank(),
                      strip.text = element_text(colour = "grey20", size = 11),
                      panel.grid.minor = element_blank(),  
                      panel.grid.major = element_blank(),
                      plot.title = element_text(colour = "grey20", size = 11))


# read in data
dat_meta <- read_csv(here("data",
                          "meta_analysis_legacies.csv"))


dat_meta %>% 
  count(biotic_unit)


dat_meta %>% 
  count(kingdom)

str_count(dat_meta$kingdom, "animalia") %>% sum(na.rm = TRUE)

str_count(dat_meta$kingdom, "plantae") %>% sum(na.rm = TRUE)

str_count(dat_meta$kingdom, "fungi") %>% sum(na.rm = TRUE)


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
  scale_fill_manual(values = alpha(c("grey60", "#de970bff"), 0.8)) +
  scale_colour_manual(values = c("grey60", "#de970bff")) +
  theme(axis.text.y = element_blank(), 
        legend.position = "none")




dat_meta %>% 
  drop_na(past_climate, quantitative) %>% 
  count(past_climate, quantitative)



testi <- glm(past_climate == "yes" ~ publish_year, data = dat_meta, family = binomial)


plot2 <- tibble(publish_year = seq(min(dat_meta$publish_year), 
                         max(dat_meta$publish_year),
                         by = 1)) %>% 
  bind_cols(predict(testi, ., se.fit = TRUE)) %>% 
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
    scale = 0.4,
    data = dat_meta %>% drop_na(past_climate)
  ) +
  stat_lineribbon(
    aes(ydist = p_legacy), 
    alpha = 0.35, 
    fill = "#de970bff"
  ) +
  labs(
    x = "Date published",
    y = "Pr( Climate Legacy )"
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_fill_manual(values = alpha(c(rep("grey10", 3), 
                               "grey60", "#de970bff"), 0.8)) +
  scale_colour_manual(values = c("grey60", "#de970bff")) +
  theme(legend.position = "none")
  

confint(testi)
coef(testi)

tibble(publish_year = 1980) %>% 
  bind_cols(predict(testi, ., se.fit = TRUE, type = "response"))
