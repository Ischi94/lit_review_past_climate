library(esc)
library(tidyverse)
library(here)

 
# convert effect sizes ----------------------------------------------------


# go through effect estimates of individual studies and convert to cohen's d
# Butler et al.
# unstandardized regression coefficient
study_1 <- tibble(beta_coef = c(1.068, 0.011, -0.011),
       sd = c(0.313, 0.009, 0.011), 
       sample_size = rep(55, 3)) %>% 
  mutate(cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                             grp1n = sample_size,
                             grp2n = sample_size)[[3]]) %>% 
  add_column(study = "Butler et al. 2017") %>% 
  select(study, cohens_d, cohens_d_var)


# Imperio et al. 
# unstandardized regression coefficient
study_2 <- tibble(beta_coef = c(rep(-0.19, 6),
                     rep(-0.18, 6), 
                     -0.20, -0.14, 
                     rep(-0.17, 3)),
       sd = c(rep(0.04, 16),
              0.05),
       sample_size = rep(15, 17)) %>% 
  mutate(cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                              grp1n = sample_size,
                              grp2n = sample_size)[[3]]) %>% 
  add_column(study = "Imperio et al. 2013") %>% 
  select(study, cohens_d, cohens_d_var)


# Varela et al.
# students t-test
study_3 <- tibble(cohens_d = esc_t(8.324, 0.0001, 273)[[1]], 
       cohens_d_var = esc_t(8.324, 0.0001, 273)[[3]]) %>% 
  add_column(study = "Varela et al. 2015", .before = 1)


# Mayhew et al.
# unstandardized correlation coefficient (pearson)
study_4 <- tibble(beta_coef = c(-0.512, -0.376, -0.166,
                     -0.427, -0.474, -0.414,
                     -0.009, -0.182, -0.190, 
                     -0.497, -0.361, -0.497,
                     -0.361),
       sd = sd(beta_coef)) %>% 
  mutate(sample_size = nrow(.), 
         cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                              grp1n = sample_size,
                              grp2n = sample_size)[[3]])  %>% 
  add_column(study = "Mayhew et al. 2008") %>% 
  select(study, cohens_d, cohens_d_var)




# Griebeler et al.
# unstandardized correlation coefficient (pearson)
study_5 <- tibble(beta_coef = c(0.98, 0.898, 0.248, 0.242, 0.271, 0.343),
       sd = sd(beta_coef)) %>% 
  mutate(sample_size = nrow(.),
         cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                              grp1n = sample_size,
                              grp2n = sample_size)[[3]]) %>% 
  add_column(study = "Griebeler & Gottschalk 2000") %>% 
  select(study, cohens_d, cohens_d_var)


# Mathes et al.
# percentage increase
dat_perc <- read_csv("https://raw.githubusercontent.com/Ischi94/pal-int-extinction/master/data/results/effect_intensity.csv")

# add sample size manually, based on https://raw.githubusercontent.com/Ischi94/pal-int-extinction/master/data/all_data_trends.csv
dat_perc <- tibble(taxon = dat_perc$taxon,
                   sample_size = c(28819, 11836, 9902,
                                   4632, 2520, 8743,
                                   2360, 3867)) %>%
  full_join(dat_perc)

# calculate sd from CI
study_6 <- dat_perc %>% 
  mutate(warm_sd = sqrt(sample_size) * (warm_high - warm_low), 
         cool_sd = sqrt(sample_size) * (cool_high - cool_low)) %>% 
  select(-contains(c("high", "low"))) %>% 
  # simulate groups with the corresponding sample size and 
  # percentage change
  mutate(ww = pmap(list(n = sample_size, 
                        mean = 1*(1 + warm), 
                        sd = warm_sd), 
                   rnorm), 
         cw = pmap(list(n = sample_size, 
                        mean = 1, 
                        sd = warm_sd), 
                   rnorm), 
         cc = pmap(list(n = sample_size, 
                        mean = 1*(1 + cool), 
                        sd = cool_sd), 
                   rnorm), 
         wc = pmap(list(n = sample_size, 
                        mean = 1, 
                        sd = cool_sd), 
                   rnorm), 
         warm = map_dbl(ww, mean), 
         warm_comp = map_dbl(cw, mean),
         cool = map_dbl(cc, mean), 
         cool_comp = map_dbl(wc, mean)) %>% 
  select(-c(ww, cw, cc, wc)) %>% 
  # calculate cohens d by hand
  mutate(cohens_d_warm = (warm - warm_comp) / sqrt((warm_sd^2 + warm_sd^2)/2), 
         cohens_d_var_warm = (warm_sd + warm_sd/
                                warm_sd*warm_sd) + 
           (cohens_d_warm^2 / 2*(warm_sd + warm_sd)), 
         cohens_d_cool = (cool - cool_comp) / sqrt((cool_sd^2 + cool_sd^2)/2), 
         cohens_d_var_cool = (cool_sd + cool_sd/
                                cool_sd*cool_sd) + 
           (cohens_d_cool^2 / 2*(cool_sd + cool_sd)), 
         across(c(cohens_d_warm, cohens_d_cool), abs))


  
# reshape
study_6 <- tibble(study = "Mathes, van Dijk, et al. 2021", 
       cohens_d = c(study_6$cohens_d_warm, 
                    study_6$cohens_d_cool), 
       cohens_d_var = c(study_6$cohens_d_var_warm,
                        study_6$cohens_d_var_cool))


# combine and save --------------------------------------------------------

bind_rows(study_1, study_2, study_3, 
          study_4, study_5, study_6) %>% 
  write_csv(here("data", 
                 "cohens_d_effect.csv"))
