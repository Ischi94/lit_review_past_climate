library(esc)
library(tidyverse)

# Butler et al.
# unstandardized regression coefficient
tibble(beta_coef = c(1.068, 0.011, -0.011),
       sd = c(0.313, 0.009, 0.011), 
       sample_size = rep(55, 3)) %>% 
  mutate(cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                             grp1n = sample_size,
                             grp2n = sample_size)[[3]]) %>% 
  mutate(cohens_d_sd = sqrt(cohens_d_var)) %>% 
  mutate(dist = list(rnorm(1e5, cohens_d, cohens_d_sd))) %>% 
  unnest(dist) %>% 
  summarise(cohens_d = mean(dist), cohens_d_sd = sd(dist)) 



# Imperio et al. 
# unstandardized regression coefficient
tibble(beta_coef = c(rep(-0.19, 6),
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
  mutate(cohens_d_sd = sqrt(cohens_d_var)) %>% 
  mutate(dist = list(rnorm(1e5, cohens_d, cohens_d_sd))) %>% 
  unnest(dist) %>% 
  summarise(cohens_d = mean(dist), cohens_d_sd = sd(dist))


# Varela et al.
# unstandardized correlation coefficient (pearson)
tibble(cohens_d = esc_t(8.324, 0.0001, 273)[[1]], 
       cohens_d_var = esc_t(8.324, 0.0001, 273)[[3]])


# Mayhew et al.
# unstandardized correlation coefficient (pearson)
tibble(pear_cor = c(-0.512, -0.376, -0.166, -0.427, -0.474, -0.414, -0.009, -0.182, -0.190)) %>% 
  summarise(beta_coef = mean(pear_cor), sd = sd(pear_cor), 
            sample_size = n()) %>% 
  mutate(cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                              grp1n = sample_size,
                              grp2n = sample_size)[[3]])


# Griebeler et al.
# unstandardized correlation coefficient (pearson)
tibble(pear_cor = c(0.98, 0.898, 0.248, 0.242, 0.271, 0.343)) %>% 
  summarise(beta_coef = mean(pear_cor), sd = sd(pear_cor), 
            sample_size = n()) %>% 
  mutate(cohens_d = esc_B(b = beta_coef, sdy = sd,
                          grp1n = sample_size,
                          grp2n = sample_size)[[1]], 
         cohens_d = abs(cohens_d), 
         cohens_d_var = esc_B(b = beta_coef, sdy = sd,
                              grp1n = sample_size,
                              grp2n = sample_size)[[3]])


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
dat_perc %>% 
  mutate(warm_sd = sqrt(sample_size) * (warm_high - warm_low), 
         cool_sd = sqrt(sample_size) * (cool_high - cool_low)) %>% 
  select(-contains(c("high", "low"))) %>% 
  # calculate cohens d
  mutate(cohens_d_warm = esc_mean_sd(grp1m = warm, grp1sd = warm_sd,
                          grp1n = sample_size,
                          grp2m = 0, grp2sd = 1,
                          grp2n = sample_size)[[1]], 
         cohens_d_cool = esc_mean_sd(grp1m = cool, grp1sd = cool_sd,
                                     grp1n = sample_size,
                                     grp2m = 0, grp2sd = 1,
                                     grp2n = sample_size)[[1]]) %>% 
  pivot_longer(c(cohens_d_warm, cohens_d_cool), 
               values_to = "cohens_d_val") %>% 
  mutate(cohens_d_val = abs(cohens_d_val)) %>% 
  summarise(cohens_d = mean(cohens_d_val), cohens_d_sd = sd(cohens_d_val))

