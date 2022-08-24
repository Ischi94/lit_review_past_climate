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

esc_B(1.068, 0.313, 55, 55)[1]

# Imperio et al. 
