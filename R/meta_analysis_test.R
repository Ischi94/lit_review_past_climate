# add variables
dat_cohens_d <- dat_cohens_d %>% 
  mutate(year = word(study, -1), 
         year = as.integer(year)) %>% 
  left_join(dat_meta %>% 
              filter(quantitative == "yes") %>% 
              mutate(year = as.integer(publish_year)) %>%  
              select(year, biotic_unit, 
                     kingdom, past_scale_in_year)) %>% 
  mutate(kingdom = as.numeric(as.factor(kingdom)), 
         biotic_unit = as.numeric(as.factor(biotic_unit)))


# model as a mixed effect model with biotic unit as moderator
model_metafor_biounit <- rma.mv(
  yi = cohens_d,
  V = cohens_d_var,
  mods = ~biotic_unit,
  data = dat_cohens_d,
  random = list( ~  1 | study, 
                 ~ 1 | id))


# model as a mixed effect model with kingdom as moderator
model_metafor_kingdom <- rma.mv(
  yi = cohens_d,
  V = cohens_d_var,
  mods = ~ kingdom,
  data = dat_cohens_d,
  random = list( ~  1 | study, 
                 ~ 1 | id))

# model as a mixed effect model with past scale as moderator
model_metafor_scale <- rma.mv(
  yi = cohens_d,
  V = cohens_d_var,
  mods = ~ past_scale_in_year,
  data = dat_cohens_d,
  random = list( ~  1 | study, 
                 ~ 1 | id))
