library (revtools)
library(here)
library(tidyverse)


# import data from working directory
dat_list <- read_bibliography(list.files(here("bibliography"), 
                                         full.names = TRUE, 
                                         pattern = "*.bib"))


# remove duplicates -------------------------------------------------------


# find duplicated DOIs within the dataset
doi_match <- find_duplicates(dat_list,
                             match_variable = "doi",
                             group_variables = NULL,
                             match_function = "exact")

# automatically extract one row per duplicate
dat_unique <- extract_unique_references(dat_list, doi_match)

# number of entries left
nrow(dat_unique) 

# additionally check titles of the remaining
# an alternative is to try fuzzy title matching

title_match <- find_duplicates(dat_unique,
                               match_variable = "title",
                               group_variables = NULL,
                               match_function = "stringdist",
                               method = "osa",
                               threshold = 5)

# extract duplicates
dat_unique <- extract_unique_references(dat_unique, title_match)

# number of entries left
nrow(dat_unique) 

# percentage of duplicates
(nrow(dat_list) - nrow(dat_unique))/nrow(dat_list) 



# manual screening --------------------------------------------------------


# manually screen titles
dat_titles <- screen_titles(dat_unique)

# save in between 
dat_titles %>% 
  write_csv(here("bibliography",
                 "cleaned",
                 "cleaned_by_title.csv"))


# extract topics
dat_topics <- dat_unique %>% 
  slice_head(n = 25) %>% 
  make_dtm(.$title) %>% 
  run_topic_model(., "lda", 
                  n_topics = 5, 100)

data_unique %>% 
  as_tibble() %>% 
  ggplot(aes(year)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45))
