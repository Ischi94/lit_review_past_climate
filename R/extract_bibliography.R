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

# remove excluded entries
dat_unique <- dat_titles %>% 
  filter(screened_titles == "selected")

# manually screen abstracts
dat_abstracts <- screen_abstracts(dat_unique)

# save in between 
dat_abstracts %>% 
  write_csv(here("bibliography",
                 "cleaned",
                 "cleaned_by_abstract.csv"))


# remove excluded entries
dat_final <- dat_abstracts %>% 
  filter(screened_abstracts == "selected")

# save final data 
dat_final %>% 
  write_csv(here("bibliography",
                 "cleaned",
                 "cleaned_final.csv"))



# quick visualisation through time
dat_final %>% 
  ggplot(aes(year)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45))



# prepare final dataset  --------------------------------------------------

# set up spreadsheet to manually enter information from papers
dat_final %>% 
  select(title, author, doi) %>% 
  add_column(biotic_unit = NA_character_, 
             clade = NA_character_, 
             scale = NA_character_, 
             methodology = NA_character_,
             past_climate = NA_character_, 
             quantitative = NA_character_, 
             past_mechanism = NA_character_, 
             past_scale = NA_character_, 
             past_effect_size = NA_integer_, 
             past_effect_unit = NA_character_) %>% 
  write_csv2(here("data", 
                  "review_spreadsheet.csv"))
