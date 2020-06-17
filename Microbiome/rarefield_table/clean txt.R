
library(tidyverse)
library(magrittr)

### step 1: clean from txt to xlsx    ----------------------------------
data <- 
  read.delim("~/Dropbox (CSU Fullerton)/Research Project/Flies/Robert -- statistical analysis/Microbiome/rarefield_table/rarefied_table_Robert.txt", 
             header = TRUE, comment.char = "#")

cols_to_clean <- names(data)[str_detect(names(data), "[_|.]", negate=TRUE)]

# data %>% 
#   select(cols_to_clean) %>% 
#   apply(., 2, function(x) str_replace(x, "k__|p__|c__|o__|f__|g__|s__","")) %>% 
#   apply(., 2, function(x) ifelse(x=="0"|x==" ",NA, x)) %>% 
#   as.data.frame() %>% 
#   bind_cols(
#     data %>% select(-cols_to_clean)
#   ) %>% 
#   write.csv(
#     file = "~/Dropbox (CSU Fullerton)/Research Project/Flies/Robert -- statistical analysis/Microbiome/rarefield_table/clean_table.csv")


### step 2: collapse table     ------------------------------------------

data <- readxl::read_excel("clean_table.xlsx", sheet = 1)
data_species <- data %>% filter(!is.na(species)) %>% select(-c(cols_to_clean[c(-1,-7)])) %>% rename(others = 2)

data_collapsed <- 
  data %>% 
  filter(is.na(species)) %>% 
  select(-c(cols_to_clean[c(-1,-2)], OTU.ID)) %>% 
  group_by(kingdom, phylum) %>% 
  summarise_all(funs(sum)) %>% 
  rename(others = 2)

bind_rows(data_collapsed, data_species) %>% write.csv(file = "dummy.csv", row.names = FALSE)













