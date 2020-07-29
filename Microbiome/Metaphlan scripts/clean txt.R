### Last Edited: 07/29/20

library(tidyverse)
library(magrittr)

data = read.delim("mpa_v20/desiccation_full_abundance.txt")


## filter out rows with 'g__' but no 's__'
tb = data %>% 
  filter(ID %>% str_detect('g__') & ID %>% str_detect('s__', negate = TRUE))


## split long labels 
tb_clean = str_split(as.character(tb$ID), "\\|", simplify = TRUE) %>% 
  apply(., 2, function(col) str_remove(col, '[kpcofg]__')) %>% 
  set_colnames(c('Kingdom','Phylum','Class','Order','Family','Genus')) %>% 
  cbind(tb %>% select(-ID))


View(tb_clean)


# write.csv(tb_clean, 'desiccation.csv')
