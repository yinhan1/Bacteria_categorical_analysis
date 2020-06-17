
library(tidyverse)
library(magrittr)
library(ggsci)

data <- read.delim("rarefied_table_Robert.txt", header = TRUE, comment.char = "#")
colNames <- read.delim("Rarefied table colNames.txt", header = TRUE, comment.char = "#")
names(data) <- names(colNames)


#### -------------------   step 1: collapse others   ------------------- ####

### separate identified OUT.ID and others
rows_identified <- grepl('^[^0-9]+$',data$OTU.ID)
cols_to_take <- c(TRUE, !grepl('^[^0-9]+$', names(data)[-1]))
tb_identified <- data[rows_identified, cols_to_take]
tb_others <- data[!rows_identified, cols_to_take]

### group repeat names
tb_identified <- 
  tb_identified %>% 
  mutate(OTU.ID = str_to_title(OTU.ID)) %>% 
  group_by(OTU.ID) %>% 
  summarise_each(sum)


### collapse others and recombine
tb_full <- rbind(
  data.frame(OTU.ID = "Others", t(colSums(tb_others[,-1]))),
  tb_identified
) %>% 
  set_rownames(NULL)



#### -------------------   step 2: average columns   ------------------- ####

names(tb_full)[-1] <- str_replace(names(tb_identified[,-1]), "\\.[0-9]", "")
tb_average <- 
  data.frame(
    OTU.ID = tb_full[,1],
    {t(rowsum(t(tb_full[,-1]), names(tb_full)[-1])/c(table(names(tb_full)[-1])))}
)
levels(tb_average$OTU.ID) <- c(
  levels(tb_average$OTU.ID)[-1], "Others"
)

View(tb_average)
  

#### ------------------------   step 3: plot   ------------------------ ####

tb_average %>%
  reshape2::melt(id.vars = "OTU.ID") %>% 
  ggplot(aes(x = variable, y = value, fill = OTU.ID)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_jco() +
  theme_minimal() +
  labs(x = "", 
       y = "Relative abundance",
       fill = "")
  















