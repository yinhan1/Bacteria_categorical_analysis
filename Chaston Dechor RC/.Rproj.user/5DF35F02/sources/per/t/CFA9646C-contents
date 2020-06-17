# Last Edited: 06/17/2020

library(tidyverse)
library(magrittr)
library(ggsci)

data <- 
  readxl::read_excel('Chaston Dechor RC.xlsx') %>% 
  select(
    'Total CFU Proteobacteria',
    'Total CFU Firmicutes',
    'Total CFU',
    'Population',
    'replicate',
    'sex',
    'genotype'
  ) %>% 
  set_colnames(
    c(
      'Proteobacteria',
      'Firmicutes',
      'Total',
      'Population',
      'Replicate',
      'Sex',
      'Genotype'
    )
  ) %>% 
  mutate(
    Sex = recode(Sex,
                 'f' = 'Female',
                 'm' = 'Male')
  )

data_melt <- 
  data %>% 
  reshape2::melt(
    variable.names = c(
      'Proteobacteria',
      'Firmicutes',
      'Total'
    ),
    id.vars = c(
      'Population',
      'Replicate',
      'Sex',
      'Genotype'
    ),
    variable.name = 'Bacteria',
    value.name = 'Count'
  ) %>% 
  mutate_at(c('Population','Replicate','Sex'), as.factor)


#### ----------------------   Step 0: EDA   ---------------------- ####

set_fills <- c(
  'ACO Female' = "#FF7F00",
  'CO Female' = "#FF7F00",
  'ACO Male' = "#1874CD",
  'CO Male' = "#1874CD"
)

#### by replicates #### 

data_melt %>% 
  # filter(Bacteria == "Total") %>% 
  group_by(Bacteria) %>% 
  ggplot(aes(x = as.factor(Replicate),
             y = Count,
             fill = Sex)) +
  geom_boxplot(alpha = 0.5) +
  geom_dotplot(dotsize = 1,
               binaxis = "y", 
               stackdir = "center", 
               position = "dodge") + 
  scale_fill_jco() +
  facet_grid(Sex~Population) +
  theme_bw() +
  labs(x = "Replicate",
       fill = "")


#### replicates combined ####

data_melt %>%  
  filter(Bacteria == "Total") %>% 
  group_by(Bacteria) %>% 
  ggplot(aes(x = Population, 
             y = Count, 
             fill = Sex)) +
  geom_boxplot(alpha = 0.7) +
  geom_dotplot(
    dotsize = 1,
    binaxis = "y", 
    stackdir = "center", 
    position = "dodge"
  ) + 
  scale_fill_jco() +
  facet_wrap(~Sex, ncol = 1) +
  theme_bw() +
  theme(legend.position="top")


#### -----------------   Step 0.2: Stats Tables  ----------------- ####

data_melt %>% 
  filter(Bacteria == "Total") %>% 
  group_by(Population, Sex) %>% 
  summarise(n = length(Count),
            Mean = mean(Count),
            Median = median(Count),
            SD = sd(Count)) %>% 
  arrange(Sex, Population)


#### ------------   Step 1: Possion Regression Selection  ------------ ####

library(MASS)

#### Proteobacteria (model 3) #### 

model_prot <- glm(Count ~ Population + Sex, 
                  family = poisson(link = "log"), 
                  data = data_melt %>% filter(Bacteria == "Proteobacteria"))
model2_prot <- glm(Count ~ Population + Sex + Population:Sex, 
                   family = poisson(link = "log"),
                   data = data_melt %>% filter(Bacteria == "Proteobacteria"))
model3_prot <- glm.nb(Count ~ Population + Sex,
                      data = data_melt %>% filter(Bacteria == "Proteobacteria"))
model4_prot <- glm.nb(Count ~ Population + Sex + Population:Sex, 
                      data = data_melt %>% filter(Bacteria == "Proteobacteria"))


summary(model_prot)
summary(model2_prot)
summary(model3_prot)
summary(model4_prot)

AIC(model_prot)
AIC(model2_prot)
AIC(model3_prot)
AIC(model4_prot)

BIC(model_prot)
BIC(model2_prot)
BIC(model3_prot)
BIC(model4_prot)


#### Firmicutes (model 3) #### 

model_firm <- glm(Count ~ Population + Sex, 
                  family = poisson(link = "log"), 
                  data = data_melt %>% filter(Bacteria == "Firmicutes"))
model2_firm <- glm(Count ~ Population + Sex + Population:Sex, 
                   family = poisson(link = "log"),
                   data = data_melt %>% filter(Bacteria == "Firmicutes"))
model3_firm <- glm.nb(Count ~ Population + Sex,
                      data = data_melt %>% filter(Bacteria == "Firmicutes"))
model4_firm <- glm.nb(Count ~ Population + Sex + Population:Sex, 
                      data = data_melt %>% filter(Bacteria == "Firmicutes"))

summary(model_firm)
summary(model2_firm)
summary(model3_firm)
summary(model4_firm)

AIC(model_firm)
AIC(model2_firm)
AIC(model3_firm)
AIC(model4_firm)

BIC(model_firm)
BIC(model2_firm)
BIC(model3_firm)
BIC(model4_firm)


#### Total (model 3) #### 
model_total <- glm(Count ~ Population + Sex, 
                  family = poisson(link = "log"), 
                  data = data_melt %>% filter(Bacteria == "Total"))
model2_total <- glm(Count ~ Population + Sex + Population:Sex, 
                   family = poisson(link = "log"),
                   data = data_melt %>% filter(Bacteria == "Total"))
model3_total <- glm.nb(Count ~ Population + Sex,
                      data = data_melt %>% filter(Bacteria == "Total"))
model4_total <- glm.nb(Count ~ Population + Sex + Population:Sex, 
                      data = data_melt %>% filter(Bacteria == "Total"))


summary(model_total)
summary(model2_total)
summary(model3_total)
summary(model4_total)

AIC(model_total)
AIC(model2_total)
AIC(model3_total)
AIC(model4_total)

BIC(model_total)
BIC(model2_total)
BIC(model3_total)
BIC(model4_total)


#### ------------   Step 2: Possion Regression Conclusion  ------------ ####

#### functions ####
new_data <- 
  data_melt %>% 
  dplyr::select(Population, Sex) %>% 
  unique()

get_CI_tb <- function(model){
  est <- predict(model, newdata = new_data, type = "link", se.fit=TRUE)
  est2 <- within(est, {
    estCount <- exp(fit)
    LL <- exp(fit - 1.96 * se.fit)
    UL <- exp(fit + 1.96 * se.fit)
  })
  cbind(new_data, data.frame(est2)) %>% return()
}

plot_CI_tb <- function(model, bacteria){
  ggplot(data = NULL) +
    geom_errorbar(data = get_CI_tb(model),
                  aes(x = Population, 
                      y = estCount, 
                      ymin = LL, 
                      ymax = UL, 
                      color = Sex),
                  width = 0.1,
                  size = 1) +
    geom_point(data = get_CI_tb(model),
               aes(x = Population, 
                   y = estCount, 
                   color = Sex),
               size = 2) +
    geom_dotplot(data = data_melt %>% filter(Bacteria == bacteria),
                 aes(x = Population, 
                     y = Count, 
                     fill = Sex),
                 dotsize = 0.6,
                 binaxis = "y", 
                 stackdir = "center", 
                 position = "dodge",
                 alpha = 0.08) +
    scale_color_jco() +
    scale_fill_jco() +
    facet_wrap(~Sex, ncol = 1) +
    theme_bw() +
    theme(legend.position="top") +
    labs(y = "Estimates")
}

#### conclusions #### 

plot_CI_tb(model3_prot, "Proteobacteria")
plot_CI_tb(model3_firm, "Firmicutes")
plot_CI_tb(model3_total, "Total")

summary(model3_prot)
summary(model3_firm)
summary(model3_total)








