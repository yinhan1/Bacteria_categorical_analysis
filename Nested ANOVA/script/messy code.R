
library(tidyverse)
library(magrittr)
library(ggsci)

####  ---------------- step 0: clean data ----------------  ####

#### load data file ####
data = readxl::read_excel('data/Per Fly Data.xlsx')

#### extract columns needed ####
cols_to_take = c('Pop','Sex','Population Replicate','CFU/mg Firm','CFU/mg Prot')
df = data %>% 
  select(cols_to_take) %>% 
  set_colnames(c('population','sex','replicate','firm','prot')) %>% 
  mutate(sex = recode(sex, 'F' = 'Female', 'M' = 'Male'),
         replicate = factor(replicate))

#### check NA #### 
df %>% 
  select(everything()) %>% 
  summarise_all(funs(sum(is.na(.))))



####  ------------------- step 1: EDA -------------------  ####

df %>% 
  reshape2::melt(id.vars = c('population','sex','replicate')) %>% 
  na.omit() %>% 
  mutate(variable = recode(variable, 'firm' = 'Firmictes', 'prot' = 'Proteobacteria')) %>% 
  ggplot(aes(x = replicate, y = value, fill = sex)) +
  geom_boxplot() +
  facet_grid(population ~ variable) +
  labs(x = 'Population replicates', y = 'Counts per mg weight', fill = '') +
  theme_bw()



####  --------------- step 2.1: Model on Firmicutes cnt  ---------------  ####

library(lme4)
library(lmerTest)

#### models setup ####
df_temp = df %>% filter(!is.na(firm)) 
model1 = lmer(firm ~ population + (1|replicate), df_temp, REML = FALSE)
model2 = lmer(firm ~ sex + (1|replicate), df_temp, REML = FALSE)
model3 = lmer(firm ~ population + sex + (1|replicate), df_temp, REML = FALSE)
model4 = lmer(firm ~ population + sex + population:sex + (1|replicate), df_temp, REML = FALSE)
model5 = lmer(firm ~ population + sex + population:sex + (population|replicate), data = df_temp, REML = FALSE)

#### model selection #### 
anova(model1, model3)
anova(model2, model3)
anova(model3, model4)
anova(model4, model5)  
rand(model5)
# best model 5


#### assumptions check #### 
model = model5
{
  par(mfrow=c(2,2))
  plot(fitted(model),resid(model))
  qqnorm(resid(model))
  qqnorm(as.data.frame(ranef(model))[,4])
}


#### model summary ####
summary(model)


####  --------------- step 2.2: Model on Proteobacteria cnt  ---------------  ####

#### model selection #### 
df_temp = df %>% filter(!is.na(prot)) 
model1 = lmer(prot ~ population + (1|replicate), df_temp, REML = FALSE)
model2 = lmer(prot ~ sex + (1|replicate), df_temp, REML = FALSE)
model3 = lmer(prot ~ population + sex + (1|replicate), df_temp, REML = FALSE)
model4 = lmer(prot ~ population + sex + population:sex + (1|replicate), df_temp, REML = FALSE)
# overfitting when (population|replicate)


#### model selection #### 
anova(model1, model3)
anova(model2, model3)
anova(model3, model4)
rand(model3)
# best model 3

#### assumptions check #### 
model = model3
{
  par(mfrow=c(2,2))
  plot(fitted(model),resid(model))
  qqnorm(resid(model))
  qqnorm(as.data.frame(ranef(model))[,4])
}


#### model summary ####
summary(model)



####  --------------- step 3: Model on ratios of cnt  ---------------  ####

#### model selection #### 
df_temp = df %>% na.omit() %>% mutate(ratio = firm / prot)
model1 = lmer(ratio ~ population + (1|replicate), df_temp, REML = FALSE)
model2 = lmer(ratio ~ sex + (1|replicate), df_temp, REML = FALSE)
model3 = lmer(ratio ~ population + sex + (1|replicate), df_temp, REML = FALSE)
model4 = lmer(ratio ~ population + sex + population:sex + (1|replicate), df_temp, REML = FALSE)

model3.0 = lm(ratio ~ population + sex, df_temp)
model4.0 = lm(ratio ~ population + sex + population:sex, df_temp)
                
#### model selection #### 
anova(model1, model3)
anova(model2, model3)
anova(model3, model4)
rand(model3)
anova(model3, model3.0)
anova(model3.0, model4.0)


#### assumptions check #### 
model = model3.0
{
  par(mfrow=c(2,2))
  plot(model)
}

#### model summary ####
summary(model)
anova(model)

library(boot)

boot.fn = function(data, index){
  return( coef(glm(default ~ income + balance, data = Default, family = binomial, subset = index )) )
}

set.seed(1)
boot(data = Default, boot.fn, 1000)
set.seed(1)
boot(data = Auto, boot.fn, 1000)




