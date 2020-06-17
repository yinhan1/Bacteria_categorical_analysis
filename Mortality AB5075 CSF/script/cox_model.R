
# First edited: 04/22/2020

library(tidyverse)
library(magrittr)
library(ggsci)
library(dplyr)

#### Preparing data ------------------------------- ####
#### load and clean data file ####
mortality_data <- 
  readxl::read_excel('Mortality Data - Ramirez Collaboration.xlsx') %>% 
  select(-Sex) %>% 
  filter(!(Treatment %>% str_detect('A42|PF'))) %>%
  mutate_at(vars(Treatment, Replicate), funs(factor)) %>% 
  na.omit() %>% 
  # rename("Remaining"= "Flies Remaining") %>% 
  select(-`Survival %`) %>% 
  droplevels() %>% 
  mutate(Treatment = fct_relevel(Treatment, "PBS"))

names(mortality_data)[5] = "Remaining"

#### create data frame for Day 0 ####
day_zero <- 
  mortality_data %>% 
  filter(Day == 1) %>% 
  mutate(Remaining = Remaining + Mortality, 
         Mortality = 0,
         Day = 0)

#### combine data sets for day 0 to the end ####
mortality_data <- 
  rbind(day_zero, mortality_data) %>% 
  group_by(Treatment, Replicate) %>% 
  mutate(Initial = max(Remaining)) 

#### convert mortality to death status ####
survival_data <- 
  mortality_data[rep(1:nrow(mortality_data), mortality_data$Mortality),] %>% 
  mutate(Death = 1) %>% 
  select(-c(Mortality, Initial))
  
#### create data frame for alive flies #### 
day_10 <-
  mortality_data %>%
  filter(Day == 10) %>%
  mutate(Death = 0)

day_10 <-
  day_10[rep(1:nrow(day_10), day_10$Remaining),] %>%
  select(-c(Mortality, Remaining, Initial))

#### combine dead flies and survival flies for survival_data
survival_data <- rbind(survival_data %>% select(-Remaining), day_10) %>% as_tibble()

#### Raw data plot: survival percent ------------------ ####

#### overall ####
mortality_data %>% 
  group_by(Treatment, Day) %>% 
  summarise(Mortality = sum(Mortality),
            Initial = sum(Initial)) %>% 
  group_by(Treatment) %>% 
  mutate(Mortality = cumsum(Mortality)) %>% 
  ggplot(aes(x = Day, y = 100*(1-(Mortality/Initial)), color = Treatment)) +
  geom_line() +
  geom_point() +
  scale_color_jco() +
  ylim(85,100) +
  labs(y = 'Survival percent (%)') + 
  theme_minimal()

# #### replicate 1 ####
# mortality_data %>%
#   filter(Replicate == 1) %>%
#   group_by(Treatment, Day) %>%
#   summarise(Mortality = sum(Mortality),
#             Initial = sum(Initial)) %>%
#   group_by(Treatment) %>%
#   mutate(Mortality = cumsum(Mortality)) %>%
#   ggplot(aes(x = Day, y = 100*(1-(Mortality/Initial)), color = Treatment)) +
#   geom_line() +
#   geom_point() +
#   scale_color_jco() +
#   labs(y = 'Survival percent (%)') +
#   ylim(80, 100) +
#   theme_minimal()
# 
# #### replicate 2 ####
# mortality_data %>%
#   filter(Replicate == 2) %>%
#   group_by(Treatment, Day) %>%
#   summarise(Mortality = sum(Mortality),
#             Initial = sum(Initial)) %>%
#   group_by(Treatment) %>%
#   mutate(Mortality = cumsum(Mortality)) %>%
#   ggplot(aes(x = Day, y = 100*(1-(Mortality/Initial)), color = Treatment)) +
#   geom_line() +
#   geom_point() +
#   scale_color_jco() +
#   labs(y = 'Survival percent (%)') +
#   ylim(80, 100) +
#   theme_minimal()
# 
# #### replicate 3 ####
# mortality_data %>%
#   filter(Replicate == 3) %>%
#   group_by(Treatment, Day) %>%
#   summarise(Mortality = sum(Mortality),
#             Initial = sum(Initial)) %>%
#   group_by(Treatment) %>%
#   mutate(Mortality = cumsum(Mortality)) %>%
#   ggplot(aes(x = Day, y = 100*(1-(Mortality/Initial)), color = Treatment)) +
#   geom_line() +
#   geom_point() +
#   scale_color_jco() +
#   labs(y = 'Survival percent (%)') +
#   ylim(80, 100) +
#   theme_minimal()


#### Test on replicates ----------------------- ####

library(survival)
library(survminer)

fit <- survfit(Surv(Day, Death) ~ Replicate, data = survival_data)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)
pairwise_survdiff(Surv(Day, Death == 1) ~ Replicate, 
                  p.adjust.method = "holm", data = survival_data)


#### Test proportional hazard model ------------------ ####

model = coxph(Surv(Day, Death) ~ Treatment, data = survival_data)
cox.zph(model)
cox.zph(model, transform=log)

summary(model)


#### Power analysis ------------------ ####

day_initial = day_zero %>% select(Treatment,Remaining) %>% 
  group_by(Treatment) %>% summarise(Initial = sum(Remaining))

day_end =  
  mortality_data %>% 
  filter(Day == 10) %>% 
  select(Treatment,Remaining) %>% 
  group_by(Treatment) %>% 
  summarise(Remaining = sum(Remaining))

tb = full_join(day_initial, day_end, by = "Treatment") %>% 
  mutate(Dead = Initial - Remaining) 
  
tb$pi = tb$Initial / sum(tb$Initial)
tb$ratio = tb$Dead/tb$Initial


top = (qnorm(0.975)+qnorm(0.9))^2
hr = seq(1.3,3,length.out = 20) 
hr = 2.9418


dummy = data.frame()
for(i in 1:3){
  bottom = (log(hr)^2)*tb$pi[4]*tb$pi[i]
  dummy_i = data.frame(
    tag = tb$Treatment[i], 
    hazard = hr, 
    death = top/bottom) %>% 
    mutate(size = death/(tb$ratio[4]+tb$ratio[i]))
  dummy = bind_rows(dummy,dummy_i)
}
dummy$tag = factor(dummy$tag, levels=levels(tb$Treatment))

dummy %>% 
  ggplot(aes(x=size, y=hazard, color=tag)) +
  geom_line(alpha=0.6, size=1, position=position_jitter(w=0.02, h=0.01)) +
  geom_point(size=1) +
  geom_vline(xintercept = 34464.690, linetype="dashed", size=0.5, color="grey60") +
  annotate("text", x=2500, y = 1.4, label="1,000", size=3, color="grey60") +
  scale_colour_jama(drop=FALSE) +
  theme_bw() +
  labs(x="Sample Size", y="Hazard Ratio", color="")






packageVersion("rlang")


