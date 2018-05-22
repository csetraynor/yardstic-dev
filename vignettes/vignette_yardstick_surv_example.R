library(survival)
library(purrr)
library(dplyr)
str(lung)

############ Extending the rsample vignette of topepo https://topepo.github.io/rsample/articles/Applications/Survival_Analysis.html

# Author: csetraynor
data(lung)

lung_data <- lung %>%
  rename(time = time,
         status = status) %>%
  mutate(status = status == 2)

library(rsample)
set.seed(9666)
mc_samp <- mc_cv(lung_data, strata = "status", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))


################ Create formulas

three_fact <- as.formula(Surv(time, status) ~ ph.ecog + age + strata(sex))
rm_ph.ecog <- as.formula(Surv(time, status) ~           age + strata(sex))
rm_age     <- as.formula(Surv(time, status) ~ ph.ecog +       strata(sex))
rm_sex     <- as.formula(Surv(time, status) ~ ph.ecog + age              )

############### Fitting function
mod_fit <- function(x, form, ...) {
  coxph(form, data = x, ...)
}

############### Create models

mc_samp$mod_full    <- map(mc_samp$splits, mod_fit, form = three_fact)
mc_samp$mod_ph.ecog <- map(mc_samp$splits, mod_fit, form = rm_ph.ecog)
mc_samp$mod_age     <- map(mc_samp$splits, mod_fit, form = rm_age)
mc_samp$mod_sex     <- map(mc_samp$splits, mod_fit, form = rm_sex)

############### Get Brier
mc_samp$brier_full <- pmap(list(mc_samp$splits, mc_samp$mod_full),
                           function(data, model){
                             get_tdbrier(data = data,
                                         mod = model)
                           })
mc_samp$brier_ph.ecog <- pmap(list(mc_samp$splits, mc_samp$mod_ph.ecog),
                           function(data, model){
                             get_tdbrier(data = data,
                                         mod = model)
                           })
mc_samp$brier_age <- pmap(list(mc_samp$splits, mc_samp$mod_age),
                           function(data, model){
                             get_tdbrier(data = data,
                                         mod = model)
                           })
mc_samp$brier_sex <- pmap(list(mc_samp$splits, mc_samp$mod_sex),
                           function(data, model){
                             get_tdbrier(data = data,
                                         mod = model)
                           })
###Get integrate
mc_samp$ibrier_full <- map_dbl(mc_samp$brier_full, integrate.tdbrier)
mc_samp$ibrier_ph.ecog <- map_dbl(mc_samp$brier_ph.ecog, integrate.tdbrier)
mc_samp$ibrier_age <- map_dbl(mc_samp$brier_age, integrate.tdbrier)
mc_samp$ibrier_sex <- map_dbl(mc_samp$brier_sex, integrate.tdbrier)


mc_samp$brier_full$`12`$AppErr

mc_samp$mod_full$`12`
assessment(mc_samp$splits$`12`)
