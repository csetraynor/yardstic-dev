#my_trials
#by Carlos S Traynor

library(yardstick)
library(dplyr)

head(two_class_example)

metrics(two_class_example, truth, predicted)

two_class_example %>% roc_auc(truth, Class1)

lvl <- levels(two_class_example$truth)

two_class_example %>% mnLogLoss(truth, !! lvl)

#### For Simulation and time-dependent covariates
library(PermAlgo)

#### Move on with the classic example

library(yardstick)
library(survival)
library(dplyr)

mc_samp <- readRDS("C:/RFactory/data_test.RDS")
mc_samp <- readRDS("/home/mtr/rfactory/data_test.RDS")

survdata <- readRDS("my_musings/mc_samp.RDS")
survdata$splits <- mc_samp$splits

survtest <- as.data.frame(survdata$mod_lasso$`1`)
survsplit <- as.data.frame(survdata$splits$`1`$data)[survdata$splits$`1`$in_id,]

surv_example <- survsplit %>%
  select(time = os_months, status = os_deceased)

surv_example <- tibble::as.tibble(surv_example)

surv_example$risk <-  yardstick::pred_lp(obj = survtest,
                                       test_data = survsplit %>%
                      dplyr::select(-patient_id) ) %>% unlist

surv_example <- dplyr::arrange(surv_example, time)

surv_example$bh <-  yardstick::base_haz(obj = survtest,
                                         test_data = survsplit %>%
                                          dplyr::select(-patient_id) )


##### Example of ROC

surv_example %>% roc_surv(risk)
surv_example %>% roc_surv(risk, integrated = TRUE)



##Put covariates in place
surv_example <- tibble::rownames_to_column(surv_example, var = "id")

covariates <- survsplit %>% select(age_std, npi, contains("feature"))
surv_example$covariates <- stick_covar(survdata = surv_example, covars = covariates)










