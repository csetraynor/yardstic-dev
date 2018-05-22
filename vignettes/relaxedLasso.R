library(glmnet)
library(dplyr)
enet <- readRDS("C:/RFactory/predsurv/performance_results/enet_model_iclust2_filter.RDS")
brca <- readRDS("/home/mtr/rfactory/brca_data.RDS")
#or
brca <- readRDS("C:/RFactory/Rdata_brca/brca_data.RDS")


iclust2 <- brca[brca$intclust == 2, ]

plot(enet)
optimal <- as.matrix(coef(enet, s = "lambda.min"))
optimal <- as.data.frame(optimal)
colnames(optimal) <- "mod"
optimal$Hugo_Symbol <- rownames(optimal)
optimal <- optimal %>% filter(mod!=0)
#### Perform Relaxed Lasso
require(survival)
clinical <- coxph(Surv(os_months, os_deceased) ~ age_std + npi , data = iclust2 )

X <- iclust2[,c("age_std", "npi",optimal$Hugo_Symbol)]

genomic <- coxph(Surv(time = iclust2 %>%
                        dplyr::select(os_months) %>%
                        unlist,
                      event = iclust2 %>%
                        dplyr::select(os_deceased) %>%
                        unlist)~.,
                 data = iclust2[,optimal$Hugo_Symbol])

clinico_genomic <- coxph(Surv(iclust2$os_months, iclust2$os_deceased) ~ . , data = iclust2[,c("age_std", "npi",optimal$Hugo_Symbol)])

#

library(rsample)
set.seed(9666)
mc_samp <- mc_cv(iclust2, strata = "os_deceased", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$os_deceased == 1)
summary(map_dbl(mc_samp$splits, cens_rate))


paste0(optimal$Hugo_Symbol, collapse = "+")

clinical <- as.formula(Surv(os_months, os_deceased) ~ age_std+npi)
genomic <- as.formula(Surv(os_months, os_deceased) ~ AGA+C8orf80+CB242622+C9orf95+MAP1B+NDUFA4L2+NGF+AA831838+CAMKK1+ZFPM1)
clinico_genomic <- as.formula(Surv(os_months, os_deceased) ~ AGA+C8orf80+CB242622+C9orf95+MAP1B+NDUFA4L2+NGF+AA831838+CAMKK1+ZFPM1+age_std+npi)


#The model fitting function will take the formula as an argument:

mod_fit <- function(x, form, ...) {
  coxph(form, data = analysis(x), ...)
}

mc_samp$mod_clinical <- map(mc_samp$splits, mod_fit, form = clinical)


mc_samp$clinical_brier <- map2(mc_samp$splits, mc_samp$mod_clinical, tdbrier.model.list)

mc_samp$clinical_ibrier_model <- map_dbl(mc_samp$clinical_brier, tdbrier.int.matrix)

mc_samp$clinical_tdROC <- map2(mc_samp$splits, mc_samp$mod_clinical, tdroc.model.list)

mc_samp$clinical_ci <- map2_dbl(mc_samp$splits, mc_samp$mod_clinical, get_cindex)

mc_samp$clinical_dev <- map_dbl(mc_samp$mod_clinical, get_deviance)

