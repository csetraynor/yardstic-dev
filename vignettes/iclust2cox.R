###### Load data

brca <- readRDS("/home/mtr/rfactory/brca_data.RDS")
brca <- readRDS("C:/RFactory/Rdata_brca/brca_data.RDS")
cvfit <- readRDS("C:/RFactory/bymetabric_files/ridge_model_brca_pooled.RDS")

iclust2 <- brca[brca$intclust == 2, ]
#rm(brca)


library(glmnet)
library(dplyr)
optimal <- as.matrix(coef(cvfit, s = "lambda.min"))
optimal <- as.data.frame(optimal)
colnames(optimal) <- "mod"
optimal$Hugo_Symbol <- rownames(optimal)
optimal <- optimal %>% filter(mod!=0)

iclust2 <- iclust2[,c( "os_months", "os_deceased" , optimal$Hugo_Symbol)]
assertthat::assert_that(all.equal.character(colnames(iclust2)[3:24372], optimal$Hugo_Symbol))


library(glmnet)
enet <- readRDS("C:/RFactory/predsurv/performance_results/enet_model_iclust2_filter.RDS")

plot(enet)

options(expressions = 5e5)
memory.limit(5e10)

#### Perform Cox model
require(survival)
clinical <- coxph(Surv(os_months, os_deceased) ~ age_std+npi , data = iclust2, init = optimal$mod[match(c("age_std", "npi"), optimal$Hugo_Symbol)], iter = 0 )

init = optimal$mod[-match(c("age_std", "npi"), optimal$Hugo_Symbol)]
X <- iclust2[,optimal$Hugo_Symbol[-match(c("age_std", "npi"),optimal$Hugo_Symbol) ]]

library(dplyr)

genomic <- coxph(Surv(time = iclust2 %>%
  dplyr::select(os_months) %>%
  unlist,
event = iclust2 %>%
  dplyr::select(os_deceased) %>%
  unlist)~.,
init = inits, iter = 0,
data = X)

clinico_genomic_fit <- coxph(Surv(iclust2$os_months, iclust2$os_deceased) ~ . , data = iclust2[,c(optimal$Hugo_Symbol)], init = optimal$mod, iter = 0)

### Resample

#To resample these data, it would be a good idea to try to maintain the same censoring rate across the splits. To do this, stratified resampling can be used where each analysis/assessment split is conducted within each value of the status indicator. To demonstrate, Monte Carlo resampling is used where 75% of the data are in the analysis set. A total of 100 splits are created.

library(rsample)
set.seed(9666)
mc_samp <- mc_cv(iclust2, strata = "os_deceased", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$os_deceased == 1)
summary(map_dbl(mc_samp$splits, cens_rate))

#To demonstrate the use of resampling with censored data, the parametric model shown above will be fit with different variable sets to characterize how important each predictor is to the outcome. To do this, a set of formulas are created for the different variable sets:

paste0(optimal$Hugo_Symbol, collapse = "+")

clinical <- as.formula(Surv(os_months, os_deceased) ~ age_std+npi)
genomic <- as.formula(Surv(os_months, os_deceased) ~ AGA+C8orf80+CB242622+C9orf95+MAP1B+NDUFA4L2+NGF+AA831838+CAMKK1+ZFPM1)
clinico_genomic <- as.formula(Surv(os_months, os_deceased) ~ AGA+C8orf80+CB242622+C9orf95+MAP1B+NDUFA4L2+NGF+AA831838+CAMKK1+ZFPM1+age_std+npi)

#The model fitting function will take the formula as an argument:

mod_fit <- function(x, form, ...)
  coxph(form, data = analysis(x), ...)

#To calculate the efficacy of the model, the concordance statistic is used (see ?survConcordance):

get_concord <- function(split, mod, ...) {
  pred_dat <- assessment(split)
  pred_dat$pred <- predict(mod, newdata = pred_dat)
  survConcordance(Surv(os_months, os_deceased) ~ pred, pred_dat, ...)$concordance
}


#With these functions, a series of models are created for each variable set.

mc_samp$mod_clinical    <- map(mc_samp$splits, mod_fit, form = clinical)
mc_samp$mod_genomic <- map(mc_samp$splits, mod_fit, form = genomic )
mc_samp$mod_clincogenomic <- map(mc_samp$splits, mod_fit, form = clinico_genomic)

# Similarly, the concordance values are computed for each model:

mc_samp$clinical_ci <- map2_dbl(mc_samp$splits, mc_samp$mod_clinical, get_concord)
mc_samp$genomic_ci     <- map2_dbl(mc_samp$splits, mc_samp$mod_genomic, get_concord)
mc_samp$clincogenomic_ci     <- map2_dbl(mc_samp$splits, mc_samp$mod_clincogenomic, get_concord)



library(dplyr)
concord_est <- mc_samp %>%
  select(-matches("^mod"))

library(tidyr)
library(ggplot2)
concord_est %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

library(tidyposterior)
concord_est <- perf_mod(concord_est, seed = 6507, iter = 5000, transform = logit_trans)
ggplot(tidy(concord_est)) +
  theme_bw()
comparisons <- contrast_models(
  concord_est,
  list_1 = rep("clincogenomic_ci", 2),
  list_2 = c("clinical_ci", "genomic_ci"),
  seed = 4654
)
ggplot(comparisons, size = 0.05) +
  theme_bw()

summary(comparisons, size = 0.05) %>%
  select(contrast, starts_with("pract"))


######## Now consider hallmarks

library(qusage)

hallmark_all <- qusage::read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/h.all.v6.1.symbols.gmt")
c6 <- qusage::read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/c6.all.v6.1.symbols.gmt")

onco_signature <- ontology_search(c6, optimal, Predictor = Hugo_Symbol, coef = mod)
ontology_search(hallmark_all, optimal, Predictor = Hugo_Symbol, coef = mod)

capture.output(print(onco_signature), file = "onco_signature.txt")



names(onco_signature)

paste0(optimal$Hugo_Symbol, collapse = "+")

lapply(onco_signature, function(x) as.character(x$features))

cox_onco_signature <- list(rep(NA, length(onco_signature)))
for(i in seq_along(onco_signature)) {
  coef <- optimal$Hugo_Symbol[!(optimal$Hugo_Symbol %in% as.character(onco_signature[[i]]$features))]

  survform <- as.formula(paste0("Surv(os_months, os_deceased)",
                                "~", paste0(coef, collapse = "+"), "+age_std+npi"))

  cox_onco_signature[[i]] <- coxph(survform, data = iclust2)
}

names(cox_onco_signature) <- names(onco_signature)


