

mc_samp <- readRDS("/home/mtr/rfactory/data_test.RDS")


splits <- readRDS("C:/RFactory/data_test.RDS")


mc_samp$brier_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                       function(splits, mod){yardstick::fun_test3(
                                         obj = mod,
                                         test_data = splits,
                                         fit = "Lasso",
                                         pred = "Brier",
                                         subject = patient_id
                                       )
                                       })


mc_samp$brier_reference <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                       function(splits, mod){yardstick::fun_test3(
                                         obj = mod,
                                         test_data = splits,
                                         fit = "Lasso",
                                         pred = "Brier",
                                         subject = patient_id,
                                         reference = TRUE
                                       )
                                       })

mc_samp$splits <- NULL
saveRDS(mc_samp, "mc_samp.RDS")

library(dplyr)
library(tidyr)
library(ggplot2)
##Plot Brier

brier_dens <- mc_samp %>%
   dplyr::select(-splits, - mod_lasso) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "bottom")
brier_dens <- brier_dens +
  labs(x = "Integrated Brier Score",
       title = "Density of iBrier") +
  geom_vline(xintercept =  0.25 , linetype = "dotted" )

library(tidyposterior)
mc_samp_brier <- tidyposterior::perf_mod(mc_samp %>%
                                           dplyr::select(-splits, - mod_lasso),
                                         seed = 6507, iter = 5000, transform = logit_trans,
                                         hetero_var = FALSE)

mbri_tab <- summary(tidy(mc_samp_brier))
mbri_tab <- as.data.frame(mbri_tab)

star = stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,digits = 3, digit.separator = ".",
                 title = "Bayesian analysis of resampling AUC")


stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          title = "Bayesian analysis of resampling AUC")


posterior_brier <- ggplot(tidy(mc_samp_brier)) +
  theme_bw()+
  labs(
    title = "Posterior probability for integrated Brier Score")
posterior_brier <- posterior_brier +   labs(
  title = "Posterior probability of iBrier")
posterior_brier

comparisons_brier <- contrast_models(
  mc_samp_brier,
  list_1 = rep("brier_reference", 1),
  list_2 = "brier_lasso",
  seed = 4654
)

compare_brier <- ggplot(comparisons_brier, size = 0.05) +
  theme_bw()+
  labs(
    title = "Posterior probability of iBrier.",
    subtitle ="Benchmark: Reference")

compare_brier <- compare_brier +   labs(
  title = "Posterior probability for iBrier",
  subtitle ="Benchmark: Reference")
compare_brier


