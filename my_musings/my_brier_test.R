lasso <- readRDS("/home/mtr/rfactory/yardstick/my_musings/mod_lasso2.RDS")

# library(predsurv)
mc_samp$mod_lasso <- lasso

mc_samp$brier_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                       function(splits, mod){fun_test3(
                                         obj = mod,
                                         test_data = splits,
                                         fit = "Lasso",
                                         pred = "Brier",
                                         integrated = TRUE,
                                         
                                       )
                                       })
