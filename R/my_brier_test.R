#' Test model
#'
#' This function takes a survival fit and following the Cox model (or random forest), estimates the hazard for individual from the model coefficient, while the baseline hazard is estimated with the Nelson-Aalen method. Gives prediction error measures Brier, c index, loglik, auc
#'
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' pred : prediction error Brier, ROC or C-Index
#' adapted : in Lasso allows to swap to "adapated Lasso", otherwise will take usual Lasso \cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
#' @import survival
fun_test3 <- function(obj, test_data, fit = "Lasso", subject = subject, time = os_months, status = os_deceased, event_type = 1,  pred = "Brier",  all = FALSE, integrated = TRUE, reference = FALSE, noboot = 0, mc = FALSE, ...){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  subject <- dplyr::enquo(subject)


  #transform splits to train data
  test_data <- as.data.frame(test_data$data)[test_data$in_id,] %>%
    dplyr::select(-!!subject)


  #this is the most important bit, gets the fitted model and extracts the selected features
  if(fit == "Univariate" | fit == "Lasso" | fit == "Adaptive Lasso" | fit == "Ridge regression" | fit == "Elastic net" ){
    # find lambda for which dev.ratio is max
    selectedBeta <- rownames(obj)
  }
  #Now prepares for getting performance measures, random forest at the bottom
  if(fit != "Random forest"){
    ##### No variables selected complete "shrinkage" tipycal in Lasso
    if(length(selectedBeta) == 0 ){
      mod <-  survival::coxph(Surv(time, status)~1, iter = 0,
                              data = train_data %>%
                                dplyr::mutate(time = !!time, status = !!status))
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      ndata <- test_data %>% dplyr::mutate(time = !!time, status =  !!status)
    }
    #Create training and test set and formula
    if(length(selectedBeta) > 0 ){
      # Create train X take only covariates for which beta is not zero
      X <- test_data %>% dplyr::select(-!!time, -!!status)
      X   <- as.data.frame(X[,colnames(X) %in% selectedBeta])
      # prepare for coxph model

    mod <-  survival::coxph(survival::Surv(time = test_data %>%
                                       dplyr::select(!!time) %>%
                                       unlist,
                                     event = test_data %>%
                                       dplyr::select(!!status) %>%
                                       unlist)~.,
                                init = as.vector(unlist(obj)), iter = 0,
                                data = X)

    }
    #Create grid of equidistant time points for testing
    timepoints <-  seq(0, max(test_data %>%
                                dplyr::filter(!!status == event_type) %>%
                                dplyr::select(!!time) %>% unlist
    ) , length.out = 100L)
    ######################################################
    #### Prediction Error
    if(pred == "Brier" | all){
      if(length(selectedBeta) == 0){
        probs <- matrix(runif(nrow(test_data)*100), ncol = 100)
        brier <- pec::pec(probs, survival::Surv(time, status) ~ 1,
                          data = test_data %>%
                            dplyr::mutate( time   = !!time, status = !!status),
                          maxtime = max(timepoints),
                          exact = F,
                          exactness = 99L)
        out <- brier
        if(integrated){
          out <- yarstick::fun_ibrier_score_reference(out)
        }
      }else{
        #Calculate probs
        probs <- pec::predictSurvProb(mod,
                                      newdata = X,
                                      times = timepoints)

        #Calculate brier score
        brier <- pec::pec(probs, Surv(time, status) ~ 1,
                          data = test_data %>% dplyr::mutate( time   = !!time, status = !!status),
                          maxtime = max(timepoints),
                          exact = F,
                          exactness = 99L)
        out <- brier
        if(integrated){
          if(reference){
            out <- yardstick::fun_ibrier_score_reference(out)
          }else{
            out <- yardstick::fun_ibrier_score(out)
          }
        }
      }
    }
    if(pred == "ROC" | all ){
      if(length(selectedBeta) == 0){
        probs <- rep(0, nrow(testX))
        names(probs) <- rownames(testX)
      }else{
        if(fit == "Ridge regression"){
          probs <- as.matrix(testX, ncol = length(selectedBeta)) %*% as.vector(unlist(obj))
        }else{
          probs <- predict(mod, newdata = testX, type = "lp")
        }
      }
      roc <- tdROC::tdROC(X = probs,
                          Y = test_data %>% dplyr::select(time = !!time)%>% unlist,
                          delta = test_data %>% dplyr::select(status = !!status)%>% unlist,
                          tau = quantile(test_data %>% dplyr::select(time = !!time)%>% unlist, .73), nboot = noboot, alpha = 0.05, n.grid = 1000,  type = "uniform"
      )
      out <- roc
      if(integrated){
        out <- out$AUC[1] %>% unlist
      }
    }
    if(pred  == "Deviance" | all  ){
      logl <- -2*(mod$loglik)
      out <- logl
      if(integrated){
        out <- out[2]
      }
    }
    if(pred == "c_index" | all ){
      if(length(selectedBeta) == 0){
        out <- 0.5
      }else{
        ###Create your survival estimates
        ci <-  pec::cindex(mod, formula = Surv(time, status) ~ 1,
                           data = ndata %>% dplyr::mutate(
                             time   = !!time, status = !!status))
        out <- ci

        if(integrated){
          out <-  out$AppCindex$coxph
        }
      }
    }
  }else{
    if(fit == "Random forest"){
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), testX)
      #Create grid of equidistant time points for testing
      timepoints <-  seq(0,
                         max(test_data %>%
                               dplyr::filter(!!status == event_type) %>%
                               dplyr::select(!!time) %>% unlist), length.out = 100L)
      #Calculate probs
      mod <- obj;
      probs <- pec::predictSurvProb(mod,
                                    newdata = ndata, times = timepoints)
      #Calculate brier score
      out <- pec::pec(probs, Surv(time, status) ~ 1,
                      data = ndata %>% dplyr::mutate( time   = !!time, status = !!status),
                      maxtime = max(timepoints),
                      exact = F,
                      exactness = 99L)
      if(integrated){
        out <- yardstick::fun_ibrier_score(out)
      }
    }else{
      print("Not accepted fit")
    }
  }

  if(all & fit != "Random forest"){
    attr(out, 'brier_pred') <- brier
    attr(out, 'ci_pred') <- ci
    attr(out, 'roc_pred') <- roc
    attr(out, 'dev_pred') <- logl
  }
  # attr(out, 'number.of.individuals.cohort') <- nrow(train_data) + nrow(test_data)
  # attr(mod, 'predction.of.model') <- fit
  return(out)
}

#' Integrated Brier Score
#'
#' Predict integration Brier Score in a new dataset.
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_ibrier_score <- function(brier){

  ibrier <- pec::crps(brier, models = "matrix")[1]

  return(ibrier)

}

#' Integrated Brier Score
#'
#' Predict integration Brier Score in a new dataset.
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_ibrier_score_reference <- function(brier){

  ibrier <- pec::crps(brier, models = "Reference")[1]

  return(ibrier)

}


#' Calculate Linear Predictor
#'
#' Predict linear components or hazard ratio.
#' @param
#' obj : coefficient predicted for a survival model e.g. with glmnet \cr
#' test_data: a test holdout dataframe, is recommended to test the model in a separate dataset from which was used for training\cr
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import survival
pred_lp <- function(obj, test_data, time = os_months, status = os_deceased, iter = 0, ...){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  X <- test_data %>%
    dplyr::select(rownames(obj))

  #Maye be more efficient for a large number of covariates with a trained model
  #out <- as.matrix(X) %*% as.vector(unlist(obj))

  mod <-  survival::coxph(survival::Surv(time = test_data %>%
                                           dplyr::select(!!time) %>%
                                           unlist,
                                         event = test_data %>%
                                           dplyr::select(!!status) %>%
                                           unlist)~.,
                          init = as.vector(unlist(obj)),
                          control = coxph.control(iter.max = iter),
                          data = X )
  out <- predict(mod, newdata = X , type = "lp")
  return(out)
}

#' Extract Baseline Hazard
#'
#' Compute baseline hazard for a Cox model
#' @param
#' obj : coefficient predicted for a survival model e.g. with glmnet \cr
#' test_data: a test holdout dataframe, is recommended to test the model in a separate dataset from which was used for training\cr
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import survival
base_haz <- function(obj, test_data, time = os_months, status = os_deceased, iter = 0, ...){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  test_data <- dplyr::arrange(test_data, !!time)

  X <- test_data %>%
    dplyr::select(rownames(obj))

  #Maye be more efficient for a large number of covariates with a trained model
  #out <- as.matrix(X) %*% as.vector(unlist(obj))

  mod <-  survival::coxph(survival::Surv(time = test_data %>%
                                           dplyr::select(!!time) %>%
                                           unlist,
                                         event = test_data %>%
                                           dplyr::select(!!status) %>%
                                           unlist)~.,
                          init = as.vector(unlist(obj)),
                          control = coxph.control(iter.max = iter),
                          data = X )
  bh <- survival::basehaz(mod)

  bhdata <- data_frame( 'cumhaz' = bh %>%
                          dplyr::select(dplyr::contains("haz")) %>%
                          unlist,
                        'time' = bh %>%
                          dplyr::select(dplyr::contains("time")) %>%
                          unlist
                        )

  out <- dplyr::arrange(bhdata, time) %>%
    dplyr::select(cumhaz) %>% unlist

  return(out)
}


#' Survival ROC
#'
#' Estimate survival ROC, option to have the integrated survival ROC.
#' @param
#' obj : coefficient predicted for a survival model e.g. with glmnet \cr
#' test_data: a test holdout dataframe, is recommended to test the model in a separate dataset from which was used for training\cr
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

roc_surv <- function(data, ..., integrated = FALSE){

  risk <- dplyr::quo(...)
  pred <- data %>% dplyr::select(!!risk) %>% unlist

  out <- tdROC::tdROC(X = pred,
                      Y = data %>% dplyr::select(time ) %>% unlist,
                      delta = data %>% dplyr::select(status)%>% unlist,
                      tau = quantile(data %>%
                                       dplyr::select(time) %>% unlist, .5),
                      nboot = 0,
                      alpha = 0.05, n.grid = 100
  )

  if(integrated){
    out <- as.numeric(out$AUC[1] %>% unlist)
  }

  return(out)

}

#' Calculate Brier Skill Score
#'
#' Estimate Brier Skill Score usually from reference
#' @param
#' obj : coefficient predicted for a survival model e.g. with glmnet \cr
#' test_data: a test holdout dataframe, is recommended to test the model in a separate dataset from which was used for training\cr
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

bss <- function(bs_forecast, bs_reference){

  bss <- (1 - (bs_forecast)/(bs_reference))
  return(bss)

}





#' Predict Surv Prob
#'
#' Estimate survival ROC, option to have the integrated survival ROC.
#' @param
#' obj : coefficient predicted for a survival model e.g. with glmnet \cr
#' test_data: a test holdout dataframe, is recommended to test the model in a separate dataset from which was used for training\cr
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

pred_surv_prob <- function(survdata, time = time, status = status,
                           coeff = covariates, iter = 0, inits = 0 , ...){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  coeff <- dplyr::enquo(coeff)

  out <- purrr::map(seq_along(1:nrow(survdata)), function(index){
    coeff <- survdata %>%
      dplyr::select(!!coeff)
    print(index)
    coeff <- coeff[index,]
    print(coeff)
})

  # X <- coeff[rownames(inits),]

  #Maye be more efficient for a large number of covariates with a trained model
  #out <- as.matrix(X) %*% as.vector(unlist(obj))

  mod <-  survival::coxph(survival::Surv(time = test_data %>%
                                           dplyr::select(!!time) %>%
                                           unlist,
                                         event = test_data %>%
                                           dplyr::select(!!status) %>%
                                           unlist)~.,
                          init = as.vector(unlist(obj)),
                          control = coxph.control(iter.max = iter),
                          data = X )
  out <- predict(mod, newdata = X , type = "lp")
}









