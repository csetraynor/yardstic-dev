#' Stick covariates to a dataframe for subject
#'
#' Sitck covars...
#' @param
#' survdata: survival dataset \cr
#' covars : covariates \cr

#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

stick_covar <- function(survdata, covars){
  out <- purrr::map(seq_along(1:nrow(survdata)), function(index){
    coeff <- covars[index,] %>% unlist
    df <- data.frame(coeff = coeff  )
    rownames(df) <- colnames(covars)
    return(df)
    })

  return(out)

}


