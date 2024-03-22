#' Impulse response functions
#'
#' Calculates impulse response functions (IRFs) from autoregressive distributive lag (ADL) models and generalized error correction models (GECMs). Note: this function only works for independent variables that include all lags from 1 to q.
#'
#' @param mod a `dynlm` object
#' @param xvars a character vector of names of the independent variables to be included
#' @param yvar a character with the name of the dependent variable
#' @param data the data frame or matrix used for the model, for calculating standard deviations
#' @param h the number of lags to calculate IRFs out to
#' @param cum logical: whether to calculate cumulative IRFs
#' @import tidyverse
#' @export


irfs = function(mod, x, y, data, h=10, cum=TRUE) {

  # Function to calculate IRFs for one variable
  irf.x = function(x) {

    # Standard deviation of X
    sd.x = sd(data[,x], na.rm=TRUE)

    # Get coefficients and such
    coefs = coef(mod)
    a = coefs[grep(y, names(coefs))]
    b = coefs[grep(x, names(coefs))]
    p = length(a)
    q = length(b)
    yt = rep(0, h+p+1)

    # Calculate IRFs
    for (t in 1:(h+1)) {
      yt[t+p] = ifelse(t <= q, b[t]*sd.x, 0) + sum(a*yt[(t+p-1):t])
    }

    # Remove leading zeros and multiply by standard deviation
    yt = yt[-c(1:p)]

    # Return cumulative or non-cumulative sum
    if (cum) {cumsum(yt)} else {yt}

  }

  # Apply to all independent variables supplied
  names(x) = x
  cbind(t=0:h, purrr::map_dfc(x, irf.x))

}




