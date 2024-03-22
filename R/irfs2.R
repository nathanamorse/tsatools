#' Impulse response functions
#'
#' Calculates impulse response functions (IRFs) from autoregressive distributive lag (ADL) models and generalized error correction models (GECMs). Note: this function only works for independent variables that include all lags from 1 to q.
#'
#' @param mod a `dynlm` object
#' @param xvars a character vector of names of the independent variables to be included
#' @param yvar a character with the name of the dependent variable
#' @param h the number of lags to calculate IRFs out to
#' @param cum logical: whether to calculate cumulative IRFs
#' @param data the data frame or matrix used for the model. If `data` is supplied, IRFs are reported in terms of standard deviations; if `data` is NULL, IRFs are reported in terms of the original units of each variable.
#' @import tidyverse
#' @export


irf.adl = function(mod, x.var, y.var, h=10, cum=TRUE, data=NULL) {

  # Standard deviations
  sd.x = ifelse(!is.null(data), sd(data[,x.var], na.rm=TRUE), 1)
  sd.y = ifelse(!is.null(data), sd(data[,y.var], na.rm=TRUE), 1)

  # Get coefficients and such
  coefs = coef(mod)
  a = coefs[grep(y.var, names(coefs))]
  b = coefs[grep(x.var, names(coefs))]
  p = length(a)
  q = length(b)
  yt = rep(0, h+p+1)

  # Calculate IRFs
  for (t in 1:(h+1)) {
    yt[t+p] = ifelse(t <= q, b[t]*sd.x, 0) + sum(a*yt[(t+p-1):t])
  }

  # Remove leading zeros and multiply by standard deviation
  yt = yt[-c(1:p)]*sd.y

  # Return cumulative or non-cumulative sum
  if (cum) {cumsum(yt)} else {yt}

}


irf.gecm = function(mod, x.var, y.var, h=10, cum=TRUE, data=NULL) {

  # Standard deviations
  sd.x = ifelse(!is.null(data), sd(data[,x.var], na.rm=TRUE), 1)
  sd.y = ifelse(!is.null(data), sd(data[,y.var], na.rm=TRUE), 1)

  # Get coefficients and such
  coefs = coef(mod)
  a1 = coefs[grepl(y.var, names(coefs)) & !grepl("^d ", names(coefs))]
  ad = coefs[grepl(y.var, names(coefs)) & grepl("^d ", names(coefs))]
  b0 = coefs[grepl(x.var, names(coefs)) & grepl("^d ", names(coefs))]
  b1 = coefs[grepl(x.var, names(coefs)) & !grepl("^d ", names(coefs))]
  p = length(ad)+1
  q = length(b0)
  yt = rep(0, h+p+1)

  # Calculate IRFs
  for (t in 1:(h+1)) {
    yt[t+p] = ifelse(t <= q, b0[t]*sd.x, 0) + ifelse(t>1, b1*sd.x, 0) +
      yt[t+p-1] + a1*yt[t+p-1] + ifelse(p>1, ad*diff(yt)[t+p-2], 0)
  }

  # Remove leading zeros and multiply by standard deviation
  yt = yt[-c(1:p)]*sd.y

  # Return cumulative or non-cumulative sum
  if (cum) {cumsum(yt)} else {yt}

}


irfs = function(mod, x.vars, y.var, type=c("adl", "gecm"), h=10, cum=TRUE, data=NULL) {

  # Match argument
  type = match.arg(type)

  # Apply to all independent variables supplied
  names(x.vars) = x.vars
  cbind(
    t = 0:h,
    switch(type,
           adl = map_dfc(x.vars, ~irf.adl(mod, .x, y.var, h, cum, data)),
           gecm = map_dfc(x.vars, ~irf.gecm(mod, .x, y.var, h, cum, data)))
  )

}


irf1 = irfs(adl2b, x.vars=c("infl", "unemp", "lei"), y.var="app", type="adl", data=df.ts)
irf2 = irfs(gecm2b, x.vars=c("infl", "unemp", "lei"), y.var="app", type="gecm", data=df.ts)
plot(ts(irf2$infl))



