#' Long-run multipliers
#'
#' Calculates long-run multipliers (LRMs) from autoregressive distributive lag (ADL) models and generalized error correction models (GECMs).
#'
#' @param mod a `dynlm` object
#' @param x a character vector of names of the independent variables to be included
#' @param y a character with the name of the dependent variable
#' @param type the type of model; either `"adl"` or `"gecm"`
#' @param vcov a character string specifying the type of heteroscedasticity-consistent covariance matrix estimation, to be passed to `sandwich::vcovHC` via `msm::deltamethod`
#' @param digits integer; the number of digits to round results to
#' @import msm sandwich tidyverse knitr
#' @export


lrms = function(mod, x, y, type=c("adl", "gecm"), vcov="const", digits=2, format="simple") {

  # LRM for one X in ADL
  lrm.adl = function(x) {

    # Coefficients
    aj = grep(y, names(mod$coef))
    bj = grep(x, names(mod$coef))
    alphas = mod$coef[aj]
    betas = mod$coef[bj]

    # Estimate LRM
    est = sum(betas)/(1-sum(alphas))

    # Estimate SE
    se = msm::deltamethod(
      formula(paste0("~(",paste0("x", bj, collapse="+"), ")/(1-",
                     paste0("x", aj, collapse="+"), ")")),
      mod$coef, sandwich::vcovHC(mod, type=vcov))

    # Significance
    tval = est/se
    pval = 1-pt(abs(tval), df=mod$df.residual)

    # Output
    c(`lrm`=est, `se`=se, `statistic`=tval, `p.value`=pval)

  }

  # LRM for one X in GECM
  lrm.gecm = function(x) {

    # Coefficients
    aj = grep(y, names(mod$coef))
    bj = which(grepl(x, names(mod$coef)) & grepl("L\\(", names(mod$coef)))
    alphas = mod$coef[aj]
    beta1s = mod$coef[bj]

    # Estimate LRM
    est = -sum(beta1s)/sum(alphas)

    # Estimate SE
    se = msm::deltamethod(
      formula(paste0("~(-",paste0("x", bj, collapse="+"), ")/(",
                     paste0("x", aj, collapse="+"), ")")),
      mod$coef, sandwich::vcovHC(mod, type=vcov))

    # Significance
    tval = est/se
    pval = 1-pt(abs(tval), df=mod$df.residual)

    # Output
    c(`lrm`=est, `se`=se, `statistic`=tval, `p.value`=pval)

  }

  # Prepare output object
  if (type=="adl") {out = purrr::map(x, lrm.adl)} else {out = purrr::map(x, lrm.gecm)}
  names(out) = x

  # Display output
  if (format!="none") {
    out1 = bind_rows(out, .id="Variable") %>%
      mutate(across(lrm:statistic, ~sprintf(paste0("%.", digits, "f"), round(.x,digits))),
             LRM = paste0(lrm, ifelse(p.value<.05, "*", " ")),
             SE = paste0("(", se, ")"),
             `p-value` = sprintf("%.4f", round(p.value,4))) %>%
      select(Variable, LRM, SE, `Test statistic`=statistic, `p-value`) %>%
      kable(format=format, align="lcccr",
            caption=paste0("Long-run effects on '", y, "'"))
    out1 = gsub("Table: ", "", out1)
    print(out1)
  }
  invisible(out)

}

