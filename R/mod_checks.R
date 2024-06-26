#' Diagnostic checks for time series models
#'
#' Conducts tests for model misspecification with `lmtest::resettest`, empirical fluctuation with `strucchange::efp` and `strucchange::sctest`, stationarity with the model coefficients, serial correlation with `lmtest::bgtest`, and heteroscedasticity with `lmtest::bptest`.
#'
#' @param mod a `dynlm` object
#' @param y a character with the name of the dependent variable
#' @param data the data frame or matrix used in the model
#' @param orders an integer vector of maximal orders of serial correlation to be tested
#' @import tidyverse lmtest strucchange
#' @export


mod_checks = function(mods, y, data, orders=c(6,12,18,24)) {

  # Function to check for stationarity
  station = function(coefs, y) {
    a = coefs[grep(y, names(coefs))]
    (sum(a)<1) & (all(abs(a)<1))
  }

  # Diagnostic tests
  ramsey = purrr::map(mods, ~lmtest::resettest(.x, power=2, type="regressor", data=data))
  efp = purrr::map(mods, ~strucchange::sctest(strucchange::efp(.x, type="OLS-CUSUM", data=data)))
  exog = purrr::map(mods, ~station(.x$coefficients, y))
  serial = purrr::map(mods, ~lapply(orders, function(i) lmtest::bgtest(.x, order=i, type="Chisq")))
  hetero = purrr::map(mods, lmtest::bptest)

  # Print results
  rbind(
    `Ramsey's RESET test` = paste(round(unlist(purrr::map(ramsey, "statistic")), 3),
                                  purrr::map(purrr::map(ramsey, "p.value"), tsatools::psig)),
    `CUSUM test` = paste(round(unlist(purrr::map(efp, "statistic")), 3),
                         purrr::map(purrr::map(efp, "p.value"), tsatools::psig)),
    `Stationary coefficients` = unlist(exog),
    `Serial correlation present` = purrr::map(purrr::map(serial, purrr::map, "p.value"), ~any(.x<.05)),
    `B-P test for heteroscedasticity` = paste(
      round(unlist(purrr::map(hetero, "statistic")), 2),
      purrr::map(purrr::map(hetero, "p.value"), tsatools::psig))
  ) %>% as.data.frame() %>% rownames_to_column("Test")

}
