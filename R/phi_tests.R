#' ADF tests with uncertain deterministics
#'
#' Conducts augmented Dickey-Fuller (ADF) tests for unit roots.
#'
#' @param y a vector or `ts` object
#' @param pmax the maximum lag to test; if NULL, `pmax` is automatically chosen with 12*(t/100)^(.25)
#' @param lmtestorder the order for the Lagrange multiplier tests for serial correlation; if NULL, `pmax` is used
#' @param lag.strategy the strategy for choosing lags `"p"` (for smallest lags with white noise residuals), `"aic"`, `"bic"`, or `"diffs"` (first significant differences)
#' @param format argument to be passed to knitr::kable for displaying the results output. Possible values are `"latex"`, `"html"`, `"pipe"`, `"simple"`, `"rst"`, and `"none"`. If `"none"` is selected, no output is displayed in the console.
#' @import tidyverse lmtest knitr
#' @export


phi_tests = function(y, pmax=NULL, lmtestorder=NULL, lag.strategy=c("p", "aic", "bic", "diffs"),
                     format="simple") {

  # Function to select ADF tests at each lag
  select.adf = function(adf) {
    p = c(`p` = which(adf$Lags==min(adf$Lags[adf$`LM p-value`>.05])),
          `aic` = which(adf$AIC==min(adf$AIC[adf$`LM p-value`>.05])),
          `bic` = which(adf$BIC==min(adf$BIC[adf$`LM p-value`>.05])),
          `diffs` = first(which(adf$`Diff p-value`[adf$`LM p-value`>.05] < .10)))
    adf[p,] %>% mutate(`Lag strategy` = names(p), .before=1) %>% remove_rownames()
  }

  # ADF tests
  adf.trend = select.adf(tsatools::adf_tests(y, pmax, lmtestorder, type="trend", format="none"))
  adf.drift = select.adf(tsatools::adf_tests(y, pmax, lmtestorder, type="constant", format="none"))

  # Lag strategy
  p = which(lag.strategy==c("p", "aic", "bic", "diffs"))
  p.name = c("smallest lags", "AIC", "BIC", "first significant differences")[p]

  # Initialize objects
  tt = list(reject=FALSE)
  phi3 = list(reject=FALSE)
  phi2 = list(reject=FALSE)
  tu = list(reject=FALSE)
  phi1 = list(reject=FALSE)

  # Test tau_tau
  tt = list(
    lags = adf.trend$Lags[p],
    statistic = adf.trend$`DF test`[p],
    p.value = adf.trend$`DF p-value`[p],
    reject = (adf.trend$`DF p-value`[p] < .05),
    result = ifelse(adf.trend$`DF p-value`[p] < .05, "Stationary", "Unit roots")
  )

  # Test phi3
  if (!tt$reject) {
    phi3 = list(
      lags = adf.trend$Lags[p],
      statistic = adf.trend$Phi3[p],
      cv = adf.trend$`Phi3 CV`[p],
      reject = (adf.trend$Phi3[p] > adf.trend$`Phi3 CV`[p]),
      result = ifelse(adf.trend$Phi3[p] > adf.trend$`Phi3 CV`[p],
                      "Deterministic trend is present", "No deterministic trend")
    )
  }

  # Test phi2
  if (!tt$reject & phi3$reject) {
    phi2 = list(
      lags = adf.trend$Lags[p],
      statistic = adf.trend$Phi2[p],
      cv = adf.trend$`Phi2 CV`[p],
      reject = (adf.trend$Phi2[p] > adf.trend$`Phi2 CV`[p]),
      result = ifelse(adf.trend$Phi2[p] > adf.trend$`Phi2 CV`[p],
                      "Drift is present", "No drift is present")
    )
  }

  # Test tau_mu
  if (!tt$reject & !phi3$reject) {
    tu = list(
      lags = adf.drift$Lags[p],
      statistic = adf.drift$`DF test`[p],
      p.value = adf.drift$`DF p-value`[p],
      reject = (adf.drift$`DF p-value`[p] < .05),
      result = ifelse(adf.drift$`DF p-value`[p] < .05, "Stationary", "Unit roots")
    )
  }

  # Test phi1
  if (!tt$reject & !phi3$reject & !tu$reject) {
    phi1 = list(
      lags = adf.drift$Lags[p],
      statistic = adf.drift$Phi1[p],
      cv = adf.drift$`Phi1 CV`[p],
      reject = (adf.drift$Phi1[p] > adf.drift$`Phi1 CV`[p]),
      result = ifelse(adf.drift$Phi1[p] > adf.drift$`Phi1 CV`[p],
                      "Drift is present", "No drift is present")
    )
  }

  # Final result
  conc = dplyr::case_when(
    tt$reject ~ "stationary",
    !tt$reject & phi3$reject & !phi2$reject ~ "unit root, no drift or trend",
    !tt$reject & phi3$reject & phi2$reject ~ "unit root with drift and trend",
    tu$reject ~ "stationary",
    phi1$reject ~ "unit root with drift",
    TRUE ~ "unit root, no drift or trend"
  )

  # Put results together
  tests = list(`tau_tau`=tt, `phi3`=phi3, `phi2`=phi2, `tau_mu`=tu, `phi1`=phi1)

  # Prepare summary table
  out = tests %>%
    bind_rows(.id="Parameter") %>%
    mutate(Dt = c(rep("c + dt", 3), rep("c", 2)),
           H0 = c("g0 = 0", "g = d = 0", "g = c = d = 0", "g = 0", "g = c = 0"),
           Statistic = paste0(sprintf("%.2f", round(statistic, 2)),
                              ifelse(reject, "*", " "))) %>%
    select(Lags=lags, Dt, Parameter, H0, Statistic, Result=result) %>%
    na.omit()

  # Add summary info to results
  tests = keep(tests, ~length(.x)>1)
  tests$summary = data.frame(out)
  tests$result = conc

  # Display output
  if (format!="none") {
    if (format %in% c("latex", "html")) {
      names(out) = c("Lags", "$D_t$", "Parameter", "$H_0$", "Statistic", "Result")
    }
    out1 = out %>%
      kable(format=format, align="llllcl",
            caption=paste("Unit root tests selected with", p.name))
    print(gsub("Table: ", "", out1))
    cat("Conclusion: Series is", conc)
  }
  invisible(tests)

}



