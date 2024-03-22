#' Unit root tests
#'
#' Conducts multiple types of tests for unit roots and automatically selects the best fitting models to display.
#'
#' @param y a vector or `ts` object
#' @param pmax the maximum lag to test; if NULL, `pmax` is automatically chosen with 12*(t/100)^(.25)
#' @param lmtestorder the order for the Lagrange multiplier tests for serial correlation; if NULL, `pmax` is used
#' @param type the model deterministics
#' @param tests a character vector specifying which tests to conduct: "adf" for augmented Dickey-Fuller tests, "pp" for Phillips-Perron tests, "dfgls" for Dickey-Fuller tests with generalized least squares, and "kpss" for Kwiatkowski-Phillips-Schmidt-Shin tests. Can include anywhere from one to four of these tests; defaults to all four.
#' @param format aargument to be passed to knitr::kable for displaying the results output. Possible values are `"latex"`, `"html"`, `"pipe"`, `"simple"`, `"rst"`, and `"none"`. If `"none"` is selected, no output is displayed in the console.
#' @import tidyverse lmtest urca knitr
#' @export


ur_tests = function(y, pmax=NULL, lmtestorder=NULL, type = c("none", "constant", "trend"),
                    tests = c("adf", "pp", "dfgls", "kpss"), format="simple") {

  # Determine pmax
  t = length(y)
  if (is.null(pmax)) pmax = floor(12*(t/100)^(.25))
  if (is.null(lmtestorder)) lmtestorder = pmax

  # Initialize list to store data
  urs = list()

  # Determine if tests and types are compatible
  if (type=="none" & !("adf" %in% tests)) {
    stop(paste("PP, DF-GLS, and KPSS tests require a constant or trend.",
               "Try setting `type` to '\"constant\"' or '\"trend\"' or setting `tests` to '\"adf\"'."))
  }

  # ADF TESTS -------------------------------------------------------------
  if ("adf" %in% tests) {
    adf = tsatools::adf_tests(y, pmax, lmtestorder, type, format="none")

    # Various selections for p
    p = c(
      `Smallest lags` = which(adf$Lags==min(adf$Lags[adf$`LM p-value`>.05])),
      `Smallest AIC` = which(adf$AIC==min(adf$AIC[adf$`LM p-value`>.05])),
      `Smallest BIC` = which(adf$BIC==min(adf$BIC[adf$`LM p-value`>.05])),
      `First sig. diffs.` = first(which(adf$`Diff p-value`[adf$`LM p-value`>.05] < .10))
    )

    # Collect ADF tests
    urs[[length(urs)+1]] = adf[p,] %>%
      mutate(Test = "ADF",
             `Lag strategy` = names(p),
             Lags = adf$Lags[p],
             Conclusion = ifelse(`DF p-value`<.05, "Stationary", "Unit roots")) %>%
      select(Test, `Lag strategy`, Lags, Statistic=`DF test`, `p-value`=`DF p-value`, Conclusion)
  }


  # PHILLIPS-PERRON TESTS -------------------------------------------------
  if (type!="none" & "pp" %in% tests) {
    model = type
    pp_short = ur.pp(y, type="Z-tau", model=model, lags="short")
    pp_long = ur.pp(y, type="Z-tau", model=model, lags="long")

    # Collect PP tests
    pps = tibble(
      Test = "Phillips-Perron",
      `Lag strategy` = c("Short", "Long"),
      Lags = c(attr(pp_short, "lag"), attr(pp_long, "lag")),
      Statistic = c(attr(pp_short, "teststat"), attr(pp_long, "teststat")),
      `Critical value` = c(attr(pp_short, "cval")[2], attr(pp_long, "cval")[2])
    )
    pps$Conclusion = ifelse(abs(pps$Statistic)>abs(pps$`Critical value`),
                            "Stationary", "Unit roots")
    urs[[length(urs)+1]] = pps
  }


  # DF-GLS TESTS ----------------------------------------------------------
  if (type!="none" & "dfgls" %in% tests) {
    model = type
    dfgls = ur.ers(y, model=model, lag.max=pmax)
    dfglss = tibble(
      Test = "DF-GLS",
      `Lag strategy` = "Manual",
      Lags = attr(dfgls, "lag"),
      Statistic = attr(dfgls, "teststat"),
      `Critical value` = attr(dfgls, "cval")[2]
    )
    dfglss$Conclusion = ifelse(abs(dfglss$Statistic)>abs(dfglss$`Critical value`),
                               "Stationary", "Unit roots")
    urs[[length(urs)+1]] = dfglss
  }


  # KPSS TESTS ------------------------------------------------------------
  if (type!="none" & "kpss" %in% tests) {
    model = ifelse(type=="trend", "tau", "mu")
    kpss_short = ur.kpss(y, type=model, lags="short")
    kpss_long = ur.kpss(y, type=model, lags="long")

    # Collect KPSS tests
    kpsss = tibble(
      Test = "KPSS",
      `Lag strategy` = c("Short", "Long"),
      Lags = c(attr(kpss_short, "lag"), attr(kpss_long, "lag")),
      Statistic = c(attr(kpss_short, "teststat"), attr(kpss_long, "teststat")),
      `Critical value` = c(attr(kpss_short, "cval")[2], attr(kpss_long, "cval")[2])
    )
    kpsss$Conclusion = ifelse(abs(kpsss$Statistic)>abs(kpsss$`Critical value`),
                              "Unit roots", "Stationary")
    urs[[length(urs)+1]] = kpsss
  }


  # COMBINE ALL TESTS -----------------------------------------------------
  if (type=="none") out = urs[[1]] %>% remove_rownames()
  if (type!="none") {
    out = bind_rows(urs) %>%
      select(1,2,3,4,5,7,6) %>%
      remove_rownames()
  }

  # Display output
  if (format!="none") {
    out1 = out %>%
      mutate(Test = ifelse(lag(Test)==Test & !is.na(lag(Test)), "", Test))
    out2 = kable(out1, format, digits=3, caption="Unit Root Test Results")
    out2 = gsub("Table: ", "", out2)
    out2 = gsub("NA", "  ", out2)
    out2[length(out2)+1] = paste("\nObservations:", t)
    print(out2)
  }
  invisible(out)

}

