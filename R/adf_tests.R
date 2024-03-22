#' ADF tests
#'
#' Conducts augmented Dickey-Fuller (ADF) tests for unit roots.
#'
#' @param y a vector or `ts` object
#' @param pmax the maximum lag to test; if NULL, `pmax` is automatically chosen with 12*(t/100)^(.25)
#' @param lmtestorder the order for the Lagrange multiplier tests for serial correlation; if NULL, `pmax` is used
#' @param type the model deterministics: "none" for no drift or trend (Dt=0), "constant" for drift (Dt=c), or "trend" for drift and trend (Dt=c+dt)
#' @param format argument to be passed to knitr::kable for displaying the results output. Possible values are `"latex"`, `"html"`, `"pipe"`, `"simple"`, `"rst"`, and `"none"`. If `"none"` is selected, no output is displayed in the console.
#' @import tidyverse lmtest knitr
#' @export


adf_tests <- function(x, pmax=NULL, lmtestorder=NULL, type = c("none", "constant", "trend"),
                      format="simple") {

  # PRELIMINARY STEPS -----------------------------------------------------

  # Make sure arguments work
  type = match.arg(type)
  format = match.arg(format, c("latex", "html", "pipe", "simple", "rst", "none"))

  # Determine pmax
  t = length(x)
  if (is.null(pmax)) pmax = floor(12*(t/100)^(.25))
  if (is.null(lmtestorder)) lmtestorder = pmax
  if (pmax < 0) {
    warning("\nlag.max must be greater or equal to 1 and integer; setting lag.max=4")
    pmax <- 4
  }

  # Prepare data
  testlags <- pmax+1
  x <- as.vector(x)
  z <- diff(x, 1)
  n <- length(z)
  x1 <- embed(z, testlags) # creates matrix of z with no missing values
  z.lag.1 <- dplyr::lag(x,1) # lags y and also adds tsp (time) attribute
  length.y <- seq(1:testlags)
  z.lag.1 <- z.lag.1[-length.y]  # gets rid of first pmax+1 rows of z.lag.1 (drops missing values due to lagging)
  TestObs <- length(z.lag.1)  # gives number of observations for all tests
  adf.data <- cbind(x1, z.lag.1) # creates matrix data for all tests
  t <- 1:length(z.lag.1)


  # DICKEY FULLER TESTS ---------------------------------------------------

  # Estimate dickey fuller regression with constant for the number of test lags down to 2
  if (type == "none")
    my_lms <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ z.lag.1 - 1 + adf.data[,2:x]))
  if (type == "constant")
    my_lms <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ z.lag.1 + 1 + adf.data[,2:x]))
  if (type == "trend")
    my_lms <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ z.lag.1 + 1 + t + adf.data[,2:x]))

  # Pull out the summary for each regression
  lm_sums  <- purrr::map(my_lms, summary)

  # Pull out the t value for the DF test, which is the t-value on z.lag.1 from each regression
  tau.tval <- purrr::map(lm_sums, function(item) item$coefficients["z.lag.1","t value"])

  # Round t-value to 3 decimal points.
  tau.tval <- purrr::map(tau.tval, round, digits = 3)

  # Get the p-value associated with each t-value, assuming a constant in the test regression
  tau.pval <- purrr::map(tau.tval, punitroot, N = TestObs, trend = "c", statistic = "t", na.rm = FALSE)

  # Round the p-value to 3 decimals.
  tau.pval <- purrr::map(tau.pval, round, digits = 3)

  # Pull out AIC from each regression
  AICresult <- purrr::map(my_lms, AIC)

  # Get number of lags in each test regression and make a list for use later
  plags <- as.list(seq(purrr::pmax, 1, -1))

  # Pull out BIC from each regression
  BICresult <- purrr::map(my_lms, BIC)

  # Breusch Godfrey test p-values for residual serial correlation
  lmresult <- purrr::map(my_lms, bgtest, order=lmtestorder, type="Chisq", fill = NA)
  lmstat <- purrr::map(lmresult, "statistic")
  lmstat <- purrr::map(lmstat, round, digits=3)
  lmpvalue <- purrr::map(lmresult, "p.value")
  lmpvalue <- purrr::map(lmpvalue, round, digits=3)


  # NO CONSTANT OR TREND --------------------------------------------------
  if (type == "none") {

    # Get t statistic on highest lag of Delta z
    highestlag <- as.list(seq(pmax,1)+1)
    diff.t <- purrr::pmap(list(lm_sums, highestlag), function(item, item2) item$coefficients[item2,"t value"])

    # Prepare information
    mylist <- list(tau.tval, tau.pval, AICresult, BICresult, diff.t, lmpvalue)
    eachel <- purrr::map(mylist, unlist)
    eachel2 <- as.vector(unlist(eachel))

    # Prepare output
    foroutput <- matrix(eachel2, nrow=pmax, ncol=6, byrow = FALSE)
    Nobsoutput <- rep(TestObs, pmax)
    foroutput3 <- cbind(pmax:1,foroutput, Nobsoutput)
    colnames(foroutput3) <- c("Lags","Test Statistic", "p-value", "AIC", "BIC",
                              "Diff t-value", "LM p-value", "N")
    table.title <- "No Constant or Trend in Test Regression"

  }

  # CONSTANT --------------------------------------------------------------
  if (type == "constant") {

    # Get t statistic on highest lag of Delta z
    highestlag <- as.list(seq(pmax,1)+2)
    diff.t <- purrr::pmap(list(lm_sums, highestlag), function(item, item2) item$coefficients[item2,"t value"])

    # Prepare information
    phi.reg <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ -1+ adf.data[,2:x]))
    Ftestresult <- purrr::pmap(list(my_lms, phi.reg), function(item, item2) anova(item, item2)$F[2])
    mylist <- list(tau.tval, tau.pval, AICresult, BICresult, diff.t, lmpvalue, Ftestresult)
    eachel <- purrr::map(mylist, unlist)
    eachel2 <- as.vector(unlist(eachel))

    # Prepare output
    foroutput <- matrix(eachel2, nrow=pmax, ncol=, byrow = FALSE)
    Nobsoutput <- rep(TestObs, pmax)
    foroutput2 <- cbind(pmax:1,foroutput, Nobsoutput)
    colnames(foroutput2) <- c("Lags","Test Statistic", "p-value", "AIC", "BIC",
                              "Diff t-value", "LM p-value", "Phi1", "N")

    # Conduct tests for p=0
    result <- lm(adf.data[,1] ~ z.lag.1 + 1)
    tau.tval <- round(coef(summary(result))[2, 3], digits = 3)
    tau <- round(coef(summary(result))[2, 1], digits = 3)
    tausq <- tau*tau
    teststat <- as.matrix(t(tau.tval))
    colnames(teststat) <- c("tau")
    phi.reg <- lm(adf.data[,1] ~ -1)
    Ftestresult <- anova(result, phi.reg)$F[2]
    AICresult <- AIC(result)
    BICresult <- BIC(result)
    lm0result <- bgtest(result, order=lmtestorder, type = "Chisq", fill = NA)
    lm0pvalue <- round(unlist(lm0result["p.value"]), digits =3)
    tau.pval <- round(punitroot(q=tau.tval, N = TestObs, trend = "c", statistic = "t", na.rm = FALSE),
                      digits = 3)
    teststats <- as.matrix(t(c(0, tau.tval, tau.pval, AICresult, BICresult,
                               NA, lm0pvalue, Ftestresult, TestObs)))
    colnames(teststats) <- c("Lags", "Test Statistic", "p-value",  "AIC", "BIC",
                             "Diff t-value", "LM p-value", "Phi1", "N")
    foroutput3 <- rbind(foroutput2, teststats)
    table.title <- "Constant in Test Regression"

  }

  # TREND -----------------------------------------------------------------
  if (type == "trend") {

    # Get t statistic on highest lag of Delta z
    highestlag <- as.list(seq(pmax,1)+3)
    diff.t <- purrr::pmap(list(lm_sums, highestlag), function(item, item2) item$coefficients[item2,"t value"])

    # Prepare information
    phi2.reg <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ -1+ adf.data[,2:x]))
    Ftestresult2 <- purrr::pmap(list(my_lms, phi2.reg), function(item, item2) anova(item, item2)$F[2])
    phi3.reg <- lapply(testlags:2, function(x) lm(adf.data[,1] ~ adf.data[,2:x]))
    Ftestresult3 <- purrr::pmap(list(my_lms, phi3.reg), function(item, item2) anova(item, item2)$F[2])
    mylist <- list(tau.tval, tau.pval, AICresult, BICresult, diff.t, lmpvalue, Ftestresult2, Ftestresult3)
    eachel <- purrr::map(mylist, unlist)
    eachel2 <- as.vector(unlist(eachel))

    # Prepare output
    foroutput <- matrix(eachel2, nrow=pmax, ncol=8, byrow = FALSE)
    Nobsoutput <- rep(TestObs, pmax)
    foroutput2 <- cbind(pmax:1,foroutput, Nobsoutput)
    colnames(foroutput2) <- c("Lags","Test Statistic", "p-value", "AIC", "BIC",
                              "Diff t-value", "LM p-value", "Phi2", "Phi3", "N")

    # Conduct tests for p=0
    result <- lm(adf.data[,1] ~ z.lag.1 + 1 + t)
    tau.tval <- round(coef(summary(result))[2, 3], digits = 3)
    tau <- round(coef(summary(result))[2, 1], digits = 3)
    tausq <- tau*tau
    teststat <- as.matrix(t(tau.tval))
    colnames(teststat) <- c("tau")
    phi2.reg <- lm(adf.data[,1] ~ -1)
    phi3.reg <- lm(adf.data[,1] ~ 1)
    Ftestresult2 <- anova(result, phi2.reg)$F[2]
    Ftestresult3 <- anova(result, phi3.reg)$F[2]
    AICresult <- AIC(result)
    BICresult <- BIC(result)
    lm0result <- bgtest(result, order=lmtestorder, type = "Chisq", fill = NA)
    lm0pvalue <- round(unlist(lm0result["p.value"]), digits =3)
    tau.pval <- round(punitroot(q=tau.tval, N = TestObs, trend = "ct", statistic = "t", na.rm = FALSE),
                      digits = 3)
    teststats <- as.matrix(t(c(0, tau.tval, tau.pval, AICresult, BICresult,
                               NA, lm0pvalue, Ftestresult2, Ftestresult3, TestObs)))
    colnames(teststats) <- c("Lags", "Test Statistic", "p-value",  "AIC",
                             "BIC", "Diff t-value", "LM p-value", "Phi2"," Phi3", "TestObs")
    foroutput3 <- rbind(foroutput2, teststats)
    table.title <- "Trend in Test Regression"

  }


  # OUTPUT ----------------------------------------------------------------
  df = as.data.frame(foroutput3)
  df$`Diff p-value` = 1-pt(abs(df$`Diff t-value`), df=TestObs)
  if (type=="none") df$`LM test` = round(unlist(lmstat),3)
  if (type!="none") df$`LM test` = round(c(unlist(lmstat), lm0result$statistic), 3)
  if (type=="constant") df$`Phi1 CV` = tsatools::phi_cv(phi=1, n=TestObs, alpha=.05)
  if (type=="trend") df$`Phi2 CV` = tsatools::phi_cv(phi=2, n=TestObs, alpha=.05)
  if (type=="trend") df$`Phi3 CV` = tsatools::phi_cv(phi=3, n=TestObs, alpha=.05)
  df = select(df, Lags, `DF test`=`Test Statistic`, `DF p-value`=`p-value`,
              AIC, BIC, `Diff t-value`, `Diff p-value`, `LM test`, `LM p-value`,
              contains("Phi1"), contains("Phi2"), contains("Phi3"))
  df = round(df, 3)

  # Function to combine test statistics and stars and align along the decimal
  aligner = function(t, p) {
    sprintf("%.3f", t) %>%
      str_replace("NA", "") %>%
      str_pad(max(nchar(.)), side="left", pad="#") %>%
      paste0(p) %>%
      str_pad(max(nchar(.)), side="right", pad="#")
  }

  # Prepare final table for display
  out = df
  out$`DF test` = aligner(df$`DF test`, tsatools::psig(df$`DF p-value`))
  out$`Diff t-value` = aligner(df$`Diff t-value`, tsatools::psig(df$`Diff p-value`))
  out$`LM test` = aligner(df$`LM test`, tsatools::psig(df$`LM p-value`))
  if (type=="constant") out$Phi1 = aligner(df$Phi1, tsatools::phi_sig(df$Phi1, phi=1, n=TestObs))
  if (type=="trend") out$Phi2 = aligner(df$Phi2, tsatools::phi_sig(df$Phi2, phi=2, n=TestObs))
  if (type=="trend") out$Phi3 = aligner(df$Phi3, tsatools::phi_sig(df$Phi3, phi=3, n=TestObs))
  out = select(out, -contains("p-value"), -contains("CV"))

  # Display output
  if (format!="none") {
    out1 = kable(out, format=format, digits=3, align="rcrrcccc",
                 caption=paste("Dickey Fuller Test Results --", table.title))
    out1 = gsub("#", " ", out1)
    out1 = gsub("Table: ", "", out1)
    out1[length(out1)+1] = paste("\nSig. codes: 0 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",
                                 "\nObservations:", TestObs)
    print(out1)
  }
  invisible(df)

}

