#' Serial correlation tests
#'
#' Conducts Breusch-Godfrey tests for serial correlation with multiple orders at once. Calls lmtest::bgtest.
#'
#' @param mod a fitted 'lm' object or formula to be passed to lmtest::bgtest
#' @param orders an integer vector of maximal orders of serial correlation to be tested
#' @param type the type of test statistic to be used; must be either "Chisq" or "F"
#' @param format argument to be passed to knitr::kable for displaying the results output. Possible values are `"latex"`, `"html"`, `"pipe"`, `"simple"`, `"rst"`, and `"none"`. If `"none"` is selected, no output is displayed in the console.
#' @import lmtest tibble knitr
#' @export


bg_tests = function(mod, orders, type="Chisq", format="simple") {

  # Make sure arguments work
  type = match.arg(type, choices=c("Chisq", "F"))

  # Conduct BG tests
  tests = lapply(orders, function(i) bgtest(mod, order=i, type=type))
  out = tibble(Order=orders,
               `LM test statistic`=unlist(map(tests, "statistic")),
               `p-value`=unlist(map(tests, "p.value")))
  out$Conclusion = ifelse(out$`p-value`<.05, "Serial correlation present",
                          "No serial correlation detected")

  # Display output
  if (format!="none") {
    out1 = kable(out, format=format, digits=3,
                 caption="Breusch-Godfrey tests for serial correlation")
    out1 = gsub("Table: ", "", out1)
    print(out1)
  }
  invisible(out)

}
