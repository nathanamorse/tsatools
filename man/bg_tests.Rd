% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bg_tests.R
\name{bg_tests}
\alias{bg_tests}
\title{Serial correlation tests}
\usage{
bg_tests(mod, orders, type = "Chisq", format = "simple")
}
\arguments{
\item{mod}{a fitted 'lm' object or formula to be passed to lmtest::bgtest}

\item{orders}{an integer vector of maximal orders of serial correlation to be tested}

\item{type}{the type of test statistic to be used; must be either "Chisq" or "F"}

\item{format}{argument to be passed to knitr::kable for displaying the results output. Possible values are `"latex"`, `"html"`, `"pipe"`, `"simple"`, `"rst"`, and `"none"`. If `"none"` is selected, no output is displayed in the console.}
}
\description{
Conducts Breusch-Godfrey tests for serial correlation with multiple orders at once. Calls lmtest::bgtest.
}
