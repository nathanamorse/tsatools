% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_cv.R
\name{phi_cv}
\alias{phi_cv}
\title{Critical values for phi tests}
\usage{
phi_cv(phi = c(1, 2, 3), n, alpha = c(0.01, 0.05, 0.1))
}
\arguments{
\item{phi}{an integer specifying which phi is being tested; must be 1 for \eqn{\phi_1}, 2 for \eqn{\phi_2}, or 3 for for \eqn{\phi_3}}

\item{n}{sample size (number of time points)}

\item{alpha}{confidence level; must be .01, .05, or .1}

\item{tstat}{the test statistic for the phi test}
}
\description{
Helps interpret phi tests from ADF procedures for unit roots.
}
\examples{
phi_cv(phi=2, n=150, alpha=.05)
}
