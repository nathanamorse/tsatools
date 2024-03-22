#' Critical values for phi tests
#'
#' Helps interpret phi tests from ADF procedures for unit roots.
#'
#' @param phi an integer specifying which phi is being tested; must be 1 for \eqn{\phi_1}, 2 for \eqn{\phi_2}, or 3 for for \eqn{\phi_3}
#' @param tstat the test statistic for the phi test
#' @param n sample size (number of time points)
#' @param alpha confidence level; must be .01, .05, or .1
#' @import dplyr
#' @export
#' @examples
#' phi_cv(phi=2, n=150, alpha=.05)


phi_cv = function(phi=c(1,2,3), n, alpha=c(.01, .05, .1)) {

  # Store critical values
  cvals = list(
    phi1 = rbind(c(7.88, 5.18, 4.12), c(7.06, 4.86, 3.94), c(6.7, 4.71, 3.86),
                 c(6.52, 4.63, 3.81), c(6.47, 4.61, 3.79), c(6.43, 4.59,3.78)),
    phi2 = rbind(c(8.21, 5.68, 4.67), c(7.02, 5.13, 4.31), c(6.5, 4.88, 4.16),
                 c(6.22, 4.75, 4.07), c(6.15, 4.71, 4.05), c(6.09, 4.68, 4.03)),
    phi3 = rbind(c(10.61, 7.24, 5.91), c(9.31, 6.73, 5.61), c(8.73, 6.49, 5.47),
                 c(8.43, 6.49, 5.47), c(8.34, 6.3, 5.36), c(8.27, 6.25, 5.34))
  )

  # Get row number based on sample size
  i = dplyr::case_when(n<25 ~ 1, n<50 ~ 2, n<100 ~ 3, n<250 ~ 4, n<500 ~ 5, TRUE ~ 6)

  # Get column number based on alpha level
  j = dplyr::case_when(alpha==.01 ~ 1, alpha==.05 ~ 2, alpha==.1 ~ 3)

  # Return critical value
  cvals[[phi]][i,j]

}



