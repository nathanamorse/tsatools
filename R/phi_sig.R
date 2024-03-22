#' Significance levels for phi tests
#'
#' Helps interpret phi tests from ADF procedures for unit roots.
#'
#' @param tstat the test statistic for the phi test
#' @param phi an integer specifying which phi is being tested; must be 1 for \eqn{\phi_1}, 2 for \eqn{\phi_2}, or 3 for for \eqn{\phi_3}
#' @param n sample size (number of time points)
#' @import dplyr
#' @export
#' @examples
#' phi_sig(tstat=3.654, phi=2, n=150)


phi_sig <- function(tstat, phi=c(1,2,3), n) {
  dplyr::case_when(
    tstat < tsatools::phi_cv(phi=phi, n=n, alpha=.1) ~ "",
    tstat < tsatools::phi_cv(phi=phi, n=n, alpha=.05) ~ " .",
    tstat < tsatools::phi_cv(phi=phi, n=n, alpha=.01) ~ " *",
    TRUE ~ " **"
  )
}
