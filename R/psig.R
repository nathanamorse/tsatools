#' Asterisk generator for p-values
#'
#' Returns asterisks based on level of significance
#'
#' @param p a numeric vector of probabilities
#' @param stars an integer indicating the maximum number of stars to show
#' @param dot logical: if true, returns "." for significance at .05 < p < 0.10
#' @param blank logical: if true, returns "" for insignificant values; if false, returns NA for insignificant values
#' @import dplyr
#' @export
#' @examples
#' psig(c(.03, .0015, .6, .0002, .09))
#' psig(c(.03, .0015, .6, .0002, .09), dot=FALSE)
#' psig(c(.03, .0015, .6, .0002, .09), blank=FALSE)


psig = function(p, stars=2, dot=TRUE, blank=TRUE) {
  dplyr::case_when(
    p < 0.001 & stars>2 ~ "***",
    p < 0.010 & stars>1 ~ "**",
    p < 0.050 ~ "*",
    p < 0.100 & dot ~ ".",
    TRUE ~ ifelse(blank, "", NA_character_)
  )
}
