#' Rename covariates from ADL and GECM models
#'
#' Makes coefficient names more readable.
#'
#' @param mod a `dynlm` object
#' @param key a named character vector for renaming coefficients, where the name of each element is the name of a variable used in the code for the model and the value of the element is the display name to replace the code name. (For example, `c(infl="Inflation rate")` replaces `infl` with `Inflation rate`.)
#' @param latex, logical: whether to include latex notation ($\\\\Delta$ instead of `d` and $t-1$ instead of `t-1`). Recommended if the output format is html or latex.
#' @import stringr
#' @export


rename_coefs = function(mod, key=NULL, latex=FALSE) {

  # Store variable info in data frame
  x = data.frame(var = names(mod$coefficients)) %>%
    mutate(var = str_remove(var, "zoo\\(coredata\\(x\\), tt\\)\\."),
           lag = str_detect(var, "L\\("),
           diff = str_detect(var, "d\\("),
           name = str_replace_all(var, ".*\\(|,.*|\\).*", "") %>%
             str_replace("Intercept", "\\(Intercept\\)")) %>%

    # Keep track of which lag is which
    mutate(lags = ifelse(grepl(",", var), str_replace_all(var, ".*\\, +|c\\(|\\).*", ""), NA),
           i = as.numeric(ifelse(grepl(",", var), str_extract(var, "[0-9]+$"), NA)),
           i0 = as.numeric(purrr::map(lags, ~first(unlist(str_split(.x, "\\:"))))),
           i1 = as.numeric(purrr::map(lags, ~last(unlist(str_split(.x, "\\:")))))) %>%
    rowwise() %>%
    mutate(h = ifelse(!is.na(lags), c(i0:i1)[i], NA),
           h = ifelse(lag & is.na(h), 1, h)) %>%

    # Paste together information into final variable names
    mutate(final = name,
           final = ifelse(diff, paste("d", final), final),
           final = ifelse(lag, paste0(final, " (t-", h, ")"), final),
           final = ifelse(lag & h==0, str_remove(final, "-0"), final),
           final = ifelse(final!="(Intercept)" & !grepl("\\)", final), paste(final, "(t)"), final))

  # Store variable names as vector
  x1 = x$final

  # Rename variables
  if(!is.null(key)) {
    for (i in 1:length(key)) x1 = str_replace(x1, names(key)[i], key[i])
  }

  # Prepare for latex
  if (latex) {
    x1 = ifelse(grepl("^d ", x1), gsub("^d ", "$\\\\Delta$ ", x1), x1) %>%
      str_replace("\\(t", "$\\(t") %>%
      str_replace("(?<=[0-9])\\)|(?<=\\(t)\\)", "\\)$")
  }

  # Return model
  names(mod$coefficients) = x1
  mod

}






