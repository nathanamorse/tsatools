#' Standard errors for VECMs
#'
#' Calculates standard errors for coefficients from vector error correction models (VECMs) fitted with `urca::ca.jo`
#'
#' @param cajo a `ca.jo` object
#' @param rls a `cajorls` object
#' @import urca
#' @export


vecm_se = function(cajo, rls) {

  # Get coefficients and info
  alpha <- coef(rls$rlm)[1,]
  beta <- rls$beta
  resids <- resid(rls$rlm)
  N <- nrow(resids)

  # Calculate SEs for alphas
  sigma <- crossprod(resids)/N
  alpha.se <- sqrt(solve(crossprod(cbind(cajo@ZK %*% beta,
                                         cajo@Z1)))[1, 1] * diag(sigma))
  alpha.t <- alpha/alpha.se
  alphas = data.frame(alpha, alpha.se, alpha.t)

  # Calculate SEs for betas
  beta.se <- c(NA, sqrt(diag(kronecker(solve(crossprod(cajo@RK[,-1])),
                                       solve(t(alpha) %*% solve(sigma) %*% alpha)))))
  beta.t <- suppressWarnings(beta[-1]/beta.se)
  betas = data.frame(beta, beta.se, beta.t)

  # Output
  list(`alphas`=alphas, `betas`=betas)

}

