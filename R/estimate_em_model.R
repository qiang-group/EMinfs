#' Fit a Gamma Frailty Model for Interval Censored Data
#'
#' @param d A numeric vector for the tumor status indicator (1 if tumor present, 0 otherwise).
#' @param s A numeric vector for the informative censoring indicator (1 if informative, 0 otherwise).
#' @param Ci A numeric vector of the observed censoring times.
#' @param Xp A numeric matrix of covariates.
#' @param n_knots An integer specifying the number of interior knots for the splines. Default is 5.
#' @param spline_order An integer specifying the order of the I-splines and M-splines. Default is 2.
#' @param tol The convergence tolerance for the EM algorithm. Default is 1e-5.
#' @param max_iter The maximum number of iterations for the EM algorithm. Default is 10000.
#' @return A list containing the final model estimates and convergence details.
#' @export
#'
estimate_em_model <- function(d, s, Ci, Xp, n_knots = 5, spline_order = 2, tol = 1e-5, max_iter = 10000) {

  # --- Setup Splines ---
  ti <- Ci
  ti_max <- max(ti) + 1e-05
  ti_min <- min(ti) / 2
  knots <- seq(ti_min, ti_max, length.out = (n_knots + 2))

  k <- length(knots) - 2 + spline_order

  bCi <- t(Ispline(x = Ci, order = spline_order, knots = knots))
  bdervativeCi <- t(Mspline(x = Ci, order = spline_order, knots = knots))

  # --- Initial Values ---
  Cg0 <- Tg0 <- rep(0.5, k)
  Cb0 <- Tb0 <- rep(0.5, ncol(Xp))
  tau0 <- 0.5

  i <- 0
  diff0 <- 1

  # --- EM Algorithm Loop ---
  while (diff0 > tol && i < max_iter) {
    i <- i + 1
    Cxb0 <- Xp %*% Cb0
    Txb0 <- Xp %*% Tb0

    Cg0[Cg0 < 1e-5] <- 1e-5
    Tg0[Tg0 < 1e-5] <- 1e-5

    bi <- tau0 + (bCi %*% Cg0) * exp(Cxb0) + (bCi %*% Tg0) * exp(Txb0)
    ci <- tau0 + (bCi %*% Cg0) * exp(Cxb0)

    Eetai <- (1 - d) * (tau0 + s) / bi + d * ((tau0 + s) / ci) * ((1 - (ci/bi)^(tau0 + s + 1)) / (1 - (ci/bi)^(tau0 + s)))
    Elogetai <- (1 - d) * (digamma(tau0 + s) - log(bi)) +
      d * (digamma(tau0 + s) - (log(ci) - (ci/bi)^(tau0 + s) * log(bi)) / (1 - (ci/bi)^(tau0 + s)))

    EZi <- d * ((bCi %*% Tg0) * exp(Txb0)) * (tau0 + s) / ci * bi^(tau0 + s) / (bi^(tau0 + s) - ci^(tau0 + s))
    EZil <- t(t(bCi) * Tg0) / as.vector(bCi %*% Tg0) * as.vector(EZi)
    dv <- bdervativeCi %*% Cg0
    EVil <- t(t(bdervativeCi) * Cg0) / as.vector(dv)

    # --- M-Step ---
    H1 <- function(tau) {
      -length(d) * log(gamma(tau)) + length(d) * tau * log(tau) + tau * sum(Elogetai - Eetai)
    }
    tau1_fit <- try(maxLik::maxLik(H1, start = tau0), silent = TRUE)
    if (inherits(tau1_fit, "try-error")) {
      warning("maxLik failed for tau; using previous value.")
      tau1 <- tau0
    } else {
      tau1 <- summary(tau1_fit)$estimate[, 1]
    }

    H2b <- function(cb1) {
      Cxb1 <- Xp %*% cb1
      rc <- (s %*% EVil) / (as.vector(exp(Cxb1) * Eetai) %*% bCi)
      colSums(as.vector(-bCi %*% as.vector(rc) * exp(Cxb1) * Eetai) * Xp + as.vector(s) * Xp)
    }
    Cb1 <- rootSolve::multiroot(f = H2b, start = Cb0)$root

    H3b <- function(tb1) {
      txb1 <- Xp %*% tb1
      rt <- colSums(EZil) / (as.vector(exp(txb1) * Eetai) %*% bCi)
      colSums(as.vector(-bCi %*% as.vector(rt) * exp(txb1) * Eetai) * Xp + as.vector(EZi) * Xp)
    }
    Tb1 <- rootSolve::multiroot(f = H3b, start = Tb0)$root

    Cg1 <- (s %*% EVil) / (as.vector(exp(Xp %*% Cb1) * Eetai) %*% bCi)
    Tg1 <- colSums(EZil) / (as.vector(exp(Xp %*% Tb1) * Eetai) %*% bCi)

    # --- Check for Convergence ---
    current_params <- c(Cb1, Cg1, Tb1, Tg1, tau1)
    prev_params <- c(Cb0, Cg0, Tb0, Tg0, tau0)
    diff0 <- max(abs(current_params - prev_params))

    Cg0 <- as.vector(Cg1); Tg0 <- as.vector(Tg1)
    Cb0 <- as.vector(Cb1); Tb0 <- as.vector(Tb1)
    tau0 <- tau1
  }

  # --- Final Calculations (Strictly following original script) ---
  variance_matrix <- Louis(tau0, Cb0, Cg0, Tb0, Tg0, bdervativeCi, bCi, s, d, Xp)

  se_list <- tryCatch({
    variances <- diag(solve(variance_matrix))
    all_se <- sqrt(variances)
    num_covariates <- ncol(Xp)
    se_c <- all_se[2:(1 + num_covariates)]
    se_t <- all_se[(2 + num_covariates):(1 + 2 * num_covariates)]
    list(se_c = se_c, se_t = se_t)
  }, error = function(e) {
    warning("Could not compute standard errors: variance-covariance matrix may be singular.")
    num_covariates <- ncol(Xp)
    list(se_c = rep(NA, num_covariates), se_t = rep(NA, num_covariates))
  })

  # --- Return a structured list of results ---
  results <- list(
    coefficients_c = Cb0,
    std_err_c = se_list$se_c,
    coefficients_t = Tb0,
    std_err_t = se_list$se_t,
    spline_coeffs_c = Cg0,
    spline_coeffs_t = Tg0,
    tau = tau0,
    variance_matrix = variance_matrix,
    iterations = i,
    convergence_diff = diff0
  )

  return(results)
}
