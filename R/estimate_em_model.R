estimate_em_model <- function(d, s, Ci, Xp, n_knots = 5, spline_order = 2, tol = 1e-5, max_iter = 1500) {

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
    # Using try() to handle potential errors in maxLik
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

    Cg0 <- as.vector(Cg1)
    Tg0 <- as.vector(Tg1)
    Cb0 <- as.vector(Cb1)
    Tb0 <- as.vector(Tb1)
    tau0 <- tau1

    # Optional: print progress
    # print(paste("Iteration:", i, "Difference:", diff0))
  }

  # --- Final Calculations (Variance) ---
  variance_matrix <- Louis(tau0, Cb0, Cg0, Tb0, Tg0, bdervativeCi, bCi, s, d, Xp)

  # --- Return a structured list of results ---
  results <- list(
    coefficients_c = Cb0,
    spline_coeffs_c = Cg0,
    coefficients_t = Tb0,
    spline_coeffs_t = Tg0,
    tau = tau0,
    variance_matrix = variance_matrix,
    iterations = i,
    convergence_diff = diff0
  )

  return(results)
}
