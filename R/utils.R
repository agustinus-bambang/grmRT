# Internal helpers -----------------------------------------------------------

# Gauss-Hermite nodes/weights (physicists' convention, weight exp(-x^2)),
# via the Golub-Welsch algorithm.
.gauss_hermite <- function(n) {
  if (n == 1) return(list(x = 0, w = sqrt(pi)))
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2)
  J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = sqrt(pi) * e$vectors[1, o]^2)
}

# Smooth distance h(u) = sqrt(u^2 + kappa) - sqrt(kappa); h(0) = 0, h'(0) = 0.
.hdist <- function(u, kappa) sqrt(u^2 + kappa) - sqrt(kappa)

# Distance term entering the log-RT mean: gamma * h(theta - delta).
# For the asymmetric model, gamma = c(lower, upper).
.dterm <- function(theta, delta, gamma, kappa) {
  u <- theta - delta
  h <- .hdist(u, kappa)
  if (length(gamma) == 1L) return(gamma * h)
  ifelse(u < 0, gamma[1L] * h, gamma[2L] * h)
}

# GRM category probabilities: returns length(theta) x K matrix.
.grm_probs <- function(a, b, theta) {
  cum <- stats::plogis(outer(theta, b, function(t, bb) a * (t - bb)))
  cum <- cbind(1, cum, 0)
  p <- cum[, -ncol(cum), drop = FALSE] - cum[, -1L, drop = FALSE]
  pmax(p, 1e-300)
}

.logsumexp_rows <- function(M) {
  m <- apply(M, 1L, max)
  m + log(rowSums(exp(M - m)))
}

# Log-density of the RT error on residual r.
.ldens_err <- function(r, sigma, nu = Inf) {
  if (is.infinite(nu)) {
    -log(sigma) - 0.5 * log(2 * pi) - 0.5 * (r / sigma)^2
  } else {
    lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi) - log(sigma) -
      (nu + 1) / 2 * log1p((r / sigma)^2 / nu)
  }
}
