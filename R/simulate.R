#' Simulate responses and response times from the joint model
#'
#' Generates Likert-type responses from the graded response model and log
#' response times from the distance-difficulty model used by [grmrt_fit()].
#' Optionally, a proportion of response times is contaminated to mimic
#' implausible timings (for example, inter-click intervals in single-page web
#' questionnaires).
#'
#' @param n_persons Number of respondents.
#' @param n_items Number of items.
#' @param n_cat Number of response categories.
#' @param a,delta Item discriminations and locations. Defaults are drawn from
#'   \code{U(1, 2.5)} and \code{U(-1.5, 1)}.
#' @param b Optional \code{n_items x (n_cat - 1)} threshold matrix. If omitted,
#'   thresholds are spread symmetrically around \code{delta}.
#' @param lambda,sigma Item time intensities and residual SDs of log time.
#' @param gamma Distance slopes: a vector (symmetric) or a two-column matrix of
#'   lower and upper slopes (asymmetric).
#' @param sd_tau,rho Speed SD and trait-speed correlation.
#' @param error \code{"normal"} or \code{"t"} errors.
#' @param nu Degrees of freedom of the \code{"t"} errors.
#' @param contamination Proportion of response times replaced by contaminated
#'   values (log time shifted by \code{contam_shift} times a standard normal).
#' @param contam_shift SD (in log seconds) of the contaminating shift.
#' @param kappa Smoothing constant of the distance function.
#' @return A list with \code{responses}, \code{logrt}, \code{theta}, \code{tau},
#'   and \code{truth} (a data frame of item parameters).
#' @examples
#' sim <- simulate_grmrt(n_persons = 200, n_items = 6)
#' str(sim$responses)
#' @export
simulate_grmrt <- function(n_persons = 1000, n_items = 10, n_cat = 5,
                           a = NULL, delta = NULL, b = NULL, lambda = NULL,
                           gamma = 0.15, sigma = NULL, sd_tau = 0.4, rho = -0.1,
                           error = c("normal", "t"), nu = 5, contamination = 0,
                           contam_shift = 1.5, kappa = 0.01) {
  error <- match.arg(error)
  J <- n_items; N <- n_persons
  if (is.null(a)) a <- stats::runif(J, 1, 2.5)
  if (is.null(delta)) delta <- stats::runif(J, -1.5, 1)
  if (is.null(b)) {
    off <- seq(-1.5, 1.5, length.out = n_cat - 1)
    b <- t(vapply(seq_len(J), function(j) delta[j] + off * stats::runif(1, 0.8, 1.2),
                  numeric(n_cat - 1)))
  }
  b <- matrix(b, nrow = J)
  delta <- rowMeans(b)
  if (is.null(lambda)) lambda <- stats::runif(J, 1, 1.6)
  if (is.null(sigma)) sigma <- stats::runif(J, 0.4, 0.6)
  G <- if (is.matrix(gamma)) gamma else matrix(rep(gamma, length.out = J), J, 1)
  theta <- stats::rnorm(N)
  tau <- sd_tau * (rho * theta + sqrt(1 - rho^2) * stats::rnorm(N))
  X <- matrix(NA_integer_, N, J); Y <- matrix(NA_real_, N, J)
  for (j in seq_len(J)) {
    P <- .grm_probs(a[j], b[j, ], theta)
    U <- stats::runif(N)
    X[, j] <- 1L + rowSums(U > t(apply(P, 1, cumsum))[, -n_cat, drop = FALSE])
    g <- if (ncol(G) == 1L) G[j, 1] else G[j, ]
    mu <- lambda[j] - tau - .dterm(theta, delta[j], g, kappa)
    e <- if (error == "normal") stats::rnorm(N, 0, sigma[j]) else sigma[j] * stats::rt(N, nu)
    Y[, j] <- mu + e
    if (contamination > 0) {
      k <- stats::runif(N) < contamination
      Y[k, j] <- Y[k, j] + contam_shift * stats::rnorm(sum(k))
    }
  }
  nm <- paste0("item", seq_len(J))
  colnames(X) <- colnames(Y) <- nm
  truth <- data.frame(item = nm, a = a, delta = delta, lambda = lambda,
                      sigma = sigma)
  if (ncol(G) == 1L) truth$gamma <- G[, 1] else {
    truth$gamma_lower <- G[, 1]; truth$gamma_upper <- G[, 2]
  }
  truth <- cbind(truth, stats::setNames(as.data.frame(b), paste0("b", seq_len(n_cat - 1))))
  list(responses = X, logrt = Y, theta = theta, tau = tau, truth = truth,
       sd_tau = sd_tau, rho = rho)
}
