#' Person scores using responses and response times
#'
#' Computes expected a posteriori (EAP) estimates of the trait (\eqn{\theta})
#' and speed (\eqn{\tau}) from the fitted joint model, together with the
#' trait EAP that uses the item responses alone. Because the distance term
#' links response times to \eqn{\theta}, response times act as collateral
#' information, and the posterior SD of \eqn{\theta} shrinks relative to the
#' responses-only score. The column \code{se_reduction} reports this
#' gain as \eqn{1 - SD_{joint}/SD_{resp}}.
#'
#' @param fit A \code{grmrt} object.
#' @return A data frame with one row per respondent.
#' @examples
#' set.seed(3)
#' sim <- simulate_grmrt(n_persons = 300, n_items = 8)
#' fit <- grmrt_fit(sim$responses, sim$logrt, se = FALSE)
#' head(score_persons(fit))
#' @export
score_persons <- function(fit) {
  stopifnot(inherits(fit, "grmrt"))
  E <- .estep(fit$items, fit$struct, fit$data, fit$ctrl, fit$spec)
  th <- fit$ctrl$theta
  wt <- rowSums(E$w, dims = 2L)                     # N x Qt
  m_th <- as.vector(wt %*% th)
  sd_th <- sqrt(pmax(as.vector(wt %*% th^2) - m_th^2, 0))
  m_tau <- apply(E$w * E$Tau, 1L, sum)
  sd_tau <- sqrt(pmax(apply(E$w * E$Tau^2, 1L, sum) - m_tau^2, 0))
  lp <- sweep(E$LR, 2L, fit$ctrl$lprior, "+")
  pr <- exp(lp - .logsumexp_rows(lp))
  m_r <- as.vector(pr %*% th)
  sd_r <- sqrt(pmax(as.vector(pr %*% th^2) - m_r^2, 0))
  data.frame(theta = m_th, theta_sd = sd_th, theta_resp = m_r,
             theta_resp_sd = sd_r, se_reduction = 1 - sd_th / sd_r,
             tau = m_tau, tau_sd = sd_tau)
}
