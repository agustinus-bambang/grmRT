#' Item-level tests of the distance-difficulty effect
#'
#' Tests, for each item, whether response times follow the inverted-U
#' (distance-difficulty) pattern, that is, \eqn{H_0: \gamma_j = 0}. Items whose
#' response times do not slow down near the item location (non-significant or
#' negative slopes) are flagged as candidates for review, because their timing
#' does not reflect the trait-relevant deliberation expected for a
#' well-functioning item (for example, ambiguous wording, extreme or socially
#' desirable content, or careless processing). For asymmetric fits, the
#' relative asymmetry (upper minus lower slope, relative to the item average)
#' is also tested.
#'
#' @param fit A \code{grmrt} object.
#' @param method \code{"wald"} (fast, uses the fitted covariance matrix) or
#'   \code{"lr"} (likelihood-ratio tests from refits with the item's slope fixed
#'   at zero, warm-started from \code{fit}; symmetric fits only).
#' @param p_adjust Multiplicity adjustment passed to [stats::p.adjust()].
#' @param alpha Significance level for flagging.
#' @return A data frame with one row per item.
#' @examples
#' set.seed(2)
#' sim <- simulate_grmrt(n_persons = 400, n_items = 5, gamma = c(0.3, 0.3, 0.3, 0.3, 0))
#' fit <- grmrt_fit(sim$responses, sim$logrt)
#' item_diagnostics(fit)
#' @export
item_diagnostics <- function(fit, method = c("wald", "lr"), p_adjust = "holm",
                             alpha = 0.05) {
  stopifnot(inherits(fit, "grmrt"))
  method <- match.arg(method)
  est <- fit$estimates; J <- length(fit$items)
  asym <- fit$spec$asymmetric
  gpar <- if (asym) c("gamma_lower", "gamma_upper") else "gamma"
  V <- fit$vcov_natural
  out <- data.frame(item = fit$item_names, stringsAsFactors = FALSE)
  g <- se <- stat <- p <- rep(NA_real_, J)
  for (j in seq_len(J)) {
    if (fit$spec$gamma_zero[j]) next
    rows <- which(est$item == fit$item_names[j] & est$parameter %in% gpar)
    g[j] <- mean(est$estimate[rows])
    if (!is.null(V)) {
      w <- rep(1 / length(rows), length(rows))
      se[j] <- sqrt(drop(t(w) %*% V[rows, rows, drop = FALSE] %*% w))
    }
  }
  out$gamma <- g; out$se <- se
  if (method == "wald") {
    if (is.null(V)) stop("Wald tests need standard errors; refit with se = TRUE or use method = 'lr'.")
    stat <- (g / se)^2; df <- 1
    p <- stats::pchisq(stat, df, lower.tail = FALSE)
    out$z <- g / se
  } else {
    if (asym) stop("Likelihood-ratio tests of gamma = 0 are available for symmetric fits only; ",
                   "in the asymmetric model, use method = 'wald' (the mean slope is identified).")
    df <- 1
    for (j in seq_len(J)) {
      if (fit$spec$gamma_zero[j]) next
      args <- fit$settings
      args$responses <- fit$data$X; args$logrt <- fit$data$Y
      args$gamma_zero <- unique(c(which(fit$spec$gamma_zero), j))
      args$se <- FALSE; args$start <- fit
      f0 <- do.call(grmrt_fit, args)
      stat[j] <- max(2 * (fit$logLik - f0$logLik), 0)
    }
    p <- stats::pchisq(stat, df, lower.tail = FALSE)
    out$LR <- stat
  }
  out$p <- p
  out$p_adj <- stats::p.adjust(p, method = p_adjust)
  out$flag <- ifelse(is.na(out$p_adj), "fixed",
                ifelse(out$p_adj < alpha & out$gamma > 0, "inverted-U",
                  ifelse(out$p_adj < alpha & out$gamma < 0, "reversed", "no evidence")))
  if (asym && !is.null(V)) {
    asy <- asy_se <- rep(NA_real_, J)
    for (j in seq_len(J)) {
      if (fit$spec$gamma_zero[j]) next
      rl <- which(est$item == fit$item_names[j] & est$parameter == "gamma_lower")
      ru <- which(est$item == fit$item_names[j] & est$parameter == "gamma_upper")
      asy[j] <- est$estimate[ru] - est$estimate[rl]
      Vs <- V[c(ru, rl), c(ru, rl)]
      asy_se[j] <- sqrt(max(Vs[1, 1] + Vs[2, 2] - 2 * Vs[1, 2], 0))
    }
    out$asymmetry <- asy; out$asymmetry_se <- asy_se
    out$asymmetry_p <- 2 * stats::pnorm(-abs(asy / asy_se))
  }
  attr(out, "method") <- method
  out
}
