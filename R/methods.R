#' Methods for grmrt objects
#'
#' Print, summarize, and extract results from a fitted [grmrt_fit()] model.
#' \code{anova()} performs a likelihood-ratio test between two nested fits to
#' the same data (for example, symmetric versus asymmetric, or normal versus
#' Student-\eqn{t} errors; the latter is a boundary test and is conservative).
#'
#' @param x,object A \code{grmrt} object.
#' @param ... For \code{anova()}, a second \code{grmrt} object; otherwise unused.
#' @return \code{summary()} returns an object of class \code{summary.grmrt};
#'   \code{coef()} a data frame of estimates; \code{logLik()} a \code{logLik}
#'   object; \code{anova()} a data frame.
#' @name grmrt-methods
NULL

#' @rdname grmrt-methods
#' @export
print.grmrt <- function(x, ...) {
  cat("Joint GRM and distance-difficulty response-time model\n")
  cat(sprintf("  %d respondents, %d items, %d observed response times\n",
              x$n, length(x$items), x$n_obs_rt))
  cat(sprintf("  Errors: %s; distance slopes: %s; item location: %s\n",
              x$spec$error, if (x$spec$asymmetric) "asymmetric" else "symmetric",
              if (x$spec$delta_free) "free" else "tied to GRM thresholds"))
  cat(sprintf("  logLik = %.2f (df = %d), %s after %d EM iterations\n",
              x$logLik, x$n_par, if (x$converged) "converged" else "NOT converged",
              x$iterations))
  e <- x$estimates
  cat(sprintf("  SD(speed) = %.3f, cor(trait, speed) = %.3f\n",
              e$estimate[e$parameter == "sd_tau"], e$estimate[e$parameter == "rho"]))
  invisible(x)
}

#' @rdname grmrt-methods
#' @export
summary.grmrt <- function(object, ...) {
  e <- object$estimates
  wide <- function(p) {
    r <- e[e$parameter == p & !is.na(e$item), ]
    stats::setNames(r$estimate, r$item)
  }
  tab <- data.frame(item = object$item_names, a = wide("a"), delta = wide("delta"),
                    lambda = wide("lambda"), sigma = wide("sigma"), row.names = NULL)
  if (object$spec$asymmetric) {
    tab$gamma_lower <- wide("gamma_lower"); tab$gamma_upper <- wide("gamma_upper")
  } else tab$gamma <- wide("gamma")
  if (object$spec$error == "t") tab$nu <- wide("nu")
  diag <- tryCatch(item_diagnostics(object), error = function(err) NULL)
  structure(list(fit = object, items = tab, diagnostics = diag,
                 person = e[is.na(e$item), ]), class = "summary.grmrt")
}

#' @rdname grmrt-methods
#' @export
print.summary.grmrt <- function(x, ...) {
  print(x$fit)
  cat("\nItem parameters:\n")
  print(format(x$items, digits = 3), row.names = FALSE)
  if (!is.null(x$diagnostics)) {
    cat("\nDistance-difficulty tests (Wald, Holm-adjusted):\n")
    d <- x$diagnostics[, c("item", "gamma", "se", "z", "p_adj", "flag")]
    print(format(d, digits = 3), row.names = FALSE)
  }
  invisible(x)
}

#' @rdname grmrt-methods
#' @export
coef.grmrt <- function(object, ...) object$estimates

#' @rdname grmrt-methods
#' @export
vcov.grmrt <- function(object, ...) object$vcov_natural

#' @rdname grmrt-methods
#' @export
logLik.grmrt <- function(object, ...) {
  structure(object$logLik, df = object$n_par, nobs = object$n, class = "logLik")
}

#' @rdname grmrt-methods
#' @export
anova.grmrt <- function(object, ...) {
  other <- list(...)
  if (length(other) != 1L || !inherits(other[[1L]], "grmrt"))
    stop("Supply exactly one other grmrt object.")
  m <- list(object, other[[1L]])
  m <- m[order(vapply(m, function(z) z$n_par, 1))]
  if (m[[1]]$n != m[[2]]$n) stop("Models were fitted to different data.")
  lr <- max(2 * (m[[2]]$logLik - m[[1]]$logLik), 0)
  df <- m[[2]]$n_par - m[[1]]$n_par
  data.frame(model = c("restricted", "general"),
             logLik = c(m[[1]]$logLik, m[[2]]$logLik),
             df = c(m[[1]]$n_par, m[[2]]$n_par),
             AIC = c(-2 * m[[1]]$logLik + 2 * m[[1]]$n_par, -2 * m[[2]]$logLik + 2 * m[[2]]$n_par),
             BIC = c(-2 * m[[1]]$logLik + log(m[[1]]$n) * m[[1]]$n_par,
                     -2 * m[[2]]$logLik + log(m[[2]]$n) * m[[2]]$n_par),
             LR = c(NA, lr), p = c(NA, stats::pchisq(lr, df, lower.tail = FALSE)))
}
