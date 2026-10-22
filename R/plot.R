#' Plot response time against person-item distance
#'
#' Displays the distance-difficulty pattern. For every observed response time,
#' the speed-adjusted log time \eqn{\log T_{ij} + \hat\tau_i - \hat\lambda_j}
#' is plotted against the signed distance \eqn{\hat\theta_i - \hat\delta_j}
#' using EAP person estimates, binned, and overlaid with the model-implied
#' curve \eqn{-\hat\gamma_j h(\theta - \delta_j)}.
#'
#' @param fit A \code{grmrt} object.
#' @param items Items to display; defaults to all items pooled. Use
#'   \code{by_item = TRUE} for one panel per item.
#' @param by_item Logical; draw one panel per item.
#' @param bin_width Width of distance bins.
#' @param min_n Minimum number of observations per displayed bin.
#' @param main Plot title (pooled plot only).
#' @param ... Further arguments passed to [graphics::plot()].
#' @return Invisibly, a data frame of the binned means.
#' @examples
#' set.seed(4)
#' sim <- simulate_grmrt(n_persons = 400, n_items = 6, gamma = 0.3)
#' fit <- grmrt_fit(sim$responses, sim$logrt, se = FALSE)
#' plot_time_distance(fit)
#' @export
plot_time_distance <- function(fit, items = NULL, by_item = FALSE,
                               bin_width = 0.5, min_n = 30,
                               main = "Distance-difficulty pattern", ...) {
  stopifnot(inherits(fit, "grmrt"))
  sc <- score_persons(fit)
  J <- length(fit$items)
  idx <- if (is.null(items)) seq_len(J) else
    if (is.character(items)) match(items, fit$item_names) else items
  idx <- idx[colSums(fit$data$obsY)[idx] > 0]
  if (!length(idx)) stop("None of the selected items has observed response times.")
  rows <- lapply(idx, function(j) {
    it <- fit$items[[j]]; y <- fit$data$Y[, j]; ok <- !is.na(y)
    dl <- .item_delta(it, fit$spec)
    data.frame(item = fit$item_names[j], j = j,
               dist = sc$theta[ok] - dl,
               adj = y[ok] + sc$tau[ok] - it$lambda)
  })
  d <- do.call(rbind, rows)
  brk <- seq(floor(min(d$dist) / bin_width) * bin_width,
             ceiling(max(d$dist) / bin_width) * bin_width + bin_width, by = bin_width)
  d$bin <- cut(d$dist, brk, include.lowest = TRUE)
  mids <- (brk[-1] + brk[-length(brk)]) / 2
  agg <- function(dd) {
    m <- tapply(dd$adj, dd$bin, mean); n <- tapply(dd$adj, dd$bin, length)
    s <- tapply(dd$adj, dd$bin, stats::sd)
    out <- data.frame(mid = mids, mean = as.vector(m), n = as.vector(n),
                      se = as.vector(s) / sqrt(as.vector(n)))
    out[!is.na(out$n) & out$n >= min_n, ]
  }
  curve_fun <- function(js, u) {
    rowMeans(vapply(js, function(j) {
      -.dterm(u, 0, fit$items[[j]]$gamma, fit$ctrl$kappa)
    }, numeric(length(u))))
  }
  draw <- function(dd, js, main) {
    a <- agg(dd)
    if (!nrow(a)) return(invisible(a))
    u <- seq(min(a$mid), max(a$mid), length.out = 200)
    cv <- curve_fun(js, u)
    off <- stats::weighted.mean(a$mean - curve_fun(js, a$mid), a$n)
    yl <- range(c(a$mean - 2 * a$se, a$mean + 2 * a$se, cv + off))
    graphics::plot(a$mid, a$mean, pch = 19, ylim = yl, main = main,
                   xlab = expression(hat(theta) - hat(delta)[j]),
                   ylab = "Speed-adjusted log time", ...)
    graphics::segments(a$mid, a$mean - 2 * a$se, a$mid, a$mean + 2 * a$se)
    graphics::lines(u, cv + off, lwd = 2)
    graphics::abline(v = 0, lty = 3)
    invisible(a)
  }
  if (by_item) {
    k <- length(idx); nc <- ceiling(sqrt(k))
    op <- graphics::par(mfrow = c(ceiling(k / nc), nc), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op))
    res <- lapply(seq_along(idx), function(q) {
      a <- draw(d[d$j == idx[q], ], idx[q], fit$item_names[idx[q]])
      if (nrow(a)) a$item <- fit$item_names[idx[q]]
      a
    })
    return(invisible(do.call(rbind, res)))
  }
  draw(d, idx, main)
}
