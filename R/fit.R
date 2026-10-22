#' Fit the joint GRM and distance-difficulty response-time model
#'
#' Estimates a joint model for Likert-type responses and item response times
#' by marginal maximum likelihood (EM algorithm). Responses follow Samejima's
#' graded response model (GRM). Log response times follow a distance-difficulty
#' (inverted-U) model in which responding is slowest when the respondent's trait
#' level is close to the item location:
#' \deqn{\log T_{ij} = \lambda_j - \tau_i - \gamma_j h(\theta_i - \delta_j) + \varepsilon_{ij},}
#' where \eqn{h(u) = \sqrt{u^2 + \kappa} - \sqrt{\kappa}} is a smooth version of
#' \eqn{|u|}, \eqn{\tau_i} is the person's speed, \eqn{\lambda_j} the item time
#' intensity, and \eqn{\gamma_j} the distance-difficulty slope. A positive
#' \eqn{\gamma_j} indicates the expected inverted-U pattern. By default the item
#' location \eqn{\delta_j} is tied to the mean of the GRM thresholds, so that
#' response times carry collateral information about \eqn{\theta}. Errors are
#' normal or Student-\eqn{t}; the latter down-weights implausible timings such as
#' inter-click intervals from single-page web questionnaires. With
#' \code{asymmetric = TRUE}, separate slopes are estimated below and above the
#' item location.
#'
#' Person parameters are \eqn{\theta_i \sim N(0, 1)} and
#' \eqn{\tau_i \mid \theta_i \sim N(\beta\theta_i, \omega)}, reported as the
#' speed standard deviation \eqn{\sigma_\tau} and the trait-speed correlation
#' \eqn{\rho}. The trait dimension is integrated on a rectangular grid; for each
#' person and grid node, speed is integrated by adaptive Gauss-Hermite
#' quadrature centred on its Gaussian-conjugate posterior, which is exact for
#' normal errors. Standard errors use the cross-product of person-level scores.
#'
#' @param responses A matrix or data frame of item responses (persons in rows),
#'   or a \code{grmrt_data} object from [prep_rt()]. Each item's observed
#'   categories are recoded to \code{1, ..., K}. Reverse-keyed items should be
#'   recoded beforehand (see [prep_rt()]).
#' @param logrt A matrix of log response times with the same dimensions as
#'   \code{responses}. Ignored if \code{responses} is a \code{grmrt_data} object.
#'   Missing values are allowed in both matrices.
#' @param error Error distribution of log response times: \code{"normal"} or
#'   \code{"t"}.
#' @param asymmetric Logical; estimate separate distance slopes below and above
#'   the item location.
#' @param delta \code{"tied"} (default) sets \eqn{\delta_j} to the mean GRM
#'   threshold; \code{"free"} estimates it as a separate item parameter.
#' @param gamma_zero Optional item names or indices whose distance slope is
#'   fixed at zero (used for likelihood-ratio tests).
#' @param asym_ref Item (name or index) whose slope is constrained to be
#'   symmetric during estimation in the asymmetric model; defaults to the first
#'   item with a free slope. Item asymmetries are identified only up to a common
#'   shift that is confounded with the trait-speed regression, so they are
#'   reported relative to the item average (sum-to-zero), which makes the
#'   reported estimates invariant to the choice of \code{asym_ref}.
#' @param n_theta,theta_range Number and range of trait grid nodes.
#' @param n_tau Number of adaptive Gauss-Hermite nodes for speed.
#' @param kappa Smoothing constant of the distance function.
#' @param se Logical; compute standard errors.
#' @param max_iter,tol Maximum EM iterations and convergence tolerance on the
#'   relative change in the marginal log-likelihood.
#' @param start Optional fitted \code{grmrt} object used for starting values.
#' @param verbose Logical; print iteration history.
#' @return An object of class \code{grmrt}.
#' @examples
#' set.seed(1)
#' sim <- simulate_grmrt(n_persons = 300, n_items = 5)
#' fit <- grmrt_fit(sim$responses, sim$logrt, se = FALSE)
#' fit
#' @export
grmrt_fit <- function(responses, logrt = NULL, error = c("normal", "t"),
                      asymmetric = FALSE, delta = c("tied", "free"),
                      gamma_zero = NULL, asym_ref = NULL, n_theta = 31, theta_range = c(-5, 5),
                      n_tau = 7, kappa = 0.01, se = TRUE, max_iter = 500,
                      tol = 1e-7, start = NULL, verbose = FALSE) {
  cl <- match.call()
  error <- match.arg(error); delta <- match.arg(delta)
  if (inherits(responses, "grmrt_data")) {
    logrt <- responses$logrt; responses <- responses$responses
  }
  dat <- .make_data(responses, logrt)
  J <- ncol(dat$X)
  gz <- rep(FALSE, J)
  if (!is.null(gamma_zero)) {
    idx <- if (is.character(gamma_zero)) match(gamma_zero, dat$items) else gamma_zero
    if (anyNA(idx)) stop("Unknown item(s) in 'gamma_zero'.")
    gz[idx] <- TRUE
  }
  no_rt <- colSums(dat$obsY) == 0L
  if (any(no_rt)) {
    message("Item(s) without observed response times (response model only): ",
            paste(dat$items[no_rt], collapse = ", "))
    gz[no_rt] <- TRUE
  }
  ref <- 1L
  if (asymmetric) {
    ref <- if (is.null(asym_ref)) which(!gz)[1L] else
      if (is.character(asym_ref)) match(asym_ref, dat$items) else asym_ref
    if (is.na(ref) || gz[ref]) stop("'asym_ref' must be an item with a free distance slope.")
  }
  spec <- list(error = error, asymmetric = asymmetric,
               delta_free = delta == "free", gamma_zero = gz, asym_ref = ref)
  theta <- seq(theta_range[1], theta_range[2], length.out = n_theta)
  lp <- stats::dnorm(theta, log = TRUE); lp <- lp - log(sum(exp(lp)))
  ctrl <- list(theta = theta, lprior = lp, gh = .gauss_hermite(n_tau),
               kappa = kappa, inflate = if (error == "t") 1.5 else 1,
               mstep_maxit = 50)

  sv <- .start_values(dat, spec, ctrl)
  if (!is.null(start)) sv <- .adopt_start(sv, start, spec)
  items <- sv$items; struct <- sv$struct

  ll_old <- -Inf; hist <- numeric(0); converged <- FALSE
  for (iter in seq_len(max_iter)) {
    E <- .estep(items, struct, dat, ctrl, spec)
    hist <- c(hist, E$ll)
    if (verbose) cat(sprintf("iter %3d  logLik = %.4f\n", iter, E$ll))
    if (is.finite(ll_old) && abs(E$ll - ll_old) < tol * (abs(E$ll) + 1)) {
      converged <- TRUE; break
    }
    ll_old <- E$ll
    for (j in seq_len(J)) items[[j]] <- .mstep_item(items[[j]], j, E, dat, ctrl, spec)
    struct <- .mstep_struct(E, ctrl, nrow(dat$X))
  }
  if (!converged) warning("EM did not converge in ", max_iter, " iterations.")
  E <- .estep(items, struct, dat, ctrl, spec)

  fit <- list(call = cl, items = items, struct = struct, spec = spec,
              ctrl = ctrl, item_names = dat$items, cat_values = dat$cat_values,
              logLik = E$ll, n = nrow(dat$X), n_obs_rt = sum(dat$obsY),
              iterations = iter, converged = converged, ll_history = hist,
              data = dat)
  fit$n_par <- length(.pack_all(items, struct, spec))
  fit$vcov <- NULL
  if (se) fit$vcov <- .vcov_xpd(fit, E)
  fit$estimates <- .natural(fit)
  fit$vcov_natural <- attr(fit$estimates, "vcov")
  attr(fit$estimates, "vcov") <- NULL
  fit$settings <- list(error = error, asymmetric = asymmetric, delta = delta,
                       n_theta = n_theta, theta_range = theta_range,
                       n_tau = n_tau, kappa = kappa, max_iter = max_iter,
                       tol = tol)
  class(fit) <- "grmrt"
  fit
}

# Data preparation --------------------------------------------------------------
.make_data <- function(responses, logrt) {
  X <- as.matrix(responses)
  if (is.null(logrt)) stop("'logrt' is required.")
  Y <- as.matrix(logrt)
  if (!all(dim(X) == dim(Y))) stop("'responses' and 'logrt' must have equal dimensions.")
  storage.mode(Y) <- "double"
  if (any(!is.finite(Y[!is.na(Y)]))) stop("'logrt' contains non-finite values.")
  items <- colnames(X); if (is.null(items)) items <- paste0("item", seq_len(ncol(X)))
  cat_values <- vector("list", ncol(X)); Xr <- matrix(NA_integer_, nrow(X), ncol(X))
  for (j in seq_len(ncol(X))) {
    v <- sort(unique(X[!is.na(X[, j]), j]))
    if (length(v) < 2L) stop("Item ", items[j], " has fewer than two observed categories.")
    cat_values[[j]] <- v
    Xr[, j] <- match(X[, j], v)
  }
  colnames(Xr) <- items
  list(X = Xr, Y = Y, obsY = !is.na(Y), K = vapply(cat_values, length, 1L),
       items = items, cat_values = cat_values)
}

.adopt_start <- function(sv, start, spec) {
  for (j in seq_along(sv$items)) {
    s <- start$items[[j]]; it <- sv$items[[j]]
    it$a <- s$a; it$b <- s$b; it$lambda <- s$lambda; it$sigma <- s$sigma
    if (isTRUE(it$delta_free)) it$delta <- if (!is.null(s$delta)) s$delta else mean(s$b)
    g <- s$gamma
    if (spec$asymmetric && length(g) == 1L) g <- c(g, g)
    if (!spec$asymmetric && length(g) == 2L) g <- mean(g)
    if (spec$asymmetric && j == spec$asym_ref) g <- rep(mean(g), 2L)
    it$gamma <- if (spec$gamma_zero[j]) 0 * g else g
    if (spec$error == "t") it$nu <- if (is.finite(s$nu)) s$nu else 10
    sv$items[[j]] <- it
  }
  sv$struct <- start$struct
  sv
}

# Full parameter vector ---------------------------------------------------------
.pack_all <- function(items, struct, spec) {
  v <- unlist(lapply(seq_along(items), function(j) .pack_item(items[[j]], spec, j)))
  if (spec$error == "t") v <- c(v, log(vapply(items[.has_rt(items)], function(it) it$nu - 2, 1)))
  c(v, struct$beta, log(struct$omega))
}

.has_rt <- function(items) vapply(items, function(it) isTRUE(it$has_rt), TRUE)

.item_par_len <- function(it, spec, j) length(.pack_item(it, spec, j))

.unpack_all <- function(v, items, struct, spec) {
  pos <- 1L
  for (j in seq_along(items)) {
    L <- .item_par_len(items[[j]], spec, j)
    items[[j]] <- .unpack_item(v[pos:(pos + L - 1L)], items[[j]], spec, j)
    pos <- pos + L
  }
  if (spec$error == "t") {
    hr <- which(.has_rt(items))
    for (q in seq_along(hr)) items[[hr[q]]]$nu <- 2 + exp(v[pos + q - 1L])
    pos <- pos + length(hr)
  }
  struct <- list(beta = v[pos], omega = exp(v[pos + 1L]))
  list(items = items, struct = struct)
}

# Natural-scale estimates (with delta-method SEs) --------------------------------
.natural_vec <- function(items, struct, spec) {
  out <- c(); nm <- c()
  shift <- 0
  if (spec$asymmetric) {
    free <- which(!spec$gamma_zero)
    dj <- vapply(items[free], function(it) (it$gamma[2L] - it$gamma[1L]) / 2, 1)
    shift <- mean(dj)
    struct$beta <- struct$beta + shift
  }
  for (j in seq_along(items)) {
    it <- items[[j]]
    g <- it$gamma
    if (spec$asymmetric && !spec$gamma_zero[j]) g <- g + c(shift, -shift)
    if (!it$has_rt) { g[] <- NA; it$lambda <- NA; it$sigma <- NA; it$nu <- NA }
    gn <- if (length(g) == 1L) "gamma" else c("gamma_lower", "gamma_upper")
    vals <- c(it$a, it$b, .item_delta(it, spec), it$lambda, g, it$sigma)
    nms <- c("a", paste0("b", seq_along(it$b)), "delta", "lambda", gn, "sigma")
    if (spec$error == "t") { vals <- c(vals, it$nu); nms <- c(nms, "nu") }
    out <- c(out, vals); nm <- c(nm, paste0(j, ":", nms))
  }
  sd_tau <- sqrt(struct$beta^2 + struct$omega)
  out <- c(out, sd_tau, struct$beta / sd_tau)
  names(out) <- c(nm, "sd_tau", "rho")
  out
}

.natural <- function(fit) {
  est <- .natural_vec(fit$items, fit$struct, fit$spec)
  se <- rep(NA_real_, length(est)); V <- NULL
  if (!is.null(fit$vcov)) {
    eta <- .pack_all(fit$items, fit$struct, fit$spec)
    f <- function(e) {
      u <- .unpack_all(e, fit$items, fit$struct, fit$spec)
      .natural_vec(u$items, u$struct, fit$spec)
    }
    Jm <- .num_jacobian(f, eta)
    V <- Jm %*% fit$vcov %*% t(Jm)
    se <- sqrt(pmax(diag(V), 0))
    se[!is.finite(diag(V)) | abs(diag(V)) < 1e-14] <- NA  # fixed parameters
  }
  item_id <- suppressWarnings(as.integer(sub(":.*", "", names(est))))
  par <- sub("^[0-9]+:", "", names(est))
  out <- data.frame(item = ifelse(is.na(item_id), NA, fit$item_names[item_id]),
                    parameter = par, estimate = unname(est), se = se,
                    stringsAsFactors = FALSE)
  attr(out, "vcov") <- V
  out
}

.num_jacobian <- function(f, x, h = 1e-6) {
  f0 <- f(x)
  Jm <- matrix(0, length(f0), length(x))
  for (p in seq_along(x)) {
    hp <- h * max(1, abs(x[p]))
    xp <- x; xm <- x; xp[p] <- x[p] + hp; xm[p] <- x[p] - hp
    Jm[, p] <- (f(xp) - f(xm)) / (2 * hp)
  }
  Jm
}
