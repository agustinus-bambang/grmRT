# Estimation engine (internal) ------------------------------------------------
#
# Model for person i and item j:
#   P(X_ij >= k | theta_i) = plogis(a_j (theta_i - b_jk))                 (GRM)
#   log T_ij = lambda_j - tau_i - gamma_j h(theta_i - delta_j) + e_ij      (RT)
#   e_ij ~ N(0, sigma_j^2)  or  sigma_j * t_{nu_j}
#   theta_i ~ N(0, 1);  tau_i | theta_i ~ N(beta * theta_i, omega)
# delta_j is tied to the GRM location mean(b_j) by default.
#
# Marginal likelihood: rectangular grid over theta; for each person and theta
# node, tau is integrated with adaptive Gauss-Hermite quadrature centred on the
# Gaussian-conjugate posterior of tau (exact for normal errors).

.item_delta <- function(it, spec) if (isTRUE(it$delta_free)) it$delta else mean(it$b)

# number of free distance-slope parameters for item j
.n_gamma <- function(spec, j) {
  if (spec$gamma_zero[j]) return(0L)
  if (spec$asymmetric && j != spec$asym_ref) 2L else 1L
}

.pack_item <- function(it, spec, j) {
  v <- c(log(it$a), it$b[1L])
  if (length(it$b) > 1L) v <- c(v, log(diff(it$b)))
  if (isTRUE(it$delta_free)) v <- c(v, it$delta)
  if (!it$has_rt) return(v)
  v <- c(v, it$lambda)
  ng <- .n_gamma(spec, j)
  if (ng == 1L) v <- c(v, it$gamma[1L]) else if (ng == 2L) v <- c(v, it$gamma)
  c(v, log(it$sigma))
}

.unpack_item <- function(v, it, spec, j) {
  K1 <- length(it$b)
  pos <- 1L
  it$a <- exp(v[pos]); pos <- pos + 1L
  b <- v[pos]; pos <- pos + 1L
  if (K1 > 1L) {
    b <- cumsum(c(b, exp(v[pos:(pos + K1 - 2L)])))
    pos <- pos + K1 - 1L
  }
  it$b <- b
  if (isTRUE(it$delta_free)) { it$delta <- v[pos]; pos <- pos + 1L }
  if (!it$has_rt) return(it)
  it$lambda <- v[pos]; pos <- pos + 1L
  ng <- .n_gamma(spec, j)
  len <- if (spec$asymmetric) 2L else 1L
  if (ng == 0L) it$gamma <- rep(0, len)
  else if (ng == 1L) it$gamma <- rep(v[pos], len)
  else it$gamma <- v[pos:(pos + 1L)]
  pos <- pos + ng
  it$sigma <- exp(v[pos])
  it
}

# Response (GRM) log-likelihood per item at theta nodes: N x Qt
.resp_ll_item <- function(it, x, theta) {
  P <- .grm_probs(it$a, it$b, theta)
  out <- matrix(0, length(x), length(theta))
  ok <- !is.na(x)
  out[ok, ] <- t(log(P[, x[ok], drop = FALSE]))
  out
}

# RT residual array r[i, t, k] for item j given nodes Tau (N x Qt x Kq)
.rt_resid <- function(it, y, Tau, theta, spec, kappa) {
  ct <- -it$lambda + .dterm(theta, .item_delta(it, spec), it$gamma, kappa)
  r <- sweep(Tau, 2L, ct, "+")
  y0 <- y; y0[is.na(y0)] <- 0
  r + y0
}

# Proposal (Gaussian-conjugate posterior of tau) for each person and theta node
.proposal <- function(items, struct, dat, ctrl, spec) {
  J <- length(items); theta <- ctrl$theta
  sig2 <- vapply(items, function(it) it$sigma^2, 1)
  if (spec$error == "t") {
    nu <- vapply(items, function(it) it$nu, 1)
    sig2 <- sig2 * nu / (nu - 2)   # match the t error variance
  }
  lam <- vapply(items, function(it) it$lambda, 1)
  O <- dat$obsY * 1
  Y0 <- dat$Y; Y0[is.na(Y0)] <- 0
  A <- 1 / struct$omega + as.vector(O %*% (1 / sig2))
  base <- as.vector(((Y0 - rep(lam, each = nrow(Y0))) * O) %*% (1 / sig2))
  G <- vapply(seq_len(J), function(j) {
    it <- items[[j]]
    .dterm(theta, .item_delta(it, spec), it$gamma, ctrl$kappa) / sig2[j]
  }, numeric(length(theta)))            # Qt x J
  G <- matrix(G, nrow = length(theta))
  Bpart <- base + O %*% t(G)            # N x Qt  (sum_j e_ijt / s_j^2)
  B <- Bpart - outer(rep(1, nrow(O)), struct$beta * theta / struct$omega)
  list(mu = -B / A, s = ctrl$inflate / sqrt(A))
}

# E-step: marginal log-likelihood and posterior over (theta node, tau node)
.estep <- function(items, struct, dat, ctrl, spec, prop = NULL,
                   keep_item_ll = FALSE) {
  N <- nrow(dat$Y); J <- length(items)
  theta <- ctrl$theta; Qt <- length(theta); Kq <- length(ctrl$gh$x)
  if (is.null(prop)) prop <- .proposal(items, struct, dat, ctrl, spec)
  # tau nodes
  Tau <- array(prop$mu, c(N, Qt, Kq)) +
    array(rep(sqrt(2) * prop$s, Qt * Kq) * rep(ctrl$gh$x, each = N * Qt),
          c(N, Qt, Kq))
  # log prior of tau | theta minus log proposal, plus GH weight
  mprior <- sweep(Tau, 2L, struct$beta * theta, "-")
  lw <- -0.5 * log(2 * pi * struct$omega) - 0.5 * mprior^2 / struct$omega
  lq <- -log(prop$s) - 0.5 * log(2 * pi)  # N vector (s depends on i only)
  lw <- lw - lq
  lw <- sweep(lw, 3L, log(ctrl$gh$w / sqrt(pi)) + ctrl$gh$x^2, "+")
  # response part
  LR <- matrix(0, N, Qt)
  item_resp <- item_rt <- NULL
  if (keep_item_ll) { item_resp <- vector("list", J); item_rt <- vector("list", J) }
  for (j in seq_len(J)) {
    lr <- .resp_ll_item(items[[j]], dat$X[, j], theta)
    LR <- LR + lr
    if (keep_item_ll) item_resp[[j]] <- lr
  }
  # RT part
  LT <- array(0, c(N, Qt, Kq))
  nu_inf <- spec$error == "normal"
  for (j in seq_len(J)) {
    it <- items[[j]]
    r <- .rt_resid(it, dat$Y[, j], Tau, theta, spec, ctrl$kappa)
    ld <- .ldens_err(r, it$sigma, if (nu_inf) Inf else it$nu)
    ld[!dat$obsY[, j], , ] <- 0
    LT <- LT + ld
    if (keep_item_ll) item_rt[[j]] <- ld
  }
  lj <- sweep(LT + lw, c(1L, 2L), LR, "+")
  lj <- sweep(lj, 2L, ctrl$lprior, "+")
  dim(lj) <- c(N, Qt * Kq)
  lli <- .logsumexp_rows(lj)
  w <- exp(lj - lli)
  dim(w) <- c(N, Qt, Kq)
  list(ll = sum(lli), lli = lli, w = w, Tau = Tau, prop = prop, LR = LR,
       lw = lw, LT = LT, item_resp = item_resp, item_rt = item_rt)
}

# Expected complete-data objective for one item (negative, for optim)
.item_obj <- function(v, it, spec, j, stats, ctrl) {
  it <- .unpack_item(v, it, spec, j)
  P <- .grm_probs(it$a, it$b, ctrl$theta)
  f <- sum(stats$R * log(P))
  if (stats$nobs > 0) {
    ct <- -it$lambda + .dterm(ctrl$theta, .item_delta(it, spec), it$gamma,
                              ctrl$kappa)
    q <- sum(stats$S2 + 2 * ct * stats$S1 + ct^2 * stats$S0)
    f <- f - stats$nobs * log(it$sigma) - q / (2 * it$sigma^2)
  }
  -f
}

.mstep_item <- function(it, j, E, dat, ctrl, spec) {
  x <- dat$X[, j]; y <- dat$Y[, j]
  Wt <- rowSums(E$w, dims = 2L)                          # N x Qt
  okx <- !is.na(x)
  R <- crossprod(Wt[okx, , drop = FALSE],
                 outer(x[okx], seq_len(it$K), "==") * 1)  # Qt x K
  oky <- dat$obsY[, j]
  stats <- list(R = R, nobs = sum(oky))
  nu_new <- it$nu
  if (stats$nobs > 0) {
    w <- E$w[oky, , , drop = FALSE]
    z <- E$Tau[oky, , , drop = FALSE] + y[oky]
    if (spec$error == "t") {
      r <- .rt_resid(it, y[oky], E$Tau[oky, , , drop = FALSE], ctrl$theta,
                     spec, ctrl$kappa)
      d2 <- (r / it$sigma)^2
      u <- (it$nu + 1) / (it$nu + d2)
      Elogu <- digamma((it$nu + 1) / 2) - log((it$nu + d2) / 2)
      Snu <- sum(w * (Elogu - u))
      fnu <- function(ln) {
        nu <- 2 + exp(ln)
        -(stats$nobs * ((nu / 2) * log(nu / 2) - lgamma(nu / 2)) + (nu / 2) * Snu)
      }
      o <- stats::optimize(fnu, c(log(0.05), log(500)))
      nu_new <- 2 + exp(o$minimum)
      psi <- w * u
    } else psi <- w
    stats$S0 <- rowSums(colSums(psi, dims = 1L))
    stats$S1 <- rowSums(colSums(psi * z, dims = 1L))
    stats$S2 <- rowSums(colSums(psi * z^2, dims = 1L))
  }
  v0 <- .pack_item(it, spec, j)
  o <- stats::optim(v0, .item_obj, it = it, spec = spec, j = j, stats = stats,
                    ctrl = ctrl, method = "BFGS",
                    control = list(maxit = ctrl$mstep_maxit, reltol = 1e-10))
  it <- .unpack_item(o$par, it, spec, j)
  it$nu <- nu_new
  if (spec$error == "t" && stats$nobs > 0) {
    # ECME step: maximize the expected t log-density directly over (sigma, nu),
    # integrating over (theta, tau) only. Much faster than the scale-mixture
    # update when nu is large (near-normal errors).
    r <- .rt_resid(it, y[oky], E$Tau[oky, , , drop = FALSE], ctrl$theta,
                   spec, ctrl$kappa)
    f2 <- function(p) -sum(w * .ldens_err(r, exp(p[1]), 2 + exp(p[2])))
    p0 <- c(log(it$sigma), log(it$nu - 2))
    o2 <- stats::optim(p0, f2, method = "L-BFGS-B",
                       lower = c(log(1e-3), log(0.05)), upper = c(log(1e3), log(500)))
    if (o2$value <= f2(p0)) { it$sigma <- exp(o2$par[1]); it$nu <- 2 + exp(o2$par[2]) }
  }
  it
}

.mstep_struct <- function(E, ctrl, N) {
  th <- array(rep(ctrl$theta, each = nrow(E$w)), dim(E$w))
  sw2 <- sum(E$w * th^2)
  beta <- sum(E$w * th * E$Tau) / sw2
  omega <- sum(E$w * (E$Tau - beta * th)^2) / N
  list(beta = beta, omega = max(omega, 1e-8))
}

.start_values <- function(dat, spec, ctrl) {
  X <- dat$X; Y <- dat$Y; J <- ncol(X); N <- nrow(X)
  ss <- rowMeans(sapply(seq_len(J), function(j) {
    x <- X[, j]; (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }), na.rm = TRUE)
  th0 <- stats::qnorm((rank(ss) - 0.5) / N)
  lam <- colMeans(Y, na.rm = TRUE)
  tau0 <- rowMeans(-(Y - rep(lam, each = N)), na.rm = TRUE)
  tau0[is.na(tau0)] <- 0
  items <- lapply(seq_len(J), function(j) {
    K <- dat$K[j]; x <- X[, j]
    pk <- vapply(2:K, function(k) mean(x >= k, na.rm = TRUE), 1)
    pk <- pmin(pmax(pk, 0.01), 0.99)
    b <- sort(-stats::qlogis(pk) * 1.3)
    if (any(diff(b) <= 0)) b <- b + seq(0, 0.01 * (length(b) - 1), length.out = length(b))
    has_rt <- any(dat$obsY[, j])
    res <- Y[, j] - lam[j] + tau0
    it <- list(K = K, a = 1.2, b = b, lambda = if (has_rt) lam[j] else 0,
               has_rt = has_rt,
               gamma = if (spec$gamma_zero[j]) rep(0, 1 + spec$asymmetric)
                       else rep(0.05, 1 + spec$asymmetric),
               sigma = if (sum(dat$obsY[, j]) > 1L) max(stats::sd(res, na.rm = TRUE), 0.05) else 1,
               nu = if (spec$error == "t") 10 else Inf)
    it$delta_free <- spec$delta_free && !spec$gamma_zero[j]
    if (it$delta_free) it$delta <- mean(b)
    it
  })
  beta <- stats::cov(th0, tau0)
  omega <- max(stats::var(tau0) - beta^2, 0.01)
  list(items = items, struct = list(beta = beta, omega = omega))
}
