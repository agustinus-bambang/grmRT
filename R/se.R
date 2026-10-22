# Standard errors: cross-product (XPD) of person-level score vectors ----------
# Scores are central finite differences of each person's marginal log-likelihood,
# evaluated with the adaptive quadrature nodes held fixed at the MLE.

.item_contrib <- function(it, j, dat, Tau, ctrl, spec) {
  resp <- .resp_ll_item(it, dat$X[, j], ctrl$theta)
  r <- .rt_resid(it, dat$Y[, j], Tau, ctrl$theta, spec, ctrl$kappa)
  ld <- .ldens_err(r, it$sigma, if (spec$error == "t") it$nu else Inf)
  ld[!dat$obsY[, j], , ] <- 0
  sweep(ld, c(1L, 2L), resp, "+")
}

.prior_part <- function(struct, Tau, theta) {
  m <- sweep(Tau, 2L, struct$beta * theta, "-")
  -0.5 * log(2 * pi * struct$omega) - 0.5 * m^2 / struct$omega
}

.vcov_xpd <- function(fit, E, h = 1e-4) {
  items <- fit$items; struct <- fit$struct; spec <- fit$spec
  ctrl <- fit$ctrl; dat <- fit$data
  Tau <- E$Tau; d <- dim(Tau); N <- d[1]
  lse <- function(arr) { dim(arr) <- c(N, d[2] * d[3]); .logsumexp_rows(arr) }
  base <- sweep(E$LT + E$lw, c(1L, 2L), E$LR, "+")
  base <- sweep(base, 2L, ctrl$lprior, "+")
  eta <- .pack_all(items, struct, spec)
  S <- matrix(0, N, length(eta))
  pos <- 0L
  J <- length(items)
  for (j in seq_len(J)) {
    v0 <- .pack_item(items[[j]], spec, j)
    rest <- base - .item_contrib(items[[j]], j, dat, Tau, ctrl, spec)
    for (p in seq_along(v0)) {
      hp <- h * max(1, abs(v0[p]))
      vp <- v0; vm <- v0; vp[p] <- v0[p] + hp; vm[p] <- v0[p] - hp
      ip <- .unpack_item(vp, items[[j]], spec, j)
      im <- .unpack_item(vm, items[[j]], spec, j)
      S[, pos + p] <- (lse(rest + .item_contrib(ip, j, dat, Tau, ctrl, spec)) -
                       lse(rest + .item_contrib(im, j, dat, Tau, ctrl, spec))) / (2 * hp)
    }
    pos <- pos + length(v0)
  }
  if (spec$error == "t") {
    for (j in which(.has_rt(items))) {
      rest <- base - .item_contrib(items[[j]], j, dat, Tau, ctrl, spec)
      l0 <- log(items[[j]]$nu - 2); hp <- h * max(1, abs(l0))
      ip <- items[[j]]; im <- items[[j]]
      ip$nu <- 2 + exp(l0 + hp); im$nu <- 2 + exp(l0 - hp)
      pos <- pos + 1L
      S[, pos] <- (lse(rest + .item_contrib(ip, j, dat, Tau, ctrl, spec)) -
                   lse(rest + .item_contrib(im, j, dat, Tau, ctrl, spec))) / (2 * hp)
    }
  }
  rest <- base - .prior_part(struct, Tau, ctrl$theta)
  s0 <- c(struct$beta, log(struct$omega))
  for (p in 1:2) {
    hp <- h * max(1, abs(s0[p]))
    sp <- s0; sm <- s0; sp[p] <- s0[p] + hp; sm[p] <- s0[p] - hp
    stp <- list(beta = sp[1], omega = exp(sp[2]))
    stm <- list(beta = sm[1], omega = exp(sm[2]))
    S[, pos + p] <- (lse(rest + .prior_part(stp, Tau, ctrl$theta)) -
                     lse(rest + .prior_part(stm, Tau, ctrl$theta))) / (2 * hp)
  }
  info <- crossprod(S)
  V <- tryCatch(solve(info), error = function(e) {
    warning("Information matrix is singular; standard errors unavailable.")
    NULL
  })
  V
}
