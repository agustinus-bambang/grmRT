# Usage: Rscript 01_simulation.R A 200 <cores> 500,1000,2000   (Study 1)
#        Rscript 01_simulation.R B 100 <cores> 0,0.05,0.10       (Study 2)
# Simulation study for grmRT ---------------------------------------------------
# Study A: parameter recovery, SE calibration, Type I error and power of the
#          item-level Wald test, and SE reduction of joint trait scores.
#          N in {500, 1000, 2000}; J = 10; gamma = (0, 0, .05 x3, .10 x3, .15 x2).
# Study B: robustness to contaminated timings. N = 1000; contamination in
#          {0, .05, .10} (log-time shift SD = 2); normal vs Student-t errors.
suppressMessages(library(grmRT))
library(parallel)
args <- commandArgs(TRUE)
study <- args[1]; R <- as.integer(args[2]); cores <- as.integer(args[3])
gam <- c(0, 0, rep(.05, 3), rep(.10, 3), rep(.15, 2))
J <- 10

gen <- function(N, contamination = 0) {
  a <- runif(J, 1, 2.5); dl <- runif(J, -1.5, 1)
  simulate_grmrt(n_persons = N, n_items = J, a = a, delta = dl, gamma = gam,
                 sd_tau = .4, rho = -.1, contamination = contamination,
                 contam_shift = 2)
}

pull <- function(fit, sim, label) {
  e <- fit$estimates
  it <- e[!is.na(e$item), ]
  tr <- sim$truth
  truth <- c(a = list(tr$a), gamma = list(tr$gamma), lambda = list(tr$lambda),
             sigma = list(tr$sigma), delta = list(tr$delta))
  out <- do.call(rbind, lapply(names(truth), function(p) {
    r <- it[it$parameter == p, ]
    data.frame(par = p, item = seq_len(J), true = truth[[p]], est = r$estimate, se = r$se)
  }))
  pp <- e[is.na(e$item), ]
  out <- rbind(out, data.frame(par = c("sd_tau", "rho"), item = NA,
                               true = c(sim$sd_tau, sim$rho),
                               est = pp$estimate, se = pp$se))
  d <- item_diagnostics(fit)
  out$p <- NA; out$p[out$par == "gamma"] <- d$p
  out$label <- label; out$converged <- fit$converged
  out
}

one_A <- function(r, N) {
  set.seed(1e5 + 100 * N + r)
  sim <- gen(N)
  t0 <- proc.time()[3]
  fit <- grmrt_fit(sim$responses, sim$logrt)
  tm <- proc.time()[3] - t0
  res <- pull(fit, sim, "normal")
  sc <- score_persons(fit)
  res$N <- N; res$rep <- r; res$time <- tm
  attr(res, "score") <- data.frame(N = N, rep = r,
    se_red = mean(sc$se_reduction),
    rmse_joint = sqrt(mean((sc$theta - sim$theta)^2)),
    rmse_resp = sqrt(mean((sc$theta_resp - sim$theta)^2)),
    cor_joint = cor(sc$theta, sim$theta), cor_resp = cor(sc$theta_resp, sim$theta))
  res
}

one_B <- function(r, cont) {
  set.seed(2e5 + 1000 * cont + r)
  sim <- gen(1000, contamination = cont)
  fn <- grmrt_fit(sim$responses, sim$logrt)
  ft <- grmrt_fit(sim$responses, sim$logrt, error = "t", start = fn)
  res <- rbind(pull(fn, sim, "normal"), pull(ft, sim, "t"))
  res$cont <- cont; res$rep <- r
  res
}

if (study == "A") {
  Ns <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else c(500, 1000, 2000)
  for (N in Ns) {
    out <- mclapply(seq_len(R), one_A, N = N, mc.cores = cores)
    ok <- !vapply(out, inherits, TRUE, "try-error")
    saveRDS(list(res = do.call(rbind, out[ok]),
                 score = do.call(rbind, lapply(out[ok], attr, "score")),
                 failed = sum(!ok)), sprintf("simA_N%d.rds", N))
    cat("done N =", N, "failed", sum(!ok), "\n")
  }
} else {
  Cs <- if (length(args) >= 4) as.numeric(strsplit(args[4], ",")[[1]]) else c(0, .05, .10)
  for (cont in Cs) {
    out <- mclapply(seq_len(R), one_B, cont = cont, mc.cores = cores)
    ok <- !vapply(out, inherits, TRUE, "try-error")
    saveRDS(list(res = do.call(rbind, out[ok]), failed = sum(!ok)),
            sprintf("simB_c%02d.rds", round(100 * cont)))
    cat("done cont =", cont, "failed", sum(!ok), "\n")
  }
}
