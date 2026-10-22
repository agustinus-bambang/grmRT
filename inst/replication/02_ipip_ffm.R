# Empirical application: IPIP-FFM (Open-Source Psychometrics Project) ---------
suppressMessages(library(grmRT))
# Data: IPIP-FFM-data-8Nov2018.zip from https://openpsychometrics.org/_rawdata/
# The analysis in the article used the first 20,000 rows of data-final.csv.
d <- read.delim("data-final.csv", nrows = 20000)
fac <- c("EXT", "AGR", "CSN", "EST", "OPN")
items <- as.vector(t(outer(fac, 1:10, paste0)))
times <- paste0(items, "_E")
pres <- as.vector(outer(fac, 1:10, paste0))       # EXT1, AGR1, CSN1, EST1, OPN1, EXT2, ...
rev <- c(paste0("EXT", c(2, 4, 6, 8, 10)), paste0("EST", c(2, 4)),
         paste0("AGR", c(1, 3, 5, 7)), paste0("CSN", c(2, 4, 6, 8)),
         paste0("OPN", c(2, 4, 6)))

n0 <- nrow(d)
d <- d[d$IPC == 1, ]; n1 <- nrow(d)
d <- d[rowSums(d[, items] == 0 | is.na(d[, items])) == 0, ]; n2 <- nrow(d)
d <- d[d$testelapse >= 60 & d$testelapse <= 3600, ]; n3 <- nrow(d)
dat <- prep_rt(d, items = items, times = times, order = pres, time_unit = "ms",
               reverse = rev, scale_range = c(1, 5), min_time = 0.5, max_time = 60,
               drop_first = TRUE)
print(dat)
flow <- c(raw = n0, ipc1 = n1, complete = n2, time_window = n3)
print(flow)

res <- list(flow = flow, screening = dat$screening)
for (f in fac) {
  cols <- grep(paste0("^", f), colnames(dat$responses))
  X <- dat$responses[, cols]; Y <- dat$logrt[, cols]
  cat("\n=====", f, "=====\n")
  t0 <- proc.time()[3]
  fn <- grmrt_fit(X, Y, se = TRUE)
  tn <- proc.time()[3] - t0
  t0 <- proc.time()[3]
  ft <- grmrt_fit(X, Y, error = "t", se = TRUE, start = fn)
  tt <- proc.time()[3] - t0
  fa <- grmrt_fit(X, Y, error = "t", asymmetric = TRUE, se = TRUE, start = ft)
  f0 <- grmrt_fit(X, Y, error = "t", gamma_zero = colnames(X), se = FALSE, start = ft)
  print(ft)
  dg <- item_diagnostics(ft)
  da <- item_diagnostics(fa)
  sc <- score_persons(ft)
  print(dg)
  cat("normal vs t:\n"); print(anova(fn, ft))
  cat("sym vs asym (t):\n"); print(anova(ft, fa))
  cat("all gamma = 0 vs free (t):\n"); print(anova(f0, ft))
  cat("time (s): normal", tn, " t", tt, "\n")
  res[[f]] <- list(fn = fn[c("estimates", "logLik", "n_par", "iterations", "converged")],
                   ft = ft[c("estimates", "logLik", "n_par", "iterations", "converged")],
                   fa = fa[c("estimates", "logLik", "n_par", "iterations", "converged")],
                   f0 = f0[c("logLik", "n_par")], diag = dg, diag_asym = da,
                   score = summary(sc$se_reduction), score_mean = mean(sc$se_reduction),
                   cor_joint_resp = cor(sc$theta, sc$theta_resp),
                   time = c(normal = tn, t = tt),
                   n = ft$n, n_rt = ft$n_obs_rt)
  saveRDS(res, "apply_results.rds")
  saveRDS(ft, sprintf("fit_t_%s.rds", f))
  png(sprintf("fig_%s.png", f), 1600, 1200, res = 200)
  plot_time_distance(ft, bin_width = 0.5, min_n = 100)
  dev.off()
  assign(paste0("fit_", f), ft)
}
save(list = paste0("fit_", fac), file = "fits_t.RData")
cat("ALL DONE\n")
