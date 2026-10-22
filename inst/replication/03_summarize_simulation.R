fA <- list.files(pattern = "^simA_N.*rds$")
outA <- do.call(rbind, lapply(fA, function(f) {
  x <- readRDS(f); r <- x$res; N <- r$N[1]
  s <- do.call(rbind, lapply(split(r, r$par), function(d) data.frame(
    N = N, par = d$par[1], bias = mean(d$est - d$true), rmse = sqrt(mean((d$est - d$true)^2)),
    se_ratio = mean(d$se, na.rm = TRUE) / sd(d$est - d$true),
    coverage = mean(abs(d$est - d$true) <= 1.96 * d$se, na.rm = TRUE))))
  s
}))
print(outA, digits = 3, row.names = FALSE)
pw <- do.call(rbind, lapply(fA, function(f) {
  r <- readRDS(f)$res; g <- r[r$par == "gamma", ]
  a <- aggregate(cbind(reject = g$p < .05) ~ true, data = g, FUN = mean); a$N <- r$N[1]; a
}))
print(reshape(pw, idvar = "N", timevar = "true", direction = "wide"), digits = 3, row.names = FALSE)
sc <- do.call(rbind, lapply(fA, function(f) readRDS(f)$score))
print(aggregate(cbind(se_red, rmse_joint, rmse_resp, cor_joint, cor_resp) ~ N, data = sc, FUN = mean), digits = 4)
tm <- do.call(rbind, lapply(fA, function(f) { r <- readRDS(f)$res; data.frame(N = r$N[1], time = median(r$time[!duplicated(r$rep)]), conv = mean(r$converged[!duplicated(r$rep)])) }))
print(tm)
fB <- list.files(pattern = "^simB_c.*rds$")
if (length(fB)) {
  outB <- do.call(rbind, lapply(fB, function(f) {
    r <- readRDS(f)$res; g <- r[r$par == "gamma", ]
    do.call(rbind, lapply(split(g, g$label), function(d) data.frame(
      cont = d$cont[1], model = d$label[1], R = length(unique(d$rep)),
      bias = mean(d$est - d$true), rmse = sqrt(mean((d$est - d$true)^2)),
      coverage = mean(abs(d$est - d$true) <= 1.96 * d$se),
      typeI = mean(d$p[d$true == 0] < .05),
      power05 = mean(d$p[d$true == .05] < .05), power10 = mean(d$p[d$true == .10] < .05),
      power15 = mean(d$p[d$true == .15] < .05))))
  }))
  print(outB, digits = 3, row.names = FALSE)
  saveRDS(outB, "summaryB.rds")
}
saveRDS(list(A = outA, power = pw, score = sc, time = tm), "summaryA.rds")
