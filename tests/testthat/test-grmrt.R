test_that("simulation returns consistent objects", {
  set.seed(1)
  s <- simulate_grmrt(n_persons = 50, n_items = 4)
  expect_equal(dim(s$responses), c(50, 4))
  expect_equal(dim(s$logrt), c(50, 4))
  expect_true(all(s$responses %in% 1:5))
})

test_that("estimation recovers the distance slope", {
  set.seed(10)
  s <- simulate_grmrt(n_persons = 1500, n_items = 6, gamma = 0.25)
  f <- grmrt_fit(s$responses, s$logrt)
  expect_true(f$converged)
  g <- f$estimates$estimate[f$estimates$parameter == "gamma"]
  expect_lt(abs(mean(g) - 0.25), 0.05)
  expect_true(all(is.finite(f$estimates$se[f$estimates$parameter == "gamma"])))
})

test_that("diagnostics flag an item without the effect", {
  set.seed(11)
  s <- simulate_grmrt(n_persons = 1500, n_items = 6,
                      gamma = c(0.3, 0.3, 0.3, 0.3, 0.3, 0))
  f <- grmrt_fit(s$responses, s$logrt)
  d <- item_diagnostics(f)
  expect_equal(nrow(d), 6)
  expect_true(all(d$flag[1:5] == "inverted-U"))
  expect_false(d$flag[6] == "inverted-U")
})

test_that("joint scoring reduces posterior SD", {
  set.seed(12)
  s <- simulate_grmrt(n_persons = 400, n_items = 6, gamma = 0.3)
  f <- grmrt_fit(s$responses, s$logrt, se = FALSE)
  sc <- score_persons(f)
  expect_equal(nrow(sc), 400)
  expect_gt(mean(sc$se_reduction), 0)
})

test_that("prep_rt reorders, reverses and screens", {
  df <- data.frame(a = c(1, 5, 3), b = c(2, 4, 5), a_ms = c(2000, 100, 3000),
                   b_ms = c(1500, 70000, 2500))
  d <- prep_rt(df, items = c("a", "b"), times = c("a_ms", "b_ms"),
               order = c("b", "a"), reverse = "a", drop_first = TRUE)
  expect_equal(colnames(d$responses), c("b", "a"))
  expect_equal(unname(d$responses[, "a"]), c(5, 1, 3))
  expect_true(all(is.na(d$logrt[, "b"])))
  expect_true(is.na(d$logrt[2, "a"]))
})

test_that("asymmetric and t-error fits run and anova works", {
  set.seed(13)
  s <- simulate_grmrt(n_persons = 600, n_items = 5, error = "t", nu = 4)
  f1 <- grmrt_fit(s$responses, s$logrt, se = FALSE)
  f2 <- grmrt_fit(s$responses, s$logrt, error = "t", se = FALSE, start = f1)
  a <- anova(f1, f2)
  expect_gt(a$LR[2], 0)
  f3 <- grmrt_fit(s$responses, s$logrt, asymmetric = TRUE, max_iter = 50)
  expect_true(all(c("gamma_lower", "gamma_upper") %in% f3$estimates$parameter))
  d <- item_diagnostics(f3)
  expect_true("asymmetry" %in% names(d))
})
