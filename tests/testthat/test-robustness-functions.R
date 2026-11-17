test_that("the standardized shift moves exactly delta to A1 and keeps the row summing to 1", {
  for (i in seq_len(t_dim)) {
    out <- harden_one_row(P_baseline_raw, i, delta = 0.05)
    expect_equal(sum(out$P[i, ]), 1, tolerance = 1e-12)
    expect_equal(out$P[i, t_dim + 1] - P_baseline_raw[i, t_dim + 1], 0.05, tolerance = 1e-12)
    expect_equal(sum(out$P != P_baseline_raw), 2)
  }
})

test_that("the proportional change raises A1 by the stated share and keeps zeros at zero", {
  for (i in seq_len(t_dim)) {
    P2 <- harden_one_row_prop(P_baseline_raw, i, rel = 0.10)
    expect_equal(sum(P2[i, ]), 1, tolerance = 1e-12)
    expect_equal(P2[i, t_dim + 1], 1.10 * P_baseline_raw[i, t_dim + 1], tolerance = 1e-12)
    expect_true(all(P2[i, P_baseline_raw[i, ] == 0] == 0))
  }
})

test_that("Dirichlet draws keep structural zeros and valid rows", {
  set.seed(1)
  P <- draw_P(P_baseline_raw, alpha = 50)
  expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-12)
  expect_true(all(P[P_baseline_raw == 0] == 0))
  expect_equal(P[absorbing_names, ], P_baseline_raw[absorbing_names, ])
})

test_that("a higher concentration keeps draws closer to the elicited matrix", {
  set.seed(1)
  spread <- \(a) mean(map_dbl(1:200, \(k) max(abs(draw_P(P_baseline_raw, a) - P_baseline_raw))))
  expect_lt(spread(200), spread(10))
})

test_that("paired draws use common random numbers", {
  set.seed(1)
  pr <- draw_pair(P_baseline_raw, P_baseline_raw, alpha = 50)
  expect_equal(pr$b, pr$h)
})

test_that("containment from S1 rises when any row is shifted toward A1", {
  base <- outcomes_from_s1(P_baseline_raw)
  for (i in seq_len(t_dim)) {
    expect_gt(shift_cont_from_s1(P_baseline_raw, i), 0)
  }
  expect_equal(unname(base["cont"] + base["cat"]), 1, tolerance = 1e-12)
})

test_that("a simulated walk always ends in an absorbing outcome", {
  set.seed(1)
  walks <- map(1:500, \(k) simulate_walk(P_baseline_raw, 1))
  expect_true(all(map_dbl(walks, "outcome") %in% 1:3))
  expect_true(all(map_dbl(walks, "steps") >= 1))
})
