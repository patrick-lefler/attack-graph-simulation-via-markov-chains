# A two-stage chain small enough to solve by hand:
#   T1: stay 0.2, move to T2 0.3, contained 0.5
#   T2: stay 0.5, contained 0.25, exfiltration 0.25
# N = (I - Q)^-1 = [[1.25, 0.75], [0, 2]]
toy_Q <- matrix(c(0.2, 0.3,
                  0.0, 0.5), 2, byrow = TRUE, dimnames = list(c("T1", "T2"), c("T1", "T2")))
toy_R <- matrix(c(0.50, 0.00, 0,
                  0.25, 0.25, 0), 2, byrow = TRUE)

test_that("the solver reproduces hand-computed N, B and expected steps", {
  sol <- solve_absorbing_markov(toy_Q, toy_R, h = c(1, 1))
  expect_equal(unname(sol$N), matrix(c(1.25, 0.75, 0, 2), 2, byrow = TRUE))
  expect_equal(unname(sol$B[, 1:2]), matrix(c(0.8125, 0.1875, 0.5, 0.5), 2, byrow = TRUE))
  expect_equal(unname(sol$steps), c(2, 2))
})

test_that("conditional time to catastrophe matches the hand calculation", {
  sol <- solve_absorbing_markov(toy_Q, toy_R, h = c(1, 1))
  # From T1: (1.25 * 0.1875 + 0.75 * 0.5) / 0.1875 = 3.25; from T2 the holding time is geometric, mean 2
  expect_equal(unname(sol$time_given_catastrophe), c(3.25, 2))
})

test_that("hours scale with the step durations h", {
  sol <- solve_absorbing_markov(toy_Q, toy_R, h = c(4, 4))
  expect_equal(unname(sol$hours), 4 * unname(sol$steps))
})

test_that("the solver refuses a chain with a trapped transient state", {
  Q <- matrix(c(0.5, 0.5,
                0.0, 1.0), 2, byrow = TRUE)
  R <- matrix(c(0, 0, 0,
                0, 0, 0), 2, byrow = TRUE)
  expect_false(all_reach_absorbing(Q, R))
  expect_error(solve_absorbing_markov(Q, R, h = c(1, 1)))
})
