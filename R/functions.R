# Every function used by the analysis. Sourced by R/00_setup.R and by the tests.
# Functions read t_dim, all_state_names and step_hours from the calling
# environment, which R/01_model.R defines. The "## ---- label ----" markers let
# pages display each section with knitr::read_chunk().

## ---- fn-params ----
# Build a 9x9 transition matrix from the rows listed in params.yaml
matrix_from_params <- function(rows, transient, absorbing) {
  all_states <- c(transient, absorbing)
  stopifnot(length(rows) == length(transient), all(lengths(rows) == length(all_states)))
  P <- matrix(0, length(all_states), length(all_states), dimnames = list(all_states, all_states))
  P[transient, ] <- do.call(rbind, map(rows, unlist))
  P[absorbing, absorbing] <- diag(length(absorbing))
  P
}

## ---- fn-solver ----
# TRUE if every transient state has a path to an absorbing state
all_reach_absorbing <- function(Q, R) {
  reach <- rowSums(R) > 0
  repeat {
    updated <- reach | as.vector((Q > 0) %*% reach) > 0
    if (identical(updated, reach)) break
    reach <- updated
  }
  all(reach)
}

solve_absorbing_markov <- function(Q, R, h) {
  stopifnot(all_reach_absorbing(Q, R))

  N <- solve(diag(nrow(Q)) - Q, tol = 1e-14)
  dimnames(N) <- dimnames(Q)

  B <- N %*% R
  stopifnot(all(abs(rowSums(B) - 1) < 1e-8))

  # Conditional expected time to each outcome: N (h * B_k) / B_k
  time_given_outcome <- apply(B, 2, \(b) as.vector(N %*% (h * b)) / b)
  rownames(time_given_outcome) <- rownames(Q)

  cat_prob <- B[, 2] + B[, 3]
  time_given_catastrophe <- as.vector(N %*% (h * cat_prob)) / cat_prob

  list(
    N = N,
    B = B,
    steps = rowSums(N),
    hours = as.vector(N %*% h) |> setNames(rownames(Q)),
    time_given_outcome = time_given_outcome,
    time_given_catastrophe = setNames(time_given_catastrophe, rownames(Q))
  )
}

## ---- fn-jump-chain ----
# B depends only on the embedded jump chain: remove self-loops and renormalize
jump_chain <- function(P) {
  for (i in seq_len(t_dim)) {
    s <- P[i, i]
    P[i, i] <- 0
    P[i, ] <- P[i, ] / (1 - s)
  }
  P
}

## ---- fn-time-distribution ----
# Probability of absorption into each outcome at each step, by matrix powers
outcome_pmf <- function(Q, R, start = 1, max_steps = 200) {
  x <- replace(numeric(nrow(Q)), start, 1)
  out <- matrix(0, max_steps, ncol(R), dimnames = list(NULL, colnames(R)))
  for (s in seq_len(max_steps)) {
    out[s, ] <- as.vector(x %*% R)
    x <- as.vector(x %*% Q)
  }
  out
}

## ---- fn-one-row ----
harden_one_row <- function(P, i, delta = 0.05) {
  donor <- unname(which.max(replace(P[i, ], c(i, 7), 0)))
  P[i, donor] <- P[i, donor] - delta
  P[i, 7] <- P[i, 7] + delta
  list(P = P, donor = all_state_names[donor])
}

outcomes_from_s1 <- function(P) {
  B <- solve_absorbing_markov(P[1:t_dim, 1:t_dim], P[1:t_dim, 7:9], step_hours)$B
  c(cont = B[1, 1], cat = B[1, 2] + B[1, 3])
}

## ---- fn-identity ----
harden_one_row_prop <- function(P, i, rel = 0.10) {
  inc <- rel * P[i, t_dim + 1]
  others <- setdiff(which(P[i, ] > 0), c(i, t_dim + 1))
  P[i, others] <- P[i, others] - inc * P[i, others] / sum(P[i, others])
  P[i, t_dim + 1] <- P[i, t_dim + 1] + inc
  P
}

## ---- fn-dirichlet ----
draw_P <- function(P, alpha) {
  for (i in seq_len(t_dim)) {
    idx <- which(P[i, ] > 0)
    g <- rgamma(length(idx), shape = alpha * P[i, idx])
    P[i, idx] <- g / sum(g)
  }
  P
}

summarise_draw <- function(P) {
  sol <- solve_absorbing_markov(P[1:t_dim, 1:t_dim], P[1:t_dim, 7:9], step_hours)
  tibble(Stage = transient_names, cont = sol$B[, 1], cat = sol$B[, 2] + sol$B[, 3], hours = sol$hours)
}

## ---- fn-paired ----
# Common random numbers: the same uniforms drive both regimes' gamma draws
draw_pair <- function(Pb, Ph, alpha) {
  for (i in seq_len(t_dim)) {
    idx <- which(Pb[i, ] > 0 | Ph[i, ] > 0)
    u <- runif(length(idx))
    gb <- qgamma(u, shape = pmax(alpha * Pb[i, idx], 1e-8))
    gh <- qgamma(u, shape = pmax(alpha * Ph[i, idx], 1e-8))
    Pb[i, idx] <- gb / sum(gb)
    Ph[i, idx] <- gh / sum(gh)
  }
  list(b = Pb, h = Ph)
}

## ---- fn-ranking ----
shift_cont_from_s1 <- function(P, i, delta = 0.05) {
  donor <- unname(which.max(replace(P[i, ], c(i, t_dim + 1), 0)))
  if (P[i, donor] < delta) return(NA_real_)
  P2 <- P
  P2[i, donor] <- P2[i, donor] - delta
  P2[i, t_dim + 1] <- P2[i, t_dim + 1] + delta
  outcomes_from_s1(P2)["cont"] - outcomes_from_s1(P)["cont"]
}

## ---- fn-simulation ----
simulate_walk <- function(P, start) {
  s <- start
  steps <- 0
  while (s <= t_dim) {
    s <- sample.int(ncol(P), 1, prob = P[s, ])
    steps <- steps + 1
  }
  c(outcome = s - t_dim, steps = steps)
}
