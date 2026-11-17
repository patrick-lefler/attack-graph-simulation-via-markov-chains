# States, transition matrices, closed-form solution and timing results.
# Requires R/00_setup.R. Functions are defined in R/functions.R. The
# "## ---- label ----" markers let pages display each section with knitr::read_chunk().

## ---- model-states ----
transient_names <- params$states$transient
absorbing_names <- params$states$absorbing
all_state_names <- c(transient_names, absorbing_names)
t_dim <- length(transient_names)
r_dim <- length(absorbing_names)

# Mean step duration per transient stage (hours); the same for every stage by assumption
step_hours <- setNames(rep(params$step_hours$value, t_dim), transient_names)

# Baseline and improved-posture transition matrices P (9x9), illustrative values from params.yaml
P_baseline_raw <- matrix_from_params(params$transition_matrices$baseline, transient_names, absorbing_names)
P_hardened_raw <- matrix_from_params(params$transition_matrices$improved, transient_names, absorbing_names)

# Every row must be a probability distribution
stopifnot(all(abs(rowSums(P_baseline_raw) - 1) < 1e-8))
stopifnot(all(abs(rowSums(P_hardened_raw) - 1) < 1e-8))

Q_base <- P_baseline_raw[1:t_dim, 1:t_dim]
R_base <- P_baseline_raw[1:t_dim, (t_dim + 1):(t_dim + r_dim)]
Q_hard <- P_hardened_raw[1:t_dim, 1:t_dim]
R_hard <- P_hardened_raw[1:t_dim, (t_dim + 1):(t_dim + r_dim)]

## ---- model-solver ----
base_solution <- solve_absorbing_markov(Q_base, R_base, step_hours)
hard_solution <- solve_absorbing_markov(Q_hard, R_hard, step_hours)

## ---- model-jump-chain ----
J_base <- jump_chain(P_baseline_raw)
B_jump <- solve_absorbing_markov(
  J_base[1:t_dim, 1:t_dim], J_base[1:t_dim, (t_dim + 1):(t_dim + r_dim)], step_hours
)$B
stopifnot(max(abs(B_jump - base_solution$B)) < 1e-10)

# Depth diagnostics used in the text
holding_steps_b <- 1 / (1 - diag(Q_base))
p_leave_first_step <- 1 - diag(Q_base)
a1_share_b <- P_baseline_raw[1:t_dim, t_dim + 1] / (1 - diag(Q_base))
a1_share_h <- P_hardened_raw[1:t_dim, t_dim + 1] / (1 - diag(Q_hard))
stopifnot(round(range(a1_share_b), 2) == c(0.11, 0.29))

## ---- model-key-values ----
# One row per starting stage; prose numbers on every page are drawn from this table
res <- tibble(
  stage = transient_names,
  steps_b = base_solution$steps,
  hours_b = base_solution$hours,
  hours_h = hard_solution$hours,
  cont_b = base_solution$B[, 1],
  cont_h = hard_solution$B[, 1],
  cat_b = base_solution$B[, 2] + base_solution$B[, 3],
  cat_h = hard_solution$B[, 2] + hard_solution$B[, 3],
  uplift = cont_h - cont_b,
  t_cont_b = base_solution$time_given_outcome[, 1],
  t_cat_b = base_solution$time_given_catastrophe
)

# First stage where catastrophic outcomes are more likely than containment
cross_stage <- which(res$cat_b > res$cont_b)[1]

# Guards for numbers quoted as static text
stopifnot(
  round(res$hours_b[1], 1) == 14.2, round(res$hours_b[6], 1) == 4.2,
  round(100 * res$cat_b[3], 1) == 52.9, round(100 * res$cont_b[3], 1) == 47.1,
  round(100 * res$cont_b[1], 1) == 72.5, round(100 * res$cont_h[1], 1) == 97.3,
  round(100 * res$uplift[1], 1) == 24.8, cross_stage == 3,
  round(100 * res$cat_b[1], 1) == 27.5
)

# Guards for static numbers quoted in Parameterization and Assumptions
stopifnot(
  isTRUE(all.equal(range(diag(Q_base)), c(0.05, 0.15))),
  isTRUE(all.equal(range(1 - diag(Q_base)), c(0.85, 0.95))),
  round(P_baseline_raw[1, 2] / (1 - P_baseline_raw[1, 1]), 2) == 0.72
)

## ---- model-time-distribution ----
stopifnot(length(unique(step_hours)) == 1)

pmf_b <- outcome_pmf(Q_base, R_base, start = 1)
pmf_h <- outcome_pmf(Q_hard, R_hard, start = 1)
hours_grid <- seq_len(nrow(pmf_b)) * unname(step_hours[1])

cat_pmf_b <- pmf_b[, 2] + pmf_b[, 3]
cond_cat_b <- cat_pmf_b / sum(cat_pmf_b)
cond_cont_b <- pmf_b[, 1] / sum(pmf_b[, 1])
q_at <- function(cond, p) hours_grid[which(cumsum(cond) >= p - 1e-12)[1]]
q_cat_b <- map_dbl(c(0.10, 0.50, 0.90, 0.99), \(p) q_at(cond_cat_b, p))
q_cont_b <- map_dbl(c(0.50, 0.90), \(p) q_at(cond_cont_b, p))
floor_cat_h <- hours_grid[which(cat_pmf_b > 0)[1]]
mean_cat_b <- sum(hours_grid * cond_cat_b)
sd_cat_b <- sqrt(sum((hours_grid - mean_cat_b)^2 * cond_cat_b))
cov_cat_b <- sd_cat_b / mean_cat_b
p_cat_by_24h <- sum(cat_pmf_b[hours_grid <= 24])

# Consistency with the closed form, and guards for numbers quoted as static text
stopifnot(abs(mean_cat_b - res$t_cat_b[1]) < 1e-6)
stopifnot(abs(sum(pmf_b) - 1) < 1e-8)
stopifnot(floor_cat_h == 16, round(res$t_cat_b[1], 1) == 25.0)
