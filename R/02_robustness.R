# One-row sensitivity, Dirichlet parameter uncertainty and the simulation check.
# Requires R/00_setup.R and R/01_model.R. Functions are defined in R/functions.R.
# Every random block sets its own seed, so results do not depend on the order in
# which sections run.

## ---- rob-one-row ----
shift_delta <- params$sensitivity$shift_delta$value
prop_rel <- params$sensitivity$prop_rel$value
seed <- params$seed$value

base_s1 <- outcomes_from_s1(P_baseline_raw)

sens_df <- map(seq_len(t_dim), \(i) {
  shifted <- harden_one_row(P_baseline_raw, i, delta = shift_delta)
  full_row <- P_baseline_raw
  full_row[i, ] <- P_hardened_raw[i, ]

  tibble(
    Stage = transient_names[i],
    Donor_Edge = str_extract(shifted$donor, "^[SA][0-9]"),
    Shift_Cont = outcomes_from_s1(shifted$P)["cont"] - base_s1["cont"],
    Shift_Cat = outcomes_from_s1(shifted$P)["cat"] - base_s1["cat"],
    Row_Cont = outcomes_from_s1(full_row)["cont"] - base_s1["cont"],
    Row_Cat = outcomes_from_s1(full_row)["cat"] - base_s1["cat"]
  )
}) |>
  list_rbind()

best_row <- which.max(sens_df$Shift_Cont)

## ---- rob-identity ----
B_full <- rbind(base_solution$B, diag(r_dim))
rownames(B_full) <- all_state_names
donor_idx <- map_int(seq_len(t_dim), \(i) unname(which.max(replace(P_baseline_raw[i, ], c(i, t_dim + 1), 0))))

identity_df <- tibble(
  Stage = transient_names,
  Donor_Edge = sens_df$Donor_Edge,
  Visits = base_solution$N[1, ],
  Donor_Cat_Share = 1 - B_full[donor_idx, 1],
  Identity_Pred = shift_delta * Visits * Donor_Cat_Share,
  Shift_Cont = sens_df$Shift_Cont,
  Prop_Cat = map_dbl(seq_len(t_dim), \(i) {
    outcomes_from_s1(harden_one_row_prop(P_baseline_raw, i, rel = prop_rel))["cat"] - base_s1["cat"]
  })
)
stopifnot(max(abs(identity_df$Identity_Pred - identity_df$Shift_Cont)) < 5e-4)

ratio_shift <- max(identity_df$Shift_Cont) / identity_df$Shift_Cont[t_dim]
ratio_prop <- min(identity_df$Prop_Cat) / identity_df$Prop_Cat[t_dim]

## ---- rob-dirichlet ----
alpha_conc <- params$dirichlet$alpha_conc$value
n_draws <- params$dirichlet$n_draws$value

set.seed(seed)
draws <- map(seq_len(n_draws), \(d) {
  b <- summarise_draw(draw_P(P_baseline_raw, alpha_conc))
  h <- summarise_draw(draw_P(P_hardened_raw, alpha_conc))
  tibble(draw = d, Stage = b$Stage,
         cont_b = b$cont, cat_b = b$cat, hours_b = b$hours,
         cont_h = h$cont, cat_h = h$cat, hours_h = h$hours)
}) |>
  list_rbind() |>
  mutate(Stage = factor(Stage, levels = transient_names))

q05 <- \(x) quantile(x, 0.05)
q95 <- \(x) quantile(x, 0.95)

unc_summary <- draws |>
  summarise(
    cont_b_lo = q05(cont_b), cont_b_hi = q95(cont_b),
    cat_b_lo = q05(cat_b), cat_b_hi = q95(cat_b),
    hours_b_lo = q05(hours_b), hours_b_hi = q95(hours_b),
    p_improved_better = mean(cont_h > cont_b),
    .by = Stage
  ) |>
  left_join(res |> select(Stage = stage, cont_b, cat_b, hours_b) |> mutate(Stage = factor(Stage, levels = transient_names)),
            by = "Stage") |>
  arrange(Stage)

crossing_dist <- draws |>
  summarise(
    cross = which(cat_b > cont_b)[1],
    .by = draw
  ) |>
  count(cross) |>
  mutate(share = n / sum(n))

p_cross_s3 <- crossing_dist |> filter(cross == 3) |> pull(share)

## ---- rob-alpha-sweep ----
alpha_grid <- as.numeric(unlist(params$dirichlet$alpha_grid$value))
n_sweep <- params$dirichlet$n_sweep$value
set.seed(seed)
alpha_sweep <- map(alpha_grid, \(a) {
  map(seq_len(n_sweep), \(d) {
    B <- summarise_draw(draw_P(P_baseline_raw, a))
    tibble(alpha = a, cross = coalesce(which(B$cat > B$cont)[1], 0L), cat_s3 = B$cat[3])
  }) |> list_rbind()
}) |>
  list_rbind() |>
  summarise(
    p_cross_s3 = mean(cross == 3),
    cat_s3_lo = quantile(cat_s3, 0.05),
    cat_s3_hi = quantile(cat_s3, 0.95),
    .by = alpha
  )

## ---- rob-shared-factor ----
shared_sd_log <- params$dirichlet$shared_sd_log$value
set.seed(seed)
shared_cat_s1 <- map_dbl(seq_len(n_draws), \(d) {
  P <- draw_P(P_baseline_raw, alpha_conc)
  f <- exp(rnorm(1, 0, shared_sd_log))
  for (i in seq_len(t_dim)) {
    a_new <- min(P[i, t_dim + 1] * f, 0.95)
    others <- setdiff(seq_len(ncol(P)), t_dim + 1)
    P[i, others] <- P[i, others] * (1 - a_new) / (1 - P[i, t_dim + 1])
    P[i, t_dim + 1] <- a_new
  }
  stopifnot(abs(sum(P[1:t_dim, ]) - t_dim) < 1e-8)
  outcomes_from_s1(P)["cat"]
})
indep_cat_s1 <- draws |> filter(Stage == transient_names[1]) |> pull(cat_b)
shared_int <- quantile(shared_cat_s1, c(0.05, 0.95))
indep_int <- quantile(indep_cat_s1, c(0.05, 0.95))

## ---- rob-paired ----
set.seed(seed)
paired <- map(seq_len(n_draws), \(d) {
  pr <- draw_pair(P_baseline_raw, P_hardened_raw, alpha_conc)
  tibble(draw = d, Stage = transient_names,
         uplift = summarise_draw(pr$h)$cont - summarise_draw(pr$b)$cont)
}) |>
  list_rbind() |>
  mutate(Stage = factor(Stage, levels = transient_names))

paired_summary <- paired |>
  summarise(med = median(uplift), lo = quantile(uplift, 0.05), hi = quantile(uplift, 0.95), .by = Stage) |>
  arrange(Stage)

## ---- rob-ranking ----
n_rank <- params$dirichlet$n_rank$value
set.seed(seed)
rank_draws <- map(seq_len(n_rank), \(d) {
  P <- draw_P(P_baseline_raw, alpha_conc)
  g <- map_dbl(seq_len(t_dim), \(i) shift_cont_from_s1(P, i, delta = shift_delta))
  tibble(draw = d, best = which.max(g), s2_beats_s1 = g[2] > g[1])
}) |> list_rbind()

best_row_share <- rank_draws |> count(best) |> mutate(share = n / sum(n))
p_best_s2 <- sum(best_row_share$share[best_row_share$best == 2])
p_best_s1 <- sum(best_row_share$share[best_row_share$best == 1])
p_s2_beats_s1 <- mean(rank_draws$s2_beats_s1, na.rm = TRUE)

## ---- rob-simulation ----
n_sim <- params$simulation$n_sim$value
set.seed(seed)
sims <- map(seq_len(n_sim), \(k) simulate_walk(P_baseline_raw, 1)) |>
  map(\(v) tibble(outcome = v[["outcome"]], steps = v[["steps"]])) |>
  list_rbind() |>
  mutate(hours = steps * unname(step_hours[1]))

sim_props <- as.vector(prop.table(table(factor(sims$outcome, levels = 1:3))))
sim_check <- tibble(
  Quantity = c(absorbing_names, "Expected steps from S1",
               "E[hours | contained], from S1", "E[hours | catastrophic], from S1",
               "P(catastrophic within 24 h), from S1"),
  Closed_form = c(base_solution$B[1, ], base_solution$steps[1],
                  res$t_cont_b[1], res$t_cat_b[1], p_cat_by_24h),
  Simulated = c(sim_props, mean(sims$steps),
                mean(sims$hours[sims$outcome == 1]), mean(sims$hours[sims$outcome > 1]),
                mean(sims$outcome > 1 & sims$hours <= 24)),
  MC_SE = c(sqrt(base_solution$B[1, ] * (1 - base_solution$B[1, ]) / n_sim),
            sd(sims$steps) / sqrt(n_sim),
            sd(sims$hours[sims$outcome == 1]) / sqrt(sum(sims$outcome == 1)),
            sd(sims$hours[sims$outcome > 1]) / sqrt(sum(sims$outcome > 1)),
            sqrt(p_cat_by_24h * (1 - p_cat_by_24h) / n_sim))
)

# Each simulated value must lie within 4 Monte Carlo standard errors of the closed form
stopifnot(all(abs(sim_check$Closed_form - sim_check$Simulated) < 4 * sim_check$MC_SE))
