library(Matrix)

source(if (file.exists("paths.R")) "paths.R" else file.path("Simulation", "paths.R"))
source(sim_path("gibbs.R"))

# ---- experiment design (B1) ----
# Non-centered Gibbs sampler; focus on V-marginal mixing

ess_est <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3) {
    return(NA_real_)
  }
  if (is.null(max_lag)) max_lag <- min(1000, n - 1)
  ac <- acf(x, plot = FALSE, lag.max = max_lag)$acf[-1]
  if (length(ac) == 0) {
    return(n)
  }
  k <- which(ac <= 0)
  if (length(k) > 0) {
    ac <- ac[1:(k[1] - 1)]
  }
  tau <- 1 + 2 * sum(ac)
  n / max(tau, 1)
}

iact_est <- function(x, max_lag = NULL) {
  n <- length(x)
  ess <- ess_est(x, max_lag = max_lag)
  if (!is.finite(ess) || ess <= 0) {
    return(NA_real_)
  }
  n / ess
}

split_rhat <- function(chains) {
  chains <- lapply(chains, as.numeric)
  m <- length(chains)
  if (m < 2) {
    return(NA_real_)
  }

  split_chains <- list()
  for (i in seq_len(m)) {
    x <- chains[[i]]
    n <- length(x)
    half <- floor(n / 2)
    if (half < 2) {
      return(NA_real_)
    }
    split_chains[[length(split_chains) + 1L]] <- x[1:half]
    split_chains[[length(split_chains) + 1L]] <- x[(half + 1):(2 * half)]
  }

  m2 <- length(split_chains)
  n2 <- length(split_chains[[1]])
  if (m2 < 2 || n2 < 2) {
    return(NA_real_)
  }

  chain_means <- vapply(split_chains, mean, numeric(1))
  chain_vars <- vapply(split_chains, var, numeric(1))

  W <- mean(chain_vars)
  B <- n2 * var(chain_means)
  var_hat <- ((n2 - 1) / n2) * W + (B / n2)
  sqrt(var_hat / W)
}

summarize_chain <- function(V_samples, q = 0.25) {
  S_plus <- rowMeans(V_samples)
  S_minus <- rowMeans(V_samples^(-q))
  S_log <- rowMeans(log(V_samples))
  list(S_plus = S_plus, S_minus = S_minus, S_log = S_log)
}

run_regime <- function(label, pars,
                       n_iter, burn, thin,
                       A, Kinv, h,
                       sigma, sigma_e,
                       Y,
                       V_inits,
                       q = 0.25,
                       seed_base = 1) {
  chains_stats <- list()
  chains_time <- numeric(length(V_inits))

  for (i in seq_along(V_inits)) {
    t0 <- proc.time()[[3]]
    out <- gibbs_nc(
      n_iter = n_iter,
      Y = Y, X = NULL, beta = NULL,
      A = A, Kinv = Kinv, mu = pars$mu, h = h,
      p = pars$p, a = pars$a, b = pars$b,
      sigma = sigma, sigma_e = sigma_e,
      V_init = V_inits[[i]],
      burn = burn, thin = thin,
      seed = seed_base + i,
      progress = TRUE,
      progress_every = 1000
    )
    t1 <- proc.time()[[3]]
    chains_time[i] <- t1 - t0

    chains_stats[[i]] <- summarize_chain(out$V, q = q)
  }

  list(chains_stats = chains_stats, chains_time = chains_time)
}

summarize_regime <- function(chains_stats, chains_time) {
  stats_names <- c("S_plus", "S_minus", "S_log")
  out <- list()

  for (nm in stats_names) {
    series_list <- lapply(chains_stats, function(x) x[[nm]])

    ess_vals <- vapply(series_list, ess_est, numeric(1))
    iact_vals <- vapply(series_list, iact_est, numeric(1))
    ess_sec <- ess_vals / chains_time

    out[[nm]] <- list(
      ess = mean(ess_vals, na.rm = TRUE),
      ess_sec = mean(ess_sec, na.rm = TRUE),
      iact = mean(iact_vals, na.rm = TRUE),
      rhat = split_rhat(series_list)
    )
  }

  out
}

# ---- global setup ----
set.seed(1)

n <- 300
rho_ar1 <- 0.5
A <- Diagonal(n)
K <- ngme2::ar1(1:n, rho_ar1)$K
Kinv <- solve(K)

sigma <- 1
sigma_e <- 1
h <- rep(1, n)

# Main experiment: keep coupling but no real data
# Use Y = 0 with a fixed sigma_e so rho does not vanish
sigma_e_eff <- sigma_e
Y <- rep(0, n)

# Regime grid (A-F)
regimes <- list(
  A = list(label = "TC-1", p = -0.5, a = 1, b = 1, mu = 1),
  B = list(label = "TC-1", p = 1.0, a = 1, b = 1, mu = 1),
  C = list(label = "TC-2", p = 0.6, a = 1, b = 0, mu = 1),
  D = list(label = "DM-III", p = 0.3, a = 1, b = 0, mu = 1),
  E = list(label = "DM-I", p = -1.5, a = 0, b = 2, mu = 0),
  F = list(label = "DM-II", p = -1.5, a = 0, b = 2, mu = 1)
)

n_iter <- 50000
burn <- 5000
thin <- 1
q <- 0.25

rgig_fun <- resolve_rgig()
V_inits <- list(
  rep(1, n),
  rep(0.1, n),
  rep(10, n),
  rgig_fun(n, p = 1, a = 1, b = 1)
)

results <- list()
summary_rows <- list()

run_one_regime <- function(key) {
  pars <- regimes[[key]]
  cat("Running regime", key, "...\n")

  res <- run_regime(
    label = key,
    pars = pars,
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    A = A,
    Kinv = Kinv,
    h = h,
    sigma = sigma,
    sigma_e = sigma_e_eff,
    Y = Y,
    V_inits = V_inits,
    q = q,
    seed_base = 100
  )

  summ <- summarize_regime(res$chains_stats, res$chains_time)
  row <- data.frame(
    regime = pars$label,
    point = key,
    p = pars$p,
    a = pars$a,
    b = pars$b,
    mu = pars$mu,
    ess_sec_S_plus = summ$S_plus$ess_sec,
    ess_sec_S_minus = summ$S_minus$ess_sec,
    ess_sec_S_log = summ$S_log$ess_sec,
    iact_S_plus = summ$S_plus$iact,
    iact_S_minus = summ$S_minus$iact,
    iact_S_log = summ$S_log$iact,
    rhat_max = max(c(summ$S_plus$rhat, summ$S_minus$rhat, summ$S_log$rhat), na.rm = TRUE)
  )

  list(key = key, raw = res, summary = summ, row = row)
}

# ---- parallel settings (cross-platform) ----
use_parallel <- TRUE
workers <- max(1, parallel::detectCores() - 1)

keys <- names(regimes)
if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = workers)
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("txtprogressbar")
    out_list <- progressr::with_progress({
      p <- progressr::progressor(along = keys)
      future.apply::future_lapply(
        keys,
        function(k) {
          res <- run_one_regime(k)
          p(sprintf("regime %s", k))
          res
        },
        future.seed = TRUE
      )
    })
  } else {
    message("progressr not available; running without parallel progress.")
    out_list <- future.apply::future_lapply(keys, run_one_regime, future.seed = TRUE)
  }
} else {
  if (use_parallel) {
    message("future.apply not available; running sequentially.")
  }
  out_list <- lapply(keys, run_one_regime)
}

for (res in out_list) {
  key <- res$key
  results[[key]] <- list(raw = res$raw, summary = res$summary)
  summary_rows[[key]] <- res$row
}

summary_table <- do.call(rbind, summary_rows)

saveRDS(results, file = sim_path("experiment_results.rds"))
write.csv(summary_table, file = sim_path("experiment_summary.csv"), row.names = FALSE)

print(summary_table)
