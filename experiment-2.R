library(Matrix)
library(ngme2)

source(if (file.exists("paths.R")) "paths.R" else file.path("Simulation", "paths.R"))
source(sim_path("gibbs.R")) # must provide gibbs_nc(), resolve_rgig()

# ============================================================
# B2-scan-mu: null-smallness scan by varying mu (geometry fixed)
# ============================================================
# Fix A (thus B = A K^{-1} fixed), scan mu to change gamma_NS(mu).
# Summaries monitored (per retained iteration t):
#   S_plus(t)  = n^{-1} sum_i V_i(t)
#   S_minus(t) = n^{-1} sum_i V_i(t)^(-q)
#   S_log(t)   = n^{-1} sum_i log V_i(t)
#   T_null(t)  = <K^{-1} u0, M(t) - mu h>,  where u0 spans Null(B) (unit)
#
# X-axis:
#   gamma_NS(mu) = ||proj_{Null(B)} mbar(mu)|| / (sigma * sqrt(tilde_a(mu)))
#   tilde_a(mu) = a + (mu/sigma)^2
#
# NOTE:
# - In parallel, inner `gibbs_nc(progress=TRUE)` usually won't show nicely.
#   Use progressr to track outer tasks (mu x chains).
# - This script runs a fixed grid and DOES NOT auto-stop on Rhat.

# -----------------------------
# Numerics
# -----------------------------
tiny <- 1e-300 # avoid log(0) if numerical underflow happens
vmax <- 1e150 # used only for diagnostic statistics, to avoid Inf; you can also use 1e100/1e200
clip_V <- function(V, tiny = 1e-300, vmax = 1e150) {
  V <- pmax(V, tiny)
  V <- pmin(V, vmax)
  V[!is.finite(V)] <- vmax
  V
}

# -----------------------------
# Diagnostics (batch means + fallback)
# -----------------------------
ess_batch <- function(x, batch = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 200) {
    return(NA_real_)
  }
  if (is.null(batch)) batch <- floor(sqrt(n))
  batch <- max(10, batch)
  b <- floor(n / batch)
  if (b < 5) {
    return(NA_real_)
  }
  m <- b * batch
  x <- x[1:m]
  xb <- colMeans(matrix(x, nrow = batch, ncol = b))
  s2 <- var(xb) * batch # spectral var at 0
  if (!is.finite(s2) || s2 <= 0) {
    return(NA_real_)
  }
  n * var(x) / s2
}

iact_batch <- function(x, batch = NULL) {
  n <- length(x)
  ess <- ess_batch(x, batch = batch)
  if (!is.finite(ess) || ess <= 0) {
    return(NA_real_)
  }
  n / ess
}

# ACF-based fallback (works for short series; more volatile)
ess_acf <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 20) {
    return(NA_real_)
  }
  if (is.null(max_lag)) max_lag <- min(1000, n - 1)
  ac <- acf(x, plot = FALSE, lag.max = max_lag)$acf[-1]
  if (length(ac) == 0) {
    return(n)
  }
  k <- which(ac <= 0)
  if (length(k) > 0) ac <- ac[1:(k[1] - 1)]
  tau <- 1 + 2 * sum(ac)
  n / max(tau, 1)
}

iact_acf <- function(x, max_lag = NULL) {
  n <- length(x)
  ess <- ess_acf(x, max_lag = max_lag)
  if (!is.finite(ess) || ess <= 0) {
    return(NA_real_)
  }
  n / ess
}

# automatic: batch if possible, else ACF fallback
ess_auto <- function(x) {
  out <- ess_batch(x)
  if (is.finite(out)) {
    return(out)
  }
  ess_acf(x)
}

iact_auto <- function(x) {
  out <- iact_batch(x)
  if (is.finite(out)) {
    return(out)
  }
  iact_acf(x)
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
    if (half < 50) {
      return(NA_real_)
    }
    split_chains[[length(split_chains) + 1L]] <- x[1:half]
    split_chains[[length(split_chains) + 1L]] <- x[(half + 1):(2 * half)]
  }

  m2 <- length(split_chains)
  n2 <- length(split_chains[[1]])
  if (m2 < 2 || n2 < 50) {
    return(NA_real_)
  }

  chain_means <- vapply(split_chains, mean, numeric(1))
  chain_vars <- vapply(split_chains, var, numeric(1))

  W <- mean(chain_vars)
  B <- n2 * var(chain_means)

  var_hat <- ((n2 - 1) / n2) * W + (B / n2)
  sqrt(var_hat / W)
}

# -----------------------------
# Geometry helpers
# -----------------------------
build_A_projection <- function(u) {
  u <- as.numeric(u)
  u <- u / sqrt(sum(u^2))
  n <- length(u)
  I_n <- Diagonal(n)
  uuT <- tcrossprod(Matrix(u, ncol = 1)) # rank-1
  I_n - uuT
}

# -----------------------------
# Null-smallness full constant for FIXED B (vary mu)
# -----------------------------
compute_mbar <- function(B, Y, Xbeta, mu, h, sigma, sigma_e) {
  # mbar = rho * B' (Y - Xbeta + mu B h) + mu * 1
  rho <- (sigma^2) / (sigma_e^2)
  ytilde <- as.numeric(Y - Xbeta + mu * as.numeric(B %*% h))
  as.numeric(rho * Matrix::crossprod(B, ytilde) + mu * rep(1, ncol(B)))
}

gamma_ns_from_B <- function(B, u0, Y, Xbeta, mu, h, sigma, sigma_e, a) {
  mbar <- compute_mbar(B, Y, Xbeta, mu, h, sigma, sigma_e)
  m0_norm <- abs(sum(u0 * mbar)) # ||proj_{Null(B)} mbar|| since u0 unit
  tilde_a <- a + (mu / sigma)^2
  gamma_ns <- m0_norm / (sigma * sqrt(tilde_a)) # your "full constant"
  list(gamma_ns = gamma_ns, m0_norm = m0_norm, tilde_a = tilde_a)
}

# ============================================================
# Setup
# ============================================================
set.seed(2)

n <- 300
phi_ar1 <- 0.5

sigma <- 1
sigma_e <- 1
h <- rep(1, n)
Y <- rep(0, n)
Xbeta <- rep(0, n)

# ---- Choose regime (fix p,a,b; scan mu) ----
# DM-III (boundary bad): a>0, b=0, 0<p<=1/2
pars <- list(p = 0.5, a = 1, b = 0)
# DM-II (drifted heavy tail): a=0, b>0, p<0
# pars <- list(p = -1.5, a = 0, b = 2)

q <- 0.25

# K: AR(1)
K <- as.matrix(ngme2::ar1(1:n, phi_ar1)$K)
Kinv <- solve(K)

# ---- FIX A and thus FIX B ----
# Pick a single fixed null direction u (feel free to change; keep fixed for "scientific" scan)
u <- rep(1, n)
u <- u / sqrt(sum(u^2))

A <- build_A_projection(u)
B <- as.matrix(A %*% Kinv)

# Null(B) = span(K u) (for B = A K^{-1}), use unit u0 for projection
u0 <- as.numeric(K %*% u)
u0 <- u0 / sqrt(sum(u0^2))

# For T_null(t) = <K^{-1} u0, M - mu h>
ku <- as.numeric(Kinv %*% u0)

# -----------------------------
# Scan grid of mu (edit this!)
# -----------------------------
mu_grid <- c(0, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0) # example
# mu_grid <- seq(0, 2, by = 0.2)

gamma_grid <- vapply(mu_grid, function(mu) {
  gamma_ns_from_B(B, u0, Y, Xbeta, mu, h, sigma, sigma_e, pars$a)$gamma_ns
}, numeric(1))

cat("mu and gamma_NS:\n")
print(data.frame(mu = mu_grid, gamma_ns = gamma_grid))

# -----------------------------
# Gibbs settings
# -----------------------------
n_iter <- 50000
burn <- 10000
thin <- 1

rgig_fun <- resolve_rgig()
V_inits <- list(
  rep(1, n),
  rep(0.1, n),
  rep(10, n),
  rgig_fun(n, p = 1, a = 1, b = 1)
)

# ============================================================
# Flatten tasks: (mu_idx x chain_idx)
# ============================================================
task_grid <- expand.grid(
  mu_idx = seq_along(mu_grid),
  chain_idx = seq_along(V_inits),
  stringsAsFactors = FALSE
)

cat(sprintf(
  "Total tasks: %d (mu points: %d, chains/mu: %d)\n",
  nrow(task_grid), length(mu_grid), length(V_inits)
))

run_single_chain <- function(mu_idx, chain_idx) {
  mu <- mu_grid[mu_idx]
  gamma_ns <- gamma_grid[mu_idx]

  seed_ij <- 20260000 + 1000L * mu_idx + 10L * chain_idx

  t0 <- proc.time()[[3]]
  out <- gibbs_nc(
    n_iter = n_iter,
    Y = Y, X = NULL, beta = NULL,
    A = A, Kinv = Kinv, mu = mu, h = h,
    p = pars$p, a = pars$a, b = pars$b,
    sigma = sigma, sigma_e = sigma_e,
    V_init = V_inits[[chain_idx]],
    burn = burn, thin = thin,
    seed = seed_ij,
    progress = FALSE,
    rgig_tol = 1e-12
  )
  t1 <- proc.time()[[3]]

  V_raw <- out$V
  M <- out$M

  # Numerical safety: avoid passing Inf/NaN to summaries
  V <- clip_V(V_raw, tiny = tiny, vmax = vmax)

  # (Optional) record how many clips/non-finite values occurred, to help you judge if "numerical explosion" happened
  bad_rate <- mean(!is.finite(V_raw) | V_raw > vmax | V_raw < tiny)

  # summaries on V (all guaranteed to be finite)
  sp_vec <- rowMeans(V)
  sm_vec <- rowMeans(V^(-q))
  slog_vec <- rowMeans(log(V))

  # null-direction summary on M
  M_centered <- sweep(M, 2, mu * h, "-")
  tn_vec <- as.numeric(M_centered %*% ku)

  list(
    bad_rate = bad_rate,
    mu_idx = mu_idx,
    chain_idx = chain_idx,
    mu = mu,
    gamma_ns = gamma_ns,
    time = t1 - t0,
    sp_vec = sp_vec,
    sm_vec = sm_vec,
    slog_vec = slog_vec,
    tn_vec = tn_vec
  )
}

# ============================================================
# Parallel + progress
# ============================================================
use_parallel <- TRUE
workers <- max(1, parallel::detectCores() - 1)

if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("txtprogressbar")
    raw_results <- progressr::with_progress({
      p <- progressr::progressor(along = seq_len(nrow(task_grid)))
      future.apply::future_lapply(
        seq_len(nrow(task_grid)),
        function(k) {
          r <- run_single_chain(task_grid$mu_idx[k], task_grid$chain_idx[k])
          p(sprintf("task %d/%d", k, nrow(task_grid)))
          r
        },
        future.seed = TRUE
      )
    })
  } else {
    message(sprintf("Running parallel with %d workers (no progressr)...", workers))
    raw_results <- future.apply::future_lapply(
      seq_len(nrow(task_grid)),
      function(k) run_single_chain(task_grid$mu_idx[k], task_grid$chain_idx[k]),
      future.seed = TRUE
    )
  }
} else {
  message("Running sequentially...")
  raw_results <- lapply(
    seq_len(nrow(task_grid)),
    function(k) run_single_chain(task_grid$mu_idx[k], task_grid$chain_idx[k])
  )
}

# ============================================================
# Aggregate by mu_idx
# ============================================================
results_by_mu <- split(raw_results, vapply(raw_results, function(x) x$mu_idx, numeric(1)))

final_rows <- list()
final_full <- list()

for (idx_str in names(results_by_mu)) {
  chain_list <- results_by_mu[[idx_str]]
  chain_list <- chain_list[order(vapply(chain_list, function(z) z$chain_idx, numeric(1)))]

  mu <- chain_list[[1]]$mu
  gamma_ns <- chain_list[[1]]$gamma_ns

  sp_list <- lapply(chain_list, function(x) x$sp_vec)
  sm_list <- lapply(chain_list, function(x) x$sm_vec)
  slog_list <- lapply(chain_list, function(x) x$slog_vec)
  tn_list <- lapply(chain_list, function(x) x$tn_vec)
  times <- vapply(chain_list, function(x) x$time, numeric(1))

  # IACT/ESS (auto: batch if long enough, else ACF fallback)
  iact_sp <- vapply(sp_list, iact_auto, numeric(1))
  iact_sm <- vapply(sm_list, iact_auto, numeric(1))
  iact_slog <- vapply(slog_list, iact_auto, numeric(1))
  iact_tn <- vapply(tn_list, iact_auto, numeric(1))

  ess_sp <- vapply(sp_list, ess_auto, numeric(1))
  ess_sm <- vapply(sm_list, ess_auto, numeric(1))
  ess_slog <- vapply(slog_list, ess_auto, numeric(1))
  ess_tn <- vapply(tn_list, ess_auto, numeric(1))

  esssec_sp <- ess_sp / times
  esssec_sm <- ess_sm / times
  esssec_slog <- ess_slog / times
  esssec_tn <- ess_tn / times

  rhat_sp <- split_rhat(sp_list)
  rhat_sm <- split_rhat(sm_list)
  rhat_slog <- split_rhat(slog_list)
  rhat_tn <- split_rhat(tn_list)

  row <- data.frame(
    mu = mu,
    gamma_ns = gamma_ns,
    q = q,
    n_iter = n_iter,
    burn = burn,
    chains = length(chain_list),
    ess_sec_S_plus = mean(esssec_sp, na.rm = TRUE),
    iact_S_plus = mean(iact_sp, na.rm = TRUE),
    rhat_S_plus = rhat_sp,
    ess_sec_S_minus = mean(esssec_sm, na.rm = TRUE),
    iact_S_minus = mean(iact_sm, na.rm = TRUE),
    rhat_S_minus = rhat_sm,
    ess_sec_S_log = mean(esssec_slog, na.rm = TRUE),
    iact_S_log = mean(iact_slog, na.rm = TRUE),
    rhat_S_log = rhat_slog,
    ess_sec_T_null = mean(esssec_tn, na.rm = TRUE),
    iact_T_null = mean(iact_tn, na.rm = TRUE),
    rhat_T_null = rhat_tn,
    rhat_max = max(c(rhat_sp, rhat_sm, rhat_slog, rhat_tn), na.rm = TRUE)
  )

  final_rows[[idx_str]] <- row
  final_full[[as.integer(idx_str)]] <- list(
    mu = mu,
    gamma_ns = gamma_ns,
    sp_list = sp_list,
    sm_list = sm_list,
    slog_list = slog_list,
    tn_list = tn_list,
    chain_time = times
  )
}

summary_table <- do.call(rbind, final_rows)
summary_table <- summary_table[order(summary_table$gamma_ns), ]

saveRDS(final_full, file = sim_path("experiment2_results_scanmu_4stats.rds"))
write.csv(summary_table, file = sim_path("experiment2_summary_scanmu_4stats.csv"), row.names = FALSE)

print(summary_table)
