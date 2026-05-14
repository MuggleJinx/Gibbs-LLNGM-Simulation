# Gibbs sampler for non-centered parameterization
# Model:
#   Y | M ~ N(X beta + B (M - mu h), sigma_e^2 I),  B = A K^{-1}
#   M | V ~ N(mu V, diag(sigma^2 V))
#   V_i ~ GIG(p, a, b)
#
# Conditional updates:
#   M | V, Y ~ N( Qbar^{-1} mbar, sigma^2 Qbar^{-1} )
#     Qbar = rho * B'B + D_V^{-1},  rho = sigma^2 / sigma_e^2
#     mbar = rho * B' (Y - X beta + mu B h) + mu * 1
#   V_i | M_i ~ GIG(p - 1/2, a + (mu/sigma)^2, b + (M_i/sigma)^2)
#
# GIG parameterization (rgig): density ∝ x^{p-1} exp(-(a x + b/x)/2), x > 0

resolve_rgig <- function(rgig_fun = NULL) {
  if (!is.null(rgig_fun)) {
    return(rgig_fun)
  }
  if (requireNamespace("ngme2", quietly = TRUE)) {
    return(ngme2::rgig)
  }
  stop("rgig() not found. Install ngme2, or pass rgig_fun.")
}

rgig_std <- function(n, p, a, b, rgig_core, tol = 0) {
  # Boundary-safe GIG sampler for:
  # density ∝ x^{p-1} exp(-(a x + b/x)/2)
  # Gamma: b == 0, a > 0, p > 0  -> Gamma(shape=p, rate=a/2)
  # Inv-Gamma: a == 0, b > 0, p < 0 -> Inv-Gamma(alpha=-p, beta=b/2)
  p <- rep(p, length.out = n)
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)

  out <- numeric(n)
  a0 <- abs(a) <= tol
  b0 <- abs(b) <= tol

  idx_gamma <- b0 & (a > 0) & (p > 0)
  idx_invgamma <- a0 & (b > 0) & (p < 0)
  idx_gig <- !(idx_gamma | idx_invgamma)

  if (any(idx_gamma)) {
    out[idx_gamma] <- rgamma(sum(idx_gamma), shape = p[idx_gamma], rate = a[idx_gamma] / 2)
  }
  if (any(idx_invgamma)) {
    out[idx_invgamma] <- 1 / rgamma(sum(idx_invgamma), shape = -p[idx_invgamma], rate = b[idx_invgamma] / 2)
  }
  if (any(idx_gig)) {
    out[idx_gig] <- rgig_core(sum(idx_gig), p[idx_gig], a[idx_gig], b[idx_gig])
  }

  out
}

sample_M_given_VY <- function(V, mbar, BtB, rho, sigma, jitter = 0) {
  n <- length(V)
  Qbar <- rho * BtB + Matrix::Diagonal(n, 1 / V)
  if (jitter > 0) {
    Qbar <- Qbar + Matrix::Diagonal(n, jitter)
  }
  Qbar_mat <- as.matrix(Qbar)
  U <- chol(Qbar_mat)

  mean_M <- backsolve(U, forwardsolve(t(U), mbar))
  z <- rnorm(n)
  M <- as.numeric(mean_M + sigma * backsolve(U, z))
  M
}

sample_V_given_M <- function(M, p, a, b, mu, sigma, rgig_fun, rgig_tol = 0) {
  n <- length(M)
  p2 <- p - 0.5
  a2 <- a + (mu / sigma)^2
  b2 <- b + (M / sigma)^2
  rgig_std(n, p2, a2, b2, rgig_core = rgig_fun, tol = rgig_tol)
}

gibbs_nc <- function(n_iter,
                     Y, X = NULL, beta = NULL,
                     A, Kinv, mu, h,
                     p, a, b,
                     sigma, sigma_e,
                     V_init = NULL,
                     burn = 0, thin = 1,
                     seed = NULL,
                     rgig_fun = NULL,
                     jitter = 0,
                     progress = FALSE,
                     progress_every = 1000,
                     rgig_tol = 0) {
  if (!is.null(seed)) set.seed(seed)
  rgig_fun <- resolve_rgig(rgig_fun)

  Y <- as.numeric(Y)
  n_obs <- length(Y)

  if (is.null(V_init)) {
    V <- NULL
  } else {
    V <- as.numeric(V_init)
  }

  B <- as.matrix(A %*% Kinv)
  n_latent <- ncol(B)
  if (nrow(B) != n_obs) {
    stop("Dimension mismatch: nrow(B) != length(Y).")
  }

  if (is.null(V)) {
    V <- rgig_std(n_latent, p, a, b, rgig_core = rgig_fun, tol = rgig_tol)
  } else if (length(V) != n_latent) {
    stop("V_init dimension mismatch: length(V_init) != ncol(B).")
  }

  BtB <- Matrix::crossprod(B)
  Bh <- as.numeric(B %*% h)

  Xbeta <- rep(0, n_obs)
  if (!is.null(X) && !is.null(beta)) {
    Xbeta <- as.numeric(X %*% beta)
    if (length(Xbeta) != n_obs) {
      stop("Dimension mismatch: X %*% beta does not match length(Y).")
    }
  }

  rho <- (sigma^2) / (sigma_e^2)
  ytilde <- Y - Xbeta + mu * Bh
  mbar <- as.numeric(rho * Matrix::crossprod(B, ytilde) + mu * rep(1, n_latent))

  Nsave <- floor((n_iter - burn) / thin)
  M_samp <- matrix(NA_real_, nrow = Nsave, ncol = n_latent)
  V_samp <- matrix(NA_real_, nrow = Nsave, ncol = n_latent)

  idx <- 0L
  for (t in seq_len(n_iter)) {
    if (progress && (t == 1L || (t %% progress_every == 0L))) {
      message(sprintf("gibbs_nc: iter %d / %d", t, n_iter))
    }
    M <- sample_M_given_VY(V, mbar, BtB, rho, sigma, jitter = jitter)
    V <- sample_V_given_M(M, p, a, b, mu, sigma, rgig_fun, rgig_tol = rgig_tol)

    if (t > burn && ((t - burn) %% thin == 0)) {
      idx <- idx + 1L
      M_samp[idx, ] <- M
      V_samp[idx, ] <- V
    }
  }

  list(
    M = M_samp, V = V_samp,
    meta = list(n_iter = n_iter, burn = burn, thin = thin)
  )
}
