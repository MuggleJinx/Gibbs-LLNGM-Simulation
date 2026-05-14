library(Matrix)

source(if (file.exists("paths.R")) "paths.R" else file.path("Simulation", "paths.R"))
source(sim_path("gibbs.R"))

# ---- minimal sanity test for non-centered Gibbs sampler ----
set.seed(123)

# dimensions
n <- 8

# parameters
p <- 1.2
 a <- 1.1
b <- 0.9
mu <- 0.7
sigma <- 1.0
sigma_e <- 0.5

# design / operators
A <- diag(n)
Kinv <- diag(n)
B <- A %*% Kinv

X <- matrix(1, n, 1)
beta <- 0.3
h <- rep(1, n)

# simulate a small synthetic dataset
rgig_fun <- resolve_rgig()
V_true <- rgig_fun(n, p, a, b)
M_true <- mu * V_true + rnorm(n, 0, sigma * sqrt(V_true))
Y <- as.numeric(X %*% beta + B %*% (M_true - mu * h) + rnorm(n, 0, sigma_e))

# run Gibbs sampler
out <- gibbs_nc(
  n_iter = 2000,
  Y = Y, X = X, beta = beta,
  A = A, Kinv = Kinv, mu = mu, h = h,
  p = p, a = a, b = b,
  sigma = sigma, sigma_e = sigma_e,
  burn = 500, thin = 5,
  seed = 999
)

cat("M samples dim:", dim(out$M), "\n")
cat("V samples dim:", dim(out$V), "\n")

cat("Posterior mean M (first 5):", round(colMeans(out$M)[1:5], 3), "\n")
cat("Posterior mean V (first 5):", round(colMeans(out$V)[1:5], 3), "\n")

cat("True M (first 5):", round(M_true[1:5], 3), "\n")
cat("True V (first 5):", round(V_true[1:5], 3), "\n")
