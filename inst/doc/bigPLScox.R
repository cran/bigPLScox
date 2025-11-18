## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/fast-big--",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")

## ----simulate-bigmatrix, cache=TRUE, eval=LOCAL-------------------------------
library(bigPLScox)
set.seed(2024)

n_obs  <- 5000
n_pred <- 100
k_true <- 3

Sigma <- diag(n_pred)
for (b in 0:2) {
  idx <- (b * 10 + 1):(b * 10 + 10)
  Sigma[idx, idx] <- 0.7
  diag(Sigma[idx, idx]) <- 1
}

L <- chol(Sigma)
Z <- matrix(rnorm(n_obs * n_pred), nrow = n_obs)
X_dense <- Z %*% L

w1 <- numeric(n_pred); w1[1:4]   <- c( 1.0,  0.8,  0.6, -0.5)
w2 <- numeric(n_pred); w2[5:8]   <- c(-0.7,  0.9, -0.4,  0.5)
w3 <- numeric(n_pred); w3[9:12]  <- c( 0.6, -0.5,  0.7,  0.8)
w1 <- w1 / sqrt(sum(w1^2))
w2 <- w2 / sqrt(sum(w2^2))
w3 <- w3 / sqrt(sum(w3^2))

t1 <- as.numeric(scale(drop(X_dense %*% w1), center = TRUE, scale = TRUE))
t2 <- as.numeric(scale(drop(X_dense %*% w2), center = TRUE, scale = TRUE))
t3 <- as.numeric(scale(drop(X_dense %*% w3), center = TRUE, scale = TRUE))

beta_true <- c(1.0, 0.6, 0.3)
eta       <- beta_true[1]*t1 + beta_true[2]*t2 + beta_true[3]*t3

lambda0 <- 0.05
u <- runif(n_obs)
time_event <- -log(u) / (lambda0 * exp(eta))

target_event <- 0.65

f <- function(lc) {
  mean(time_event <= rexp(n_obs, rate = lc)) - target_event
}
lambda_c <- uniroot(f, c(1e-4, 1), tol = 1e-4)$root
time_cens <- rexp(n_obs, rate = lambda_c)

time   <- pmin(time_event, time_cens)
status <- as.integer(time_event <= time_cens)

## ----bigmatrixbigmatrix, cache=TRUE, eval=LOCAL-------------------------------
if (requireNamespace("bigmemory", quietly = TRUE)) {
  library(bigmemory)
  big_dir <- tempfile("bigPLScox-")
  dir.create(big_dir)
  if(file.exists(paste(big_dir,"X.bin",sep="/"))){unlink(paste(big_dir,"X.bin",sep="/"))}
  if(file.exists(paste(big_dir,"X.desc",sep="/"))){unlink(paste(big_dir,"X.desc",sep="/"))}
  X_big <- bigmemory::as.big.matrix(
    X_dense,
    backingpath = big_dir,
    backingfile = "X.bin",
    descriptorfile = "X.desc"
  )
  X_big[1:6, 1:6]
}

## ----dense-solvers, cache=TRUE, eval=LOCAL------------------------------------
fit_legacy <- big_pls_cox(
  X = X_dense,
  time = time,
  status = status,
  ncomp = k_true
)
fit_fast_dense <- big_pls_cox_fast(
  X = X_dense,
  time = time,
  status = status,
  ncomp = k_true
)

list(
  legacy = head(fit_legacy$scores),
  fast = head(fit_fast_dense$scores)
)

list(
  legacy = apply(fit_legacy$scores, 2, var),
  fast = apply(fit_fast_dense$scores, 2, var)
)

## ----dense-summary, cache=TRUE, eval=LOCAL------------------------------------
summary(fit_fast_dense)

## ----fast-big, cache=TRUE, eval=LOCAL-----------------------------------------
if (exists("X_big")) {
  fit_fast_big <- big_pls_cox_fast(
    X = X_big,
    time = time,
    status = status,
    ncomp = k_true
  )
  summary(fit_fast_big)
}

## ----gd-fit, cache=TRUE, eval=LOCAL-------------------------------------------
if (exists("X_big")) {
  fit_gd <- big_pls_cox_gd(
    X = X_big,
    time = time,
    status = status,
    ncomp = k_true,
    max_iter = 2000,
    learning_rate = 0.05,
    tol = 1e-8
  )
  str(fit_gd)
}

## ----compare-solvers, cache=TRUE, eval=LOCAL----------------------------------
if (exists("fit_fast_dense") && exists("fit_gd")) {
  eta_fast_dense <- drop(fit_fast_dense$scores %*% fit_fast_dense$cox_fit$coefficients)
  eta_fast_big <- drop(fit_fast_big$scores %*% fit_fast_big$cox_fit$coefficients)
  eta_gd <- drop(fit_gd$scores %*% fit_gd$cox_fit$coefficients)

  list(
    correlation_fast_dense_big = cor(eta_fast_dense, eta_fast_big),
    correlation_fast_dense_gd = cor(eta_fast_dense, eta_gd),
    concordance = c(
      fast_dense = survival::concordance(survival::Surv(time, status) ~ eta_fast_dense)$concordance,
      fast_big = survival::concordance(survival::Surv(time, status) ~ eta_fast_big)$concordance,
      gd = survival::concordance(survival::Surv(time, status) ~ eta_gd)$concordance
    )
  )
}

## ----predict-new, cache=TRUE, eval=LOCAL--------------------------------------
X_new <- matrix(rnorm(10 * n_pred), nrow = 10)

scores_new <- predict(
  fit_fast_dense,
  newdata = X_new,
  type = "components"
)

risk_new <- predict(
  fit_fast_dense,
  newdata = X_new,
  type = "risk"
)

list(scores = scores_new, risk = risk_new)

## ----timing, cache=TRUE, eval=LOCAL-------------------------------------------
if (requireNamespace("bench", quietly = TRUE) && exists("X_big")) {
  bench_res <- bench::mark(
    big_pls = big_pls_cox(X_dense, time, status, ncomp = k_true),
    fast_dense = big_pls_cox_fast(X_dense, time, status, ncomp = k_true),
    fast_big = big_pls_cox_fast(X_big, time, status, ncomp = k_true),
    gd = big_pls_cox_gd(X_big, time, status, ncomp = k_true, max_iter = 1500),
    iterations = 30,
    check = FALSE
  )
  bench_res$expression <- names(bench_res$expression)
  bench_res[, c("expression", "median", "itr/sec", "mem_alloc")]
}

## ----timing-plot, cache=TRUE, eval=LOCAL--------------------------------------
if (exists("bench_res")) {
  plot(bench_res, type = "jitter")
}

## ----cleanup, cache=TRUE, eval=LOCAL------------------------------------------
if (exists("X_big")) {
  rm(X_big)
  file.remove(file.path(big_dir, "X.bin"))
  file.remove(file.path(big_dir, "X.desc"))
  unlink(big_dir, recursive = TRUE)
}

