## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/benchmarking-",
  fig.width = 7,
  fig.height = 5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")

## ----packages, include=FALSE--------------------------------------------------
needed <- c("bigmemory", "bench", "survival", "plsRcox", "bigPLScox")
has    <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)

if (!all(has)) {
  missing <- paste(needed[!has], collapse = ", ")
  knitr::opts_chunk$set(eval = FALSE)
  message("Note: skipping code execution because these packages are missing: ", missing)
}

## ----packages_lib-------------------------------------------------------------
library(bigPLScox)
library(survival)
library(bench)
library(bigmemory)
library(plsRcox)

## ----simulate-data------------------------------------------------------------
set.seed(2024)
sim_design <- dataCox(
  n = 2000,
  lambda = 2,
  rho = 1.5,
  x = matrix(rnorm(2000 * 50), ncol = 50),
  beta = c(1, 3, rep(0, 48)),
  cens.rate = 5
)

cox_data <- list(
  x = as.matrix(sim_design[, -(1:3)]),
  time = sim_design$time,
  status = sim_design$status
)

X_big <- bigmemory::as.big.matrix(cox_data$x)

## ----dense-benchmark----------------------------------------------------------
bench_dense <- bench::mark(
  coxgpls = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = 5
  ),
  big_pls = big_pls_cox(cox_data$x, cox_data$time, cox_data$status, ncomp = 5),
  big_pls_fast = big_pls_cox_fast(cox_data$x, cox_data$time, cox_data$status, ncomp = 5),
  big_pls_bb = big_pls_cox_gd(X_big, cox_data$time, cox_data$status, ncomp = 5),
  big_pls_gd = big_pls_cox_gd(X_big, cox_data$time, cox_data$status, ncomp = 5, method = "gd"),
  coxph = coxph(Surv(cox_data$time, cox_data$status) ~ cox_data$x, ties = "breslow"),
  iterations = 75,
  check = FALSE
)
bench_dense$expression <- names(bench_dense$expression)
bench_dense[, c("expression", "median", "itr/sec", "mem_alloc")]

## ----big-benchmark------------------------------------------------------------
bench_big <- bench::mark(
  big_pls_fast = big_pls_cox_fast(X_big, cox_data$time, cox_data$status, ncomp = 5),
  big_pls = big_pls_cox(cox_data$x, cox_data$time, cox_data$status, ncomp = 5),
  big_pls_gd = big_pls_cox_gd(X_big, cox_data$time, cox_data$status, ncomp = 5, max_iter = 300),
  iterations = 75,
  check = FALSE
)
bench_big$expression <- names(bench_big$expression)
bench_big[, c("expression", "median", "itr/sec", "mem_alloc")]

## ----dense-plot---------------------------------------------------------------
plot(bench_dense, type = "jitter")

## ----big-plot-----------------------------------------------------------------
plot(bench_big, type = "jitter")

## ----export-benchmark, eval = FALSE-------------------------------------------
# if (!dir.exists("inst/benchmarks/results")) {
#   dir.create("inst/benchmarks/results", recursive = TRUE)
# }
# write.csv(bench_dense[, 1:9], file = "inst/benchmarks/results/benchmark-dense.csv", row.names = FALSE)
# write.csv(bench_big[, 1:9], file = "inst/benchmarks/results/benchmark-big.csv", row.names = FALSE)

