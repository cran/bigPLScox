## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/getting-started-",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")

## ----load-data, cache=TRUE, eval=LOCAL----------------------------------------
library(bigPLScox)

data(micro.censure)
data(Xmicro.censure_compl_imp)

Y_all <- micro.censure$survyear[1:80]
status_all <- micro.censure$DC[1:80]
X_all <- apply(
  as.matrix(Xmicro.censure_compl_imp),
  MARGIN = 2,
  FUN = as.numeric
)[1:80, ]

set.seed(2024)
train_id <- 1:60
test_id <- 61:80

X_train <- X_all[train_id, ]
X_test <- X_all[test_id, ]
Y_train <- Y_all[train_id]
Y_test <- Y_all[test_id]
status_train <- status_all[train_id]
status_test <- status_all[test_id]

## ----original-design, cache=TRUE, eval=LOCAL----------------------------------
X_train_raw <- Xmicro.censure_compl_imp[train_id, ]
X_test_raw <- Xmicro.censure_compl_imp[test_id, ]

## ----deviance-residuals, cache=TRUE, eval=LOCAL-------------------------------
residuals_overview <- computeDR(Y_train, status_train, plot = TRUE)
eta_null <- rep(0, length(Y_train))
head(residuals_overview)

if (requireNamespace("bench", quietly = TRUE)) {
  benchmark_dr <- bench::mark(
    survival = computeDR(Y_train, status_train, engine = "survival"),
    cpp = computeDR(Y_train, status_train, engine = "cpp", eta = eta_null),
    iterations = 10,
    check = FALSE
  )
  benchmark_dr
}

all.equal(
  as.numeric(computeDR(Y_train, status_train, engine = "survival")),
  as.numeric(computeDR(Y_train, status_train, engine = "cpp", eta = eta_null)),
  tolerance = 1e-7
)

## ----fit-coxgpls, cache=TRUE, eval=LOCAL--------------------------------------
set.seed(123)
cox_pls_fit <- coxgpls(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 20)
)
cox_pls_fit

## ----fit-formula, cache=TRUE, eval=LOCAL--------------------------------------
cox_pls_fit_formula <- coxgpls(
  ~ ., Y_train, status_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 20),
  dataXplan = data.frame(X_train_raw)
)
cox_pls_fit_formula

## ----cv-coxgpls, cache=TRUE, eval=LOCAL---------------------------------------
set.seed(123456)
cv_results <- suppressWarnings(cv.coxgpls(
  list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20)
))
cv_results

## ----fit-coxgplsdr, cache=TRUE, eval=LOCAL------------------------------------
set.seed(123456)
cox_pls_dr <- coxgplsDR(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = cv_results$nt,
  ind.block.x = c(3, 10, 20)
)
cox_pls_dr

## ----fast-dense, cache=TRUE, eval=LOCAL---------------------------------------
fast_fit_dense <- big_pls_cox_fast(
  X = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 4
)
summary(fast_fit_dense)

## ----fast-dense-predict, cache=TRUE, eval=LOCAL-------------------------------
linear_predictor_fast <- predict(fast_fit_dense, newdata = X_test, type = "link")
head(linear_predictor_fast)

## ----classic-dense, cache=TRUE, eval=LOCAL------------------------------------
legacy_fit_dense <- big_pls_cox(
  X = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 4
)
legacy_fit_dense$cox_fit

## ----fastvsgd, cache=TRUE, eval=LOCAL-----------------------------------------
set.seed(2024)
# Exact fast PLS Cox (dense)
fast_fit_dense <- big_pls_cox_fast(
  X     = X_train,
  time  = Y_train,
  status = status_train,
  ncomp = 4
)

# Exact fast PLS Cox (big.matrix)
if (requireNamespace("bigmemory", quietly = TRUE)) {
  library(bigmemory)
  X_big_train <- bigmemory::as.big.matrix(X_train)
  X_big_test  <- bigmemory::as.big.matrix(X_test)

  fast_fit_big <- big_pls_cox_fast(
    X      = X_big_train,
    time   = Y_train,
    status = status_train,
    ncomp  = 4
  )

  # Gradient based fit in the same latent space
  gd_fit <- big_pls_cox_gd(
    X        = X_big_train,
    time     = Y_train,
    status   = status_train,
    ncomp    = 4,
    max_iter = 2000
  )

  risk_table <- data.frame(
    subject   = seq_along(test_id),
    fast_dense = predict(fast_fit_dense, newdata = X_test, type = "link"),
    fast_big   = predict(fast_fit_big,   newdata = X_big_test, type = "link"),
    gd         = predict(gd_fit,         newdata = X_big_test, type = "link")
  )

  concordances <- apply(
    risk_table[-1],
    2,
    function(lp) {
      survival::concordance(
        survival::Surv(Y_test, status_test) ~ lp
      )$concordance
    }
  )

  concordances
}

## ----prediction-components, cache=TRUE, eval=LOCAL----------------------------
if (exists("fast_fit_dense")) {
  component_scores <- predict(fast_fit_dense, newdata = X_test, type = "components")
  head(component_scores)
}

## ----concordance, cache=TRUE, eval=LOCAL--------------------------------------
if (exists("fast_fit_dense")) {
  concordance_fast <- survival::concordance(
    survival::Surv(Y_test, status_test) ~ linear_predictor_fast
  )$concordance
  concordance_fast
}

## ----dk-splines, cache=TRUE, eval=LOCAL---------------------------------------
cox_DKsplsDR_fit <- coxDKgplsDR(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
  validation = "CV",
  ind.block.x = c(3, 10, 20),
  verbose = FALSE
)
cox_DKsplsDR_fit

## ----cv-dk-splines, cache=TRUE, eval=LOCAL------------------------------------
set.seed(2468)
cv_coxDKgplsDR_res <- suppressWarnings(cv.coxDKgplsDR(
  list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20)
))
cv_coxDKgplsDR_res

