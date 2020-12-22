#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function tests nonstationary of coefficient function.
#' @param data
#' @param mfit.full estimated full model
#' @param test.idx index of the test covariate in nonstationary coefficient
#' function.
#' @param nBoot number of bootstrap samples
#' @param ncore number of cores used in the parallel computing.
gof.test <- function(data, mfit.full, test.idx, nBoot = 100, ncore = 1) {

  boot.all <- matrix(NA, ncol = 1, nrow = nBoot)
  data.reduce <- data
  data.reduce$Z <- cbind(data.reduce$Z, data.reduce$X[, test.idx])
  data.reduce$X <- data.reduce$X[, -test.idx]
  if (ncol(data.reduce$X) == 0) data.reduce$X <- NULL

  mfit.reduce <- star.fit(
    data = data.reduce,
    matrix(mfit.full$best.lambda, ncol = 2),
    mfit.full$pars$V, mfit.full$pars$Tri,
    mfit.full$pars$d, mfit.full$pars$r,
    mfit.full$pars$time.knots, mfit.full$pars$rho,
    mfit.full$pars$time.bound, return.se = FALSE
  )

  test.stat <- mfit.full$MSE.Y / mfit.reduce$MSE.Y


  Y <- data$Y
  Y.hat <- mfit.full$Y.hat
  n <- length(Y)
  epsilon <- Y - Y.hat - mean(Y - Y.hat)

  registerDoParallel(cores = ncore)

  foreach(iter = 1:nBoot, .combine = comb, .multicombine = TRUE) %do% {
    set.seed(iter)

    # generate bootstrap sample
    epsilon.iter <- sample(epsilon, length(epsilon))
    mu.iter <- mfit.reduce$mu.hat + epsilon.iter
    alpha0 <- mfit.reduce$alpha.hat
    Y.iter <- solve(diag(n) - alpha0 * data$W, mu.iter)

    # fit the full model
    data.boot.full <- data
    data.boot.full$Y <- Y.iter
    mfit.full.iter <- star.fit(
      data = data.boot.full,
      matrix(mfit.full$best.lambda, ncol = 2),
      mfit.full$pars$V, mfit.full$pars$Tri,
      mfit.full$pars$d, mfit.full$pars$r,
      mfit.full$pars$time.knots, mfit.full$pars$rho,
      mfit.full$pars$time.bound, return.se = FALSE
    )

    data.boot.reduce <- data
    data.boot.reduce$Y <- Y.iter
    data.boot.reduce$Z <- cbind(data.boot.reduce$Z, data.boot.reduce$X[, test.idx])
    data.boot.reduce$X <- data.boot.reduce$X[, -test.idx]
    if (ncol(data.boot.reduce$X) == 0) data.boot.reduce$X <- NULL

    mfit.reduce.iter <- star.fit(
      data = data.boot.reduce,
      matrix(mfit.full$best.lambda, ncol = 2),
      mfit.full$pars$V, mfit.full$pars$Tri,
      mfit.full$pars$d, mfit.full$pars$r,
      mfit.full$pars$time.knots, mfit.full$pars$rho,
      mfit.full$pars$time.bound, return.se = FALSE
    )

    test.boot <- mfit.full.iter$MSE.Y / mfit.reduce.iter$MSE.Y
    boot.all[iter, ] <- test.boot
  }

  p.value <- mean(boot.all <= test.stat)

  list(test.stat = test.stat, boot.stat = boot.all, p.value = p.value)
}
