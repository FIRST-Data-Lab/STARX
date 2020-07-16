#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function evaluates the performance of estimated STAR model by
#' comparing with the true model.
#'
#' @importFrom TPST basis.tensor
#' @param fitted Estimated STAR model.
#' @param beta.true True beta functions.
#' @param eta.true True linear coefficients \eqn{\eta}.
#' @param alpha.true True \eqn{\alpha}.
#' @param data.grid Grid points to evaluate beta functions.
#' @return
#' \item{mse.alpha}{squared error of \eqn{\alpha}.}
#' \item{mse.eta}{squared error of \eqn{\eta}.}
#' \item{mise.beta}{integrated squared error of \eqn{\beta(s1, s2, t)}.}
#' \item{beta.hat.grid}{Estimated \eqn{\beta(s1, s2, t)} at grid points
#' \code{data.grid}.}
#' \item{beta.func.grid}{True \eqn{\beta(s1, s2, t)} at grid points
#' \code{data.grid}.}
#' @export
eval.star <- function(fitted, beta.true = NULL, eta.true = NULL,
                      alpha.true = NULL, data.grid = NULL){

  V <- fitted$pars$V
  Tri <- fitted$pars$Tri
  d <- fitted$pars$d
  r <- fitted$pars$r
  time.knots <- fitted$pars$time.knots
  rho <- fitted$pars$rho
  time.bound <- fitted$pars$time.bound
  n.Z <- fitted$n.Z
  n.X <- fitted$n.X

  mse.alpha = NULL
  mse.eta = NULL
  mise.beta = NULL
  beta.hat.grid = NULL
  beta.func.grid = NULL

  if(!is.null(alpha.true)) mse.alpha <- (fitted$alpha.hat - alpha.true)^2
  if(!is.null(eta.true)) mse.eta <- (fitted$theta.hat[1:n.Z] - eta.true)^2

  beta.func.grid <- c()
  if(!is.null(beta.true)) {
    for (iter in 1:n.X) {
      beta.func.grid <-
        cbind(beta.func.grid,
              beta.true[[iter]](data.grid[, 1], data.grid[, 2], data.grid[, 3]))
    }
    basis.grid <- TPST::basis.tensor(data.grid[, 1:2], data.grid[, 3],
                                     V = V, Tri = Tri, d = d, r = r,
                                     time.knots = time.knots,
                                     rho = rho, time.bound = time.bound)

    if(n.Z != 0) {
      theta.hat <- matrix(fitted$theta.hat[-c(1:n.Z)], ncol = n.X)
    } else {
      theta.hat <- matrix(fitted$theta.hat, ncol = n.X)
    }
    beta.hat.grid <- basis.grid$Psi.Q2 %*% theta.hat

    mise.beta <- colMeans((beta.hat.grid - beta.func.grid)^2, na.rm = TRUE)
  }

  list(mse.alpha = mse.alpha, mse.eta = mse.eta, mise.beta = mise.beta,
       beta.hat.grid = beta.hat.grid, beta.func.grid = beta.func.grid)
}
