#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function fits a spatiotemporal varying coefficient model.
#'
#' @import TPST Matrix stats MGLM rgeos
#' @param data A list of five containing the data to fit model. \code{Y} is the
#' response variable. \code{Z} is exogenous variables with constant linear
#' coefficients. \code{X} is exogenous variables with varying linear
#' coefficients. \code{location} is the location of data points.
#' @param Lambda A matrix of tuning parameters of the smoothness penalty.
#' @param V The \code{N} by two matrix of verities of a triangulation,
#' where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' @param Tri The triangulation matrix of dimension \code{nTr} by three,
#' where \code{nTr} is the number of triangles in the triangulation.
#' Each row is the indices of vertices in \code{V}.
#' @param d The degree of piecewise polynomials -- default is 2,
#' and usually \code{d} is greater than one. -1 represents piecewise constant.
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' @param time.knots The vector of interior time.knots for univariate spline.
#' @param time.bound The vector of two. The boundary of univariate spline.
#' @param rho The order of univariate spline.
#' @param intercept.X A logical number indicating whether to include varying intercept term
#' in the model -- default is \code{TRUE}.
#' @param intercept.Z A logical number indicating whether to include constant intercept term
#' in the model -- default is \code{FALSE}.
#' @param return.se A logical number indicating whether to calculate standard deviation
#' of the linear estimators.
#' @return
#' \item{n.Z}{number of linear coefficients.}
#' \item{n.X}{number of varying coefficient functions.}
#' \item{best.lambda}{the selected smoothing tuning parameters.}
#' \item{theta.hat}{Estimated coefficents.}
#' \item{gamma.hat}{Estimated tensor spline coefficients.}
#' \item{beta.hat}{Estimated tensor spline coefficients with matrix \eqn{Q_2}.}
#' \item{se.eta}{Standard deviation of \eqn{\alpha, \sigma^2} and linear coefficients.}
#' \item{Y.hat}{Estimated response variable.}
#' \item{mse}{Estimated \eqn{\sigma^2}.}
#' \item{mle}{MLE of fitted model.}
#' \item{MSE.Y}{Mean squared error of response variable.}
#' @export
stvcm.fit <- function(data, Lambda, V, Tri, d, r, time.knots, rho,
                     time.bound, intercept.X = TRUE, intercept.Z = FALSE,
                     return.se = FALSE) {

  Y <- data$Y
  Z <- data$Z
  X <- data$X

  coords <- data$location

  if (intercept.Z == TRUE) {
    Z <- as.matrix(cbind(rep(1, length(Y)), Z))
  } else if(!is.null(Z)){
    Z <- as.matrix(Z)
  }

  if (intercept.X == TRUE) {
    X <- as.matrix(cbind(rep(1, length(Y)), X))
  } else {
    X <- as.matrix(X)
  }

  n.X = 0
  n.Z = 0
  n <- length(Y)
  n.X <- ncol(X)
  if(!is.null(Z)) n.Z <- ncol(Z)
  coords <- as.matrix(coords)

  if(!is.null(X)) {
    X <- as.matrix(X)
    Basis <- TPST::basis.tensor(coords[, 1:2], coords[, 3], V, Tri, d, r, time.knots,
                                rho = rho, time.bound = time.bound)
    inside <- which(!is.na(Basis$Psi[, 1]))
    X <- matrix(X[inside, ], ncol = n.X)
    Z <- Z[inside, ]
    Y <- Y[inside]
    PsiX0 <- kr(X, Basis$Psi.Q2[inside, ])
    Q2 <- Basis$Q2.all

    # generate penalty function
    energQ21 <- as.matrix(Basis$D1)
    energQ22 <- as.matrix(Basis$D2)

    ncol.Psi <- ncol(PsiX0)
    EnergQ21 <- matrix(0, (n.Z + ncol.Psi), (n.Z + ncol.Psi))
    EnergQ21[(n.Z + 1):(n.Z + ncol.Psi), (n.Z + 1):(n.Z + ncol.Psi)] <-
      kronecker(diag(ncol(X)), energQ21)

    EnergQ22 <- matrix(0, (n.Z + ncol.Psi), (n.Z + ncol.Psi))
    EnergQ22[(n.Z + 1):(n.Z + ncol.Psi), (n.Z + 1):(n.Z + ncol.Psi)] <-
      kronecker(diag(ncol(X)), energQ22)

    W <- as.matrix(cbind(Z, PsiX0))
    PPsiX <- crossprod(W)
  } else {
    W <- Z
  }


  if (!is.null(X)) {
    mfit <- stvcm.cv(Y, coords[, 1:2], coords[, 3], Lambda,
                     PsiX0 = W, EnergQ21 = EnergQ21, EnergQ22 = EnergQ21)
    MSE <- mfit$MSE
    theta.hat <- mfit$THETA
    best.lambda <- mfit$best.lambda
  } else {
    Y <- as.matrix(Y)
    mfit <- stats::lm(Y ~ W + 0)
    MSE <- mean(mfit$residuals^2)
    theta.hat <- mfit$coefficients
    mfit$THETA <- mfit$coefficients
    best.lambda <- NULL
  }

  MLE <- -n / 2 * (log(2 * pi) + 1) - n * log(sqrt(MSE))
  Y.hat <- W %*% theta.hat
  MSE.Y <- Matrix::mean((Y - Y.hat)^2)
  #### calculate standard deviation ####
  se.eta <- NULL
  if(return.se == TRUE & n.Z != 0){
    if(!is.null(X)){
      Z = matrix(c(Z), ncol = n.Z)
      EnergQ21 <- Matrix::kronecker(diag(n.X), energQ21)
      EnergQ22 <- Matrix::kronecker(diag(n.X), energQ22)
      tmp <- Matrix::solve(Matrix::crossprod(PsiX0) +
                             mfit$best.lambda[1] * EnergQ21 +
                             mfit$best.lambda[2] * EnergQ22, Matrix::t(PsiX0))
      P.Z <- Z - PsiX0 %*% (tmp %*% Z)
      Sigma.all <-
        as.matrix(Matrix::crossprod(P.Z) / MSE)
      Sigma.all.inv <- solve(Sigma.all)
      se.eta <- sqrt(diag(Sigma.all.inv))
     } else {
      PsiX0 <- NULL
      Z = matrix(c(Z), ncol = n.Z)
      Sigma.all <- as.matrix(Matrix::crossprod(Z) / MSE)
      Sigma.all.inv <- solve(Sigma.all)
      se.eta <- sqrt(diag(Sigma.all.inv))
    }
  }


  if(!is.null(X)) {
    if(n.Z != 0){
      gamma.hat <- Q2 %*% matrix(theta.hat[-(1:n.Z)], ncol = n.X)
      beta.hat <- Basis$Psi.Q2 %*% matrix(theta.hat[-c(1:n.Z)], ncol = n.X)
    } else {
      gamma.hat <- Q2 %*% matrix(theta.hat, ncol = n.X)
      beta.hat <- Basis$Psi.Q2 %*% matrix(theta.hat, ncol = n.X)
    }

    pars <- list(V = V, Tri = Tri, d = d, r = r,
                 time.knots = time.knots, rho = rho,
                 time.bound = time.bound,
                 intercept.X = intercept.X, intercept.Z = intercept.Z)
  } else {
    pars <- list(intercept.X = intercept.X, intercept.Z = intercept.Z)
    gamma.hat <- NULL
    beta.hat <- NULL
  }

  list(n.Z = n.Z, n.X = n.X, best.lambda  = best.lambda,
       theta.hat = theta.hat,
       gamma.hat = gamma.hat, beta.hat = beta.hat,
       se.eta = se.eta,
       Y.hat = Y.hat, mse = MSE, mle = MLE, MSE.Y = MSE.Y,
       pars = pars
  )
}
