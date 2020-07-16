#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function fits a spatiotemporal autoregressive (STAR) model
#' with exogenous variables.
#'
#' @import TPST Matrix MGLM splines2
#' @param data A list of five containing the data to fit model. \code{Y} is the
#' response variable. \code{Z} is exogenous variables with constant linear
#' coefficients. \code{X} is exogenous variables with varying linear
#' coefficients. \code{location} is the location of data points. \code{W} is the
#' weight matrix in STAR.
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
#' @param intercept A logical number indicating whether to include intercept term
#' in the model -- default is \code{TRUE}.
#' @param det.func A function used to calculate determinate of a matrix --
#' default is \code{det} in the R package \code{Matrix}.
#' @param return.se A logical number indicating whether to calculate standard deviation
#' of the linear estimators.
#' @param Weight.type The type of weight matrix. "\code{normal}" means
#' the weight matrix does not have special structure. "\code{upper}" means
#' the weight matrix is upper triangle matrix.
#' @return
#' \item{n.Z}{number of linear coefficients.}
#' \item{n.X}{number of varying coefficient functions.}
#' \item{best.lambda}{the selected smoothing tuning parameters.}
#' \item{theta.hat}{Estimated coefficents.}
#' \item{alpha.hat}{Estimated alpha in STAR.}
#' \item{gamma.hat}{Estimated tensor spline coefficients.}
#' \item{beta.hat}{Estimated tensor spline coefficients with matrix \eqn{Q_2}.}
#' \item{se.eta}{Standard deviation of \eqn{\alpha, \sigma^2} and linear coefficients.}
#' \item{Y.hat}{Estimated response variable.}
#' \item{mse}{Estimated \eqn{\sigma^2}.}
#' \item{mle}{MLE of fitted model.}
#' \item{MSE.Y}{Mean squared error of response variable.}
#' @export
star.fit <- function(data, Lambda, V, Tri, d, r, time.knots, rho,
                     time.bound, intercept = TRUE, det.func = Matrix::det,
                     return.se = FALSE, Weight.type = "normal") {

  Y <- data$Y
  Z <- data$Z
  X <- data$X

  coords <- data$location
  Weight <- data$W

  if (intercept == TRUE) {
    Z <- as.matrix(cbind(rep(1, length(Y)), Z))
  } else {
    Z <- as.matrix(Z)
  }


  n <- length(Y)
  n.X <- ncol(X)
  n.Z <- ncol(Z)
  coords <- as.matrix(coords)

  if(!is.null(X)) {
    X <- as.matrix(X)
    Basis <- TPST::basis.tensor(coords[, 1:2], coords[, 3], V, Tri, d, r, time.knots,
                          rho = rho, time.bound = time.bound)
    inside <- which(!is.na(Basis$Psi[, 1]))
    X <- matrix(X[inside, ], ncol = n.X)
    Z <- Z[inside, ]
    Y <- Y[inside]
    Weight <- Weight[inside, inside]
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


  mle_alpha <- c()
  mse_alpha <- c()

  Det.all <- c()
  start <- -0.99
  end <- 0.99
  change.start <- 1
  change.end <- 1
  while (abs(start - end) > 0.002) {
    if (change.start == 1) {

      Ai.start <- -start * Weight
      Matrix::diag(Ai.start) <- Matrix::diag(Ai.start) + 1
      Ystari <- Ai.start %*% Y
      if (!is.null(X)) {
        FittedModel <- stvcm.cv(Ystari, coords.S = coords[, 1:2],
                                coords.T = coords[, 3], Lambda,
                                PsiX0 = W, EnergQ21 = EnergQ21,
                                EnergQ22 = EnergQ21)
        MSE <- FittedModel$MSE
      } else {
        Ystari <- as.matrix(Ystari)
        FittedModel <- stats::lm(Ystari ~ W + 0)
        MSE <- mean(FittedModel$residuals^2)
      }
      MLE.start <- -n / 2 * (log(2 * pi) + 1) - n * log(sqrt(MSE)) +
        log(abs(det.func(Ai.start)))

    }

    if (change.end == 1) {
      Ai.end <- -end * Weight
      Matrix::diag(Ai.end) <- Matrix::diag(Ai.end) + 1
      Ystari <- Ai.end %*% Y
      if (!is.null(X)) {
        FittedModel <- stvcm.cv(Ystari, coords[, 1:2], coords[, 3],
          Lambda, PsiX0 = W, EnergQ21 = EnergQ21, EnergQ22 = EnergQ21)
        MSE <- FittedModel$MSE
      } else {
        Ystari <- as.matrix(Ystari)
        FittedModel <- stats::lm(Ystari ~ W + 0)
        MSE <- mean(FittedModel$residuals^2)
      }
      MLE.end <- -n / 2 * (log(2 * pi) + 1) - n * log(sqrt(MSE)) +
        log(abs(det.func(Ai.end)))
    }


    if (MLE.end >= MLE.start) {
      start <- start + (end - start) / 2
      change.start <- 1
      change.end <- 0
    } else {
      end <- end - (end - start) / 2
      change.start <- 0
      change.end <- 1
    }

  }
  alpha.hat <- end

  Ai <- -alpha.hat * Weight
  Matrix::diag(Ai) <- Matrix::diag(Ai) + 1
  Ystar <- Ai %*% Y

  if (!is.null(X)) {
    mfit <- stvcm.cv(Ystar, coords[, 1:2], coords[, 3], Lambda,
                     PsiX0 = W, EnergQ21 = EnergQ21, EnergQ22 = EnergQ21)
    MSE <- mfit$MSE
    theta.hat <- mfit$THETA
    best.lambda <- mfit$best.lambda
  } else {
    Ystar <- as.matrix(Ystar)
    mfit <- stats::lm(Ystar ~ W + 0)
    MSE <- mean(mfit$residuals^2)
    theta.hat <- mfit$coefficients
    mfit$THETA <- mfit$coefficients
    best.lambda <- NULL
  }

  MLE <- -n / 2 * (log(2 * pi) + 1) - n * log(sqrt(MSE)) + log(abs(det.func(Ai)))
  Y.hat <- W %*% theta.hat + alpha.hat * Weight %*% Y
  MSE.Y <- Matrix::mean((Y - Y.hat)^2)
  #### calculate standard deviation ####
  se.eta <- NULL
  time.start <- proc.time()
  # save(file = paste0('application/data_for_se_',date.est + 1, '_', est.h,'_d', d, '_rho', rho, '_V', VT.infec, '_phi', phi,'.rda'),
  #      Ystar, X, Z, W, energQ21, energQ22, Weight, Ai, PsiX0, mfit)
  if(return.se == TRUE){
     if(!is.null(X)){
       Z = matrix(c(Z), ncol = n.Z)
       se.eta <- se.star.eta(Ystar, X, Z, W,
                             energQ21, energQ22, Weight, Ai,
                             PsiX0, mfit, Weight.type = Weight.type)
     } else {
       PsiX0 <- NULL
       Z = matrix(c(Z), ncol = n.Z)
       se.eta <- se.star.eta(Ystar, X, Z, W,
                             energQ21, energQ22, Weight, Ai,
                             PsiX0, mfit, Weight.type = Weight.type)
    }
  }


  if(!is.null(X)) {
    gamma.hat <- Q2 %*% matrix(theta.hat[-(1:n.Z)], ncol = n.X)
    beta.hat <- Basis$Psi.Q2 %*% matrix(theta.hat[-c(1:n.Z)], ncol = n.X)
    pars <- list(V = V, Tri = Tri, d = d, r = r,
                 time.knots = time.knots, rho = rho,
                 time.bound = time.bound, intercept = intercept)
  } else {
    pars <- list(intercept = intercept)
    gamma.hat <- NULL
    beta.hat <- NULL
  }

  list(n.Z = n.Z, n.X = n.X, best.lambda  = best.lambda,
    theta.hat = theta.hat, alpha.hat = alpha.hat,
    gamma.hat = gamma.hat, beta.hat = beta.hat,
    se.eta = se.eta,
    Y.hat = Y.hat, mse = MSE, mle = MLE, MSE.Y = MSE.Y,
    pars = pars
  )
}
