#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates simulation example for STAR with partially linear
#' varying coefficient model.
#'
#' @importFrom mgcv fs.boundary
#' @importFrom sp point.in.polygon
#' @param nS Number of points randomly generated from spatial domain.
#' @param nT Number of times that spatial points are observed across the time.
#' @param sigma Noise level.
#' @param alpha \eqn{\alpha} in STAR model.
#' @param k The weight matrix is generated based on k-nearest neighbor.
#' The default value in this simulation example is 10.
#' @param boundary The boundary of the spatial domain.
#' @return
#' A list of five containing the data to fit model.
#' \item{Y}{response variable.}
#' \item{Z}{exogenous variables with constant linear
#' coefficients}
#' \item{X}{exogenous variables with varying linear
#' coefficients.}
#' \item{location}{A \code{nS * nT} by three matrix with locations of data points.}
#' \item{W}{weight matrix in STAR.}
#' @param type distribution type of data points.
#' @export
simu2.data.generating <- function(nS, nT, sigma, alpha, k = 10, boundary,
                                  type = "uniform") {

  if (type == "uniform"){
    coord.x <- runif(3 * nS) * 4.5 - 1
    coord.y <- runif(3 * nS) * 2 - 1

    # generate time grids
    t.all <- rep((1:nT) / nT, each = nS)
  }

  if (type == "space_uneven") {
    coord.x <- rnorm(5 * nS, sd = 2) + 3
    coord.y <- rnorm(5 * nS)

    # generate time grids
    t.all <- rep((1:nT) / nT, each = nS)
  }

  if (type == "time_uneven") {
    coord.x <- rnorm(5 * nS, sd = 2) + 3
    coord.y <- rnorm(5 * nS)

    # generate time grids
    nT0 <- floor(nT * 1.5)
    nT.sample <- sample(1:nT0, nT)
    t.all <- rep(nT.sample / nT0, each = nS)
  }
  coord <- cbind(coord.x, coord.y)
  fsb <- list(mgcv::fs.boundary())[[1]]
  names(fsb) <- c("coord.x", "coord.y")
  coord.x <- coord[, 1]
  coord.y <- coord[, 2]
  ind <- which(sp::point.in.polygon(coord.x, coord.y,
                                boundary[, 1], boundary[, 2]) == 1)
  S0 <- coord[ind, ]
  S0 <- S0[1:nS, ]



  # locations
  loc.Sample <- cbind(rep(S0[, 1], nT), rep(S0[, 2], nT), t.all)

  eta0 <- 5
  eta1 <- 1
  eta2 <- -1
  #### evaluate at beta functions ####
  # beta1
  beta1 <- beta1.func(loc.Sample[, 1], loc.Sample[, 2], loc.Sample[, 3])
  # beta2
  beta2 <- beta2.func(loc.Sample[, 1], loc.Sample[, 2], loc.Sample[, 3])

  n <- nS * nT
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- rnorm(n)
  eps <- rnorm(n, 0, sigma)

  Ystar <- eta0 + eta1 * z1 + eta2 * z2 + x1 * beta1 + eps

  # generate W.matrix
  W <- w.knn(loc.Sample, k = k)

  # Response Variable
  A <- diag(n) - alpha * W
  Y <- solve(as.matrix(A), Ystar)

  # generate data to fit the model
  dat <- list()
  dat$Y <- Y
  dat$Z <- cbind(z1,z2)
  dat$X <- matrix(x1, ncol = 1)
  dat$location <- loc.Sample
  dat$W <- W

  return(dat)
}
