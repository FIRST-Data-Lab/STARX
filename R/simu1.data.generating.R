#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates simulation example for STVCM.
#'
#' @importFrom mgcv fs.boundary inSide
#' @param nS Number of points randomly generated from spatial domain.
#' @param nT Number of times that spatial points are observed across the time.
#' @param sigma Noise level.
#' @return
#' A list of five containing the data to fit model.
#' \item{Y}{response variable.}
#' \item{Z}{exogenous variables with constant linear
#' coefficients}
#' \item{X}{exogenous variables with varying linear
#' coefficients.}
#' \item{location}{A \code{nS * nT} by three matrix with locations of data points.}
#' @export
simu1.data.generating <- function(nS, nT, sigma) {

  #### generate lattice ####
  coord.x <- seq(-0.89, 3.39, length.out = 5 * floor(sqrt(nS)))
  coord.y <- seq(-0.89, 0.89, length.out = 2 * floor(sqrt(nS)))

  coord.x0 <- 1:(5 * floor(sqrt(nS)))
  coord.y0 <- 1:(2 * floor(sqrt(nS)))

  coord <- expand.grid(coord.x, coord.y)
  coord0 <- expand.grid(coord.x0, coord.y0)

  fsb <- list(mgcv::fs.boundary())[[1]]
  names(fsb) <- c("coord.x", "coord.y")
  coord.x <- coord[, 1]
  coord.y <- coord[, 2]
  ind <- mgcv::inSide(fsb, x = coord.x, y = coord.y)

  # generate spatial locations
  ind.sample <- sample(which(ind), nS)
  S0 <- coord[ind.sample, ]

  # generate time grids
  T0 <- 1:nT

  # locations
  loc.Sample <- cbind(rep(S0[, 1], nT), rep(S0[, 2], nT), rep((1:nT) / nT, each = nS))

  #### evaluate at beta functions ####
  # beta0
  beta0 <- beta1.func(loc.Sample[, 1], loc.Sample[, 2], loc.Sample[, 3])
  # beta1
  beta1 <- beta2.func(loc.Sample[, 1], loc.Sample[, 2], loc.Sample[, 3])
  # beta2
  beta2 <- beta3.func(loc.Sample[, 1], loc.Sample[, 2], loc.Sample[, 3])

  # Covariates X
  n <- nS * nT
  X1 <- rnorm(n)
  X2 <- rnorm(n)

  # Response Y
  Y <- beta0 + beta1 * X1 + beta2 * X2 + rnorm(n, sd = sigma)

  dat <- list()
  dat$Y <- Y
  dat$Z <- NULL
  dat$X <- cbind(X1, X2)
  dat$location <- loc.Sample

  dat
}
