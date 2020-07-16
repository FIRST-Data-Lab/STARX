#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function is the coefficient function in simulation example.
#' @importFrom mgcv fs.test
#' @param x x coordinate
#' @param y y coordinate
#' @param z z coordinate
#' @export
beta1.func <- Vectorize(function(x, y, z) {
  2 * mgcv::fs.test(x, y, b = 1) * (z - 0.5)^2
})

