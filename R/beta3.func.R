#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function is the coefficient function in simulation example.
#'
#' @param x x coordinate
#' @param y y coordinate
#' @param z z coordinate
#' @export
beta3.func <- Vectorize(function(x, y, z) {
  2 * sin(pi * y * (z - 0.5))
})
