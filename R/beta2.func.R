#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function is the coefficient function in simulation example.
#'
#' @param x x coordinate
#' @param y y coordinate
#' @param z z coordinate
#' @export
beta2.func <- Vectorize(function(x, y, z) {
    2 * cos((0.5 * x + y^2)) * z
})
