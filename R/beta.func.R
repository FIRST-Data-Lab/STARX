#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function is the coefficient function in simulation example.
#' @import mgcv
#' @param x x coordinate
#' @param y y coordinate
#' @param t t coordinate

beta1.func <- function(x, y, z) {
  2 * fs.test(x, y, b = 1) * (z - 0.5)^2
}
beta1.func <- Vectorize(beta1.func)

beta2.func <- function(x, y, z) {
  2 * cos((0.5 * x + y^2)) * z
}
beta2.func <- Vectorize(beta2.func)

# beta2_func=function(x,y,z){
#   2*sin(pi*y*(z-0.5))
# }
# beta2_func=Vectorize(beta2_func)
#
# sigma.simu=function(t,a){
#   0.5/(1+exp(t-2.5))+a
# }
# sigma.simu=Vectorize(sigma.simu)
