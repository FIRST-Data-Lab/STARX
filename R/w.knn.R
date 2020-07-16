#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates weight matrix based on k nearest neighbor.
#' @import FNN
w.knn <- function(coords, k = 10) {
  coords <- as.matrix(coords)
  n.coords <- nrow(coords)
  i.vector <- rep(1:n.coords, each = k)
  j.vector <- c(FNN::get.knn(coords, k = k)$nn.index)
  W <- sparseMatrix(
    i = i.vector,
    j = j.vector, x = 1 / k,
    dims = c(n.coords, n.coords)
  )
  W
}
