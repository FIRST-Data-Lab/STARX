#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates testing sets for blockwise cross-validation
#' in spatiotemporal data analysis.
#'
#' @import raster
#' @import sp
#' @param locations A matrix with three columns stores the
#' locations of data points. The first two columns are spatial locations.
#' @param nblock.x Number of bins in x-dimension.
#' @param nblock.y Number of bins in y-dimension.
#' @param nblock.t Number of bins in t-dimension. The total number of
#' blocks is equal to \code{nblock.x * nblock.y * nblock.t}
#' @param nfold Number of folds in cross-validation.
#' @return
#' \item{Fold}{A list of length \code{nfold}. Each list contains the
#' index of data points in the testing set.}
#' @export

block.cv <- function(locations, nblock.x, nblock.y, nblock.t, nfold) {

  locations <- data.frame(locations)
  locations.s <- unique(locations[, 1:2])
  names(locations) <- c("xx", "yy", "tt")
  names(locations.s) <- c("xx", "yy")
  sp::coordinates(locations) <- ~ xx + yy
  sp::coordinates(locations.s) <- ~ xx + yy

  net <- rasterNet(locations.s, xbin = nblock.x, ybin = nblock.y)
  # plot(net)
  subBlocks <- raster::intersect(net, locations.s)
  # print(plot(subBlocks))
  nrowBlocks <- length(subBlocks)

  # identify spatial block for each point
  locations$block.s <- 1
  for (iter in 1:nrowBlocks) {
    sp.over <- sp::over(locations, subBlocks[iter, ])
    locations$block.s[which(!is.na(sp.over[, 1]))] <- iter
  }

  # identify temporal block for each point
  locations$block.t <- 1
  locations$block.t <- as.numeric(cut(locations$tt, breaks = seq(0 - 0.001, max(locations$tt) + 0.001,
    length.out = nblock.t + 1
  ), labels = 1:nblock.t))

  # identify spatiotemporal block for each point
  locations$block <- (locations$block.t - 1) * nrowBlocks + locations$block.s

  # identify cross-validation fold for each point
  locations$fold <- 1
  Folds <- createFolds(1:max(locations$block), k = nfold)
  Fold <- list()
  for (iter in 1:nfold) {
    Fold[[iter]] <- which(locations$block %in% Folds[[iter]])
  }
  return(Fold)
}
