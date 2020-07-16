#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates sequence of spatial plots of estimated coefficient functions.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom sp point.in.polygon
#' @import TPST
#' @param fitted Estimated STAR model.
#' @param ngrid.x Number of grid points in x axis.
#' @param ngrid.y Number of grid points in y axis.
#' @param ngrid.t Number of grid points in t axis.
#' @param boundary Boundary of spatial domain.
#' @param boundary.t Boundary of time domain.
#' @param beta.true.ind A logical number incidating whether to draw the true coefficient functions.
#' @param beta.true A list of true coefficient functions.
#' @param coord.ratio Aspect ratio, expressed as y / x.
#' @return A list of plots
#' \item{beta.true.all}{A list of length \code{n.X} containing sequence of spatial plots
#' of each true coefficient functions.}
#' \item{beta.est.all}{A list of length \code{n.X} containing sequence of spatial plots
#' of each estiamted coefficient functions.}
#' @export
beta.func.plot <- function(fitted, ngrid.x, ngrid.y, ngrid.t, boundary, boundary.t,
                           beta.true.ind = TRUE, beta.true = NULL, coord.ratio = 1) {
  # the regression points
  xx <- seq(range(boundary[, 1])[1], range(boundary[, 1])[2], length.out = ngrid.x)
  yy <- seq(range(boundary[, 2])[1], range(boundary[, 2])[2], length.out = ngrid.y)
  tt <- seq(boundary.t[1], boundary.t[2], length.out = ngrid.t)

  ss <- expand.grid(xx, yy)

  data.grid <- data.frame(
    x = rep(ss[, 1], ngrid.t), y = rep(ss[, 2], ngrid.t),
    t = rep(tt, each = dim(ss)[1])
  )

  V <- fitted$pars$V
  Tri <- fitted$pars$Tri
  d <- fitted$pars$d
  r <- fitted$pars$r
  time.knots <- fitted$pars$time.knots
  rho <- fitted$pars$rho
  time.bound <- fitted$pars$time.bound
  n.Z <- fitted$n.Z
  n.X <- fitted$n.X

  beta.func.grid <- c()
  if(beta.true.ind == TRUE) {
    for (iter in 1:n.X) {
      beta.func.grid <-
        cbind(beta.func.grid,
              beta.true[[iter]](data.grid[, 1], data.grid[, 2], data.grid[, 3]))
    }
    ind <- which(sp::point.in.polygon(data.grid[, 1], data.grid[, 2],
                                      boundary[, 1], boundary[, 2]) == 1)
    beta.func.grid[-ind, ] <- NA
  }

  basis.grid <- TPST::basis.tensor(data.grid[, 1:2], data.grid[, 3],
                                   V = V, Tri = Tri, d = d, r = r,
                                   time.knots = time.knots,
                                   rho = rho, time.bound = time.bound)
  if(n.Z != 0) {
    theta.hat <- matrix(fitted$theta.hat[-c(1:n.Z)], ncol = n.X)
  } else {
    theta.hat <- matrix(fitted$theta.hat, ncol = n.X)
  }

  beta.hat.grid <- basis.grid$Psi.Q2 %*% theta.hat

  beta.true.all <- list()
  beta.est.all <- list()
  for (iter.X in 1:n.X) {
    if (beta.true.ind == TRUE) {
      data.true <- data.frame(x = data.grid$x, y = data.grid$y,
                              t = data.grid$t,
                              beta.true = beta.func.grid[, iter.X])
    }

    data.est <- data.frame(x = data.grid$x, y = data.grid$y,
                           t = data.grid$t,
                           beta.est = as.matrix(beta.hat.grid[, iter.X]))
    names(data.true) <- c("x", "y", "t", "value")
    names(data.est) <- c("x", "y", "t", "value")
    beta.true.plot <- list()
    beta.est.plot <- list()
    if (beta.true.ind == TRUE) {
      value.range.max <- max(c(data.true$value, data.est$value), na.rm = TRUE)
      value.range.min <- min(c(data.true$value, data.est$value), na.rm = TRUE)
    } else {
      value.range.max <- max(c(data.est$value), na.rm = TRUE)
      value.range.min <- min(c(data.est$value), na.rm = TRUE)
    }

    for(iter in 1:length(tt)) {
      if (beta.true.ind == TRUE) {
        beta.true.plot[[iter]] <-
          raster.plot(data.true[data.true$t == tt[iter], ],
                      c(value.range.min, value.range.max),
                      coord.ratio = coord.ratio, boundary)
      }
      beta.est.plot[[iter]] <-
        raster.plot(data.est[data.est$t == tt[iter], ],
                    c(value.range.min, value.range.max),
                    coord.ratio = coord.ratio, boundary)
    }
    beta.true <- NULL
    if (beta.true.ind == TRUE) {
    beta.true <- ggpubr::ggarrange(plotlist = beta.true.plot,
                                   labels = paste("t =", tt),
                                   ncol = 3, nrow = ceiling(length(tt)/3),
                                   common.legend = TRUE, legend = "right")
    }
    beta.est <- ggpubr::ggarrange(plotlist = beta.est.plot,
                                  labels = paste("t =", tt),
                                  ncol = 3, nrow = ceiling(length(tt)/3),
                                  common.legend = TRUE, legend = "right")
    beta.true.all[[iter.X]] <- beta.true
    beta.est.all[[iter.X]] <- beta.est
  }


  list(beta.true.all = beta.true.all, beta.est.all = beta.est.all)
}
