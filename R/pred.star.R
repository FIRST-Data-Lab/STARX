pred.star <- function(mfit.star, newdata, n.Z, n.X, basis.newdata = NULL, pars_only = TRUE){
  pars <- mfit.star$pars
  if (is.null(basis.newdata)){
    basis.newdata <- basis.tensor(newdata$location[, 1:2], newdata$location[, 3],
                                  V = pars$V, Tri = pars$Tri, d = pars$d, r = pars$r,
                                  time.knots = pars$time.knots,
                                  rho = pars$rho, time.bound = pars$time.bound)
  }
  inside <- which(!is.na(basis.newdata$Psi[, 1]))
  eta.hat <- mfit.star$theta.hat[1:n.Z]
  theta.hat <- matrix(mfit.star$theta.hat[-c(1:n.Z)], ncol = n.X)
  beta.pred <- basis.newdata$Psi.Q2[inside,] %*% theta.hat
  
  if (pars_only){
    return(list(loc.inside = inside, eta.hat = eta.hat, beta.pred = beta.pred))
  }else{
    new.X <- newdata$X
    new.Z <- newdata$Z
    new.Y <- newdata$Y
    new.Weight <- newdata$W
    if (pars$intercept){
      new.Z <- as.matrix(cbind(1, new.Z))
    }
    n.X <- ncol(new.X)
    n.Z <- ncol(new.Z)
    
    Y.pred <- alpha.hat * new.Weight %*% new.Y + new.Z %*% eta.hat + new.X * beta.pred
    
    return(list(loc.inside = inside, eta.hat = eta.hat, beta.pred = beta.pred, Y.pred = Y.pred))
  }

}
