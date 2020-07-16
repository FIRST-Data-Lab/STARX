#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function calculates standard deviation of parameters in STAR models.
#'
#' @import Matrix
se.star.eta <- function(Ystar, X, Z, W,
                        energQ21, energQ22, Weight, Ai,
                        PsiX0, mfit, Weight.type = "normal") {
  n.Z <- ncol(Z)
  n <- nrow(Z)
  if (Weight.type == "normal") {
    G <- Weight %*% Matrix::solve(Ai)
  } else if (Weight.type == "upper"){
    Ai[lower.tri(Ai)]  <- 0
    G <- Weight %*% bdsmatrix::backsolve(Ai, diag(dim(Ai)[1]))
  }

  if(!is.null(X)) {
    EnergQ21 <- Matrix::kronecker(diag(ncol(X)), energQ21)
    EnergQ22 <- Matrix::kronecker(diag(ncol(X)), energQ22)
    tmp <- Matrix::solve(Matrix::crossprod(PsiX0) +
                   mfit$best.lambda[1] * EnergQ21 + mfit$best.lambda[2] * EnergQ22, Matrix::t(PsiX0))

    G.diag <- Matrix::diag(G)
    tmp.G <- G %*% (W %*% mfit$THETA)
    P.mu <- tmp.G - PsiX0 %*% (tmp %*% tmp.G)
    P.Z <- Z - PsiX0 %*% (tmp %*% Z)
  } else{
    G.diag <- Matrix::diag(G)
    P.mu <- G %*% (W %*% mfit$THETA)
    P.Z <- Z
  }

  tr.G <- sum(G.diag)
  tr.G2 <- sum(c(as.matrix(G))*c(t(as.matrix(G))))
  MSE <- Matrix::mean((Ystar - W %*% mfit$THETA)^2)

  # calculate sigma matrix
  Sigma.all <- matrix(0, ncol = n.Z + 2, nrow = n.Z + 2)
  Sigma.all[1, 1] <- as.matrix(tr.G2 + sum(G^2) + Matrix::crossprod(P.mu) / MSE)
  Sigma.all[2, 2] <- n / (2 * MSE^2)
  Sigma.all[1, 2] <- Sigma.all[2, 1] <- tr.G / MSE
  Sigma.all[1, 3:(n.Z + 2)] <- Sigma.all[3:(n.Z + 2), 1] <- as.vector(Matrix::crossprod(P.Z, P.mu) / MSE)
  Sigma.all[3:(n.Z + 2), 3:(n.Z + 2)] <-
    as.matrix(Matrix::crossprod(P.Z) / MSE)
  Sigma.all.inv <- solve(Sigma.all)


  # calculate omega matrix
  mu3 <- mean(as.matrix(Ystar - W %*% mfit$THETA)^3)
  mu4 <- mean(as.matrix(Ystar - W %*% mfit$THETA)^4) - 3 * (MSE^2)
  Omega.all <- matrix(0, ncol = n.Z + 2, nrow = n.Z + 2)
  Omega.all[1, 1] <- mu4 * sum(G.diag^2) / (MSE^2)  + 2 * mu3 / (MSE^2) * sum(G.diag * P.mu)
  Omega.all[2, 1] <- Omega.all[1, 2] <- mu4 / (2 * MSE^3) * tr.G + mu3 * sum(P.mu)
  Omega.all[2, 2] <- mu4 / (4 * MSE^4)
  Omega.all[1, 3:(n.Z + 2)] <- Omega.all[3:(n.Z + 2), 1] <- mu3 / (2 * MSE^2) * Matrix::colSums(G.diag * P.Z)
  Omega.all[2, 3:(n.Z + 2)] <- Omega.all[3:(n.Z + 2), 2] <- mu3 / (2 * MSE^3) * Matrix::colSums(P.Z)


  # calculate standard error
  se.all <- sqrt(diag(Sigma.all.inv %*% Omega.all %*% Sigma.all.inv + Sigma.all.inv))
  se.all
}

