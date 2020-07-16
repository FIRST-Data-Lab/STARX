stvcm.cv <- function(Ystari, coords.S, coords.T, Lambda, PsiX0, EnergQ21, EnergQ22) {

  n <- length(Ystari)
  set.seed(50)
  cv.folds <- block.cv(cbind(coords.S, coords.T), 20, 10, 10, 5)
  MSPE <- c()
  if (nrow(Lambda) > 1) {
    for (iter in 1:length(cv.folds)) {
      test <- cv.folds[[iter]]
      train <- (1:n)[-test]
      PsiX0.train <- PsiX0[train, ]
      Y.train <- Ystari[train]

      PsiX0.test <- PsiX0[test, ]
      Y.test <- Ystari[test]
      PP.train <- Matrix::crossprod(PsiX0.train)
      PsiX0Y.train <- Matrix::crossprod(PsiX0.train, Y.train)
      mspe <- apply(Lambda, 1, function(x) {
        sum((Y.test - PsiX0.test %*% Matrix::solve(PP.train +
                                             as.numeric(x[1]) * EnergQ21 + as.numeric(x[2]) * EnergQ22, PsiX0Y.train))^2)
      })
      MSPE <- rbind(MSPE, mspe)
    }
    index <- which.min(colSums(MSPE))

    best.lambda <- as.numeric(Lambda[index, ])
  } else {
    best.lambda <- as.numeric(Lambda[1, ])
    print(best.lambda)
  }

  PP <- Matrix::crossprod(PsiX0)
  PY <- Matrix::crossprod(PsiX0, Ystari)
  THETA <- Matrix::solve(PP + best.lambda[1] * EnergQ21 + best.lambda[2] * EnergQ22, PY)
  MSE <- Matrix::mean((Ystari - PsiX0 %*% THETA)^2)
  list(THETA = THETA, MSE = MSE, best.lambda = best.lambda)
}
