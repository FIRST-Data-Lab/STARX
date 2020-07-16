## ----Packages-----------------------------------------------------------------
rm(list = ls())
library(Triangulation)
library(TPST)
library(STARX)

## ---- warning = FALSE---------------------------------------------------------
# set up
nS <- 100
nT <- 30
sigma <- 0.5
alpha <- 0.5

# boundary of space and time
data("horseshoe.b")
plot(horseshoe.b$V1, horseshoe.b$V2, type = "l", xlab = '', ylab = '')
time.bound <- c(0, 1)

# generate simulation data
data.simu <- simu.data.generating(nS, nT, sigma, alpha, k = 10, horseshoe.b)

## -----------------------------------------------------------------------------
# tensor product spline degree and order
N <- 4 # Number of interior knots.
d <- 2 # Degree of bivariate spline.
r <- 1 # smoothness of bivariate spline.
rho <- 3 # order of univariate spline.

# load triangulation
data("Tr.hs.1")
data("V.hs.1")
TriPlot(V.hs.1, Tr.hs.1)

# knots
probs <- seq(time.bound[1], time.bound[2], length.out = N + 2)
time.knots <- quantile(data.simu$location[, 3], probs = probs)[-c(1, N + 2)]

time.knots

## ----model fit, warning = FALSE-----------------------------------------------
# fit model -------------------------------------------------------------
Lambda1 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda2 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda <- expand.grid(Lambda1, Lambda2)
# 
# mfit1 <- star.fit(data = data.simu, Lambda, V.hs.1, Tr.hs.1, d, r,
#                   time.knots, rho, time.bound, return.se = TRUE)
# save(mfit1, file = "/Users/shanyu/Dropbox/Research/STVCM2018/Packages/STARX/vignettesmfit1.rda")

## ---- warning = FALSE, fig.width = 4, fig.height = 3, message = FALSE---------
load("/Users/shanyu/Dropbox/Research/STVCM2018/Packages/STARX/vignettesmfit1.rda")
# Estimated parameters in STAR-PLVCM.
c(mfit1$alpha.hat, mfit1$mse, mfit1$theta.hat[1:mfit1$n.Z])
# Estimated standard deviation
mfit1$se.eta

# sequence of spatial plots of estimated coefficient functions.
ngrid.x <- 40
ngrid.y <- 20
ngrid.t <- 6

fitted.beta <- beta.func.plot(fitted = mfit1, ngrid.x, ngrid.y, ngrid.t, 
               boundary = horseshoe.b, boundary.t = time.bound,
               beta.true.ind = TRUE, beta.true = list(beta1.func), coord.ratio = 2)

fitted.beta$beta.true.all[[1]]
fitted.beta$beta.est.all[[1]]

## ---- warning = FALSE---------------------------------------------------------
# model evaluation ------------------------------------------------------
ngrid.x <- 40
ngrid.y <- 20
ngrid.t <- 6

# the regression points
xx <- seq(-0.89, 3.39, length.out = ngrid.x)
yy <- seq(-0.89, 0.89, length.out = ngrid.y)
ss <- expand.grid(xx, yy)
tt <- (0:(ngrid.t - 1)) / (ngrid.t - 1)

data.grid <- data.frame(
  x = rep(ss[, 1], ngrid.t), y = rep(ss[, 2], ngrid.t),
  t = rep(tt, each = dim(ss)[1]), z1 = 1, x1 = 1, x2 = 1
)

eta.true <- c(5, 1, -1)
beta.true <- list(beta1.func)
alpha.true <- 0.5

result <- eval.star(fitted = mfit1, beta.true, eta.true,
                    alpha.true, data.grid)

# MSE of alpha
result$mse.alpha
# MSE of eta
result$mse.eta
# MISE of coefficient functions
result$mise.beta

## ---- warning = FALSE---------------------------------------------------------
data.simu.lm <- data.simu
data.simu.lm$X <- NULL

mfit2 <- star.fit(data = data.simu.lm, Lambda, V.hs.1, Tr.hs.1, d, r,
                  time.knots, rho, time.bound, return.se = TRUE)

## -----------------------------------------------------------------------------
# Estimated parameters in STAR-PLVCM.
c(mfit2$alpha.hat, mfit2$mse, mfit2$theta.hat[1:mfit2$n.Z])
# Estimated standard deviation
mfit2$se.eta

## ----source package, warning = FALSE------------------------------------------
data.simu.vcm <- data.simu
data.simu.vcm$Z <- NULL

# mfit3 <- star.fit(data = data.simu.vcm, Lambda, V.hs.1, Tr.hs.1, d, r,
#                   time.knots, rho, time.bound, return.se = TRUE)
# 
# save(mfit3, file = "/Users/shanyu/Dropbox/Research/STVCM2018/Packages/STARX/mfit3.rda")

## ---- warning = FALSE, fig.width = 4, fig.height = 3, message = FALSE---------
load("/Users/shanyu/Dropbox/Research/STVCM2018/Packages/STARX/mfit3.rda")
# Estimated parameters in STAR-PLVCM.
c(mfit3$alpha.hat, mfit3$mse, mfit3$theta.hat[1:mfit3$n.Z])
# Estimated standard deviation
mfit3$se.eta

# sequence of spatial plots of estimated coefficient functions.
ngrid.x <- 40
ngrid.y <- 20
ngrid.t <- 6

fitted.beta3 <- beta.func.plot(fitted = mfit3, ngrid.x, ngrid.y, ngrid.t, 
               boundary = horseshoe.b, boundary.t = time.bound,
               beta.true.ind = TRUE, beta.true = list(beta1.func), coord.ratio = 2)

fitted.beta3$beta.true.all[[1]]
fitted.beta3$beta.est.all[[1]]

## ---- warning = FALSE---------------------------------------------------------
data.simu.stvcm <- data.simu
data.simu.stvcm$W <- NULL
mfit4 <- stvcm.fit(data = data.simu.stvcm, Lambda, V.hs.1, Tr.hs.1, d, r,
                  time.knots, rho, time.bound, return.se = TRUE)

## ---- warning = FALSE, fig.width = 4, fig.height = 3, message = FALSE---------
# Estimated parameters in STAR-PLVCM.
c(mfit4$theta.hat[1:mfit4$n.Z])
# Estimated standard deviation
mfit4$se.eta

# sequence of spatial plots of estimated coefficient functions.
ngrid.x <- 40
ngrid.y <- 20
ngrid.t <- 6

fitted.beta4 <- beta.func.plot(fitted = mfit4, ngrid.x, ngrid.y, ngrid.t, 
               boundary = horseshoe.b, boundary.t = time.bound,
               beta.true.ind = TRUE, beta.true = list(beta1.func), coord.ratio = 2)

fitted.beta4$beta.true.all[[1]]
fitted.beta4$beta.est.all[[1]]

