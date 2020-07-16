rm(list = ls())

# library
library(devtools)
install_github("funstatpackages/BPST")
install_github("funstatpackages/Triangulation")
install_github("funstatpackages/TPST")
install_github("funstatpackages/STARX")
library(Triangulation)
library(STARX)

nS <- 100
nT <- 30
sigma <- 0.5
alpha <- 0.5
N <- 4 # Number of interior knots.
d <- 2 # Degree of bivariate spline.
r <- 1 # smoothness of bivariate spline.
rho <- 3 # order of univariate spline.
time.bound <- c(0, 1)
# boundary of spatial domain
data("horseshoe.b")
plot(horseshoe.b$V1, horseshoe.b$V2, type = "l", xlab = '', ylab = '')


# load triangulation
data("Tr.hs.1")
data("V.hs.1")
TriPlot(V.hs.1, Tr.hs.1)

# generate simulation data
data.simu <- simu2.data.generating(nS, nT, sigma, alpha, k = 10, horseshoe.b)

# knots
probs <- seq(time.bound[1], time.bound[2], length.out = N + 2)
time.knots <- quantile(data.simu$location[, 3], probs = probs)[-c(1, N + 2)]

# fit model -------------------------------------------------------------
Lambda1 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda2 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda <- expand.grid(Lambda1, Lambda2)

mfit1 <- star.fit(data = data.simu, Lambda, V.hs.1, Tr.hs.1, d, r, time.knots, rho,
                 time.bound, return.se = TRUE)

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
                              beta.true.ind = TRUE, beta.true = list(beta1.func),
                              coord.ratio = 2)

fitted.beta$beta.true.all[[1]]
fitted.beta$beta.est.all[[1]]

# model evaluation ------------------------------------------------------
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

result$mse.alpha
result$mse.eta
result$mise.beta
