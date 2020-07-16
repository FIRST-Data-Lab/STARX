rm(list = ls())

# library
library(devtools)
install_github("funstatpackages/BPST")
install_github("funstatpackages/Triangulation")
install_github("funstatpackages/TPST")
install_github("funstatpackages/STARX")
library(Triangulation)
library(STARX)

nS <- 200
nT <- 50
sigma <- 1
N <- 4 # number of interior knots
d <- 2
r <- 1
rho <- 3
time.bound <- c(0, 1)

# load triangulation
data("Tr.hs.1")
data("V.hs.1")
data("horseshoe.b")
plot(horseshoe.b$V1, horseshoe.b$V2, type = "l", xlab = '', ylab = '')
TriPlot(V.hs.1, Tr.hs.1)

# generate simulation data
data.simu <- simu1.data.generating(nS, nT, sigma)

# knots
probs <- seq(time.bound[1], time.bound[2], length.out = N + 2)
time.knots <- quantile(data.simu$location[, 3], probs = probs)[-c(1, N + 2)]

# fit model -------------------------------------------------------------
Lambda1 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda2 <- exp(seq(log(0.001), log(1000), length.out = 5))
Lambda <- expand.grid(Lambda1, Lambda2)

mfit1 <- stvcm.fit(data = data.simu, Lambda, V.hs.1, Tr.hs.1, d, r,
                   time.knots, rho, time.bound)
mfit1$n.X

ngrid.x <- 40
ngrid.y <- 20
ngrid.t <- 6

fitted.beta <- beta.func.plot(fitted = mfit1, ngrid.x, ngrid.y, ngrid.t,
                               boundary = horseshoe.b, boundary.t = time.bound,
                               beta.true.ind = TRUE,
                               beta.true = list(beta1.func, beta2.func, beta3.func),
                               coord.ratio = 2)

fitted.beta$beta.true.all[[1]]
fitted.beta$beta.true.all[[2]]
fitted.beta$beta.true.all[[3]]
fitted.beta$beta.est.all[[1]]
fitted.beta$beta.est.all[[2]]
fitted.beta$beta.est.all[[3]]

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

result <- eval.star(fitted = mfit1,
                    beta.true = list(beta1.func, beta2.func, beta3.func),
                    data.grid = data.grid)
result$mise.beta
