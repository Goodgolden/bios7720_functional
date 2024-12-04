library("tidyverse")
library("mgcv")

## Simulate data from the model
# y_i(s) = f_0(s) + f_1(s)x_i + b_i(s) + \epsilon_i(s)
# where f_0(s) = 0
# x_i ~ N(0,1)
# b_i ~ GP(0, \Sigma)
# \epsilon_i ~ N(0,2^2)


set.seed(19840)
# simulation settings
N <- 200 # number of functions to simulate
ns <- 100 # number of observations per function
sind <- seq(0, 1, len = ns) # functional domain of observed functions
K <- 4 # number of true eigenfunctions
lambda <- 0.5^(0:(K - 1)) # true egenfunctions
sig2 <- 2 # error variance
# set up true eigenfunctions
Phi <- sqrt(2) * cbind(
  sin(2 * pi * sind), cos(2 * pi * sind),
  sin(4 * pi * sind), cos(4 * pi * sind))

# simulate coefficients
# first, simulate standard normals, then multiply by the
# standard deviation to get correct variance
xi_raw <- matrix(rnorm(N * K), N, K)
xi <- xi_raw %*% diag(sqrt(lambda))
# simulate b_i(s) as \sum_k \xi_ik \phi_k(t)
bi <- xi %*% t(Phi)

## define f(s)
f <- function(s) sin(2 * pi * s)
## get f(s) for all i, s
## fS is an N x ns matrix with rows repeated
fS <- kronecker(matrix(f(sind), 1, ns), matrix(1, N, 1))
x <- rnorm(N)
## get f(s)*x for each individual
fX <- fS * kronecker(matrix(x, N, 1), matrix(1, 1, ns))
## simulate the outcome
y <- bi + fX + matrix(rnorm(N * ns, sd = 2), N, ns)
## combine data into a data frame for model fitting
df_fit <-
  data.frame(
    id = factor(rep(1:N, each = ns)),
    y = as.vector(t(y)),
    bi = as.vector(t(bi)),
    x = rep(x, each = ns),
    id = rep(1:N, each = ns),
    sind = rep(sind, N),
    phi1 = rep(Phi[, 1], N),
    phi2 = rep(Phi[, 2], N),
    phi3 = rep(Phi[, 3], N),
    phi4 = rep(Phi[, 4], N))

head(df_fit)

## do the naive fit (ignore within subject correlation)
fit_naive <- gam(y ~ s(sind, k = 20, bs = "cr") +
                     s(sind, by = x, k = 20),
                 data = df_fit,
                 method = "REML")

## do the "oracle" fit where we know the true b_i(s) for each person
fit_oracle <- gam(y ~ s(sind, k = 20, bs = "cr") +
                      s(sind, by = x, k = 20),
                  data = df_fit,
                  method = "REML",
                  offset = df_fit$bi)

## do second oracle fit
df_fit$ystar <- df_fit$y - df_fit$bi
fit_oracle2 <- gam(ystar ~ s(sind, k = 20, bs = "cr") +
                       s(sind, by = x, k = 20),
                   data = df_fit,
                   method = "REML")

## plot the results
par(mfrow = c(1, 3))
plot(fit_naive,
     xlab = "s",
     ylab = expression(hat(y)),
     main = "Naive Fit",
     select = 2)
lines(sind,
      f(sind),
      col = "red",
      lty = 2)
plot(fit_oracle,
     xlab = "s",
     ylab = expression(hat(y)),
     main = "Oracle Fit",
     select = 2)
lines(sind,
      f(sind),
      col = "red",
      lty = 2)
plot(fit_oracle2,
     xlab = "s",
     ylab = expression(hat(y)),
     main = "Oracle Fit (v2)",
     select = 2)
lines(sind,
      f(sind),
      col = "red",
      lty = 2)

## note that I'm using only the first eigenfunction
## for computational efficiency
## try adding in the rest of the eigenfunctions!
fit_oracle_phi <- bam(y ~ s(sind, bs = "cr", k = 20) +
                          s(sind, bs = "cr", k = 20, by = x) +
                          s(phi1, by = id, bs = "re"),
                      data = df_fit,
                      method = "fREML",
                      discrete = TRUE)

fit_oracle_phi <- bam(y ~ s(sind, bs = "cr", k = 20) +
                          s(sind, bs = "cr", k = 20, by = x) +
                          s(phi1, by = id, bs = "re") +
                          s(phi2, by = id, bs = "re") +
                          s(phi3, by = id, bs = "re"),
                      data = df_fit,
                      method = "fREML",
                      discrete = TRUE)
