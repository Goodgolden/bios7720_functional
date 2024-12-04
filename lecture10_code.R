## load some libraries
library("microbenchmark")
library("tidyverse")
library("mgcv")
library("refund")
library("gamm4")
## set data path where your NHANES data are stored --
## youll need to change this in
## order for the NHANES data to load
data_path <- here::here("NHANES_AC_processed.rds")

## Homework1 ---------------------------------------------

## slow gcv function
my_gcv1 <- function(x, y,
                    bs = "cr",
                    k = 50,
                    loglambda = seq(-3, 20, len = 100)) {
  sm <- smoothCon(s(x, bs = bs, k = k),
                  data = data.frame(x = x))
  S <- sm[[1]]$S[[1]]
  Phi <- sm[[1]]$X
  nlambda <- length(loglambda)
  gcv_vec <- rep(NA, nlambda)
  for (l in seq_along(1:nlambda)) {
    H <- Phi %*%
        solve(crossprod(Phi) + exp(loglambda[l]) * S) %*%
        t(Phi)
    trH <- sum(diag(H))
    xi_hat <- solve(crossprod(Phi) +
                        exp(loglambda[l]) * S) %*%
        t(Phi) %*% y
    y_hat <- Phi %*% xi_hat
    gcv_vec[l] <- length(y) * sum((y_hat - y)^2 /
                                      (length(y) - trH)^2)
  }
  exp(loglambda[which.min(gcv_vec)])
}


## fast gcv function
my_gcv2 <- function(x, y,
                    bs = "cr",
                    k = 50,
                    loglambda = seq(-3, 20, len = 100)) {
  sm <- smoothCon(s(x, bs = bs, k = k),
                  data = data.frame(x = x))
  S <- sm[[1]]$S[[1]]
  Phi <- sm[[1]]$X
  Phi_sq <- crossprod(Phi)
  nlambda <- length(loglambda)
  gcv_vec <- rep(NA, nlambda)
  N <- length(y)
  lambda <- exp(loglambda)
  for (l in seq_along(1:nlambda)) {
    Phisq_S_inv <- chol2inv(chol(Phi_sq + lambda[l] * S))
    trH <- sum(diag(Phi_sq %*% Phisq_S_inv))
    xi_hat <- Phisq_S_inv %*% (t(Phi) %*% y)
    y_hat <- Phi %*% xi_hat
    gcv_vec[l] <- N * sum((y_hat - y)^2 / (N - trH)^2)
  }
  lambda[which.min(gcv_vec)]
}

## comparing the two gcv functions

N <- 1000
x <- rnorm(N)
f <- function(x) sin(2 * pi * x)
y <- f(x) + rnorm(N)

my_gcv1(x = x, y = y)
my_gcv2(x = x, y = y)
microbenchmark(my_gcv2(x = x, y = y), my_gcv1(x = x, y = y), times = 5)



## code for appending data
append_ls1 <- function(K = 1000, mat) {
  ret <- vector(mode = "list", length = K)
  names(ret) <- 1:K
  for (k in 1:K) ret[[k]] <- mat
  dplyr::bind_rows(ret)
}
append_c1 <- function(K = 1000, mat) {
  ret <- c()
  for (k in 1:K) {
    ret <- rbind(ret, mat)
  }
  ret
}
append_c2 <- function(K = 1000, mat) {
  dims <- dim(mat)
  ret <- matrix(NA, dims[1] * K, dims[2])
  inx <- 1:dims[1]
  for (k in 1:K) {
    ret[inx, ] <- mat
    inx <- inx + N
  }
  ret
}

N <- 100
P <- 100
mat <- matrix(rnorm(N * P), N, P)
microbenchmark(append_ls1(K = 100, mat = mat),
  append_c1(K = 100, mat = mat),
  append_c2(K = 100, mat = mat),
  times = 5
)




################
##            ##
## FoSR       ##
##            ##
################


## simulating data
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
  sin(4 * pi * sind), cos(4 * pi * sind)
)
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

df_fit <-
  data.frame(
    id = factor(rep(1:N, each = ns)), # important that this is a factor variable!
    y = as.vector(t(y)), # stack the vectors of Y^t (stacks individual functions)
    bi = as.vector(t(bi)), # same thing with the random intercept
    x = rep(x, each = ns), # repeat the fixed covariate for each function
    sind = rep(sind, N), # include the functional domain, repeat for each function
    phi1 = rep(Phi[, 1], N), # incliude the true eigenfunctions, repeat for each function
    phi2 = rep(Phi[, 2], N),
    phi3 = rep(Phi[, 3], N),
    phi4 = rep(Phi[, 4], N)
  )
View(df_fit)



#### Fit FoSR using the iterative approach
## fit the indepedence model
fit_naive <- gam(y ~ s(sind, k = 20, bs = "cr") + s(sind, by = x, bs = "cr", k = 20),
  data = df_fit, method = "REML"
)
## extract residuals
resid_mat <- matrix(fit_naive$residuals, N, ns, byrow = TRUE)
## fit fpca
fpca_fit <- fpca.face(resid_mat, var = TRUE)
## get the estimated eigenfunctions
efuncs <- fpca_fit$efunctions
## add estimated eigenfunctions to the dataframe
for (k in 1:dim(efuncs)[2]) {
  df_fit[[paste0("Phi", k, "_hat")]] <- efuncs[, k]
}
## fit the random functional intercept model --
## just use the first 4 estimated eigenfucntions
fit_rfi <- bam(y ~ s(sind, k = 20, bs = "cr") + s(sind, by = x, bs = "cr", k = 20) +
  s(id, by = Phi1_hat, bs = "re") + s(id, by = Phi2_hat, bs = "re") +
  s(id, by = Phi3_hat, bs = "re") + s(id, by = Phi4_hat, bs = "re"),
data = df_fit, method = "fREML", discrete = TRUE
)


## get the estiamted coefficients + SEs
df_pred <- data.frame(
  "sind" = sind, "x" = 1,
  ## we need to supply all columns of the
  ## terms used in fitting the functional random
  ## intercept model. Since we're only interested in the
  ## fixed effects here, set to whatever you want
  ## note that id must be included in the dataset used for
  ## model fitting
  Phi2_hat = 0, Phi1_hat = 0, Phi3_hat = 0, Phi4_hat = 0,
  id = df_fit$id[1]
)
fhat_naive <- predict(fit_naive, newdata = df_pred, se.fit = TRUE, type = "terms")
fhat_rfi <- predict(fit_rfi, newdata = df_pred, se.fit = TRUE, type = "terms")
head(fhat_rfi$fit)


# FoSR using mixed model software
library("gamm4")
system.time({
  fit_rfit_gamm4 <- gamm4(y ~
  s(sind, bs = "cr", k = 20) +
    s(sind, by = x, bs = "cr", k = 20),
  random = ~ (Phi1_hat + 0 | id) + (Phi2_hat + 0 | id) +
    (Phi3_hat + 0 | id) + (Phi4_hat + 0 | id),
  data = df_fit, REML = TRUE
  )
})
system.time({
  fit_rfit_gamm <- gamm(y ~ s(sind, bs = "cr", k = 20) + s(sind, by = x, bs = "cr", k = 20),
    random = list(id = pdDiag(~ Phi1_hat + Phi2_hat +
      Phi3_hat + Phi4_hat + 0)),
    data = df_fit, REML = TRUE
  )
})



## NHANES application

## library("tidyverse")
df <- read_rds(here::here("data", "data_processed", "NHANES_AC_processed.rds"))
## extract the PA data
lX <- log(1 + as.matrix(df[, paste0("MIN", 1:1440)]))
lX[is.na(lX)] <- 0
N <- nrow(lX)
## bin the data into 20 minute intervals
tlen <- 20
nt <- ceiling(1440 / tlen)
inx_cols <- split(1:1440, rep(1:nt, each = tlen)[1:1440])
lX_bin <- vapply(inx_cols, function(x) rowMeans(lX[, x], na.rm = TRUE), numeric(N))
## get subject average curves
inx_rows <- split(1:N, factor(df$SEQN, levels = unique(df$SEQN)))
lX_bin_ind <- t(vapply(inx_rows, function(x) colMeans(lX_bin[x, ], na.rm = TRUE), numeric(nt)))
nid <- nrow(lX_bin_ind)
# get a data frame for model fitting
sind <- seq(0, 1, len = nt)
df_fit <-
  data.frame(
    lAC = as.vector(t(lX_bin_ind)),
    sind = rep(sind, nid),
    SEQN = rep(unique(df$SEQN), each = nt)
  ) %>%
  left_join(dplyr::select(df, SEQN, Age), by = "SEQN") %>%
  mutate(id = factor(SEQN)) %>%
  filter(!is.na(Age))

set.seed(10110)
nid_samp <- 500
id_samp <- sample(unique(df_fit$id), size = nid_samp, replace = FALSE)
df_fit_sub <- subset(df_fit, id %in% id_samp)

## fit the naive model
fit_naive <- bam(lAC ~ s(sind, bs = "cc", k = 20) + s(sind, by = Age, bs = "cc", k = 20),
  method = "fREML", data = df_fit_sub, discrete = TRUE
)
## extract the resiudals
resid_mat <- matrix(fit_naive$residuals,
  nid_samp, nt,
  byrow = TRUE
)
## fit fpca
fpca_fit <- fpca.face(resid_mat, knots = 15)
## add in eigenfunctiosn
for (k in 1:length(fpca_fit$evalues)) {
  df_fit_sub[[paste0("Phi", k)]] <- rep(fpca_fit$efunctions[, k], nid_samp)
}

## fit the fri models
system.time({
  fit_fri <- bam(lAC ~ s(sind, bs = "cc", k = 20) + s(sind, by = Age, bs = "cc", k = 20) +
    s(id, by = Phi1, bs = "re") + s(id, by = Phi2, bs = "re") +
    s(id, by = Phi3, bs = "re") + s(id, by = Phi4, bs = "re"),
  method = "fREML", data = df_fit_sub, discrete = TRUE
  )
})
