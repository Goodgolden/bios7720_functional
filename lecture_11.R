rm(list = ls())
library("dplyr")
library("mgcv")
library("refund")

## change data_path to whereever your NHANES file is located
data_path <- here::here("NHANES_AC_processed.rds")
df <- readr::read_rds(data_path)
## extract the PA data
lX <- log(1 + as.matrix(df[, paste0("MIN", 1:1440)]))
lX[is.na(lX)] <- 0
N <- nrow(lX)

## bin the data into 60 minute intervals
tlen <- 60
nt <- ceiling(1440 / tlen)
inx_cols <- split(1:1440, rep(1:nt, each = tlen)[1:1440])
lX_bin <- vapply(inx_cols,
                 function(x) rowMeans(lX[, x], na.rm = TRUE),
                 numeric(N))

## get subject average curves
inx_rows <- split(1:N,
                  factor(df$SEQN,
                         levels = unique(df$SEQN)))
lX_bin_ind <- t(vapply(inx_rows,
                       function(x) colMeans(lX_bin[x, ],
                                            na.rm = TRUE),
                       numeric(nt)))
nid <- nrow(lX_bin_ind)

# get a data frame for model fitting
sind <- seq(0, 1, len = nt)
udf <- df %>%
  dplyr::select(SEQN, Age) %>%
  group_by(SEQN) %>%
  slice(1) %>%
  ungroup()

df_fit <- data.frame(lAC = as.vector(t(lX_bin_ind)),
                     sind = rep(sind, nid),
                     SEQN = rep(unique(df$SEQN), each = nt)) %>%
  left_join(udf, by = "SEQN") %>%
  mutate(id = factor(SEQN)) %>%
  filter(!is.na(Age))


## subset the data to just 500 participants
# set.seed(10110)
set.seed(1011011)
nid_samp <- 500
id_samp <- sample(unique(df_fit$id),
                  size = nid_samp,
                  replace = FALSE)
df_fit_sub <- subset(df_fit, id %in% id_samp)

## fit the naive model
fit_naive <- bam(lAC ~ s(sind, bs = "cc", k = 20) +
                   s(sind, by = Age, bs = "cc", k = 20),
                 method = "fREML",
                 data = df_fit_sub,
                 discrete = TRUE)

## extract the resiudals
resid_mat <- matrix(fit_naive$residuals,
                    nid_samp, nt,
                    byrow = TRUE)
## fit fpca
fpca_fit <- fpca.face(resid_mat, knots = 15)
## add in eigen-function
for (k in 1:length(fpca_fit$evalues)) {
  df_fit_sub[[paste0("Phi", k)]] <- rep(fpca_fit$efunctions[, k], nid_samp)
}



system.time({
  fit_fri <- bam(lAC ~  s(sind, bs = "cc", k = 20) +
                   s(sind, by = Age, bs = "cc", k = 20) +
                   s(id, by = Phi1, bs = "re") +
                   s(id, by = Phi2, bs = "re") +
                   s(id, by = Phi3, bs = "re") +
                   s(id, by = Phi4, bs = "re"),
                 method = "fREML",
                 data = df_fit_sub,
                 discrete = TRUE)
  })


## get indices for the leave-function(s)-out CV
# get unqiue participant IDs
uid <- unique(df_fit_sub$SEQN)
nid <- length(uid)
# number of folds to split the participants
nfolds <- 5
# set the seed
set.seed(1010)
# create a vector of fold indicators
fold_vec <- sample(rep(1:nfolds,
                       ceiling(nid / nfolds))[1:nid])
# combine with IDs
df_folds <- data.frame(SEQN = uid,
                       fold = sample(1:nfolds))
# merge with the data frame for fitting
df_fit_sub <-
  left_join(df_fit_sub,
            df_folds,
            by = "SEQN")
# create a list with each element as the vector of row indices
# associated with each fold
inx_ls <- split(1:nrow(df_fit_sub),
                df_fit_sub$fold)


## set up grid search (here do a coarse search for computation time)
nlambda <- 25
loglambda <- seq(5, 20, len = nlambda)
lambda <- exp(loglambda)
## set up progress bar
pb <- txtProgressBar(min = 0,
                     max = nlambda * nfolds,
                     style = 3)
inx_pb <- 1
## loop over candidate smoothing parameters
mse <- rep(NA, nlambda)
## get smoothing parameters from the FRI model
sp_fri <- fit_fri$sp[1:2]
for (l in 1:nlambda) {
  ## loop over "folds"
  mse_l <- rep(NA, nfolds)
  for (i in 1:nfolds) {
    # get the training and test data
    inx_li <- inx_ls[[i]]
    df_train <- df_fit_sub[-inx_li, ]
    df_test <- df_fit_sub[inx_li, ]
    # fit the model using the training data
    fit_li <- bam(lAC ~ s(sind, bs = "cc", k = 20, sp = sp_fri[1]) +
                    s(sind, by = Age, bs = "cc", k = 20, sp = lambda[l]),
                  data = df_train)
    # predict on the test data
    fhat_test <- predict(fit_li,
                         newdata = df_test,
                         type = "response")
    mse_l[i] <- mean((fhat_test - df_test$lAC)^2)
    ## update progress bar
    setTxtProgressBar(pb, inx_pb)
    inx_pb <- inx_pb + 1
  }
  # average over the folds for this candidate smoothing parameter
  mse[l] <- mean(mse_l)
}


## heres a faster way to do it

## this is a function I wrote which does the calcualtions very quickly
my_LOS_cv <- function(lambda, Y, Phi, S, inx_ls) {
  lambda_mat <- lapply(1:length(lambda), function(l) lambda[l] * S[[l]])
  lambda_mat <- Matrix::bdiag(1, Matrix::bdiag(lambda_mat))
  N <- length(inx_ls)
  MSE_vec <- rep(NA, N)
  for (i in 1:N) {
    inx <- inx_ls[[i]]
    Y_i <- Y[-inx]
    Phi_i <- Phi[-inx, ]
    Phi_sq_i <- crossprod(Phi_i)
    Phi_sq_S_inv <- solve(as.matrix(Phi_sq_i + lambda_mat), tol = 1e-30)
    Yhat_i <- Phi[inx, ] %*% (Phi_sq_S_inv %*% (t(Phi_i) %*% Y_i))
    MSE_vec[i] <- mean((Y[inx] - Yhat_i)^2)
  }
  mean(MSE_vec)
}

mse_my_CV <- rep(NA, nlambda)
Phi <- predict(fit_naive, type = "lpmatrix")
pb2 <- txtProgressBar(min = 0, max = nlambda, style = 3)
S_ls <- lapply(fit_naive$smooth, function(x) x$S[[1]])
inx_pb2 <- 1
Y <- df_fit_sub$lAC
for (l in 1:nlambda) {
  mse_my_CV[l] <- my_LOS_cv(lambda = c(sp_fri[1], lambda[l]),
                            Y = Y,
                            Phi = Phi,
                            S = S_ls,
                            inx_ls = inx_ls)
  setTxtProgressBar(pb2, value = inx_pb2)
  inx_pb2 <- inx_pb2 + 1
}








