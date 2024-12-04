## remove all objects from the working directory except for
## code_path = location to source code
## results_path = locaiton to save results
## figure_path = location to save figures
rm(list = setdiff(ls(), c("code_path", "results_path", "figure_path")))

df <- read_rds(here("data", "data_processed", "NHANES_AC_processed.rds"))
min_cols <- paste0("MIN", 1:1440)
df$TAC <- rowSums(df[, min_cols], na.rm = TRUE)
df_fit <-
  df %>%
  filter(good_day == 1, !is.na(Age)) %>%
  group_by(SEQN) %>%
  mutate(TAC_mn = mean(TAC))
df_full <- df_fit
df_ind <- df_fit[!duplicated(df_fit$SEQN), ]


## 7a --------------------------------------------------------------------------
# fit the additive model the the repeated obs data
fit_full <- gam(TAC_mn ~ s(Age, bs = "cr", k = 50), data = df_full)
# fit the additive model the the single obs data
fit_ind <- gam(TAC_mn ~ s(Age, bs = "cr", k = 50), data = df_ind)
# get predictions and standard errors over a range of new age values
aind_pred <- seq(min(df_full$Age), max(df_full$Age), len = 1000)
df_pred <- data.frame(Age = aind_pred)
ests_full <- predict(fit_full, newdata = df_pred, type = "response", se.fit = TRUE)
ests_ind <- predict(fit_ind, newdata = df_pred, type = "response", se.fit = TRUE)

preds_full <- cbind(ests_full$fit, ests_full$fit + 2 * ests_full$se.fit, ests_full$fit - 2 * ests_full$se.fit)
preds_ind <- cbind(ests_ind$fit, ests_ind$fit + 2 * ests_ind$se.fit, ests_ind$fit - 2 * ests_ind$se.fit)


## do the plotting

jpeg(here(figure_path, "fig_HW1_P7a.jpeg"), height = 500, width = 500, quality = 100)
matplot(aind_pred, preds_full,
  type = "l", lty = c(1, 2, 2), xlab = "Age", ylab = expression(hat(f)(Age)), main = "Model Results", col = "black"
)
matplot(aind_pred, preds_ind,
  type = "l", lty = c(1, 2, 2), xlab = "", ylab = "", col = "red", add = TRUE
)
legend("topright", c("Repeated Obs (df_full)", "Single Obs (df_ind)"),
  col = c("black", "red"), lwd = 2, lty = 1, bty = "n"
)
dev.off()



## 7b --------------------------------------------------------------------------

## split the data into 10 folds based on participant IDs
# set the seed
set.seed(189)
# get the unique participant IDs (SEQN) and the total # of IDs
uid <- unique(df_full$SEQN)
nid <- length(uid)
# split the data into 10 training and test splits
nfolds <- 10
uid_ls <- split(uid, sample(rep(1:nfolds, ceiling(nid / nfolds))[1:nid]))
# pre-calculate indices for each fold
# not technically necessary, but offers a slight speed up of the loop below
inx_i_ls <- lapply(uid_ls, function(x) which(df_full$SEQN %in% x))

## set up smoothing parameters to loop over
nlambda <- 100
loglambda <- seq(-3, 20, len = nlambda)

## empty vector for storing MSE
MSE <- rep(NA, nlambda)

## set up progress bar to keep track of where we're at in the loop
pb <- txtProgressBar(min = 0, max = nlambda * nfolds, initial = 0, style = 3)
inx <- 0
## loop over smoothing parameters
for (l in seq_along(loglambda)) {

  # set up empty container to store MSE for each fold
  # for the current candidate smoothing parameter
  MSE_l <- rep(NA, nfolds)
  ## loop over the 10 folds
  for (i in 1:nfolds) {
    ## fit the model leaving the "ith" fold out (1/10 of participant IDs)
    ## using the current candidate smoothing parameter (sp==exp(loglambda[l]))
    df_train_li <- df_full[-inx_i_ls[[i]], ]
    fit_li <- gam(TAC_mn ~ s(Age, bs = "cr", k = 50, sp = exp(loglambda[l])), data = df_train_li)

    # get MSE for the participant-days left out
    # note that here I'm calculating one error term per participant-day
    # this will implciitly upweight participants with more "good" days of data
    df_test_li <- df_full[inx_i_ls[[i]], ]
    MSE_l[i] <- mean((predict(fit_li, newdata = df_test_li, type = "response") - df_test_li$TAC_mn)^2)

    ## update progress bar
    setTxtProgressBar(pb, inx)
    inx <- inx + 1
  }
  ## get the average over the 10 folds
  MSE[l] <- mean(MSE_l)
}

## find "best" \lambda, fit the corresponding model
lambda <- exp(loglambda[which.min(MSE)])
fit_7b <- gam(TAC_mn ~ s(Age, bs = "cr", k = 50, sp = lambda), data = df_full)
ests_7b <- predict(fit_7b, newdata = df_pred, type = "response", se.fit = TRUE)
preds_7b <- cbind(ests_7b$fit, ests_7b$fit + 2 * ests_7b$se.fit, ests_7b$fit - 2 * ests_7b$se.fit)


## compare measures of wiggliness
d2f2_full <- t(coef(fit_full)[-1]) %*% fit_full$smooth[[1]]$S[[1]] %*% coef(fit_full)[-1]
d2f2_ind <- t(coef(fit_ind)[-1]) %*% fit_ind$smooth[[1]]$S[[1]] %*% coef(fit_ind)[-1]
d2f2_7b <- t(coef(fit_7b)[-1]) %*% fit_7b$smooth[[1]]$S[[1]] %*% coef(fit_7b)[-1]


jpeg(here(figure_path, "fig_HW1_P7b.jpeg"), height = 400, width = 1200, quality = 100)
par(mfrow = c(1, 3), mar = c(7, 6, 5, 5))
matplot(aind_pred, preds_7b,
  type = "l", lty = c(1, 2, 2), xlab = "Age",
  ylab = expression(hat(f)(Age)), main = "(A) 10-fold LSS CV Results", col = "black"
)
plot(c(d2f2_full, d2f2_ind, d2f2_7b),
  xaxt = "n", xlab = "Fit", ylab = expression(hat(xi) * S * hat(xi)^t),
  main = "(B) \"Wiggliness\" of Estimated f() ", pch = 16, cex = 2
)
axis(1, at = 1:3, labels = c("df_full", "df_ind", "LSS Out CV"))
plot(aind_pred, preds_7b[, 1] - preds_ind[, 1],
  xlab = "Age", ylab = "Difference", type = "l",
  main = "(C) LSS Estimate - df_ind Estimate"
)
dev.off()
