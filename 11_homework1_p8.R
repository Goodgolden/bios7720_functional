## remove all objects from the working directory except for
## code_path = location to source code
## results_path = locaiton to save results
## figure_path = location to save figures
rm(list = setdiff(ls(), c("code_path", "results_path", "figure_path")))

df <- read_rds(here("data", "data_processed", "NHANES_AC_processed.rds"))
min_cols <- paste0("MIN", 1:1440)
df$TAC <- rowSums(df[, min_cols], na.rm = TRUE)
df_ind <-
  df %>%
  dplyr::select(-one_of(min_cols)) %>%
  filter(good_day == 1, !is.na(Age)) %>%
  group_by(SEQN) %>%
  summarize(Age = Age[1], TAC_mn = mean(TAC))

fit_GCV <- gam(TAC_mn ~ s(Age, k = 50, bs = "cr"), method = "GCV.Cp", data = df_ind)
fit_REML <- gam(TAC_mn ~ s(Age, k = 50, bs = "cr"), method = "REML", data = df_ind)

## get estimates + standard errors
aind_pred <- seq(min(df_ind$Age), max(df_ind$Age), len = 100)
yhat_GCV <- predict(fit_GCV, newdata = data.frame(Age = aind_pred), type = "response", se.fit = TRUE)
yhat_REML <- predict(fit_REML, newdata = data.frame(Age = aind_pred), type = "response", se.fit = TRUE)

preds_GCV <- cbind(yhat_GCV$fit, yhat_GCV$fit - 2 * yhat_GCV$se.fit, yhat_GCV$fit + 2 * yhat_GCV$se.fit)
preds_REML <- cbind(yhat_REML$fit, yhat_REML$fit - 2 * yhat_REML$se.fit, yhat_REML$fit + 2 * yhat_REML$se.fit)

## plot the figures for part (A)
jpeg(here(figure_path, "fig_HW1_P8a.jpeg"), height = 500, width = 1000, quality = 100)
par(mfrow = c(1, 2))
ylims <- range(c(preds_GCV, preds_REML))
matplot(aind_pred, preds_GCV, type = "l", lty = c(1, 2, 2), ylab = "TAC", main = "GCV", ylim = ylims, col = "black")
matplot(aind_pred, preds_REML, type = "l", lty = c(1, 2, 2), ylab = "TAC", main = "REML", ylim = ylims, col = "black")
dev.off()


## answer part (B)

## number of subsamples to draw
nsamp <- 1000
## range of sizes of each subsample
Ns <- c(100, 500, 1000, 2000)
## total number in the sample
Ntot <- nrow(df_ind)
## empty array for storing results
arr_wiggle <- array(NA,
  dim = c(nsamp, 2, length(Ns)),
  dimnames = list("sim" = 1:nsamp, "method" = c("GCV", "REML"), "N" = Ns)
)
## set up progress bar to keep track of where we're at in the loop
pb <- txtProgressBar(min = 0, max = nsamp * length(Ns), initial = 0, style = 3)
inx <- 0
for (n in 1:nsamp) {
  for (N in seq_along(Ns)) {
    ## do subsampling without replace (with replacement is fine)
    inx_nN <- sample(1:Ntot, size = Ns[N], replace = FALSE)
    df_nN <- df_ind[inx_nN, ]
    ## fit GCV and REML on the sub-sample
    fit_GCV_nN <- gam(TAC_mn ~ s(Age, k = 50, bs = "cr"), method = "GCV.Cp", data = df_nN)
    fit_REML_nN <- gam(TAC_mn ~ s(Age, k = 50, bs = "cr"), method = "REML", data = df_nN)
    ## calculate xi^t S xi for REML and GCV
    xi_GCV_nN <- coef(fit_GCV_nN)[-1]
    xi_REML_nN <- coef(fit_REML_nN)[-1]
    ## store results
    arr_wiggle[n, 1, N] <- t(xi_GCV_nN) %*% fit_GCV_nN$smooth[[1]]$S[[1]] %*% xi_GCV_nN
    arr_wiggle[n, 2, N] <- t(xi_REML_nN) %*% fit_REML_nN$smooth[[1]]$S[[1]] %*% xi_REML_nN
    ## update progress bar
    setTxtProgressBar(pb, inx)
    inx <- inx + 1
  }
}
## transform array to a data frame
df_wiggle <- as.data.frame.table(arr_wiggle, stringsAsFactors = FALSE, responseName = "wiggle")
plt_wiggle <-
  df_wiggle %>%
  mutate(
    log10wiggle = log10(1 + wiggle),
    N_fac = factor(N, levels = Ns, labels = paste0("N = ", Ns))
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = log10wiggle)) +
  facet_grid(~N_fac, scales = "free_y") +
  xlab("") +
  ylab(expression(log[10](1 + xi^t ~ S ~ xi))) +
  theme_classic(base_size = 18)

jpeg(here(figure_path, "fig_HW1_P8b.jpeg"), height = 400, width = 1200, quality = 100)
print(plt_wiggle)
dev.off()
