

## source the code for problem 2
## containing the grid search and optimization
## functions for GCV
source(here::here(code_path, "11_homework1_p2.R"))

## set to FALSE if you don't want to save the results
save_results <- TRUE

## set number of simulated datasets to create for each scenario
nsim <- 1000
## set up simulation scenarios
Ns <- c(100, 500, 1000)
fs <- list(
  f1 = function(x) x^2,
  f2 = function(x) sin(x),
  f3 = function(x) cos(x)^2 + x^3)
Ks <- c(5, 20, 50)

## set up empty containers to store results
# range of x values to assess simulations on, -3,3 assumes
# X ~ N(0,1), may need to modify range depedning on how you simulate X
nx_pred <- 100
xind_pred <- seq(-1, 1, len = nx_pred)
arr_MSE_pen <- arr_bias_pen <-
    arr_coverage_pen <- arr_MSE_up <-
        arr_bias_up <- arr_coverage_up <-
    array(NA, dim = c(nsim, length(Ns), length(fs), length(Ks), nx_pred),
          dimnames = list("sim" = 1:nsim, "N" = Ns,
                          "f" = paste0("f", 1:length(fs)),
                          "K" = Ks, "xind" = 1:nx_pred))

## set up a progressbar
pb <- txtProgressBar(min = 0,
                     max = nsim * length(Ns) * length(fs) * length(Ks),
                     style = 3)
inx <- 0
## loop over number of observations
for (N in seq_along(Ns)) {
  N_cur <- Ns[N]
  ## loop over association structures
  for (f in seq_along(fs)) {
    f_cur <- fs[[f]]
    ftrue <- f_cur(xind_pred)
    ## loop over number of basis functions
    for (K in seq_along(Ks)) {
      K_cur <- Ks[K]
      ## loop over simulated datasets
      for (n in 1:nsim) {
        ## simulate covariate data
        x <- runif(N_cur, min = -1, max = 1)
        ## simulate outcome data
        y <- f_cur(x) + rnorm(N_cur, mean = 0, sd = 0.5)
        ## set up basis functions, get S matrix
        sm <- smoothCon(s(x, bs = "cr", k = K_cur),
                        data = data.frame(x = x))
        Phi <- sm[[1]]$X
        S <- sm[[1]]$S[[1]]
        S_sqrt <- chol(S, pivot = TRUE)
        S_sqrt <- S_sqrt[, order(attr(S_sqrt, "pivot"))]
        ## get basis matrix for new predictions
        Phi_pred <- PredictMat(sm[[1]], data = data.frame(x = xind_pred))


        ## fit the unpenalized model
        fit_up <- lm(y ~ Phi - 1)
        ## get estimated function +/- 2SE over the range of xind
        fhat_up <- Phi_pred %*% coef(fit_up)
        se_fhat_up <- sqrt(diag(Phi_pred %*% vcov(fit_up) %*% t(Phi_pred)))
        LB_up <- fhat_up - 2 * se_fhat_up
        UB_up <- fhat_up + 2 * se_fhat_up
        ## calculate bias, MSE, coverage probability
        ## store results
        arr_bias_up[n, N, f, K, ] <- fhat_up - ftrue
        arr_coverage_up[n, N, f, K, ] <- as.numeric(ftrue >= LB_up & ftrue <= UB_up)
        arr_MSE_up[n, N, f, K, ] <- (fhat_up - ftrue)^2


        ## penalized fit
        ## find optiomal smoothing parameter
        ## given the optimal smoothing parameter, get final fit
        fit_NfK <- fn_PENSSE_grid(y = y, Phi = Phi, S = S, S_sqrt = S_sqrt, loglambda = seq(-3, 20, len = 1000))


        ## get estimated function +/- 2SE over the range of xind
        fhat_pen <- Phi_pred %*% fit_NfK$xi
        se_fhat_pen <- sqrt(diag(Phi_pred %*% fit_NfK$var_xi %*% t(Phi_pred)))
        LB_pen <- fhat_pen - 2 * se_fhat_pen
        UB_pen <- fhat_pen + 2 * se_fhat_pen
        ## calculate bias, MSE, coverage probability
        ## store results
        arr_bias_pen[n, N, f, K, ] <- fhat_pen - ftrue
        arr_coverage_pen[n, N, f, K, ] <- as.numeric(ftrue >= LB_pen & ftrue <= UB_pen)
        arr_MSE_pen[n, N, f, K, ] <- (fhat_pen - ftrue)^2
        inx <- inx + 1
        setTxtProgressBar(pb, inx)
      }
    }
  }
}

if (save_results) {
  ## save the results
  HW1_P5_results <- list(
    arr_bias_pen = arr_bias_pen,
    arr_coverage_pen = arr_coverage_pen,
    arr_MSE_pen = arr_MSE_pen,
    arr_bias_up = arr_bias_up,
    arr_coverage_up = arr_coverage_up,
    arr_MSE_up = arr_MSE_up,
    xind_pred = xind_pred,
    fs = fs, Ns = Ns, Ks = Ks)
  write_rds(HW1_P5_results,
            file.path(results_path, "HW1_P5_results.rds"))
}



## summarize and plot
## convert arrays to data frames
## for plotting using as.data.frame.table
## unpenalized results
df_bias_up <- as.data.frame.table(arr_bias_up,
                                  stringsAsFactors = FALSE,
                                  responseName = "bias")
df_coverage_up <- as.data.frame.table(arr_coverage_up,
                                      stringsAsFactors = FALSE,
                                      responseName = "coverage")
df_MSE_up <- as.data.frame.table(arr_MSE_up,
                                 stringsAsFactors = FALSE,
                                 responseName = "MSE")
df_results_up <- df_bias_up %>%
  left_join(df_coverage_up,
            by = c("sim", "N", "f", "K", "xind")) %>%
  left_join(df_MSE_up,
            by = c("sim", "N", "f", "K", "xind")) %>%
  mutate(penalized = "unpenalized")

## penalized results
df_bias_pen <-
    as.data.frame.table(
        arr_bias_pen,
        stringsAsFactors = FALSE,
        responseName = "bias")
df_coverage_pen <-
    as.data.frame.table(
        arr_coverage_pen,
        stringsAsFactors = FALSE,
        responseName = "coverage")
df_MSE_pen <-
    as.data.frame.table(
        arr_MSE_pen,
        stringsAsFactors = FALSE,
        responseName = "MSE")
df_results_pen <- df_bias_pen %>%
    left_join(df_coverage_pen,
              by = c("sim", "N", "f", "K", "xind")) %>%
  left_join(df_MSE_pen,
            by = c("sim", "N", "f", "K", "xind")) %>%
  mutate(penalized = "penalized")
## combine into a single data frame
df_results <- bind_rows(df_results_up, df_results_pen)

## clean up the workspace
rm(df_bias_up,
   df_coverage_up,
   df_MSE_up,
   df_bias_pen,
   df_coverage_pen,
   df_MSE_pen,
   df_results_pen,
   df_results_up)

gc()

df_results_summary <-
  df_results %>%
  group_by(penalized, N, f, K, xind) %>%
  summarize(bias = mean(bias),
            coverage = mean(coverage),
            MSE = mean(MSE)) %>%
  mutate(
    f = factor(f, levels = names(fs),
               labels = paste0("f = ", names(fs))),
    N = factor(N, levels = Ns,
               labels = paste0("N = ", Ns)),
    K = factor(K, levels = Ks,
               labels = paste0("K = ", Ks)),
    xind = xind_pred[as.numeric(xind)]) %>%
  ungroup()


plt_bias <- df_results_summary %>%
  ggplot(aes(x = xind,
             y = bias,
             group = K,
             color = K)) +
  geom_line() +
  facet_grid(f ~ penalized + N) +
  theme_classic(base_size = 24) +
  scale_x_continuous() +
  ylab("Bias") +
  xlab("x") +
  ggtitle(expression("(A) Bias = " ~ hat(f) - f(x))) +
  ylim(-0.1, 0.1)


plt_MSE <- df_results_summary %>%
  ggplot(aes(x = xind,
             y = MSE,
             group = K,
             color = K)) +
  geom_line() +
  facet_grid(f ~ penalized + N) +
  theme_classic(base_size = 24) +
  scale_x_continuous() +
  ylab("MSE") +
  xlab("x") +
  ggtitle(expression("(B): MSE = " ~ (hat(f)(x) - f(x))^2)) +
  ylim(0, 0.1)

plt_coverage <-
  df_results_summary %>%
  ggplot(aes(x = xind,
             y = coverage,
             group = K,
             color = K)) +
  geom_line() +
  facet_grid(f ~ penalized + N) +
  theme_classic(base_size = 24) +
  scale_x_continuous() +
  ylab("95% CI Coverage") +
  geom_abline(slope = 0,
              intercept = 0.95,
              col = "grey",
              lty = 2) +
  xlab("x") +
  ggtitle("(C): 95% Coverage Probability")



jpeg(here(figure_path, "HW1_p5.jpeg"),
     height = 1500,
     width = 1000,
     quality = 100)
grid.arrange(plt_bias,
             plt_MSE,
             plt_coverage,
             ncol = 1)
dev.off()
