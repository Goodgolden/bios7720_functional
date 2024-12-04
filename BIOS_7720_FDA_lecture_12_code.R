library("tidyverse")
library("mgcv")
library("refund")

set.seed(9454785)
## read in the data
df <- here::here("NHANES_AC_processed.rds") %>%
  readr::read_rds() %>%
  ## subset to only "good" Saturdays,
  ## participants age 5-25
  filter(DoW %in% "Saturday",
         good_day %in% 1,
         Age <= 25) %>%
  ## subset to 200 randomly selected individuals
  sample_n(size = 200)
# View(df)


## extract activity counts (Y),
## create binary RV (Z)
## totally 1440 mins for the dataset \columns
Y <- as.matrix(df[, paste0("MIN", 1:1440)])
Y[is.na(Y)] <- 0
Z <- apply(Y >= 100, 2, as.numeric)
# View(Z)


## bin the data into 20 minute intervals
N <- nrow(Y)
tlen <- 20
nt <- ceiling(1440 / tlen)
inx_cols <- split(1:1440, rep(1:nt, each = tlen)[1:1440])
Z_bin <- vapply(inx_cols, function(x) rowSums(Z[, x, drop = FALSE]), numeric(N))
## add binned data back into our data frame
df[["Z_bin"]] <- Z_bin
# View(df)


## remove unnecesarry minute columns
df <- dplyr::select(df, Age, SEQN, BMI, Z_bin, Gender) %>%
  ## create factor variable for ID
  mutate(ID = factor(SEQN))
## create long format data frame
N <- nrow(df)
sind <- seq(0, 1, len = nt)
df_long <- data.frame(
    "Z" = as.vector(t(df$Z_bin)),
    "Age" = rep(df$Age, each = nt),
    "Gender" = rep(df$Gender, each = nt),
    "sind" = rep(sind, N),
    "ID" = rep(df$ID, each = nt))
df_long$Male <- as.numeric(df_long$Gender %in% "Male")
df_long$Z_n <- tlen - df_long$Z


## estimate the marginal and GFoSR models
fit_marginal <- bam(cbind(Z, Z_n) ~ s(sind, bs = "cc", k = 10) +
                      s(sind, by = Age, bs = "cc", k = 10) +
                      s(sind, by = Male, bs = "cc", k = 10),
                    family = "binomial",
                    data = df_long,
                    method = "fREML",
                    discrete = TRUE)
gfosr_time_st <- Sys.time()
fit_gfosr <- bam(cbind(Z, Z_n) ~
                   s(sind, bs = "cc", k = 10) +
                   s(sind, by = Age, bs = "cc", k = 10) +
                   s(sind, by = Male, bs = "cc", k = 10) +
                   ## For ti smooths you can specify which marginals
                   ## should have centering constraints applied,
                   ## by supplying 0/1 or FALSE/TRUE values
                   ## for each marginal in this vector.
                   ##
                   ## By default all marginals are constrained,
                   ## which is what is appropriate for, e.g.,
                   ## functional ANOVA models.
                   ##
                   ## Note that 'ti' only applies constraints to the marginals,
                   ## so if you turn off all marginal constraints
                   ## the term will have no identifiability constraints.
                   ## Only use this if you really understand
                   ## how marginal constraints work.
                   ti(ID, sind, bs = c("re", "cr"), mc = c(TRUE, FALSE), k = c(5, 5)),
                 family = "binomial",
                 data = df_long,
                 method = "fREML",
                 chunk.size = 10000,
                 discrete = TRUE)
gfosr_end_st <- Sys.time()
difftime(gfosr_end_st,
         gfosr_time_st,
         units = "mins")


## get estimated coefficients
df_pred <- data.frame(sind = sind,
                      Age = 1,
                      Male = 1,
                      ID = df_long$ID[1])
coefs_marginal <- predict(fit_marginal,
                          newdata = df_pred,
                          type = "terms",
                          se.fit = TRUE)
coefs_gfosr <- predict(fit_gfosr,
                       newdata = df_pred,
                       type = "terms",
                       se.fit = TRUE)


## plot the estimated coefficients
ests_marginal <- coefs_marginal$fit
ests_marginal[, 1] <- ests_marginal[, 1] + coef(fit_marginal)[1]
ests_gfosr <- coefs_gfosr$fit
ests_gfosr[, 1] <- ests_gfosr[, 1] + coef(fit_gfosr)[1]

LB_marginal <- ests_marginal - 1.96 * coefs_marginal$se.fit
UB_marginal <- ests_marginal + 1.96 * coefs_marginal$se.fit
LB_gfosr <- ests_gfosr - 1.96 * coefs_gfosr$se.fit
UB_gfosr <- ests_gfosr + 1.96 * coefs_gfosr$se.fit


par(mfrow = c(1, 3),
    oma = c(2, 2, 0, 0))
xinx <- (c(1, 6, 12, 18, 23) * 60 + 1) / 1440
xinx_lab <- c("01:00", "06:00", "12:00", "18:00", "23:00")
coef_labs <- list(expression("Intercept: " ~ f[0](s)),
                  expression("Age: " ~ f[1](s)),
                  expression("Male: " ~ f[2](s)))

for (p in 1:3) {
  ylims_p <- range(c(LB_marginal[, p],
                     UB_marginal[, p],
                     LB_gfosr[, p],
                     UB_gfosr[, p]))
  matplot(sind,
          cbind(ests_marginal[, p],
                LB_marginal[, p],
                UB_marginal[, p]),
          type = "l",
          lty = c(1, 2, 2),
          col = "black",
          ylim = ylims_p,
          xaxt = "n",
          ylab = "",
          xlab = "",
          main = coef_labs[p])

  axis(1, at = xinx,
       labels = xinx_lab)
  matplot(sind,
          cbind(ests_gfosr[, p],
                LB_gfosr[, p],
                UB_gfosr[, p]),
          type = "l",
          lty = c(1, 2, 2),
          col = "red",
          add = TRUE)
  abline(h = 0,
         col = "grey",
         lty = 2,
         lwd = 2)
  if (p == 2) {
    legend("top",
           c("Marginal", "GFoSR"),
           lty = c(1, 1),
           lwd = 2,
           bty = "n",
           col = c("black", "red"))
  }
}

mtext("Time of Day (s)", side = 1, outer = TRUE)
mtext(expression(hat(f)(s) ~ "+/-" ~ "2SE(" ~ hat(f)(s) ~ ")"),
      side = 2,
      outer = TRUE)


## obtain subject predictions
## using the returned values from mgcv::bam
ests_direct <- fit_gfosr$fitted.values
## obtain subject predictions
## using predict.bam() with type="terms"
ests_terms <- predict(fit_gfosr,
                      newdata = df_long,
                      type = "terms")
expit <- function(x) 1 / (1 + exp(-x))
ests_terms <- expit(coef(fit_gfosr)[1] + rowSums(ests_terms))
## obtain subject predictions
## using predict.bam() with type="lpmatrix"
ests_lp <- predict(fit_gfosr,
                   newdata = df_long,
                   type = "lpmatrix")
ests_lp <- expit(ests_lp %*% coef(fit_gfosr))


##  plot subject predictions
set.seed(22388)
nplt <- 20
uid <- unique(df_long$ID)

id_samp <- sample(uid,
                  size = nplt,
                  replace = FALSE)


op <- par(mfrow = c(4, 5),
          oma = c(2, 2, 0, 0),
          mar = c(2, 2, 2, 2))
xinx <- (c(1, 12, 23) * 60 + 1) / 1440
xinx_lab <- c("01:00", "12:00", "23:00")

for (i in 1:nplt) {
  inx_i <- which(df_long$ID == id_samp[i])

  plot(sind,
       df_long$Z[inx_i] / 20,
       xlab = "",
       ylab = "",
       ylim = c(0, 1),
       xaxt = "n",
       main = id_samp[i],
       pch = 16,
       cex = 0.75,
       col = rgb(0, 0, 0, 0.8),
       las = 1)
  axis(1, at = xinx,
       labels = xinx_lab)
  lines(sind,
        fit_marginal$fitted.values[inx_i],
        col = "black",
        lwd = 1.25)
  lines(sind,
        fit_gfosr$fitted.values[inx_i],
        col = "red",
        lwd = 1.25)
}
mtext("Time of Day (s)",
      side = 1,
      outer = TRUE)
mtext("Pr(Active)",
      side = 2,
      outer = TRUE)

par(op)
