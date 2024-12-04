library("tidyr")
library("here")
library("readr")
library("dplyr")
library("refund")
library("mgcv")


data_mort <- here("data_mort.rds") %>%
    read_rds()

data <-
    here("NHANES_AC_processed.rds") %>%
    read_rds() %>%
    ## only consider good days of data and indiviudals age 50 or over
    filter(good_day %in% 1, Age > 50) %>%
    ## get mortality data from the rnhanesdata package
    ## merge with our data and derive 5-year mortality indicator
    left_join(data_mort, by = "SEQN") %>%
    mutate(mort_5yr = as.numeric(permth_exm / 12 <= 5 & mortstat %in% 1),
           ## replace accidental deaths within 5 years as NA
           mort_5yr = ifelse(mort_5yr == 1 & ucod_leading %in% "004", NA, mort_5yr)) %>%
    ## drop anyone missing mortality data or who had accidental deaths within 5 years
    filter(!is.na(mort_5yr))



## extract just the activity count data
Z <- log(as.matrix(data[, paste0("MIN", 1:1440)]) + 1)
Z[is.na(Z)] <- 0
## average across days within participants (SEQN)
## unique subject identifiers
uid <- unique(data$SEQN)
## number of participants
nid <- length(uid)
## empty container to store average profiles
Zmat <- matrix(NA, nid, 1440)

## list of indices
inx_ls <- lapply(uid, function(x) which(data$SEQN %in% x))
for (i in seq_along(uid)) {
  Zmat[i, ] <- colMeans(Z[inx_ls[[i]], , drop = FALSE])
}
## do fpca on the log(1+AC)
fpca_Z <- fpca.face(Y = Zmat, knots = 50)
Zsm <- fpca_Z$Yhat


## Get a data frame for analysis
## which contains one row per participant
df <- data[!duplicated(data$SEQN), ]
## drop the activity count columns
df <-
  df %>%
  dplyr::select(-one_of(paste0("MIN", 1:1440)))
## add in the activity count matrix
## using the AsIs class via I()
## note!! be careful when working
## with dataframes which contain matrixes

df$Zsm <- I(Zsm)
df$Zraw <- I(Zmat)

## clean up the workspace a bit
rm(Zsm)
rm(Zmat)
rm(Z)


## set up the functional domain matrix
## mgcv will use this to construct the basis \phi_k^\gamma(s)
sind <- seq(0, 1, len = 1440)
smat <- matrix(sind, nrow(df), 1440, byrow = TRUE)
df$smat <- I(smat)

## set up the matrix of integration weights
df$lmat <- I(matrix(1 / 1440, nrow(df), 1440))

## multiply integration weights by the functional predictor
df$zlmat <- I(df$lmat * df$Zsm)
fit_fglm_ps <- gam(mort_5yr ~ s(smat, by = zlmat, bs = "cc", k = 30),
                   data = df,
                   method = "REML",
                   family = binomial)
