library(pompp)
library(bayesplot)
library(tidyverse)
library(GGally)

sardines <- read.csv("cleaned_south_data.txt", sep = " ")
depthMat <- as.data.frame(
  readRDS("depthMatrix.rds"))
indices <- which(
  depthMat$long > -9 & depthMat$long < -7 &
    depthMat$lat > 36.6 & depthMat$lat < 37.5
)
depthMat <- na.omit(depthMat[indices, ])
mmm <- mean(depthMat[, 3], na.rm = TRUE); sdsdsd <- sd(depthMat[, 3])
depthMat[, 3] <- (depthMat[, 3] - mmm) / sdsdsd
sardines$bathy <- (sardines$bathy - mmm) / sdsdsd

sardineMini <- as.matrix(sardines[sardines$sardine > 0, c(1, 2, 5, 13)][-23, ])

maxDist <- sqrt((min(depthMat[, 1]) - max(depthMat[, 1])) ^ 2 + (min(depthMat[, 2]) - max(depthMat[, 2])) ^ 2)


#### Running ####

time <- Sys.time()
runMCMC <- fit_pompp(
  rep(0, 2), rep(0, 3), 10,
  list(mean = rep(0, 2), covariance = 10 * diag(2)),
  list(mean = rep(0, 3), covariance = 10 * diag(3)),
  list(a = 0.001, b = 0.001),
  as.matrix(depthMat),
  1, 0, 1, list(mean = 0, variance = 100),
  list(a = 0.001, b = 0.001),
  sardineMini,
  sardineMini[, 3], sardineMini[, 1:2],
  3, 3, 4, 4,
  maxDist, 0.3693, 0.1742,
  5, 1, 2,
  200, 1, 1000, 6
)
print(Sys.time() - time)

saveRDS(runMCMC, "mcmc.rds")

# Analysis
mcmcMat <- do.call(cbind, runMCMC)
summary(do.call(cbind, runMCMC))
mcmcMat %>% as.data.frame %>% mcmc_trace()
mcmcMat %>% as.data.frame %>% ggpairs()
# 
# 

# Run more

time <- Sys.time()
runMCMC <- fit_pompp(
  c(0.9930604, 2.347191), c(-18.20512, 1.541259, 5.145346), 47135.40,
  list(mean = rep(0, 2), covariance = 10 * diag(2)),
  list(mean = rep(0, 3), covariance = 10 * diag(3)),
  list(a = 0.001, b = 0.001),
  as.matrix(depthMat),
  1, 3.861628, 2.496795, list(mean = 0, variance = 100),
  list(a = 0.001, b = 0.001),
  sardineMini,
  sardineMini[, 3], sardineMini[, 1:2],
  3, 3, 4, 4,
  maxDist, 0.3693, 0.1742,
  5, 1, 2,
  0, 1, 1000, 6
)
print(Sys.time() - time)

saveRDS(runMCMC, "mcmc2.rds")

mcmcMat2 <- do.call(cbind, runMCMC)
summary(mcmcMat2)
mcmcMat2 %>% as.data.frame %>% mcmc_trace()
mcmcMat2[-(1:500),] %>% as.data.frame %>% mcmc_dens()
mcmcMat2[-(1:500)] %>% as.data.frame %>% ggpairs()


# Run even more
runMCMC <- readRDS("mcmc2.rds")
do.call(cbind, runMCMC)[1000,]

time <- Sys.time()
runMCMC <- fit_pompp(
  c(1.170994, 2.710585), c(-19.22068, 2.101303, 5.281989), 44980.07,
  list(mean = rep(0, 2), covariance = 10 * diag(2)),
  list(mean = rep(0, 3), covariance = 10 * diag(3)),
  list(a = 0.001, b = 0.001),
  as.matrix(depthMat),
  1, 3.855177, 2.311119, list(mean = 0, variance = 100),
  list(a = 0.001, b = 0.001),
  sardineMini,
  sardineMini[, 3], sardineMini[, 1:2],
  3, 3, 4, 4,
  maxDist, 0.3693, 0.1742,
  5, 1, 2,
  0, 1, 1000, 6
)
print(Sys.time() - time)

saveRDS(runMCMC, "mcmc3.rds")

mcmcMat3 <- do.call(cbind, runMCMC)
summary(mcmcMat2)
mcmcMat3 %>% as.data.frame %>% mcmc_trace()
mcmcMat2[-(1:500),] %>% as.data.frame %>% mcmc_dens()
mcmcMat2[-(1:500)] %>% as.data.frame %>% ggpairs()

# Run even more
runMCMC <- readRDS("mcmc3.rds")
do.call(cbind, runMCMC)[1000,]

time <- Sys.time()
runMCMC <- fit_pompp(
  c(2.841046e-01, 3.421696), c(-1.638075e+01, -2.225203e-02, 5.081798),
  4.504447e+04,
  list(mean = rep(0, 2), covariance = 10 * diag(2)),
  list(mean = rep(0, 3), covariance = 10 * diag(3)),
  list(a = 0.001, b = 0.001),
  as.matrix(depthMat),
  1, 3.775866, 2.364108, list(mean = 0, variance = 100),
  list(a = 0.001, b = 0.001),
  sardineMini,
  sardineMini[, 3], sardineMini[, 1:2],
  3, 3, 4, 4,
  maxDist, 0.3693, 0.1742,
  5, 1, 2,
  0, 1, 7000, 6
)
print(Sys.time() - time)

saveRDS(runMCMC, "mcmc7.rds")

mcmcMat7 <- do.call(cbind, runMCMC)
summary(mcmcMat2)
mcmcMat7[-(1:4000),] %>% as.data.frame %>% mcmc_trace()
mcmcMat7[-(1:4000),] %>% as.data.frame %>% mcmc_dens()
mcmcMat2[-(1:4000),] %>% as.data.frame %>% ggpairs()


