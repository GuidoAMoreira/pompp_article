runMCMC <- readRDS("mcmc.rds")
runMCMC2 <- readRDS("mcmc2.rds")
runMCMC3 <- readRDS("mcmc3.rds")
runMCMC4 <- readRDS("mcmc4.rds")
runMCMC7 <- readRDS("mcmc7.rds")
library(bayesplot)
color_scheme_set("gray")
library(tidyverse)
theme_set(theme_bw())
library(GGally)
library(gridExtra)
gr <- (sqrt(5) + 1) / 2
matricize <- function(...) {
  ll <- list(...)
  do.call(rbind, lapply(ll, function(x) do.call(cbind, x)))
}

mcmcMat <- matricize(runMCMC7)
mcmcMat2 <- mcmcMat[-which(mcmcMat[, 11] < 15e6), ]
mcmcMat2 <- mcmcMat[-(1:3000), ]
summary(mcmcMat)
mcmcMat %>% as.data.frame %>% mcmc_trace()
mcmcMat2 %>% as.data.frame %>% mcmc_dens()
mcmcMat2 %>% as.data.frame %>% dplyr::select(V1, V2) %>% ggpairs()
mcmcMat2 %>% as.data.frame %>% dplyr::select(V2, V4) %>% ggpairs()

pp <- par(mfrow = c(4, 4))
for (i in 1:15) {
  acf(mcmcMat[, i], main = colnames(mcmcMat[i]))
}
par(pp)

mcmcMat[-(1:1000),] %>% as.data.frame %>% select(V2, V4) %>% ggpairs()

#### Saving graphs ####
graphs <- list()
mcmc <- as.data.frame(mcmcMat2)
parNames <- c(
  expression(beta["int, intercept"]), expression(beta["int, depth"]),
  expression(beta["obs, intercept"]), expression(beta["obs, depth"]),
  expression(gamma), expression(lambda^"*"), expression(mu),
  expression(tau^2), expression(n[U]), expression(n["X\'"]),
  "sum of unobs. marks", "var of unobs. marks"
)
for (i in seq_along(mcmc[, -c(13:15)])) {
  graphs[[i]] <- ggplot(data.frame(x = mcmc[, i]), aes(x)) +
    geom_density() + labs(x = parNames[i]) +
    theme(axis.text.y = element_blank())
  if (i == 6)
    graphs[[i]] <- graphs[[i]] +
      scale_x_continuous(breaks = as.numeric(quantile(mcmc[, i], c(0.01, 0.5, 0.99))))
  else if (i == 10)
    graphs[[i]] <- graphs[[i]] +
      scale_x_continuous(breaks = as.numeric(quantile(mcmc[, i], c(0.01, 0.6, 0.99))))
  else if (i == 11)
    graphs[[i]] <- graphs[[i]] +
      scale_x_continuous(breaks = as.numeric(quantile(mcmc[, i], c(0.01, 0.99))))
}
ggsave("application.eps", plot = arrangeGrob(grobs = graphs, ncol = 4),
       device = "eps", width = 8, height = 8 / gr)

setEPS()
postscript("pairs.eps")
mcmcMat2[, 1:4] %>% as.data.frame() %>% ggpairs(columnLabels = 
                                                  c('beta["int, intercept"]',
                                                    'beta["int, depth"]',
                                                    'beta["obs, intercept"]',
                                                    'beta["obs, depth"]'), labeller = "label_parsed")
dev.off()
