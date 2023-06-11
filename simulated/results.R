library(pompp)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(bayesplot)
library(GGally)
library(gridExtra)
library(glue)
gr <- (sqrt(5) + 1) / 2

## Need to setwd() to the folder where this file is.

#### True values ####
beta <- c(-2, -1, 2, -1.5)
delta <- c(-1, -1, -2)
lambdaStar <- 6000
gamma <- 2
mu <- 5
tau2 <- 0.5
sigma2 <- 2
phi <- 0.3

#### Getting data ####
chainToArray <- function(path) {
  chain <- readRDS(path)
  mcmc_list <- list()
  for (i in 1:1) {
    eval(parse(text = glue('mcmc_mat{i} <- do.call(cbind, chain[[{i}]])')))
    eval(parse(text = glue('colnames(mcmc_mat{i})[1:8] <- c(paste0("beta", 0:3),
                             paste0("delta", 0:2), "gamma")')))
    mcmc_list[[i]] <- eval(parse(text = glue('mcmc_mat{i}')))
  }
  array(
    do.call(rbind, mcmc_list),
    dim = c(nrow(mcmc_list[[1]]), length(mcmc_list),
            ncol(mcmc_list[[1]])),
    dimnames = list(
      iterations = NULL,
      chains = paste("chain:", 1:length(mcmc_list)),
      parameters = colnames(mcmc_mat1)
    )
  )
}

dfs_to_ggplot <- function(dfs, trueSum, trueVar, trueNu, trueNxp) {
  nn  <- length(dfs)
  b0  <- do.call(c, lapply(dfs, function(d) d[, 1]))
  b1  <- do.call(c, lapply(dfs, function(d) d[, 2]))
  b2  <- do.call(c, lapply(dfs, function(d) d[, 3]))
  b3  <- do.call(c, lapply(dfs, function(d) d[, 4]))
  d0  <- do.call(c, lapply(dfs, function(d) d[, 5]))
  d1  <- do.call(c, lapply(dfs, function(d) d[, 6]))
  d2  <- do.call(c, lapply(dfs, function(d) d[, 7]))
  g   <- do.call(c, lapply(dfs, function(d) d[, 8]))
  ls  <- do.call(c, lapply(dfs, function(d) d[, 9]))
  mu2 <- do.call(c, lapply(dfs, function(d) d[, 10]))
  t2  <- do.call(c, lapply(dfs, function(d) d[, 11]))
  dt <- data.frame(Sample = c(
    b0, b1, b2, b3, d0, d1, d2, g, ls, mu2, t2
  ), dataset = factor(rep(rep(1:nn, each = nrow(dfs[[1]])), 11)),
  parameter = rep(c(
    "beta[0]", "beta[1]", "beta[2]", "beta[3]",
    "delta[0]", "delta[1]", "delta[2]", "gamma",
    "lambda^'*'", "mu", "tau^2"
  ), each = length(b0)),
  trueV = rep(c(beta, delta, gamma, lambdaStar, mu, tau2), each = length(b0)))
  ggplot(dt, aes(dataset, Sample)) + geom_violin() + # geom_boxplot() +
    geom_hline(aes(yintercept = trueV)) +
    facet_wrap(~ parameter, labeller = label_parsed, scales = "free")
}

parseFitName <- function(path) {
  filename <- basename(path)
  as.integer(strsplit(strsplit(filename, ".", fixed = TRUE)[[1]], "fit", fixed = TRUE)[[1]][2])
}

#### Max likelihood parameters ####
estMaxLikeCorrect <- trues <- readRDS(file.path("data", "otherParams.rds")) %>%
  arrange(i) %>% select(i, sigma2, phi)
estMaxLikeIncorrect <- trues <- readRDS(file.path("misspec", "data", "otherParams.rds")) %>%
  arrange(i) %>% select(i, sigma2, phi)
estMaxLike <- rbind(
  estMaxLikeCorrect %>% mutate(generator = "Correctly specified"),
  estMaxLikeIncorrect %>% mutate(generator = "Misspecified")
)
(g <- ggplot(estMaxLike, aes(i, sigma2)) +
    geom_col(col = "black", fill = "white") +
    geom_hline(yintercept = sigma2, linewidth = 2) +
    labs(y = expression(sigma^2)) +
    scale_x_continuous(breaks = 1:30, name = "dataset") +
    facet_wrap(~ generator) +
    theme(text = element_text(size = 16),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()))
ggsave("figs/maxlike_sigma2.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)
(g <- ggplot(estMaxLike, aes(i, phi)) +
    geom_col(col = "black", fill = "white") +
    geom_hline(yintercept = phi, linewidth = 2) +
    labs(y = expression(phi)) +
    scale_x_continuous(breaks = 1:30, name = "dataset") +
    facet_wrap(~ generator) +
    theme(text = element_text(size = 16),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()))
ggsave("figs/maxlike_phi.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)

#### Correctly specified ####
fitted <- sort(sapply(list.files("chains"), parseFitName))
allArr <- vector(mode = "list", length = length(fitted))
j <- 1
for (i in fitted) {
  eval(parse(text = glue('arr{i} <- chainToArray("chains/fit{i}.rds")')))
  allArr[[j]] <- get(glue("arr{i}"))
  j <- j + 1
}
(g <- allArr %>% lapply(as.data.frame) %>% dfs_to_ggplot() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
ggsave("figs/correct_params.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)

#### Misspecified ####
mfitted <- sort(sapply(list.files("misspec/chains"), parseFitName))
allMarr <- vector(mode = "list", length = length(mfitted))
j <- 1
for (i in mfitted) {
  eval(parse(text = glue('marr{i} <- chainToArray("misspec/chains/fit{i}.rds")')))
  allMarr[[j]] <- get(glue("marr{i}"))
  j <- j + 1
}
(g <- allMarr %>% lapply(as.data.frame) %>% dfs_to_ggplot() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
ggsave("figs/incorrect_params.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)

#### Individual per dataset ####
bias <- function(dfs, trueSum, trueVar, trueNu, trueNxp) {
  nn   <- length(dfs)
  nu   <- do.call(c, lapply(1:nn, function(i) dfs[[i]][, 12] - trueNu))
  nxp  <- do.call(c, lapply(1:nn, function(i) dfs[[i]][, 13] - trueNxp))
  sumz <- do.call(c, lapply(1:nn, function(i) dfs[[i]][, 14] - trueSum))
  varz <- do.call(c, lapply(1:nn, function(i) dfs[[i]][, 15] - trueVar))
  dt <- data.frame(rb = c(
    nu, nxp, sumz, varz
  ), dataset = factor(rep(rep(1:nn, each = nrow(dfs[[1]])), 4)),
  parameter = rep(c(
    c("nU", "nXp", "Sum of Z'", "Variance of Z'")
  ), each = length(nu)),
  noBias = rep(0, each = length(nu) * 4))
  ggplot(dt, aes(dataset, rb)) + geom_violin() + # geom_boxplot() +
    geom_hline(aes(yintercept = noBias)) +
    facet_wrap(~ parameter, scales = "free") +
    ylab("Bias")
}
trues <- readRDS(file.path("data", "otherParams.rds")) %>% arrange(i)
(g <- allArr %>% lapply(as.data.frame) %>% bias(trues$totalZprime,
                                                   trues$varZprime,
                                                   trues$nU,
                                                   trues$nXp))
ggsave("figs/bias_correct.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)

trues <- readRDS(file.path("misspec", "data", "otherParams.rds")) %>% arrange(i)
(g <- allMarr %>% lapply(as.data.frame) %>% bias(trues$totalZprime,
                                                trues$varZprime,
                                                trues$nU,
                                                trues$nXp))
ggsave("figs/bias_incorrect.eps", plot = g, device = "eps", width = 8,
       height = 8 / gr)


