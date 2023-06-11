library(pompp)
suppressWarnings(suppressMessages(library(geoR)))
library(glue)

maxSims <- 30

checkLastSim <- function() {
  path <- file.path("data", "tracker.rds")
  if (!file.exists(path)) return(0)

  max(readRDS(path))
}
updateSim <- function(newSim) {
  path <- file.path("data", "tracker.rds")
  if (file.exists(path)) {
    current <- readRDS(path)
    saveRDS(c(current, newSim), path)
  } else saveRDS(newSim, path)
}
saveOtherParams <- function(newSim, totalZprime, varZprime, nU, nXp, nX, s2, phi) {
  path <- file.path("data", "otherParams.rds")
  if (file.exists(path)) {
    current <- readRDS(path)
    saveRDS(rbind(current, c(newSim, totalZprime, varZprime, nU, nXp, nX, s2, phi)), path)
  } else saveRDS(data.frame(
    i = newSim, totalZprime = totalZprime, varZprime = varZprime, nU = nU,
    nXp = nXp, nX = nX, sigma2 = s2, phi = phi
  ), path)
}

while (checkLastSim() < maxSims) {
  updateSim(i <- checkLastSim() + 1)
  cat(glue("Starting simulation {i}."), "\n")
  set.seed(i)

  beta <- c(-2, -1, 2, -1.5)
  delta <- c(-1, -1, -2)
  lambdaStar <- 6000
  gamma <- 2
  mu <- 5

  total_points <- rpois(1, lambdaStar)
  random_points <- cbind(runif(total_points), runif(total_points))
  grid_size <- 50

  # Find covariate values to explain the species occurrence.
  # We give them a Gaussian spatial structure.
  reg_grid <- as.matrix(expand.grid(seq(0, 1, len = grid_size),
                                    seq(0, 1, len = grid_size)))
  Z <- cbind(
    grf(1, grid = rbind(random_points, reg_grid), cov.model = "exponential",
        cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data,
    grf(1, grid = rbind(random_points, reg_grid), cov.model = "exponential",
        cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data,
    grf(1, grid = rbind(random_points, reg_grid), cov.model = "exponential",
        cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data
  )
  Z1 <- Z[1:total_points, ]; Z2 <- Z[-(1:total_points), ]

  # Thin the points by comparing the retaining probabilities with uniforms
  # in the log scale to find the occurrences
  u <- runif(total_points)
  # occurrences <- log(-log(1 - u)) <= beta[1] + Z1 %*% beta[-1] # cloglog link
  occurrences <- log(u) - log1p(-u) <= beta[1] + Z1 %*% beta[-1] # logit link
  n_occurrences <- sum(occurrences)
  occurrences_points <- random_points[occurrences,]
  occurrences_Z <- Z1[occurrences, ]

  # Find covariate values to explain the observation bias.
  # Additionally create a regular grid to plot the covariate later.
  W <- cbind(
    grf(1, grid = rbind(occurrences_points, reg_grid), cov.model = "exponential",
        cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data,
    grf(1, grid = rbind(occurrences_points, reg_grid), cov.model = "exponential",
        cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data
  )
  W1 <- W[1:n_occurrences, ]; W2 <- W[-(1:n_occurrences), ]
  S <- grf(1, grid = occurrences_points, cov.model = "exponential",
           cov.pars = c(2, 0.3), kappa = 0, method = "cholesky", messages = FALSE)$data

  # Find the presence-only observations.
  u <- runif(n_occurrences)
  po_sightings <- log(u) - log1p(-u) <= delta[1] + W1 %*% delta[-1] + gamma * S
  marks <- exp(mu + S + rnorm(n_occurrences, 0, sqrt(0.5)))
  n_po <- sum(po_sightings)
  po_points <- occurrences_points[po_sightings, ]
  po_Z <- occurrences_Z[po_sightings, ]
  po_W <- W1[po_sightings, ]
  po_marks <- marks[po_sightings]

  cat(glue("Simulation ended with {n_po} recorded sightings."), "\n")

  geodata <- as.geodata(cbind(po_points, po_marks))
  maxlik <- likfit(geodata, lambda = 0,
                   ini.cov.pars = c(1, 1), messages = FALSE,
                   fix.nugget = FALSE, cov.model = "exponential")

  saveOtherParams(i, sum(marks[-po_sightings]), var(marks[-po_sightings]),
                  sum(!occurrences), sum(!po_sightings), n_po, maxlik$sigmasq, maxlik$phi)

  time <- Sys.time()
  runMCMC1 <- fit_pompp(
    rep(0, 4), rep(0, 4), 10,
    list(mean = rep(0, 4), covariance = 100 * diag(4)),
    list(mean = rep(0, 4), covariance = 100 * diag(4)),
    list(a = 0.001, b = 0.001),
    cbind(Z2, W2, reg_grid),
    1, 0, 1, list(mean = 0, variance = 100),
    list(a = 0.001, b = 0.001),
    cbind(po_Z, po_W),
    po_marks, po_points,
    1:3, 4:5, 1:3, 4:5,
    sqrt(2), maxlik$sigmasq, maxlik$phi,
    20, 6, 7,
    10000, 3, 30000, 2
  )
  print(Sys.time() - time)

  saveRDS(list(runMCMC), glue("chains/fit{i}.rds"))
}

