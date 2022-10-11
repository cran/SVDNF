# Run help(DNF) for description and argument details.
DNF <- function(data, N = 50, K = 20, R = 1, mu = 0.038, kappa = 3.689, theta = 0.032,
                sigma = 0.446, rho = -0.745, omega = 5.125, delta = 0.03, alpha = -0.014,
                rho_z = -1.809, nu = 0.004, p = 0.01, phi = 0.965,
                h = 1 / 252, model = "Heston", mu_x, mu_y, sigma_x, sigma_y, jump_dist, jump_params,
                mu_x_params, mu_y_params, sigma_x_params, sigma_y_params, grids = "Default") {
  T <- length(data)
  likelihoods <- seq(from = 0, to = 0, length = T)
  log_likelihood <- 0

  # Creating grid for listed models
  if (length(grids) == 1) { # Check if users input a custom grid
    grids <- gridMaker(
      R = R, N = N, K = K, sigma = sigma, kappa = kappa,
      theta = theta, rho = rho, omega = omega, alpha = alpha,
      delta = delta, mu = mu, nu = nu, rho_z = rho_z,
      h = h, model = model, phi = phi
    )
  }

  # Defining drifts/diffusions and jump distribution for various models (ends at line 188):
  # Model from Duffie, Pan, and Singleton (2000):
  if (model == "DuffiePanSingleton") {
    # Jump compensator
    alpha_bar <- exp(alpha + 0.5 * delta^2) / (1 - rho_z * nu) - 1

    # Returns drift and diffusion
    mu_y <- function(x, mu, alpha_bar, omega, h) {
      return(h * (mu - x / 2 - alpha_bar * omega))
    }
    mu_y_params <- list(mu, alpha_bar, omega, h)
    sigma_y <- function(x, h) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(h)

    # Volatility drift and diffusion
    mu_x <- function(x, kappa, theta, h) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(kappa, theta, h)

    sigma_x <- function(x, sigma, h) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma, h)

    # Jump distribution for the DuffiePanSingleton Model
    jump_dist <- dpois
    jump_params <- c(h * omega)
  }

  # Model from Bates (1996):
  if (model == "Bates") {
    # Jump compensator
    alpha_bar <- exp(alpha + 0.5 * delta^2) - 1

    # Returns drift and diffusion
    mu_y <- function(x, mu, alpha_bar, omega, h) {
      return(h * (mu - x / 2 - alpha_bar * omega))
    }
    mu_y_params <- list(mu, alpha_bar, omega, h)

    sigma_y <- function(x, h) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(h)

    # Volatility drift and diffusion
    mu_x <- function(x, h, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(h, kappa, theta)

    sigma_x <- function(x, sigma, h) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma, h)

    # Jumps for the Bates Model
    jump_dist <- dpois
    jump_params <- c(h * omega)
  }

  # Model from Heston (1993):
  if (model == "Heston") {
    # Returns drift and diffusion
    mu_y <- function(x, h, mu) {
      return(h * (mu - x / 2))
    }
    mu_y_params <- list(h, mu)
    sigma_y <- function(x, h) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(h)

    # Volatility drift and diffusion
    mu_x <- function(x, h, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(h, kappa, theta)

    sigma_x <- function(x, sigma, h) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma, h)
  }

  # Model from Pitt, Malik, and Doucet (2014)
  if (model == "PittMalikDoucet") {

    # Returns drift and diffusion
    mu_y <- function(x, mu_y_params) { # Dummy parameter to pass to the do.call function
      return(0)
    }
    mu_y_params <- list(0)
    sigma_y <- function(x, sigma_y_params) { # Dummy parameter to pass to the do.call function
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)

    # Volatility drift and diffusion
    mu_x <- function(x, theta, phi) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(theta, phi)
    sigma_x <- function(x, sigma) {
      return(sigma)
    }
    sigma_x_params <- list(sigma)

    # Jumps for the PMD Model
    jump_dist <- dbinom
    jump_params <- c(1, p) # size = 1, p=p for dbinom
  }

  # Model of  Taylor (1986) with the leverage effect
  if (model == "TaylorWithLeverage") {

    # Returns drift and diffusion
    mu_y <- function(x, mu_y_params) { # mu_y_params is included to pass to do.call
      return(0)
    }
    mu_y_params <- list(0)
    sigma_y <- function(x, sigma_y_params) { # sigma_y_params is included to pass to do.call
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)

    # Volatility drift and diffusion
    mu_x <- function(x, theta, phi) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(theta, phi)
    sigma_x <- function(x, sigma) {
      return(sigma)
    }
    sigma_x_params <- list(sigma)
  }

  # Model of  Taylor (1986)
  if (model == "Taylor") {
    rho <- 0 # No leverage effects in this model.

    # Returns drift and diffusion
    mu_y <- function(x, mu_y_params) { # mu_y_params is included to pass to do.call
      return(0) # return a vector of zeroes with the length of x
    }
    mu_y_params <- list(0)

    sigma_y <- function(x, sigma_y_params) { # sigma_y_params is included to pass to do.call
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)
    # Volatility drift and diffusion
    mu_x <- function(x, theta, phi) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(theta, phi)
    sigma_x <- function(x, sigma) {
      return(sigma) # return a vector of sigmas with the length of x
    }
    sigma_x_params <- list(sigma)
  }

  # Get the lengths from the grids list (because the grid could be tailor-made and provided by the user)
  var_mid_points <- grids$var_mid_points
  N <- length(var_mid_points)
  R <- length(grids$j_nums)
  K <- length(grids$jump_mid_points)

  # Set a matrix to store the filtering distribution at each time step
  filter_grid <- matrix(NA, nrow = N, ncol = T + 1)
  filter_grid[(N:1), 1] <- rep(1 / N, times = N) # Uniform prior

  d_probs <- probCalculator(
    grids = grids, R = R, N = N, K = K, data = data, rho = rho, rho_z = rho_z, nu = nu, alpha = alpha, delta = delta,
    mu_x = mu_x, mu_y = mu_y, sigma_x = sigma_x, sigma_y = sigma_y, jump_dist = jump_dist,
    mu_x_params = mu_x_params, mu_y_params = mu_y_params, sigma_x_params = sigma_x_params, sigma_y_params = sigma_y_params,
    jump_params = jump_params
  ) # Function to compute measurement equation times the other probabilities.
  for (t in (1:T)) {
    # Time t Likelihood contribution:
    probs <- Cpp_prodfun(d_probs[, t], filter_grid[(N:1), t]) # C++ function to multiply the previous filtering distribution with the appropriate probabilities.
    likelihoods[t] <- sum(probs)
    l = likelihoods[t]
    # Get the posterior filtering distribution:
    if(is.na(l) | is.nan(l)){
      return(list(log_likelihood = -Inf, filter_grid = filter_grid, likelihoods = likelihoods, grids = grids))
    }
    if(l > 0){
      filter_grid[(N:1), t + 1] <- Cpp_rowSums_modN(probs, N) / l
          }
  else{
    filter_grid[, t + 1] <- filter_grid[, t]
  }
  }
  # Compute the likelihood:
  log_likelihood <- sum(log(likelihoods))

  return(list(log_likelihood = log_likelihood, filter_grid = filter_grid, likelihoods = likelihoods, grids = grids))
}
