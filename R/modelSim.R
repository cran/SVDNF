# Run help(modelSim) for description and argument details.
modelSim <- function(t, mu = 0.038, kappa = 3.689, theta = 0.032, sigma = 0.446, rho = -0.745,
                     omega = 5.125, delta = 0.03, alpha = -0.014,
                     rho_z = -1.809, nu = 0.004, p = 0.01, phi = 0.965,
                     model = "Heston", h = 1 / 252,
                     mu_x, mu_y, sigma_x, sigma_y,
                     j_x = FALSE,
                     j_y = FALSE,
                     jump_dist = rpois,
                     mu_x_params, mu_y_params, sigma_x_params, sigma_y_params,
                     jump_params = 0) {
  # Default jump vectors
  j_x_vec <- rep(0, length = t)
  j_y_vec <- rep(0, length = t)

  # Defining drifts/diffusions and jump distribution for various models:
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

    # Volatility factor drift and diffusion
    mu_x <- function(x, kappa, theta, h) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(kappa, theta, h)

    sigma_x <- function(x, sigma, h) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma, h)

    # Jump distribution for the DuffiePanSingleton Model
    jump_dist <- rpois
    jump_params <- c(h * omega)

    # Setting jumps in both returns and variance
    j_y <- TRUE
    j_x <- TRUE
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

    # Volatility factor drift and diffusion
    mu_x <- function(x, h, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(h, kappa, theta)

    sigma_x <- function(x, sigma, h) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma, h)

    # Jumps for the Bates Model
    jump_dist <- rpois
    jump_params <- c(h * omega)

    # Setting jumps in returns
    j_y <- TRUE
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

    # Volatility factor drift and diffusion
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

    # Volatility factor drift and diffusion
    mu_x <- function(x, theta, phi) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(theta, phi)
    sigma_x <- function(x, sigma) {
      return(sigma)
    }
    sigma_x_params <- list(sigma)

    # Jumps for the PMD Model
    jump_dist <- rbinom
    jump_params <- c(1, p) # size = 1, prob=p for rbinom

    # Setting jumps in returns
    j_y <- TRUE
  }

  # Model of  Taylor (1986) with leverage effect
  if (model == "TaylorWithLeverage") {

    # Returns drift and diffusion
    mu_y <- function(x, mu_y_params) { # Dummy parameter to pass to the do.call function
      return(0)
    }
    mu_y_params <- list(0)
    sigma_y <- function(x, sigma_y_params) { # Dummy parameter to pass to the do.call function
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)

    # Volatility factor drift and diffusion
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
    mu_y <- function(x, mu_y_params) { # Dummy parameter to pass to the do.call function
      return(rep(0, length = length(x))) # return a vector of zeroes with the length of x
    }
    mu_y_params <- list(0)
    sigma_y <- function(x, sigma_y_params) { # Dummy parameter to pass to the do.call function
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)

    # Volatility factor drift and diffusion
    mu_x <- function(x, theta, phi) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(theta, phi)
    sigma_x <- function(x, sigma) {
      return(rep(sigma, length = length(x))) # return a vector of sigmas with the length of x
    }
    sigma_x_params <- list(sigma)
  }

  j_nums <- do.call(jump_dist, c(list(t), jump_params))
  if (j_x == TRUE) { # Generating exponential variance jumps (if they are included in the model)
    j_x_vec <- rgamma(t, shape = j_nums, scale = nu)
  }
  if (j_y == TRUE) { # Generating exponential variance jumps (if they are included in the model)
    j_y_vec <- rnorm(t, mean = j_nums * alpha, sd = sqrt(j_nums * delta^2)) + rho_z * j_x_vec
  }

  # Draw our noise
  epsilon_x <- rnorm(t, 0, 1)
  epsilon_y <- rnorm(t, 0, 1)
  x_zero <- theta # Assume x_0 = theta
  x_vec <- c(x_zero)
  y_vec <- c()
  for (i in 1:t) {
    x_tmin1 <- x_vec[i]

    # Evaluate the drift/diffusion:
    mu_y_eval <- do.call(mu_y, c(list(x_tmin1), mu_y_params))
    mu_x_eval <- do.call(mu_x, c(list(x_tmin1), mu_x_params))

    sigma_y_eval <- do.call(sigma_y, c(list(x_tmin1), sigma_y_params))
    sigma_x_eval <- do.call(sigma_x, c(list(x_tmin1), sigma_x_params))

    # Generate returns and the volatility factor:
    y_t <- mu_y_eval + sigma_y_eval * (rho * epsilon_x[i] + sqrt(1 - rho^2) * epsilon_y[i]) + j_y_vec[i]
    x_t <- mu_x_eval + sigma_x_eval * epsilon_x[i] + j_x_vec[i]

    # Update the sequences:
    x_vec <- c(x_vec, x_t)
    y_vec <- c(y_vec, y_t)
  }
  x_vec <- x_vec[-1] # Remove x_zero = theta
  model_sim <- list(volatility_factor = x_vec, returns = y_vec)
  return(model_sim)
}
