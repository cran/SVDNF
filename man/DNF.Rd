\name{DNF}
\alias{DNF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete nonlinear filtering algorithm for stochastic volatility models
}
\description{
The \code{DNF} function applies the discrete nonlinear filter (DNF) of Kitagawa (1987) as per the implementation of Bégin & Boudreault (2020) to obtain likelihood evaluations and filtering distribution estimates for a wide class of stochastic volatility models.
}
\usage{
DNF(data, N = 50, K = 20, R = 1,
  mu = 0.038, kappa = 3.689, theta = 0.032,
  sigma = 0.446, rho = -0.745,
  omega = 5.125, delta = 0.03, alpha = -0.014,
  rho_z = -1.809, nu = 0.004,
  p = 0.01, phi = 0.965, h = 1/252,
  model = "Heston",
  mu_x, mu_y, sigma_x, sigma_y,
  jump_dist, jump_params,
  mu_x_params, mu_y_params,
  sigma_x_params, sigma_y_params, grids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A series of asset returns for which we want to run the DNF.}
\item{N}{Number of nodes in the variance grid.}
\item{K}{Number of nodes in the jump size grid.}
\item{R}{Maximum number of jumps used in the numerical integration at each timestep.} 
\item{mu}{Annual expected rate of return.}
\item{kappa}{Variance rate of mean reversion.}
\item{theta}{Unconditional mean variance.}
\item{sigma}{Volatility of the variance.}
\item{rho}{Correlation between the return and the variance noise terms.}
\item{omega}{Jump arrival intensity for models with Poisson jumps.}
\item{delta}{Standard deviation of return jumps.}
\item{alpha}{Mean of return jumps.}
\item{rho_z}{Pseudo-correlation parameter between return and variance jumps.}
\item{nu}{Mean for variance jumps.}
\item{p}{Jump probability for models with Bernoulli jumps.}
\item{phi}{Volatility persistence parameter.}
\item{h}{Time interval between observations (e.g., \code{h = 1/252} for daily data).}
\item{model}{Model used by the discrete nonlinear filter. The options are \code{"Heston"}, \code{"Bates"}, \code{"DuffiePanSingleton"}, \code{"Taylor"}, \code{"TaylorWithLeverage"}, \cr \code{"PittMalikDoucet"}, and \code{"Custom"}. If \code{model = "Custom"}, users should pass the drift functions (i.e., \code{mu_x} and \code{mu_y}), the diffusion functions (i.e., \code{sigma_x} and \code{sigma_y}), and the jump distribution, (i.e., \code{jump_dist}) as well as their parameters to the \code{DNF} function. See Examples.}
\item{mu_x}{Function for variance drift (to be used with a custom model).}
\item{mu_y}{Function for returns drift (to be used with a custom model).}
\item{sigma_x}{Function for variance diffusion (to be used with a custom model).}
\item{sigma_y}{Function for returns diffusion (to be used with a custom model).}
\item{jump_dist}{Probability mass function used to generate return or volatility jumps at each timestep (if both types of jumps are in the model, they are assumed to occur simulaneously).}
\item{jump_params}{List of parameters to be used as arguments in the jump_dist function (parameters should be listed in the order that jump_dist uses them).}
\item{mu_x_params}{List of parameters to be used as arguments in the \code{mu_x} function (parameters should be listed in the order that \code{mu_x} uses them).}
\item{mu_y_params}{List of parameters to be used as arguments in the \code{mu_y} function (parameters should be listed in the order that \code{mu_y} uses them).}
\item{sigma_x_params}{List of parameters to be used as arguments in the \code{sigma_x} function (parameters should be listed in the order that \code{sigma_x} uses them).}
\item{sigma_y_params}{List of parameters to be used as arguments in the \code{sigma_y} function (parameters should be listed in the order that \code{sigma_y} uses them).}
\item{grids}{Grids to be used for numerical integration by the \code{DNF} function. The \code{DNF} function will create grids for built-in models. However, this arguments must be used for custom models. It should contain a list of three sequences: \code{var_mid_points} (variance mid-point sequence), \code{j_nums} (sequence for the number of jumps), and \code{jump_mid_points} (jump mid-point sequence). If there are no variance jumps in the model, set \code{jump_mid_points} equal to zero. If there are no jumps in the model, both \code{j_nums} and \code{jump_mid_points} should be set to zero.}
}
\value{
\item{log_likelihood}{Log-likelihood evaluation based on the DNF.}
\item{filter_grid}{Grid of dimensions \code{N} by \code{T+1} that stores each timesteps' filtering distributions (we assume the filtering distribution is uniform at \eqn{t = 0}).}
\item{likelihoods}{Likelihood contribution at each time-step throughout the series.}
\item{grids}{List of grids used for numerical integration by the DNF.}
}
\references{
Bégin, J.F., Boudreault, M. (2020). Likelihood evaluation of jump-diffusion models using deterministic nonlinear filters. Journal of Computational and Graphical Statistics, 30(2), 452–466.

Kitagawa, G. (1987). Non-Gaussian state-space modeling of nonstationary time series.
Journal of the American Statistical Association, 82(400), 1032-1041.
}

\examples{
set.seed(1)
# Generate 250 returns from the DuffiePanSingleton model
DuffiePanSingleton_sim = modelSim(t = 200, model = "DuffiePanSingleton") 

# Run DNF on the data
dnf_filter = DNF(data = DuffiePanSingleton_sim$returns,
  model = "DuffiePanSingleton") 

# Prints log-likelihood evaluation.
dnf_filter$log_likelihood

# Using a custom model.
# Here, we define the DuffiePanSingleton model as a custom model
# to get the same log-likelihood found using the built-in option

# Daily observations
h <- 1/252

# Parameter values 
mu <- 0.038; kappa <- 3.689; theta <- 0.032
sigma <- 0.446; rho <- -0.745; omega <- 5.125
delta <- 0.03; alpha <- -0.014; rho_z <- -1.809; nu <- 0.004

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
jump_dist <- dpois
jump_params <- c(h * omega)
      
# Define the grid for DNF
N <- 50; R <- 1; K <- 20
var_mid_points <- seq(from = sqrt(0.0000001),
  to = sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), length = N)^2
  
j_nums <- seq(from = 0, to = R, by = 1)

jump_mid_points <- seq(from = 0.000001, to = (3 + log(K)) * sqrt(R) * nu, length = K)

grids <- list(var_mid_points = var_mid_points,
  j_nums = j_nums, jump_mid_points = jump_mid_points)

# Run the DNF function with the custom model
dnf_custom <- DNF(data = DuffiePanSingleton_sim$returns, grids = grid, model = 'Custom',
  mu_x = mu_x, mu_y = mu_y, sigma_x = sigma_x, sigma_y = sigma_y,
  mu_x_params = mu_x_params, mu_y_params = mu_y_params,
  sigma_x_params = sigma_x_params, sigma_y_params = sigma_y_params,
  jump_params = jump_params, jump_dist = jump_dist,
  nu = nu, rho_z = rho_z)
  
# Check if we get the same log-likelihoods
dnf_custom$log_likelihood; dnf_filter$log_likelihood 
}
