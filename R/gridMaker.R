# Function gridMaker:
#   This function creates grids for the pre-specified models in the SVDNF package.
#
# Inputs:
#    - N, K, R : Sizes for the volatility, volatility jumps, and number of jumps grids.
#   	- models parameters : List of parameters used in the
#                          pre-specified models (see help(DNF) for a description).
#    - model : Model for which we are applying the DNF.
# Output:
#   	- grids : A list that contains the sequence of volatility nodes, var_mid_points,
#              the nodes for the number of jumps, j_nums, and
#              the volatility jump size nodes, jump_mid_points.

gridMaker <- function(R = 1, N = 50, K = 20, mu = 0.038, kappa = 3.689, theta = 0.032,
                      sigma = 0.446, rho = -0.745, omega = 5.125, delta = 0.003, alpha = -0.014,
                      rho_z = -1.809, nu = 0.004, p = 0.01, phi = 0.965,
                      h = 1 / 252, model) {
  # Model from Duffie, Pan, and Singleton (2000):
  if (model == "DuffiePanSingleton") {
    var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), length = N)^2
    j_nums <- seq(from = 0, to = R, by = 1)
    jump_mid_points <- seq(from = 0.000001, to = (3 + log(K)) * sqrt(R) * nu, length = K)
  }
  # Model from Bates (1996):
  if (model == "Bates") {
    var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = max(sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), sqrt(0.15)), length = N)^2
    j_nums <- seq(from = 0, to = R, by = 1)
    jump_mid_points <- 0
  }

  # Model from Heston (1993):
  if (model == "Heston") {
    var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = max(sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), sqrt(0.15)), length = N)^2
    j_nums <- 0
    jump_mid_points <- 0
  }

  # Model from Pitt, Malik, and Doucet (2014)
  if (model == "PittMalikDoucet") {
    mean <- theta
    sd <- sqrt((sigma^2) / (1 - phi^2))
    # One node point at the mean and floor(N/2) points on each side
    var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
    # Pick the first N points (or else some grids generated this way have N+1 points)
    var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
    j_nums <- seq(from = 0, to = 1, by = 1)
    jump_mid_points <- 0
  }

  # Model of  Taylor (1986) with leverage effect
  if (model == "TaylorWithLeverage") {
    mean <- theta
    sd <- sqrt((sigma^2) / (1 - phi^2))

    # One node point at the mean and floor(N/2) points on each side
    var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
    # Pick the first N points (or else some grids generated this way have N+1 points)
    var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
    j_nums <- 0
    jump_mid_points <- 0
  }

  # Model of  Taylor (1986)
  if (model == "Taylor") {
    mean <- theta
    sd <- sqrt((sigma^2) / (1 - phi^2))

    # One node point at the mean and floor(N/2) points on each side
    var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
    # Pick the first N points (or else some grids generated this way have N+1 points)
    var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
    j_nums <- 0
    jump_mid_points <- 0
  }

  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}
