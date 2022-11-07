# Run help(DNF) for description and argument details.
DNF <- function(dynamics, ...) UseMethod("DNF") 

DNF.default <- function(dynamics, ...){
  stop("This class of state-space models is not supported by the DNF function.")
}

DNF.dynamicsSVM <- function(dynamics, data, N = 50, K = 20, R = 1, grids = "Default", ...) {
  T <- length(data)
  likelihoods <- seq(from = 0, to = 0, length = T)
  log_likelihood <- 0
  
  # Creating grid for listed models
  if (length(grids) == 1) { # Check if users input a custom grid
    grids <- gridMaker(dynamics, R = R, N = N, K = K)
  }
  # Get the lengths from the grids list (because the grid could be tailor-made and provided by the user)
  var_mid_points <- grids$var_mid_points
  N <- length(var_mid_points)
  R <- length(grids$j_nums)
  K <- length(grids$jump_mid_points)

  # Set a matrix to store the filtering distribution at each time step
  filter_grid <- matrix(NA, nrow = N, ncol = T + 1)
  filter_grid[(N:1), 1] <- rep(1 / N, times = N) # Uniform prior

  d_probs <- probCalculator(grids = grids, R = R, N = N, K = K, data = data, dynamics = dynamics)
    
    # Function to compute measurement equation times the other probabilities.
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

  SVDNF <- list(log_likelihood = log_likelihood, filter_grid = filter_grid, likelihoods = likelihoods, grids = grids)
  class(SVDNF) = "SVDNF"
  return(SVDNF)
  }
