# Run help(DNFPOptim) for description and argument details.
DNFOptim <- function(data, model = "Heston", N = 50, K = 20, R = 1, h = 1 / 252, grids = "Default",
                     rho = 0, delta = 0, alpha = 0, rho_z = 0, nu = 0,
                     jump_dist = 0, jump_params = 0,
                     mu_x, mu_y, sigma_x, sigma_y, ...) {
  # Check if custom model has repeated arguments
  if(model == 'Custom'){
    # Getting the arguments from drift/diffusion functions
    args_vec <- c(formalArgs(mu_y), formalArgs(sigma_y),
                  formalArgs(mu_x), formalArgs(sigma_x))
    # Set jump_dist to a function with no arguments
    if(typeof(jump_dist) != "closure"){
      jump_dist <- function(){
      return(0)
    }
    }
    # Adding possible parameters that aren't in the drit/diffusion
    if(rho == 'var'){
      args_vec <- c(args_vec, 'rho')
    }
    if(delta == 'var'){
      args_vec <- c(args_vec, 'delta')
    }
    if(alpha == 'var'){
      args_vec <- c(args_vec, 'alpha')
    }
    if(rho_z == 'var'){
      args_vec <- c(args_vec, 'rho_z')
    }
    if(nu == 'var'){
      args_vec <- c(args_vec, 'nu')
    }
    args_vec <- c(args_vec, jump_params)
    
    args_vec <- args_vec[-which(args_vec == "dummy")]
    args_vec <- args_vec[-which(args_vec == "x")]
  }
   MLE_f <- function(params) {
    # Model from Duffie, Pan, and Singleton (2000):
    if (model == "DuffiePanSingleton") {
      mu <- params[1]
      kappa <- params[2]
      theta <- params[3]
      sigma <- params[4]
      rho <- params[5]
      omega <- params[6]
      delta <- params[7]
      alpha <- params[8]
      rho_z <- params[9]
      nu <- params[10]
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
      mu <- params[1]
      kappa <- params[2]
      theta <- params[3]
      sigma <- params[4]
      rho <- params[5]
      omega <- params[6]
      delta <- params[7]
      alpha <- params[8]
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
      mu <- params[1]
      kappa <- params[2]
      theta <- params[3]
      sigma <- params[4]
      rho <- params[5]
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
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]
      rho <- params[4]
      p <- params[5]
      delta <- params[6]
      alpha <- params[7]
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
      mu_x <- function(x, phi, theta) {
        return(theta + phi * (x - theta))
      }
      mu_x_params <- list(phi, theta)
      sigma_x <- function(x, sigma) {
        return(sigma)
      }
      sigma_x_params <- list(sigma)

      # Jumps for the PMD Model
      jump_dist <- dbinom
      jump_params <- c(1, p) # size = 1, p=p for dbinom
    }

    # Model of  Taylor (1986) with leverage effect
    if (model == "TaylorWithLeverage") {
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]
      rho <- params[4]

      # Returns drift and diffusion
      mu_y <- function(x, dummy) { # mu_y_params is included to pass to do.call
        return(0)
      }
      mu_y_params <- list(0)
      sigma_y <- function(x, dummy) { # sigma_y_params is included to pass to do.call
        return(exp(x / 2))
      }
      sigma_y_params <- list(0)

      # Volatility drift and diffusion
      mu_x <- function(x, phi, theta) {
        return(theta + phi * (x - theta))
      }
      mu_x_params <- list(phi, theta)
      sigma_x <- function(x, sigma) {
        return(sigma)
      }
      sigma_x_params <- list(sigma)
    }

    # Model of  Taylor (1986)
    if (model == "Taylor") {
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]
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
    if (model == "Custom") {
      # Checking if arguments are repeated in the the drift/diffusion/jumps functions:
      if(any(duplicated(args_vec))){
        index_formals <- which(duplicated(args_vec, fromLast = TRUE) & !duplicated(args_vec))
        index_copied_params <- c(1:length(args_vec))[which(unique(args_vec) %in% args_vec[index_formals])]
        copied_params <- params[index_copied_params]
        
        copied_formals <- unique(args_vec[which(duplicated(args_vec, fromLast = TRUE))])
        remaining_formals <- args_vec[-index_formals]
        remaining_index <- c(1:length(args_vec))[-index_formals]
        
        unique_formals <- unique(args_vec)
        unique_indices <- remaining_index[-which(args_vec[remaining_index] %in% copied_formals)]
        unique_formals_index <- sort(c(index_formals, unique_indices))
        
        params2 <- rep(params, length = length(args_vec)) 
        params2[unique_formals_index] <- params
        for (i in (1:length(copied_params))){
          if(any(copied_formals[i] ==  remaining_formals)){
            
            for(j in (1:sum(copied_formals[i] ==  remaining_formals))){
              
              copies_index <- remaining_index[which(remaining_formals == copied_formals[i])]
              params2[copies_index[j]] <- copied_params[i]
              
            }
          }
          params <- params2
        }
        #}
      }
      # Getting number of parameters to go in each function
      # (note that we subtract one since x is an argument in each function,
      # but not a parameter that we should count).
      length_mu_y <- length(formals(mu_y)) - 1; length_sigma_y <- length(formals(sigma_y)) - 1
      length_mu_x <- length(formals(mu_x)) - 1; length_sigma_x <- length(formals(sigma_x)) - 1
      if(formalArgs(mu_y)[2] == "dummy"){
        mu_y_params <- list(0)
      }
      else{
        mu_y_params <- as.list(params[1:length_mu_y])
        params <- params[-1:-length_mu_y]
      }
      if(formalArgs(sigma_y)[2] == 'dummy'){
        sigma_y_params <- list(0)
      }
      else{
        sigma_y_params <- as.list(params[1:length_sigma_y])
        params <- params[-1:-length_sigma_y]
      }
      if(formalArgs(mu_x)[2] == "dummy"){
        mu_x_params <- list(0)
      }
      else{
        mu_x_params <- as.list(params[1:length_mu_x])
        params <- params[-1:-length_mu_x]
      }
      if(formalArgs(sigma_x)[2] == "dummy"){
        sigma_x_params <- list(0)
      }
      else{
        sigma_x_params <- as.list(params[1:length_sigma_x])
        params <- params[-1:-length_sigma_x]
      }
      if(rho == 'var'){
        rho <- params[1]
        params <- params[-1]
      }
      if(delta == 'var'){
        delta <- params[1]
        params <- params[-1]
      }
      if(alpha == 'var'){
        alpha <- params[1]
        params <- params[-1]
      }
      if(rho_z == 'var'){
        rho_z <- params[1]
        params <- params[-1]
      }
      if(nu == 'var'){
        nu <- params[1]
        params <- params[-1]
      }
      if(typeof(jump_params) != 'double'){
        jump_params <- as.list(params)
      }

    }
    # Creating grid for listed models
    if (length(grids) == 1) { # Check if users input a custom grid
      grids <- gridMaker(
        R = R, N = N, K = K, sigma = sigma, kappa = kappa,
        theta = theta, rho = rho, omega = omega, alpha = alpha,
        delta = delta, mu = mu, nu = nu, rho_z = rho_z,
        h = h, model = model, phi = phi
      )
    }

     if( delta < 0 |  rho < -1 |  rho > 1){
       LL = -9999999999999999999
     }
     else{
    LL <- DNF(
      data = data, N = N, K = K, R = R,
      mu_x = mu_x, mu_y = mu_y, sigma_x = sigma_x, sigma_y = sigma_y,
      mu_x_params = mu_x_params, mu_y_params = mu_y_params, rho = rho,
      sigma_x_params = sigma_x_params, sigma_y_params = sigma_y_params,
      jump_dist = jump_dist, jump_params = jump_params, 
      delta = delta, rho_z = rho_z, nu = nu, alpha = alpha,
      h = h, model = "Custom", grids = grids
    )$log_likelihood
     }
    if (LL == Inf | LL == -Inf) {
      # There can be issues with unusual parameter combinations at the boundary of the optimizer's exploration
      LL <- -9999999999999999999 # Need finite results
    }
    return(-LL)
  }

  optimResults <- optim(fn = MLE_f, ...)
  return(optimResults)
}
