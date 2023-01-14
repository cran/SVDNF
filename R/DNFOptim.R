# Run help(DNFPOptim) for description and argument details.
DNFOptim <- function(dynamics, data, N = 50, K = 20, R = 1, grids = "Default",
                     rho = 0, delta = 0, alpha = 0, rho_z = 0, nu = 0, jump_params_list = "dummy",
                     ...) {
  model <- dynamics$model
  # Check if custom model has repeated arguments
  if(model == 'Custom'){
    mu_y <- dynamics$mu_y; sigma_y <- dynamics$sigma_y
    mu_x <- dynamics$mu_x; sigma_x <- dynamics$sigma_x
    
    # Getting the arguments from drift/diffusion functions
    args_vec <- c(formalArgs(mu_y), formalArgs(sigma_y),
                  formalArgs(mu_x), formalArgs(sigma_x))
    # Adding possible parameters that aren't in the drift/diffusion
    # If set to "var", they will be optimized, or else they are fixed at given values.
    if(rho == 'var'){
      args_vec <- c(args_vec, 'rho')
    }
    else{
      dynamics$rho = rho
    }
    if(delta == 'var'){
      args_vec <- c(args_vec, 'delta')
    } 
    else{
      dynamics$delta = delta
    }
    if(alpha == 'var'){
      args_vec <- c(args_vec, 'alpha')
    }
    else{
      dynamics$alpha = alpha
    }
    if(rho_z == 'var'){
      args_vec <- c(args_vec, 'rho_z')
    }
    else{
      dynamics$rho_z = rho_z
    }
    if(nu == 'var'){
      args_vec <- c(args_vec, 'nu')
    }
    else{
      dynamics$nu = nu
    }
    args_vec <- c(args_vec, jump_params_list)
    
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
      
      dynamics$mu_y_params <- list(mu, alpha, delta, rho_z, nu, omega)
      dynamics$mu_x_params <- list(kappa, theta)
      dynamics$sigma_x_params <- list(sigma)
      
      # Jump distribution for the DuffiePanSingleton Model
      dynamics$jump_params <- c(dynamics$h * omega)
      dynamics$rho <- rho
      dynamics$delta <- delta
      dynamics$alpha <- alpha
      dynamics$rho_z <- rho_z
      dynamics$nu <- nu
    }
    # Model from Bates (1996):
    else if (model == "Bates") {
      mu <- params[1]
      kappa <- params[2]
      theta <- params[3]
      sigma <- params[4]
      rho <- params[5]
      omega <- params[6]
      delta <- params[7]
      alpha <- params[8]
      
      dynamics$mu_y_params <- list(mu, alpha, delta, omega)
      dynamics$mu_x_params <- list(kappa, theta)
      dynamics$sigma_x_params <- list(sigma)
      
      # Jump distribution for the Bates Model
      dynamics$jump_params <- c(dynamics$h * omega)
      dynamics$rho <- rho
      dynamics$delta <- delta
      dynamics$alpha <- alpha
    }

    # Model from Heston (1993):
    else if (model == "Heston") {
      mu <- params[1]
      kappa <- params[2]
      theta <- params[3]
      sigma <- params[4]
      rho <- params[5]
      
      dynamics$mu_y_params <- list(mu)
      dynamics$mu_x_params <- list(kappa, theta)
      dynamics$sigma_x_params <- list(sigma)
      dynamics$rho <- rho
  
      
    }

    # Model from Pitt, Malik, and Doucet (2014)
    else if (model == "PittMalikDoucet") {
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]
      rho <- params[4]
      p <- params[5]
      delta <- params[6]
      alpha <- params[7]
      
      dynamics$mu_x_params <- list(theta, phi)
      dynamics$sigma_x_params <- list(sigma)
      
      # Jumps for the PMD Model
      dynamics$jump_params <- c(1, p) # size = 1, p=p for dbinom
      dynamics$delta <- delta
      dynamics$alpha <- alpha
      dynamics$rho <- rho
      
    }

    # Model of  Taylor (1986) with leverage effect
    else if  (model == "TaylorWithLeverage") {
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]
      rho <- params[4]
      
      dynamics$mu_x_params <- list(theta, phi)
      dynamics$sigma_x_params <- list(sigma)
      dynamics$rho <- rho
    }

    # Model of  Taylor (1986)
    else if (model == "Taylor") {
      phi <- params[1]
      theta <- params[2]
      sigma <- params[3]

      dynamics$mu_x_params <- list(theta, phi)
      dynamics$sigma_x_params <- list(sigma)
    }
    else if (model == "Custom") {
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
        dynamics$mu_y_params <- as.list(params[1:length_mu_y])
        params <- params[-1:-length_mu_y]
      }
      if(formalArgs(sigma_y)[2] == 'dummy'){
        dynamics$sigma_y_params <- list(0)
      }
      else{
        dynamics$sigma_y_params <- as.list(params[1:length_sigma_y])
        params <- params[-1:-length_sigma_y]
      }
      if(formalArgs(mu_x)[2] == "dummy"){
        mu_x_params <- list(0)
      }
      else{
        dynamics$mu_x_params <- as.list(params[1:length_mu_x])
        params <- params[-1:-length_mu_x]
      }
      if(formalArgs(sigma_x)[2] == "dummy"){
        sigma_x_params <- list(0)
      }
      else{
        dynamics$sigma_x_params <- as.list(params[1:length_sigma_x])
        params <- params[-1:-length_sigma_x]
      }
      if(rho == 'var'){
        dynamics$rho <- params[1]
        params <- params[-1]
      }
      if(delta == 'var'){
        dynamics$delta <- params[1]
        params <- params[-1]
      }
      if(alpha == 'var'){
        dynamics$alpha <- params[1]
        params <- params[-1]
      }
      if(rho_z == 'var'){
        dynamics$rho_z <- params[1]
        params <- params[-1]
      }
      if(nu == 'var'){
        dynamics$nu <- params[1]
        params <- params[-1]
      }
      if(jump_params_list != 'dummy'){
        dynamics$jump_params <- as.list(params)
      }
    }

     if(dynamics$delta < 0 |  dynamics$rho <= -1 |  dynamics$rho >= 1 |  dynamics$nu < 0){
       LL <- -9999999999999999999
       return(-LL)
     }
     if(any(dynamics$model == c("DuffiePanSingleton", "Bates", "Heston"))){
       if(theta <= 0){ 
         LL <- -9999999999999999999
         return(-LL)
       }
       if(kappa <= 0){ 
         LL <- -9999999999999999999
         return(-LL)
       }
       if(sigma <= 0){ 
         LL <- -9999999999999999999
         return(-LL)
       }
       if(dynamics$model != "Heston"){
         if(omega <= 0){
           LL <- -9999999999999999999
           return(-LL)
         }
       } 
     }
     if(any(dynamics$model == c("PittMalikDoucet", "TaylorWithLeverage", "Taylor"))){
       if(sigma <= 0){ 
         LL <- -9999999999999999999
         return(-LL)
       }
       if(dynamics$model == "PittMalikDoucet"){
         if(p <= 0 || p >= 1){
           LL <- -9999999999999999999
           return(-LL)
         }
       }
       if(phi <= 0 || phi >= 1){
         LL <- -9999999999999999999
         return(-LL)
       }
     }
     # Creating grid for listed models
     if (length(grids) == 1) { # Check if users input a custom grid
       grids <- gridMaker(dynamics, R = R, N = N, K = K)
       
     }
    LL <- DNF(data = data, N = N, K = K, R = R, dynamics = dynamics, grids = grids)$log_likelihood
    if (LL == Inf | LL == -Inf) {
      # There can be issues with unusual parameter combinations
      LL <- -9999999999999999999 # Need finite results
    }
    return(-LL)
  }

  optimResults <- optim(fn = MLE_f, ...)
  return(optimResults)
}
