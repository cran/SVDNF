# Run help(plot.SVDNF) for description and argument details.
plot.SVDNF <- function(x, lower_p = 0.05, upper_p = 0.95, tlim = "default", location = 'topright', ...) {
  # Extract grids from x
  var_mid_points <- x$grids$var_mid_points
  jump_mid_points <- x$grids$jump_mid_points
  j_num <- x$grids$j_nums
  
  # Extract dynamics from x
  dynamics = x$dynamics
  # Define N, K, and R from grids
  N <- length(var_mid_points); K <- length(jump_mid_points); R <- max(j_num)
  
  # Filtering distribution from the DNF and length of the series
  filtering <- x$filter_grid
  T <- dim(filtering)[2]-1 # remove one as filtering from t = 0 ,..., t = T
  
  NNKR_grid <- expand.grid(var_mid_points,var_mid_points, jump_mid_points, j_num)
  v_t = unlist(NNKR_grid[1]); v_tmin1 <- unlist(NNKR_grid[2]); j_mt <- unlist(NNKR_grid[3]); j_nums <- unlist(NNKR_grid[4])
  
  # Transition probabilities : 
  p_v <- dnorm(v_t, do.call(dynamics$mu_x, c(list(v_tmin1), dynamics$mu_x_params)), sd = do.call(dynamics$sigma_x, c(list(v_tmin1), dynamics$sigma_x_params)))
  
  p_n <- 1 # default value of 1, if there are no jumps in the model
  p_j_vol <- 1 # default value of 1, if there are no volatility jumps in the model
  
  if (any(j_nums != 0)) { # If there are no jumps, leave jump_mat = 1.
    # Probability of having j_nums jumps
    p_n <- do.call(dynamics$jump_density, c(list(j_nums), dynamics$jump_params))
  }
  
  if (any(jump_mid_points != 0)) {
    # Probability of having jumps of certain size within the jump interval grid.
    p_j_vol <- dgamma(j_mt, shape = j_nums, scale = dynamics$nu)
  }
  
  Transition <- p_v * p_n * p_j_vol
  #Sum across jump size:
  Transition <- rowSums(matrix(Transition, nrow = N * N * (R + 1), ncol = K, byrow = F))
  #Sum across jump number:
  Transition <- rowSums(matrix(Transition, nrow = N * N, ncol = R + 1, byrow = F))
  Transition <- matrix(Transition, nrow = N, ncol = N) 
  
  # If only one number for tlim is given, plot entire density
  if((length(tlim) == 1) && tlim != 'default'){
    # Var intervals
    var_intervals <- c(var_mid_points[1] - (var_mid_points[2] - var_mid_points[1]), var_mid_points, Inf)
    var_intervals <- (var_intervals[1:(N + 1)] + var_intervals[2:(N + 2)]) / 2
    # Length of the var_intervals
    lengths <- c(var_intervals[2:(N + 1)]-var_intervals[1:N])
    lengths[N] = Inf
    
    # Adjust by lengths to get the height of the density within the intervals
    #pred <- ((filtering[N:1,tlim - 1] %*% t(Transition)) /lengths)
    pred <- (filtering[N:1,tlim -1] /lengths) %*% t(Transition)
    # Normalizing constant
    filtering_constant <- sum(filtering[N:1,tlim]/lengths)
    
    plot(x = var_mid_points,
         y =  (filtering[N:1,tlim]/lengths)/filtering_constant , type = 'l',
         ylab = "Density", xlab = "Volatility Factor", col = 'magenta',
         main = substitute(paste("Plot of the Prediction and Filtering Distributions at "
                                 *italic(t) *" = ", x), list(x = tlim)), lwd = 1.5)
    points(x = var_mid_points, y = (pred) / sum(pred), type = 'l',
           col = 'blue', lwd = 1.5, lty = 2)
    legend(x = location, legend = c('Prediction', 'Filtering'),
           lwd = c(1.5,1.5), col = c('blue', 'magenta'),
           lty = c(2, 1))
  }
  else{
  pred <- matrix(NA, nrow = N, ncol = T)
  # Compute the prediction distribution at each step 
  for(i in (1:T)){
  pred[N:1,i] <- (x$filter_grid[N:1,i]) %*% t(Transition)
  }
  
  filtering_CDF <- filtering
  pred_CDF <- pred
  
  if (any(tlim == "default")){
    tlim <- c(1:T)
  }
  else{
    tlim <- c(tlim[1]:tlim[2])
  }
  # Compute the volatility's CDF at each time point from the filtering distribution
  for (i in (2:N)) {
    filtering_CDF[i, ] <- filtering_CDF[(i - 1), ] + filtering[i, ]
    pred_CDF[i, ] <- pred_CDF[(i - 1), ] + pred[i, ]
  }
  # Normalize so the CDF sums to 1
  filtering_CDF <- filtering_CDF %*% diag(1 / filtering_CDF[N, ])
  pred_CDF <- pred_CDF %*% diag(1 / pred_CDF[N, ])
  
  # Function to find the percentiles that we want to plot
  percentiles <- function(p, CDF, nodes) {
    per_vec <- c()
    T <- dim(CDF)[2]
    N <- dim(CDF)[1]
    for (i in (1:T)) {
      per_vec <- c(per_vec, nodes[N - which((CDF[, i]) > p)[1] + 1])
    }
    return(per_vec)
  }
  
  # Obtain the median, lower and upper percentiles of the filtering distribution
  p_medians <- percentiles(p = 0.5, CDF = filtering_CDF, var_mid_points)[2:(T+1)]
  p_upper <- percentiles(p = upper_p, CDF = filtering_CDF, var_mid_points)[2:(T+1)]
  p_lower <- percentiles(p = lower_p, CDF = filtering_CDF, var_mid_points)[2:(T+1)]

  # Obtain the median, lower and upper percentiles of the prediction distribution
  p_medians_pred <- percentiles(p = 0.5, CDF = pred_CDF, var_mid_points)[1:T]
  p_upper_pred <- percentiles(p = upper_p, CDF = pred_CDF, var_mid_points)[1:T]
  p_lower_pred <- percentiles(p = lower_p, CDF = pred_CDF, var_mid_points)[1:T]
  
  # Plot the median filtering distribution and add the lower/upper percentiles
  plot(y = p_medians_pred[tlim], x = tlim,
       main = "Plot of Prediction Distribution Percentiles", ...)
  points(p_upper_pred[tlim], x = tlim, type = "l", col = "grey", lty = 2)
  points(p_lower_pred[tlim], x = tlim, type = "l", col = "grey", lty = 2)
  
  # Plot the median filtering distribution and add the lower/upper percentiles
  plot(y = p_medians[tlim], x = tlim,
       main = "Plot of Filtering Distribution Percentiles", ...)
  points(p_upper[tlim], x = tlim, type = "l", col = "grey", lty = 2)
  points(p_lower[tlim], x = tlim, type = "l", col = "grey", lty = 2)
  }
}
