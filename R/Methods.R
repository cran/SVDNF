# Run help(plot.SVDNF) for description and argument details.
plot.SVDNF <- function(x, lower_p = 0.05, upper_p = 0.95, 
                       tlim = "default", location = 'topright', ...) {
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

# Run help(plot.DNFOptim) for description and argument details.
plot.DNFOptim <- function(x, lower_p = 0.05, upper_p = 0.95, tlim = "default", location = 'topright', ...) {
dnf <- x$SVDNF
plot(dnf, lower_p = lower_p, upper_p = upper_p, tlim = tlim, location = location, ...)
}

# Method to return the log-likelihood of SVDNF object
logLik.SVDNF <- function(object, ...){
  val <- object$log_likelihood
  attr(val, "nobs") <- length(object$likelihoods)
  class(val) <- "logLik"
  val
}

# Method to return the log-likelihood of SVDNF object
logLik.DNFOptim <- function(object, ...){
  val <- object$SVDNF$log_likelihood
  attr(val, "nobs") <- length(object$SVDNF$likelihoods)
  attr(val, "df") <- length(object$optim$par)
  class(val) <- "logLik"
  val}

# Summary method for DNFOptim objects
summary.DNFOptim <- function(object, confidence = 0.95, ...){
  optim <- object$optim
  SVDNF <- object$SVDNF
  dynamics <- SVDNF$dynamics
  model <- dynamics$model
  
  rho <- object$rho
  alpha <- object$alpha
  delta <- object$delta
  nu <- object$nu
  rho_z <- object$rho_z
  
  if (model == "DuffiePanSingleton") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
    nu <- 'var'
    rho_z <- 'var'
  }
  # Model from Bates (1996):
  else if (model == "Bates") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
  }
  
  # Model from Heston (1993):
  else if (model == "Heston") {
    rho <- 'var'
  }
  
  # Model from Pitt, Malik, and Doucet (2014)
  else if (model == "PittMalikDoucet") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
    
  }
  
  # Model of  Taylor (1986) with leverage effect
  else if  (model == "TaylorWithLeverage") {
    rho <- 'var'
    
  }
  
  jump_params_list <- object$jump_params_list

 # Making a list of the parameter names in the optimization: 
  
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
  args_vec <- unique(args_vec)
  # Check if the Hessian is invertible
  possibleError <- tryCatch(
    solve(optim$hessian),
    error=function(e) e
  )
  
  if(is.null(optim$hessian) | (inherits(possibleError, "error"))){
    tab <- cbind(Estimates = optim$par)
    dimnames(tab) <- list(args_vec)
    colnames(tab) <- list("Estimate")
    if(is.null(optim$hessian)){
      print("To get standard errors and confidence intervals, you should set hessian = TRUE in the DNFOptim function.")
    }
    else if((inherits(possibleError, "error"))){
      print("There was a numerical issue in inverting the Hessian. The issue is :")
      print(possibleError)
    }
  }
  else{
    if(confidence <=0 | confidence >=1){
      stop("Confidence should be between 0 and 1.")
    }
    
    # Z value for confidence intervals
    alpha <- 1 - confidence
    upper_p <- 1 - alpha / 2 
    z <- qnorm(upper_p)
  
    upper_percentage <- paste(100*upper_p,"%")
    lower_percentage <- paste(100-100*upper_p,"%")
  # Standard deviations
  se <- sqrt(diag(solve(optim$hessian)))
  lower_bound <- optim$par - z*se
  upper_bound <- optim$par + z*se
  tab <- cbind(Estimates = optim$par, StdErr = se, lower_bound = lower_bound,
               upper_bound = upper_bound)
  dimnames(tab) <- list(args_vec)
  colnames(tab) <- list("Estimate", "Std Error", lower_percentage, upper_percentage)  
  }
  #str = paste('Model:', object$SVDNF$dynamics$model, "\n")
  #cat(str)
  res = list(model =  object$SVDNF$dynamics$model, coefficients = tab,
             logLik = logLik(object))
 # model_dynamics <- c(paste("Model:", object$SVDNF$dynamics$model)) 
  #cat("Coefficients: \n")
  class(res) <- "summary.DNFOptim"
  res
}
print.summary.DNFOptim <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # Model estimated
  cat("\nModel:\n")
  cat(x$model, "\n")
  # MLE coefficient estimates
  cat("\nCoefficients:\n")
  print(x$coefficients, digits = digits)
  cat("\n")
  cat("\nLog-Likelihood:\n")
  print(x$logLik)
  cat("\n")
  invisible(x)
}

#Predict method for DNFOptim objects. See help(predict.DNFOptim) for more details.
predict.DNFOptim <- function(object, n_ahead = 15, n_sim = 1000, confidence = 0.95, ...){
  dynamics <- object$SVDNF$dynamics
  
  if(confidence <=0 | confidence >=1){
    stop("Confidence should be between 0 and 1.")
  }
  
  # Z value for confidence intervals
  alpha <- 1 - confidence
  upper_p <- 1 - alpha / 2 
  z <- qnorm(upper_p)
  
  var_mid_points <- object$SVDNF$grids$var_mid_points
  T <- length(object$SVDNF$likelihoods)
  final_filtering <- object$SVDNF$filter_grid[,T]
  N <- length(object$SVDNF$grids$var_mid_points)
  # Sampling volatilities from the time T filtering distribution:
  inits_samp <- sample(x = var_mid_points, size = n_sim, replace = TRUE,
                           prob = final_filtering[N:1])
  pred_sim <- sapply(X = inits_samp, FUN = modelSim, t = n_ahead, dynamics = dynamics)
  
  vol_sim<- pred_sim[1,]
  vol_means <- sapply((1:n_ahead), function(i) mean(sapply(vol_sim, "[[", i)))
  vol_sd <- sapply((1:n_ahead), function(i) sd(sapply(vol_sim, "[[", i)))
  
  # 95% confidence intervals
  LB_vol <- vol_means - z * vol_sd
  UB_vol <- vol_means + z * vol_sd
  
  ret_sim <- pred_sim[2,]
  ret_means <- sapply((1:n_ahead), function(i) mean(sapply(ret_sim, "[[", i)))
  ret_sd <- sapply((1:n_ahead), function(i) sd(sapply(ret_sim, "[[", i)))
  
  # 95% confidence intervals
  LB_ret <- ret_means - z * ret_sd
  UB_ret <- ret_means + z * ret_sd
  
  volatility_pred <- list(UB_vol = UB_vol, mean_vol_pred = vol_means,
                          LB_vol = LB_vol)
  
  ret_pred <- list(UB_ret = UB_ret, mean_ret_pred = ret_means,
                          LB_ret = LB_ret)
  
  pred <- list(volatility_pred = volatility_pred, ret_pred = ret_pred,
               object = object, confidence = confidence)
  class(pred) <- "predict.DNFOptim"
  pred
}

# Print method for pred.DNFOptim
print.predict.DNFOptim <- function(x, digits = max(3, getOption("digits") - 3), limit = 5, ...){
  n_ahead <- length(x$ret_pred$mean_ret_pred)
  confidence <- x$confidence
  percent <- confidence * 100
  if(limit == n_ahead){
  # Return predictions
  cat("\nReturn Mean Prediction and", percent ,"% Confidence interval\n")
  cat("\nUpper Bound:\n")
  cat(format(x$ret_pred$UB_ret[1:limit], nsmall = digits), sep = ", ")
  cat("\nMean\n")
  cat(format(x$ret_pred$mean_ret_pred[1:limit], nsmall = digits), sep = ", ")
  cat("\nLower Bound\n")
  cat(format(x$ret_pred$LB_ret[1:limit], nsmall = digits), sep = ", ")
  cat('\n \n')
  # Volatility predictions
  cat("\nVolatility Mean Prediction and", percent ,"% Confidence interval\n")
  cat("\nUpper Bound:\n")
  cat(format(x$volatility_pred$UB_vol[1:limit], nsmall = digits), sep = ", ")
  cat("\nMean\n")
  cat(format(x$volatility_pred$mean_vol_pred[1:limit], nsmall = digits), sep = ", ")
  cat("\nLower Bound\n")
  cat(format(x$volatility_pred$LB_vol[1:limit], nsmall = digits), sep = ", ")
  
  invisible(x)
  }
  
  if((0 < limit) & (limit < n_ahead)){
    cat("\nReturn Mean Prediction and", percent ,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$ret_pred$UB_ret[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nMean\n")
    cat(format(x$ret_pred$mean_ret_pred[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nLower Bound\n")
    cat(format(x$ret_pred$LB_ret[1:limit], nsmall = digits), sep = ", ", "...")
    cat('\n \n')
    # Volatility predictions
    cat("\nVolatility Mean Prediction and", percent,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$volatility_pred$UB_vol[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nMean\n")
    cat(format(x$volatility_pred$mean_vol_pred[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nLower Bound\n")
    cat(format(x$volatility_pred$LB_vol[1:limit], nsmall = digits), sep = ", ", "...")
    
    invisible(x)
  }
  else{
    stop("The number of predictions printed (limit), should be between 1 and n_ahead (the number of predictions made in the predict function).")
  }
}

plot.predict.DNFOptim <- function(x, ...){
  rets <- x$object$SVDNF$data
  threshold <- length(rets)
  n_ahead <- length(x$ret_pred$mean_ret_pred)
  
  UB_vol <- x$volatility_pred$UB_vol
  LB_vol <- x$volatility_pred$LB_vol
  vol_means <- x$volatility_pred$mean_vol_pred
  
  UB_ret <- x$ret_pred$UB_ret
  LB_ret <- x$ret_pred$LB_ret
  ret_means <- x$ret_pred$mean_ret_pred
  
  filtering <- x$object$SVDNF$filter_grid
  filtering_CDF <- filtering

  var_mid_points <- x$object$SVDNF$grids$var_mid_points
  N <- length(var_mid_points )
  T <- length(rets)
  tlim <- c(1:length(T))
  # Compute the volatility's CDF at each time point from the filtering distribution
  for (i in (2:N)) {
    filtering_CDF[i, ] <- filtering_CDF[(i - 1), ] + filtering[i, ]
  }
  # Normalize so the CDF sums to 1
  filtering_CDF <- filtering_CDF %*% diag(1 / filtering_CDF[N, ])

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
  
  
  vol_filt <- percentiles(p = 0.5, CDF = filtering_CDF, var_mid_points)[2:(T+1)]

  x_seq = seq(from = 1, to = (threshold+n_ahead), by = 1)
  y = c(vol_filt, vol_means)
  v_ylim <- c(min(LB_vol, vol_filt), max(UB_vol, vol_filt))
  # Initialize the plot
  plot(x_seq, y = y, type = "n", xlab = "Time", main = "Plot of Filtered Volatility and Predictions",
       ylab = "Volatility Factor", ylim = v_ylim)
  # Plot the line segments before the threshold in black
  segments(x_seq[1:threshold], y[1:threshold], x_seq[2:(threshold+1)], y[2:(threshold+1)], col = "black")

  # Plot the line segments after the threshold in red
  segments(x_seq[(threshold+1):(length(x_seq)-1)], y[(threshold+1):(length(y)-1)], x_seq[(threshold+2):length(x_seq)], y[(threshold+2):length(y)], col = "magenta")
   
  # Plot the lower confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), LB_vol, col = "blue")

  # Plot the upper confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), UB_vol, col = "blue")

  y = c(rets, ret_means)
  r_ylim <- c(min(c(LB_ret, rets)), max(c(UB_ret, rets)))
  
  # Initialize the plot
  plot(x_seq, y = c(rets, ret_means), type = "n", xlab = "Time", main = "Plot of Past Returns and Predictions",
       ylab = "Returns", ylim = r_ylim)
  # Plot the line segments before the threshold in black
  segments(x_seq[1:threshold], y[1:threshold], x_seq[2:(threshold+1)], y[2:(threshold+1)], col = "black")

  # Plot the line segments after the threshold in red
  segments(x_seq[(threshold+1):(length(x_seq)-1)], y[(threshold+1):(length(y)-1)], x_seq[(threshold+2):length(x_seq)], y[(threshold+2):length(y)], col = "magenta")

  # Plot the lower confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), LB_ret, col = "blue")

  # Plot the upper confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), UB_ret, col = "blue")
}

#Predict method for DNFOptim objects. See help(predict.DNFOptim) for more details.
predict.SVDNF <- function(object, n_ahead = 15, n_sim = 1000,
                          confidence = 0.95, ...){
  dynamics <- object$dynamics
  
  if(confidence <=0 | confidence >=1){
    stop("Confidence should be between 0 and 1.")
  }

  # Z value for confidence intervals
  alpha <- 1 - confidence
  upper_p <- 1 - alpha / 2 
  lower_p <- 1 - upper_p
  z <- qnorm(upper_p)
  
  var_mid_points <- object$grids$var_mid_points
  jump_mid_points <- object$grids$jump_mid_points
  j_num <- object$grids$j_nums
  
  T <- length(object$likelihoods)
  final_filtering <- object$filter_grid[,T]
  N <- length(object$grids$var_mid_points)
  # Sampling volatilities from the time T filtering distribution:
  inits_samp <- sample(x = var_mid_points, size = n_sim, replace = TRUE,
                       prob = final_filtering[N:1])
  pred_sim <- sapply(X = inits_samp, FUN = modelSim, t = n_ahead, dynamics = dynamics)
  
  vol_sim<- pred_sim[1,]
  vol_means <- sapply((1:n_ahead), function(i) mean(sapply(vol_sim, "[[", i)))
  vol_sd <- sapply((1:n_ahead), function(i) sd(sapply(vol_sim, "[[", i)))
  
  # 95% confidence intervals
  LB_vol <- vol_means - z * vol_sd
  UB_vol <- vol_means + z * vol_sd
  
  ret_sim <- pred_sim[2,]
  ret_means <- sapply((1:n_ahead), function(i) mean(sapply(ret_sim, "[[", i)))
  ret_sd <- sapply((1:n_ahead), function(i) sd(sapply(ret_sim, "[[", i)))
  
  # 95% confidence intervals
  LB_ret <- ret_means - z * ret_sd
  UB_ret <- ret_means + z * ret_sd
  
  

  volatility_pred <- list(UB_vol = UB_vol, mean_vol_pred = vol_means,
                          LB_vol = LB_vol)
  
  ret_pred <- list(UB_ret = UB_ret, mean_ret_pred = ret_means,
                   LB_ret = LB_ret)
  
  pred <- list(volatility_pred = volatility_pred, ret_pred = ret_pred,
               object = object, confidence = confidence)
  class(pred) <- "predict.SVDNF"
  pred
}

# Print method for pred.DNFOptim
print.predict.SVDNF <- function(x, digits = max(3, getOption("digits") - 3), limit = 5, ...){
  n_ahead <- length(x$ret_pred$mean_ret_pred)
  confidence <- x$confidence
  percent <- confidence * 100
  if(limit == n_ahead){
    # Return predictions
    cat("\nReturn Mean Prediction and", percent ,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$ret_pred$UB_ret[1:limit], nsmall = digits), sep = ", ")
    cat("\nMean\n")
    cat(format(x$ret_pred$mean_ret_pred[1:limit], nsmall = digits), sep = ", ")
    cat("\nLower Bound\n")
    cat(format(x$ret_pred$LB_ret[1:limit], nsmall = digits), sep = ", ")
    cat('\n \n')
    # Volatility predictions
    cat("\nVolatility Mean Prediction and", percent ,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$volatility_pred$UB_vol[1:limit], nsmall = digits), sep = ", ")
    cat("\nMean\n")
    cat(format(x$volatility_pred$mean_vol_pred[1:limit], nsmall = digits), sep = ", ")
    cat("\nLower Bound\n")
    cat(format(x$volatility_pred$LB_vol[1:limit], nsmall = digits), sep = ", ")
    
    invisible(x)
  }
  
  if((0 < limit) & (limit < n_ahead)){
    cat("\nReturn Mean Prediction and", percent ,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$ret_pred$UB_ret[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nMean\n")
    cat(format(x$ret_pred$mean_ret_pred[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nLower Bound\n")
    cat(format(x$ret_pred$LB_ret[1:limit], nsmall = digits), sep = ", ", "...")
    cat('\n \n')
    # Volatility predictions
    cat("\nVolatility Mean Prediction and", percent,"% Confidence interval\n")
    cat("\nUpper Bound:\n")
    cat(format(x$volatility_pred$UB_vol[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nMean\n")
    cat(format(x$volatility_pred$mean_vol_pred[1:limit], nsmall = digits), sep = ", ", "...")
    cat("\nLower Bound\n")
    cat(format(x$volatility_pred$LB_vol[1:limit], nsmall = digits), sep = ", ", "...")
    
    invisible(x)
  }
  else{
    stop("The number of predictions printed (limit), should be between 1 and n_ahead (the number of predictions made in the predict function).")
  }
}

plot.predict.SVDNF <- function(x, ...){
  rets <- x$object$data
  threshold <- length(rets)
  n_ahead <- length(x$ret_pred$mean_ret_pred)
  
  UB_vol <- x$volatility_pred$UB_vol
  LB_vol <- x$volatility_pred$LB_vol
  vol_means <- x$volatility_pred$mean_vol_pred
  
  UB_ret <- x$ret_pred$UB_ret
  LB_ret <- x$ret_pred$LB_ret
  ret_means <- x$ret_pred$mean_ret_pred
  
  filtering <- x$object$filter_grid
  filtering_CDF <- filtering
  
  var_mid_points <- x$object$grids$var_mid_points
  N <- length(var_mid_points )
  T <- length(rets)
  tlim <- c(1:length(T))
  # Compute the volatility's CDF at each time point from the filtering distribution
  for (i in (2:N)) {
    filtering_CDF[i, ] <- filtering_CDF[(i - 1), ] + filtering[i, ]
  }
  # Normalize so the CDF sums to 1
  filtering_CDF <- filtering_CDF %*% diag(1 / filtering_CDF[N, ])
  
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
  
  
  vol_filt <- percentiles(p = 0.5, CDF = filtering_CDF, var_mid_points)[2:(T+1)]
  
  x_seq = seq(from = 1, to = (threshold+n_ahead), by = 1)
  y = c(vol_filt, vol_means)
  v_ylim <- c(min(c(LB_vol, vol_filt)), max(c(UB_vol, vol_filt)))
  # Initialize the plot
  plot(x_seq, y = y, type = "n", xlab = "Time", main = "Plot of Filtered Volatility with Prediction",
       ylab = "Volatility Factor", ylim = v_ylim)
  # Plot the line segments before the threshold in black
  segments(x_seq[1:threshold], y[1:threshold], x_seq[2:(threshold+1)], y[2:(threshold+1)], col = "black")
  
  # Plot the line segments after the threshold in red
  segments(x_seq[(threshold+1):(length(x_seq)-1)], y[(threshold+1):(length(y)-1)], x_seq[(threshold+2):length(x_seq)], y[(threshold+2):length(y)], col = "magenta")
  
  # Plot the lower confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), LB_vol, col = "blue")
  
  # Plot the upper confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), UB_vol, col = "blue")
  
  y = c(rets, ret_means)
  r_ylim <- c(min(c(LB_ret, rets)), max(c(UB_ret, rets)))
  
  # Initialize the plot
  plot(x_seq, y = c(rets, ret_means), type = "n", xlab = "Time", main = "Plot of Predicted Returns",
       ylab = "Returns", ylim = r_ylim)
  # Plot the line segments before the threshold in black
  segments(x_seq[1:threshold], y[1:threshold], x_seq[2:(threshold+1)], y[2:(threshold+1)], col = "black")
  
  # Plot the line segments after the threshold in red
  segments(x_seq[(threshold+1):(length(x_seq)-1)], y[(threshold+1):(length(y)-1)], x_seq[(threshold+2):length(x_seq)], y[(threshold+2):length(y)], col = "magenta")
  
  # Plot the lower confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), LB_ret, col = "blue")
  
  # Plot the upper confidence interval after the threshold in blue
  lines(threshold+(1:n_ahead), UB_ret, col = "blue")
}


