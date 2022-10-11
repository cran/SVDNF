# Run help(DNFPlot) for description and argument details.
DNFPlot <- function(dnf, lower_p = 0.05, upper_p = 0.95, ylim = 0) {
  # Compute the volatility's CDF at each time point from the filtering distribution
  filtering <- dnf$filter_grid
  N <- dim(filtering)[1]
  filtering_CDF <- filtering
  for (i in (2:N)) {
    filtering_CDF[i, ] <- filtering_CDF[(i - 1), ] + filtering[i, ]
  }
  # Normalizing so the CDF sums to 1
  filtering_CDF <- filtering_CDF %*% diag(1 / filtering_CDF[N, ])


  var_mid_points <- dnf$grids$var_mid_points
  var_intervals <- c(var_mid_points[1] - (var_mid_points[2] - var_mid_points[1]), var_mid_points, Inf)
  var_intervals <- (var_intervals[1:(N + 1)] + var_intervals[2:(N + 2)]) / 2

  nodes <- var_intervals

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
  p_medians <- percentiles(p = 0.5, CDF = filtering_CDF, nodes)[2:(length(dnf$likelihoods) + 1)]
  p_upper <- percentiles(p = upper_p, CDF = filtering_CDF, nodes)[2:(length(dnf$likelihoods) + 1)]
  p_lower <- percentiles(p = lower_p, CDF = filtering_CDF, nodes)[2:(length(dnf$likelihoods) + 1)]

  # Use grid min/max values as default for ylim
  if (sum(ylim^2) == 0) {
    ylim <- c(min(var_mid_points), max(var_mid_points))
  }

  # Plot the median filtering distribution and add the lower/upper percentiles
  x <- c(1:length(dnf$likelihoods))
  plot(y = p_medians, x = x, ylim = ylim, type = "l", main = "Plot of Filtering Distribution Percentiles", ylab = "Volatility Factor", xlab = "Time", col = "black", lty = 1)
  points(p_upper, type = "l", col = "grey", lty = 2)
  points(p_lower, type = "l", col = "grey", lty = 2)
}
