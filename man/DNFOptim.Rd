\name{DNFOptim}
\alias{DNFOptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete nonlinear filter maximum likelihood estimation function
}
\description{
The \code{DNFOptim} function finds maximum likelihood estimates for stochastic volatility models parameters using the discrete nonlinear filter.
}
\usage{
DNFOptim(data, model = "Heston", N = 50, K = 20, R = 1,
  h = 1/252, grids = 'Default',
  rho = 0, delta = 0, alpha = 0, rho_z = 0, nu = 0,
  jump_dist = 0, jump_params = 0,
  mu_x, mu_y, sigma_x, sigma_y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A series of asset returns for which we want to run the DNF.}
  \item{model}{Model used by the discrete nonlinear filter. The options are \code{"Heston"}, \code{"Bates"}, \code{"DuffiePanSingleton"}, \code{"Taylor"}, \code{"TaylorWithLeverage"}, and \cr \code{"PittMalikDoucet"}.}
  \item{N}{Number of nodes in the variance grid.}
  \item{K}{Number of nodes in the jump size grid.}
  \item{R}{Maximum number of jumps used in the numerical integration at each timestep.} 
  \item{rho, delta, alpha, rho_z, nu}{See \code{help(DNF)} for a description of each of these arguments individually. These arguments should be used only for custom models and can be fixed to a certain value. If they are estimated, they should be set to \code{'var'} (e.g., to estimate \code{rho} set \code{rho = 'var'}) and include it in the vector \code{par} to be passed to the \code{optim} function. See Note for the order in which custom models should receive parameters.}
  \item{h}{Time interval between observations (e.g., \code{h = 1/252} for daily data).}
  \item{grids}{Grids to be used for numerical integration by the \code{DNF} function. The \code{DNF} function will create grids for built-in models. However, this arguments must be used for custom models. It should contain a list of three sequences: \code{var_mid_points} (variance mid-point sequence), \code{j_nums} (sequence for the number of jumps), and \code{jump_mid_points} (jump mid-point sequence). If there are no variance jumps in the model, set \code{jump_mid_points} equal to zero. If there are no jumps in the model, both \code{j_nums} and \code{jump_mid_points} should be set to zero.}
  \item{jump_dist}{Probability mass function used to generate return or volatility jumps at each timestep (if both types of jumps are in the model, they are assumed to occur simulaneously).}
  \item{jump_params}{List of parameters to be used as arguments in the \code{jump_dist} function (parameters should be listed in order that \code{jump_dist} uses them).}
  \item{mu_x}{Function for variance drift (to be used with a custom model).}
  \item{mu_y}{Function for returns drift (to be used with a custom model).}
  \item{sigma_x}{Function for variance diffusion (to be used with a custom model).}
  \item{sigma_y}{Function for returns diffusion (to be used with a custom model).}
  \item{...}{Further arguments to be passed to the \code{optim} function. See Note.}
}
\note{
When passing the initial parameter vector \code{par} to the \code{optim} function (via \code{...}), the parameters should follow a specific order.

For the \code{PittMalikDoucet} model, the parameters should be in the following order: \code{phi}, \code{theta}, \code{sigma}, \code{rho}, \code{p}, \code{delta}, and \code{alpha}.

For the \code{DuffiePanSingleton} model, the parameters should be in the following order: \code{mu}, \code{kappa}, \code{theta}, \code{sigma}, \code{rho}, \code{omega}, \code{delta}, \code{alpha}, \code{rho_z}, and \code{nu}.

All other built-in models can be seen as being nested within these two  models (i.e., \code{Heston} and \code{Bates} models are nested in the \code{DuffiePanSingleton} model, while \code{Taylor} and \cr \code{TaylorWithLeverage} are nested in the \code{PittMalikDoucet} model). Their parameters should be passed in the same order as those in the more general models, minus the parameters not found in these nested models.

For example, the \code{Taylor} model contains neither jumps nor correlation between volatility and returns innovations. Thus, its three parameters are passed in the order: \code{phi}, \code{theta}, and \code{sigma}.

When \code{models = "Custom"}, parameters should be passed in the following order: \code{mu_y_params}, \code{sigma_y_params}, \code{mu_x_params}, \code{sigma_x_params}, \code{rho}, \code{delta}, \code{alpha}, \code{rho_z}, \code{nu}, and \code{jump_params}. If an argument is repeated (e.g., both \code{mu_y_params} and  \code{sigma_y_params} use the same parameter), write it only when it first appears in the custom model order.

}
\value{
Returns a list obtained from \proglang{R}'s\code{optim} function. See \code{help(optim)} for details about the output.
}
\references{
 \proglang{R} Core Team (2019). \proglang{R}: A language and environment for statistical computing. \proglang{R} Foundation
  for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  }
\examples{
set.seed(1)

# Generating return data
Taylor_sim <- modelSim(t = 75, model = "Taylor", phi = 0.9,
  theta = -7.36, sigma = 0.363)
plot(Taylor_sim$volatility_factor, type = 'l')
plot(Taylor_sim$returns, type = 'l')

# Initial values and optimization bounds
init_par <- c( 0.7, -5, 0.3)
lower <- c(0.01, -20, 0.01); upper <- c(0.99, 0, 1)

# Running DNFOptim to get MLEs
optim_test <- DNFOptim(data = Taylor_sim$returns,
  model = "Taylor",
  par = init_par, lower = lower, upper = upper, method = "L-BFGS-B")
                
# Parameter estimates
optim_test$par
}
