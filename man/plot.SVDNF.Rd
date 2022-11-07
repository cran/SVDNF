\name{plot.SVDNF}
\alias{plot.SVDNF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DNF Filtering Distribution Plot Function}
\description{
This function plots the median of the filtering distribution estimated from the \code{DNF} function along with user-selected upper and lower percentiles.
}
\usage{
\method{plot}{SVDNF}(x, lower_p = 0.05, upper_p = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{SVDNF} object. The \code{plot} will plot the median and percentiles from its filtering distribution.}
  \item{lower_p}{Lower percentile of the filtering distribution to plot.}
  \item{upper_p}{Upper percentile of the filtering distribution to plot.}
  \item{...}{Other parameters to be passed through to plotting functions.}
}
\value{
No return value, this function generates a plot with the median of the volatility factor obtained from the filtering distribution as well as its upper and lower percentiles from \code{lower_p} and \code{upper_p}.
}
\examples{
set.seed(1)
# Generate 500 returns from the Bates model.
Bates_mod <- dynamicsSVM(model = "Bates")
Bates_sim <- modelSim(t = 500, dynamics = Bates_mod)

# Runs DNF on the data.
dnf_filter <- DNF(data = Bates_sim$returns, dynamics = Bates_mod)
plot(dnf_filter, ylim = c(0, 0.15))
}