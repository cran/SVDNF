\name{DNFPlot}
\alias{DNFPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DNF filtering distribution plot function}
\description{
The \code{DNFPlot} function plots the median of the filtering distribution estimated from the \code{DNF} function along with user-selected upper and lower percentiles.
}
\usage{
DNFPlot(dnf, lower_p = 0.05, upper_p = 0.95, ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnf}{Result from the \code{DNF} function. \code{DNFPlot} will plot the median and percentiles from its filtering distribution.}
  \item{lower_p}{Lower percentile of the filtering distribution to plot.}
  \item{upper_p}{Upper percentile of the filtering distribution to plot.}
  \item{ylim}{Vector with two elements containing the lower and upper bound for the graph (by default, these bounds will be given by the lower and upper bounds of volatility grid used by the \code{DNF} function).}
}
\value{
No return value, this function generates a plot with the median of the volatility factor obtained from the filtering distribution as well as its upper and lower percentiles from \code{lower_p} and \code{upper_p}.
}
\examples{
set.seed(1)
# Generate 500 returns from the Bates model.
Bates_sim <- modelSim(t = 500, model = 'Bates')

# Runs DNF on the data.
dnf_filter <- DNF(data = Bates_sim$returns, model = "Bates")
DNFPlot(dnf_filter, ylim = c(0, 0.15))
}
