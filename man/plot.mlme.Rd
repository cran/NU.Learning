\name{plot.mlme}
\alias{plot.mlme}
\title{Display a Pair (or Pairs) of Histograms showing LOCAL effect-sizes for Patients "Most-Like-Me".}
\description{Display Pair(s) of Histograms of Local effect-size (LTD or LRC) Distributions for a specified
  Number (or combinations of Numbers) of "Nearest-Neighbors in X-confounder Space.
}
\usage{
  \method{plot}{mlme}(x, NN=50, breaks=50, \dots)
}
\arguments{
  \item{x}{An object output by mlme().}
  \item{NN}{Number(s) of Nearest Neighbors displayed in Bottom Histogram(s). NN can be
    a single integer like NN = 40 or a combination of integers like NN = c( 50, 250,
    2500 ).}
  \item{breaks}{Integer number of breaks in the Top Histogram for the full LTD or LRC distribution.
    Because the Bottom Histogram may include only a few Nearest Neighbors, it is always displayed
    using breaks = "Sturges".}
  \item{\dots}{Other arguments passed on to plot().}
}
\value{ NULL }
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{mlme.stats}},\code{\link{print.mlme}},\code{\link{mlme}}}
\keyword{methods}
