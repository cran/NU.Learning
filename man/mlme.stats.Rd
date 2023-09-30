\name{mlme.stats}
\alias{mlme.stats}
\title{Print Summary Statistics for One or More "Most-Like-Me" Histogram Pairs.}
\description{Print Summary Statistics for Local effect-size (LTD or LRC) Distributions
  associated with given Numbers of "Nearest-Neighbors" in X-confounder Space.
}
\usage{mlme.stats(x, NN = 50, \dots)}
\arguments{
  \item{x}{An object output by mlme.data().}
  \item{NN}{Number(s) of "Nearest Neighbors" displayed in Histogram(s). NN can be either
    a single integer like NN = 40 or a combination of integers like NN = c( 50, 250,
    2500 ).}
  \item{\dots}{Other arguments passed on to print().}
}
\value{ NULL }
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{plot.mlme}},\code{\link{print.mlme}},\code{\link{mlme}}}
\keyword{methods}
