\name{plot.ltdagg}
\alias{plot.ltdagg}
\title{Display Visualizations of an Observed LTD Distribution in NU.Learning}
\description{Display a Histogram, Box-Whisker Diagram and/or empirical Cumulative Distribution
  Function depicting the Observed Local Treatment Difference (LTD) Distribution across K Clusters.
}
\usage{
  \method{plot}{ltdagg}(x, envir, show="all", breaks="Sturges", \dots)
}
\arguments{
  \item{x}{An object output by ltdagg() for K = Number of Clusters in baseline X-covariate space.}
  \item{envir}{R environment output by a previous call to NUsetup().}
  \item{show}{Choice of "all", "seq", "hist", "boxp", or "ecdf".}
  \item{breaks}{Parameter setting for hist(); May be an integer value ...like 25 or 50.}
  \item{\dots}{Other arguments passed on to plot().}
}
\value{ NULL }
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{plot.lrcagg}}}
\keyword{methods}
