\name{plot.ivadj}
\alias{plot.ivadj}
\title{Display an Instrumental Variable (LAO) plot with Linear and smooth.spline Fits}
\description{For a given number of patient clusters, K, in baseline X-covariate space and a
  specified Y-outcome variable, display the distribution of Local Average Outcomes (LAOs)
  plotted versus Within-Cluster Propensity-like Scores: Treatment Selection Fractions
  or Relative Exposure Levels.
}
\usage{
  \method{plot}{ivadj}(x, maxsiz = 0.15, \dots)
}
\arguments{
  \item{x}{An object output by ivadj() for K given Clusters in baseline X-covariate space.}
  \item{maxsiz}{Radius of the Circle plotting symbol for the largest Cluster. Usually < 0.6}
  \item{\dots}{Other arguments passed on to plot().}
}
\value{ NULL }
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{plot.ltdagg}}, \code{\link{plot.lrcagg}}}
\keyword{methods}
