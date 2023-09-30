\name{reveal.data}
\alias{reveal.data}
\title{Create a data.frame for use in Prediction of a LTD/LRC effect-size Distribution}
\description{reveal.data() forms a data.frame by sorting and appending the LTD or LRC exposure
  effect-size measures from ltdagg() or lrcagg() -- as well as a Cluster membership-number
  variable -- to a copy of the data.frame specified in NUsetup(). In the fourth and final REVEAL
  Phase of NU.Learning, a stretch-goal is to predict variation in LTD/LRC effect-size
  distributions using the known (baseline) X-covariate characteristics of experimental units. For
  example, the data.frame output by reveal.data() is suitable for input to party::ctree() as well
  as to a number of other "less Visual" prediction methods available in \bold{R}.
}
\usage{reveal.data(x, clus.var="Clus", effe.var="eSiz")}
\arguments{
  \item{x}{An output object resulting from a call to ltdagg() or lrcagg().}
  \item{clus.var}{Quoted NAME for the Cluster-Number variable.}
  \item{effe.var}{Quoted NAME for the LTD/LRC effect-size variable.}
}
\value{ The desired data.frame:
  \item{outdf}{A data.frame containing clus.var, effe.var plus (X, trex & Y) variables.}
}
\references{
  Obenchain RL. (2023) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ltdagg}}, \code{\link{lrcagg}}, and \code{\link{NUsetup}}.}
