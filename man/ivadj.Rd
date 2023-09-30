\name{ivadj}
\alias{ivadj}
\title{Instrumental Variable LAO Fitting and Smoothing}
\description{For a given number of patient clusters in baseline X-covariate space and a
  specified Y-outcome variable, smooth the distribution of Local Average Outcomes (LAOs)
  plotted versus Within-Cluster Propensity-like Scores: the Treatment Selection Fraction
  or the Relative Exposure Level.
}
\usage{
  ivadj(x)
}
\arguments{
  \item{x}{An output object from ltdagg() or lrcagg() using K Clusters in X-covariate space.}
}
\details{Multiple invocations of ivadj(ltdagg()) or ivadj(lrgagg()) using varying numbers of
  clusters, K, can be made. Each invocation of ivadj() displays a linear lm() fit and a
  smooth.spline() fit to the scatter of LAO estimates plotted versus their within-cluster
  propensity-like score estimates.
}
\value{ An output list object of class ivadj:
  \item{hclobj}{Name of clustering object output by NUcluster().}
  \item{dframe}{Name of data.frame containing X, trtm & Y variables.}
  \item{trtm}{Name of the numeric treatment variable.}
  \item{yvar}{Name of the numeric outcome Y variable.}
  \item{K}{Number of Clusters Requested.}
  \item{actclust}{Number of Clusters actually produced.}
}
\references{
  McClellan M, McNeil BJ, Newhouse JP. (1994) Does More Intensive Treatment of
  Myocardial Infarction in the Elderly Reduce Mortality?: Analysis Using Instrumental
  Variables. \emph{JAMA} \bold{272}: 859-866.

  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192.

  Obenchain RL. (2023) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41-55.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ltdagg}}, \code{\link{lrcagg}} and \code{\link{NUcompare}}.}
\examples{
  \donttest{
  # Running takes about 7 seconds...
  data(pci15k)
  xvars = c("stent", "height", "female", "diabetic", "acutemi", "ejfract", "ves1proc")
  hclobj = NUcluster(pci15k, xvars)
  NU.env = NUsetup(hclobj, pci15k, thin, surv6mo)
  surv050 = ltdagg(50, NU.env)
  iv050 = ivadj(surv050)
  iv050
  plot(iv050)
  }
}
\keyword{nonparametric}
