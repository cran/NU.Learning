\name{lrcagg}
\alias{lrcagg}
\title{Calculate the observed Distribution of LRCs in NU.Learning}
\description{For a given number, K, of Clusters of Experimental Units in baseline X-covariate
  space, lrcagg() calculates the observed distribution of "Local Rank Correlations" (LRCs)
  across Clusters ...where each LRC = cor(trex, Y, method = "spearman") within a Cluster,
  trex is a numeric measure of Exposure, and Y is a numeric measure of Outcome.}
\usage{
  lrcagg(K, envir)
}
\arguments{
  \item{K}{Number of Clusters in baseline X-covariate space.}
  \item{envir}{R environment output by a previous call to NUsetup().}
}
\details{Multiple calls to lrcagg(K) for varying numbers of clusters, K, are typically made
  after first invoking NUcluster() to hierarchically cluster patients in X-space and then
  invoking NUsetup() to specify a Y Outcome variable and a continuous, numerical treatment
  Exposure: trex. lrcagg() computes an observed LRC Distribution, updates information stored
  in its envir object, and outputs an object that is typically saved in the user's .GlobalEnv
  to allow subsequent use by print.lrcagg(), plot.lrcagg(), confirm() or KSperm(). Uninformative
  Clusters (those containing only 1 or 2 experimental units) contribute \bold{NA}
  values to the LRCtabl$LRC and LRCdist$LRC objects within the lrcagg() output list.   
}
\value{ An output list of 12 objects, of class lrcagg:
  \item{hclobj}{Name of clustering dendrogram object created by NUcluster().}
  \item{dframe}{Name of data.frame containing X, trex & Y variables.}
  \item{trex}{Name of numerical treatment/exposure level variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{K}{Number of Clusters Requested.}
  \item{actclust}{Number of Clusters delivered.}
  \item{LRCtabl}{data.frame with 5 columns and K rows for Clusters.}
  \item{LRCtabl$c}{Cluster ID Factor, "1", "2", ..., "K".}
  \item{LRCtabl$LRC}{Numerical value of Local Treatment Difference for a Cluster.}
  \item{LRCtabl$w}{Integer value of "weight" = Cluster Size.}
  \item{LRCtabl$LAO}{Numerical value of within-cluster Local Average Outcome (Y-value).} 
  \item{LRCtabl$PS}{Numerical value of Local Relative Propensity for Exposure, 0.0 to 1.0.}   
  \item{LRCdist}{data.frame with 5 columns and same number of rows as the data: dframe.}
  \item{LRCdist$c.K}{Cluster ID Variable of the form: "c.K" }
  \item{LRCdist$ID}{Observation ID Variable for the rows of the input dframe.}
  \item{LRCdist$y}{Numerical values of Y-Outcomes for Experimental Units.}
  \item{LRCdist$t}{Numerical values of Treatment-Exposure Levels for Experimental Units.}
  \item{LRCdist$LRC}{Numerical values of the LRC for the Cluster containing each Unit.} 
  \item{infoclus}{Integer value of Number of Informative Clusters.}
  \item{infounits}{Integer value of Number of Units within Informative Clusters.}
  \item{LRCmean}{Numerical value of mean(LRCdist$LRC) = Weighted Average of
    LRCtabl$LRC values.}
  \item{LRCstde}{Numerical value of sqrt(var(LRCdist$LRC)) = Weighted Standard
    Deviation of LRCtabl$LRC values.}
}
\references{
  Obenchain RL. (2010) The Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192.

  Obenchain RL. (2019) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ivadj}}, \code{\link{ltdagg}} and \code{\link{NUcompare}}.}
\examples{
    data(radon)
    xvars = c("obesity", "over65", "cursmoke")
    hclobj = NUcluster(radon, xvars)
    e = NUsetup(hclobj, radon, lnradon, lcanmort)
    lrc050 = lrcagg(50, e)
    lrc050
    plot(lrc050, e)
}
\keyword{nonparametric}