\name{ltdagg}
\alias{ltdagg}
\title{Calculate the Observed Distribution of LTDs in NU.Learning}
\description{For a given number, K, of Clusters of Experimental Units in baseline X-covariate
  space, ltdagg() calculates the observed distribution of "Local Treatment Differences" (LTDs)
  of the form LTD = (( mean(Y) for units receiving trtm==1 ) - ( mean(Y) for units receiving
  trtm==0 )).}
\usage{
  ltdagg(K, envir)
}
\arguments{
  \item{K}{Number of Clusters in baseline X-covariate space.}
  \item{envir}{R environment output by a previous call to NUsetup().}
}
\details{Multiple calls to ltdagg(K) for varying numbers of clusters, K, are typically made
  after first invoking NUcluster() to hierarchically cluster patients in X-space and then
  invoking NUsetup() to specify a Y Outcome variable and a two-level, numerical treatment
  variable: trtm. ltdagg() computes an observed LTD Distribution, updates information stored
  in its envir object, and outputs an object that is typically saved in the user's .GlobalEnv
  to allow subsequent use by print.ltdagg(), plot.ltdagg(), confirm() or KSperm(). Uninformative
  Clusters (those containing either only trtm==1 units or else only trtm==0 units) contribute
  \bold{NA} values to the LTDtabl$LTD and LTDdist$LTD objects within the ltdagg() output list
  object.   
}
\value{ An output list of 12 objects, of class ltdagg:
  \item{hiclus}{Name of clustering object created by NUcluster().}
  \item{dframe}{Name of data.frame containing X, trtm & Y variables.}
  \item{trtm}{Name of treatment factor variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{K}{Number of Clusters Requested.}
  \item{actclust}{Number of Clusters delivered.}
  \item{LTDtabl}{data.frame with 5 columns and K rows for Clusters.}
  \item{LTDtabl$c}{Cluster ID Factor, "1", "2", ..., "K".}
  \item{LTDtabl$LTD}{Numerical value of Local Treatment Difference for a Cluster.}
  \item{LTDtabl$w}{Integer value of "weight" = Cluster Size.}
  \item{LTDtabl$LAO}{Numerical value of within-cluster Local Average Outcome (Y-value).} 
  \item{LTDtabl$PS}{Numerical value of Propensity Score = Local Fraction of Experimental
    Units receiving trtm==1; 0.0 <= PS <= 1.0.}   
  \item{LTDdist}{data.frame with 5 columns and same number of rows as the data: dframe.}
  \item{LTDdist$c.K}{Factor values within c("1", "2", ..., "K").}
  \item{LTDdist$ID}{Observation ID Variable for the rows of the input dframe.}
  \item{LTDdist$y}{Numerical value of the Y-Outcome for an Experimental Unit.}
  \item{LTDdist$t}{Numerical value of trtm (0 or 1) for an Experimental Unit.}
  \item{LTDdist$LTD}{Numerical value of the LTD for the Cluster containing each Exp. Unit.} 
  \item{infoclus}{Integer value of Number of Informative Clusters.}
  \item{infounits}{Integer value of Number of Units within Informative Clusters.}
  \item{LTDmean}{Numerical value of mean(LTDdist$LTD) = Weighted Average of
    LTDtabl$LTD values.}
  \item{LTDstde}{Numerical value of sqrt(var(LTDdist$LTD)) = Weighted Standard
    Deviation of LTDtabl$LTD values.}
}
\references{
  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192.  

  Obenchain RL. (2019) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ivadj}}, \code{\link{lrcagg}} and \code{\link{NUcompare}}.}
\examples{
  \donttest{
  # Running takes more than 7 seconds...
  data(pci15k)
  xvars = c("stent", "height", "female", "diabetic", "acutemi", "ejfract", "ves1proc")
  hclobj = NUcluster(pci15k, xvars)
  NUe = NUsetup(hclobj, pci15k, thin, surv6mo)
  surv050 = ltdagg(50, NUe)
  surv050
  plot(surv050, NUe)
  }
}
\keyword{nonparametric}