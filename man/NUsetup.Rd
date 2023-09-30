\name{NUsetup}
\alias{NUsetup}
\title{Specify KEY parameters used in NU.Learning to "design" analyses of Observational Data.
}
\description{Invoke NUsetup() to specify the name of the Hierarchical Clustering object
  output by NUcluster() and the name of the data.frame containing all desired X-covariates, the
  Treatment/Exposure variable and the Y-Outcome variable. It is ESSENTIAL to save the Environment
  output by NUsetup() as a named object within the user's .GlobalEnv space.
}
\usage{NUsetup(hclobj, dframe, trex, yvar)}
\arguments{
  \item{hclobj}{Name of a NUcluster() output object created using a cluster::diana or
     stats::hclust method.}
  \item{dframe}{Name of the data.frame containing all X-covariates, the Treatment/Exposure
     variable and the Y-Outcome variable.}
  \item{trex}{Name of the numerical Treatment/Exposure variable.}
  \item{yvar}{Name of the numerical Y-Outcome variable.}
}
\value{The environment output by NUsetup() must be saved to the user's .GlobalEnv space. It's
  contents will be automatically updated by calls to other NU.Learning functions:
  \item{aggdf}{data.frame with 4 columns and 1 row for each call to ltdagg() or lrcagg().}
  \item{aggdf$Label}{Factor value of "LTD" or "LRC".}
  \item{aggdf$Blocks}{K = integer Number of Clusters requested.}
  \item{aggdf$LTDmean or aggdf$LRCmean}{numerical value of cluster mean of LTD or LRC estimates.}
  \item{aggdf$LTDstde or aggdf$LRCstde}{numerical value of the within-cluster standard deviation.} 
  \item{boxdf}{data.frame of 2 variables ...for input to boxplot() by NUcompare().}
  \item{boxdf$NUstat}{LTD or LRC estimate for a single experimental unit from ltdagg() or lrcagg().}
  \item{boxdf$K}{Number of Cluters used in forming the LTD or LRC estimate for each Experimental Unit.}
  \item{Kmax}{Maximum Number of Clusters so that Average Size will be >= 12 experimental units.}
  \item{LTDmax or LRCmax}{Maximum Treatment Effect-Size estimate across Clusters.}  
  \item{LTDmin or LRCmin}{Minimum Treatment Effect-Size estimate across Clusters.}
  \item{NumLevels}{Integer number of distinct Levels of the Treatment/Exposure variable: trex.}
  \item{pars}{Character data.frame with 4 columns and 1 row.}   
  \item{pars[1,1]}{Name of the diana or hclust object created by NUcluster().}  
  \item{pars[1,2]}{Name of data.frame containing the X, Treatment/Exposure and Y variables.}
  \item{pars[1,3]}{Name of Treatment/Exposure variable within data.frame pars[1,2].}
  \item{pars[1,4]}{Name of Y-outcome variable within data.frame pars[1,2].}
  }
\references{
  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192. 

  Obenchain RL. (2023) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ltdagg}}, \code{\link{ivadj}} and \code{\link{lrcagg}}.}
\examples{
  \donttest{
  # Running takes about 7 seconds...
  data(pci15k)
  xvars  = c("stent", "height", "female", "diabetic", "acutemi", "ejfract", "ves1proc")
  hclobj = NUcluster(pci15k, xvars)
  NUe    = NUsetup(hclobj, pci15k, thin, surv6mo)
  ls.str(NUe)
  }
}
\keyword{design}