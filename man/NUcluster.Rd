\name{NUcluster}
\alias{NUcluster}
\title{Hierarchical Clustering of experimental units (such as patients) in X-covariate Space}
\description{Form the full, hierarchical clustering tree (dendrogram) for all units (regardless
  of Treatment/Exposure status) using Mahalonobis distances computed from specified baseline
  X-covariate characteristics.}
\usage{
  NUcluster(dframe, xvars, method="ward.D")
}
\arguments{
  \item{dframe}{Name of data.frame containing baseline X covariates.}
  \item{xvars}{List of names of X variable(s).}
  \item{method}{Hierarchical Clustering Method of "diana", "ward.D", "ward.D2", "complete",
        "average", "mcquitty", "median" or "centroid".}
}
\details{The first step in applying NU.Learning to data is to hierarchically cluster
  experimental units in baseline X-covariate space ...thereby creating "Blocks" of relatively
  well-matched units.  NUcluster first calls stats::prcomp() to calculate Mahalanobis distances
  using standardized and orthogonal Principal Coordinates. NUcluster then uses either the
  divisive cluster::diana() method or one of seven agglomerative methods from stats::hclust()
  to compute a dendrogram tree. The hclust function is based on Fortran code contributed to
  STATLIB by F. Murtagh.
}
\value{ An output list object of class NUcluster, derived from cluster::diana or stats::hclust.
  \item{dframe}{Name of data.frame containing all baseline X-covariates.}
  \item{xvars}{List of 1 or more X-variable names.}
  \item{method}{Hierarchical Clustering Method: "diana", "ward.D", "ward.D2", "complete",
        "average", "mcquitty", "median" or "centroid".}
  \item{hclobj}{Hierarchical clustering object created by the designated method.}
}
\references{
  Kaufman L, Rousseeuw PJ.  (1990) \bold{Finding Groups in Data.  An Introduction to
  Cluster Analysis}.  New York: John Wiley and Sons.

  Kereiakes DJ, Obenchain RL, Barber BL, et al. (2000) Abciximab provides
  cost effective survival advantage in high volume interventional practice.
  \emph{Am Heart J} \bold{140}: 603-610.
  
  Murtagh F. (1985) Multidimensional Clustering Algorithms. \bold{COMPSTAT Lectures} 4.

  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192.
  
  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293-298.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{NUsetup}}, \code{\link{ltdagg}} and \code{\link{lrcagg}}.}
\examples{
data(radon)
xvars  = c("obesity", "over65", "cursmoke")
hclobj = NUcluster(radon, xvars)  # ...using default method = "ward.D"
plot(hclobj)  
}
\keyword{cluster}
\keyword{design}