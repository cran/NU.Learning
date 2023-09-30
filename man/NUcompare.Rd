\name{NUcompare}
\alias{NUcompare}
\title{Display NU Sensitivity Graphic for help in choice of K = Number of Clusters}
\description{This function displays Box-Whisker diagrams that compare Treatment Effect-Size
  distributions for different values of K = Number of Clusters requested in X-covariate space.
  After an initial call to NUsetup(), the analyst typically makes multiple calls to either
  ltdagg() or lrcagg() for different values of K.  The analyst then invokes NUcompare() to
  see how choice of K changes the location, spread and/or skewness of the distribution of
  Treatment Effect-Size estimates across Clusters. Variance-Bias trade-offs occur as K
  increases; large values of K may reduce Bias, but they definitely inflate the Variance of
  LTD and LRC distributions.}

\usage{
  NUcompare(envir)
}
\arguments{
  \item{envir}{R environment output by an earlier call to NUsetup().}
}
\details{The third phase of NU.Learning is called EXPLORE and uses graphical
  Sensitivity Analyses to show how Treatment Effect-Size distributions change with choice
  of NU parameter settings. Choice of K = Number of Clusters requested is guided, primarily,
  by NUcompare() graphics. Equally important are the analyst's choices of (i) which [and how
  many] of the available baseline X-covariates to "adjust for" and (ii) which clustering
  algorithm and dissimilarity metric to use. Unfortunately, changing these latter choices
  requires the analyst to essentially "start over" ...i.e. invoking NUcluster() with changed
  arguments, followed by an invocation of NUsetup() with a different 1st argument. To change
  only one's choice of y-Outcome variable and/or the Treatment/Exposure variable, a new
  NUsetup() invocation is all that is needed.
}
\value{ NULL }
\references{
  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192.
   
  Obenchain RL. (2015) \bold{NU_Confirm_Guidelines.pdf} http://localcontrolstatistics.org 

  Obenchain RL. (2023) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293-298.
  
  Tukey JW. (1977) \bold{Exploratory Data Analysis}, \emph{New York: Addison-Wesley}, Section 2C.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ltdagg}}, \code{\link{ivadj}} and \code{\link{lrcagg}}.}
\examples{
  \donttest{
  # Running takes more than 7 seconds...
  data(pci15k)
  xvars   = c("stent", "height", "female", "diabetic", "acutemi", "ejfract", "ves1proc")
  hclobj  = NUcluster(pci15k, xvars)
  NU.env  = NUsetup(hclobj, pci15k, thin, surv6mo)
  surv050 = ltdagg( 50, NU.env)
  surv100 = ltdagg(100, NU.env)
  surv200 = ltdagg(200, NU.env)
  NUcompare(NU.env)
  }
}
