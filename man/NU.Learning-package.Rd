\name{NU.Learning-package}
\alias{NU.Learning-package}
\docType{package}
\title{NU.Learning: Nonparametric and Unsupervised Adjustment for Bias and Confounding}
\description{NU.Learning forms Local Treatment Differences (LTDs) or Local Rank Correlations
  (LRCs) within Clusters of experimental units (patients, etc.) who have been relatively well-matched
  on their baseline X-confounder characteristics.  The resulting distribution of LTD/LRC effect-size
  estimates can be interpreted much like a Bayesian posterior. Yet these distributions have been
  formed, via Nonparametric and Unsupervised Preprocessing, in purely Objective Ways.
}
\details{
  \tabular{ll}{
  Package: \tab NU.Learning\cr
  Type: \tab Package\cr
  Version: \tab 1.5\cr
  Date: \tab 2023-09-15\cr
  License: \tab GPL-2\cr
  }

  UNSUPERVISED LOCAL TREATMENT DIFFERENCES or LOCAL RANK CORRELATIONS:
  
  Multiple calls to ltdagg(K) or lrcagg(K) for varying numbers of clusters, K, are typically
  made after first invoking NUcluster() to hierarchically cluster patients in X-space and
  invoking NUsetup() to specify a numeric y-Outcome variable and a numeric treatment choice
  or exposure level measure, trex. 
  
  UNSUPERVISED INSTRUMENTAL VARIABLES = LOCAL AVERAGE y-OUTCOME EFFECTS:
    
  An OBSERVED Propensity Score (PS) is defined here to be either (i) the local (within-cluster)
  fraction of experimental units (patients) receiving trex==1 (new) rather than trex==0 (control)
  or else (ii) a measure of "relative exposure" when the numeric trex measure has (many) more
  than 2 observed levels. Multiple calls to ivadj(K) for varying numbers of clusters, K, then
  yield alternative scatters of Local Average Outcomes (LAOs) for Clusters when plotted against
  their PS estimates and, thus, different possible linear fits or smooth.splines() yielding
  potentially different inferences about across-cluster Treatment or Exposure Effects.
  
  CONFIRMATION and SENSITIVITY ANALYSES of LOCAL EFFECT-SIZE DISTRIBUTIONS:
  
  For a given value of K = Number of Clusters requested, the output object from ltdagg(K) or
  lrcagg(K) can be input to confirm() to use (nonparametric) permutation theory to display 
  visual evidence (empirical CDF comparisons) concerning the Question: "Does x-matching Truly
  Matter?" The NULL hypothesis here is that the x-Covariates used in Clustering / Matching of
  Experimental Units are actually IGNORABLE. Evidence against this hypothesis is provided when
  the observed LOCAL Effect-Size Distribution clearly deviates from the purely RANDOM, NULL
  distribution computed (to any desired precision) by randomly PERMUTING cluster ID labels
  across experimental units. Furthermore, the statistical significance of differences between
  the observed and random NULL distributions can be estimated using KSperm(), which simulates
  the random permutation distribution of the Kolmogorov-Smirnov D-statistic when many tied
  values occur in both distributions being compared. Finally, the NUcompare() function helps
  users of NU.Learning decide which Number of Clusters, K, optimizes Variance-Bias trade-offs.
  Larger values of K tend to yield smaller clusters with better matches and, thus, potentially
  reduced BIAS. On the other hand, smaller values of K usually yield local effect-size
  estimates with much lower Variability (higher Precision).
  
  "Most-Like-Me" HISTOGRAMS for DOCTOR-PATIENT discussions of PERSONALIZED MEDICINE:
  
  For a specified vector, xvec, of numerical values of the X-confounder variables used in the
  current CLUSTERING of eUnits, display histograms of observed LTD or LRC effect-sizes for
  (i) all available patients and (ii) for the specified number, NN, of "Nearest-Neighbors" in
  X-confounder space of the TARGET eUnit ...i.e. xvec defines "Me".
}
\author{Bob Obenchain <wizbob@att.net>}
\references{ 
  McClellan M, McNeil BJ, Newhouse JP. (1994) Does More Intensive Treatment of
  Myocardial Infarction in the Elderly Reduce Mortality?: Analysis Using Instrumental
  Variables. \emph{JAMA} \bold{272}: 859-866.  
  
  Obenchain RL. (2010) The Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192. 
 
  Obenchain RL, Young SS. (2013) Advancing Statistical Thinking in Observational Health
  Care Research. \emph{J. Stat. Theory and Practice}, \bold{7}: 456-469,
  \doi{10.1080/15598608.2013.772821}.

  Lopiano KK, Obenchain RL, Young SS. (2014) Fair treatment comparisons in observational
  research. \emph{Statistical Analysis and Data Mining}, \bold{7}: 376-384,
  \doi{10.1002/sam.11235}.

  Obenchain RL. NU.Learning-vignette. (2023) \bold{NU.Learning_in_R.pdf}
  http://localcontrolstatistics.org 

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41-55.

  Rosenbaum PR, Rubin RB. (1984) Reducing Bias in Observational Studies Using
  Subclassification on a Propensity Score. \emph{JASA} \bold{79}: 516-524.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293-298.
  
  Stuart EA. (2010) Matching Methods for Causal Inference: A Review and a Look Forward.
  \emph{Statistical Science} \bold{25}: 1-21.
}
\keyword{package}
