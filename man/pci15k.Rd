\name{pci15k}
\alias{pci15k}
\docType{data}
\title{Six-month Survival, Cardiac cost and Baseline Covariate data for 15,487 PCI patients.}
\description{Using observational data on 996 patients who received a Percutaneous Coronary
  Intervention (PCI) at Ohio Heart Health, Lindner Center, Christ Hospital, Cincinnati
  (Kereiakes et al, 2000), we generated this much larger dataset via "plasmode simulation."}
\usage{data(pci15k)}
\format{
  A data frame of 11 variables on 15,487 patients; no NAs.
  \describe{
    \item{patid}{Patient ID number: 1 to 15487.}
    \item{surv6mo}{Binary PCI Survival variable: 1 => Survival for at least 6 months following
       PCI, 0 => Survival for less than 6 months.}
    \item{cardcost}{Cardiac related costs incurred within 6 months of patient's initial PCI;
       numeric value in 1998 dollars; costs were truncated by death for the 404 patients with
       surv6mo == 0.}
    \item{thin}{Numeric treatment selection indicator: thin = 0 implies usual PCI care alone;
       thin = 1 implies usual PCI care augmented by either planned or rescue treatment with a
       new blood thinning agent.}
    \item{stent}{Coronary stent deployment; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{height}{Height in centimeters; numeric integer from 133 to 198.}
    \item{female}{Female gender; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{diabetic}{Diabetes mellitus diagnosis; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{acutemi}{Acute myocardial infarction within the previous 7 days; numeric, with 1
       meaning YES and 0 meaning NO.}
    \item{ejfract}{Left ejection fraction; numeric value from 17 percent to 77 percent.}
    \item{ves1proc}{Number of vessels involved in the patient's initial PCI procedure; numeric
       integer from 0 to 5.}
  }
}
\references{
  Kereiakes DJ, Obenchain RL, Barber BL, et al.  Abciximab provides cost effective survival
  advantage in high volume interventional practice. \emph{Am Heart J} 2000; \bold{140}: 603-610.
  
  Gadbury GL, Xiang Q, Yang L, Barnes S, Page GP, Allison DB. Evaluating Statistical Methods Using
  Plasmode Data Sets in the Age of Massive Public Databases: An Illustration Using False Discovery
  Rates. \emph{PLOS Genetics} 2008; \bold{4}: 1-8, e1000098 (Open Access). 

  Obenchain RL. (2023) \bold{NU.Learning_in_R.pdf} http://localcontrolstatistics.org 
}
\examples{
    data(pci15k)
    str(pci15k)
}
\keyword{datasets}
