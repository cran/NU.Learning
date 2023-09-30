\name{radon}
\alias{radon}
\docType{data}
\title{Radon exposure and lung cancer mortality data for 2,881 US counties in 46 States.}
\description{Federal EPA and state government agencies have been reporting observational
  data at the US County level since about 1980. The data given here include 5 potential
  X-confounder variables of the relationship between lung cancer mortality and radon
  exposure; they were amassed and checked by Goran Krstic, Fraser Health Authority,
  Vancouver, BC, Canada.}
\usage{data(radon)}
\format{
  A data frame of 11 variables for 2881 US counties. One Missing Value; row 778 for
  Shannon County, SD, fips == 46113, has hhincome == NA.	  
  \describe{
    \item{fips}{County FIPS code. Codes are 4 or 5 digit integers; 2881 unique values.}
    \item{state}{State Factor variable (2-character codes); 46 unique levels.}
    \item{county}{County or Parish Factor variable (character codes); 1703 unique levels.}
    \item{lcanmort}{Lung Cancer Mortality rate (deaths per 100,000 person-years), 1980-2004.}
    \item{radon}{County Radon Exposure level in picocuries per liter (pCi/L) for some
        unspecified period within 1986-1992; rounded to nearest single decimal place.}
    \item{lnradon}{Natural logarithm of County Radon Exposure level. Radon levels reported as
	    0.0 for 10 US counties are Windsorized here to ln(0.05), which is roughly -3.}
    \item{obesity}{Percentage of County Residents considered Obese (age adjusted), 2008.}
    \item{over65}{Percentage of County Residents of Age 65 and over, 2000 Census.}
    \item{cursmoke}{Percentage of County Residents who Currently Smoke, 1997-2003.}
    \item{evrsmoke}{Percentage of County Residents who Ever Smoked, 1997-2003.}
	\item{hhincome}{Average Median HouseHold Income in Thousands ($1,000), 1989-2004.}
  }
}
\references{
  Krstic G, Obenchain RL. (2016) Radon dataset documentation and downloads.
  http://localcontrolstatistics.org 

  Obenchain RL. (2018) \bold{RADON_short.pdf} http://localcontrolstatistics.org  40 PPT Slides and
  Commentary in Notes Pages format.
}
\examples{
    data(radon)
    str(radon)
}
\keyword{datasets}
