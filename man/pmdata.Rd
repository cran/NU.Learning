\name{pmdata}
\alias{pmdata}
\docType{data}
\title{Particulate Matter, Mortality and Other data for 2980 US Counties} 
\description{This data.frame combines 122 variables from the 5 sources referenced below. Several PM
  variables appear to be predictions from EPA \dQuote{CMAQ} models rather than values from validated
  measuring instruments. NU.Learning concepts are illustrated in demo(pmdata) using Clustering of 2973
  Counties and Parishes within the contiguous 48 US States and Washington, D.C.}
\usage{data(pmdata)}
\format{
  This data.frame contains 122 variables for 2,980 US counties. A total of 738 "NA"s imply that only about
  two tenths of one percent of these 363,560 values are missing.	  
  \describe{
    \item{fips}{Federal Information Processing Standard code; 4 or 5 digits; 2980 unique values}
    \item{C50}{Cluster ID Number between 1 and 50. Total of 50 unique values}
    \item{LRC50}{Local (Spearman) Rank Correlation between Bvoc and AACRmort within Cluster}
    \item{County}{County or Parish name is a Factor variable (character code)}
    \item{State}{State name is a 2-Character Factor code; 49 unique levels}
    \item{Deaths}{CDC: Total number of Deaths in the County in 2016}
    \item{Population}{CDC: Total population of County in 2016}
    \item{CRmort}{CDC: Crude Rate of Circulatory-Respiratory Mortality for the County in 2016}
    \item{CrudeL95}{CDC: Lower 95\% confidence limit for Crude Rate of CR Mortality}
    \item{CrudeU95}{CDC: Upper 95\% confidence limit for Crude Rate of CR Mortality}
    \item{CrudeSE}{CDC: Standard Error for Crude Rate of CR Mortality}
    \item{AACRmort}{CDC: Age Adjusted Rate of Circulatory-Respiratory Mortality in 2016}
    \item{AACRL95}{CDC: Lower 95\% confidence limit for Age Adjusted Rate of CR Mortality}
    \item{AACRU95}{CDC: Upper 95\% confidence limit for Age Adjusted Rate of CR Mortality}
    \item{AACR_SE}{CDC: Standard Error for Age Adjusted Rate of CR Mortality}
    \item{TotDeathPct}{CDC: Total Death Percentage - Circulatory-Respiratory Mortality in 2016}
    \item{lat}{EPA: County Latitude used in EPA "CMAQ" model calculations}
    \item{long}{EPA: County Longitude used in EPA "CMAQ" model calculations}
    \item{RHpct}{EPA: Relative Humidity Percentage for 2016}
    \item{SFCtmpC}{EPA: Surface Temperature in Degrees Centigrade for 2016}
    \item{NO2}{EPA: Nitrogen Dioxide level (NO2.ppbV) for 2016}
    \item{O3}{EPA: Ozone level (O3.ppbV) for 2016}   
    \item{pmCL}{EPA: Chlorine level in Particulate Matter (PM25_CL.ugm3) for 2016}
    \item{pmEC}{EPA: Ethylene Carbonate level in Particulate Matter (PM25_EC.ugm3) for 2016}
    \item{pmNA}{EPA: Sodium level in Particulate Matter (PM25_NA.ugm3) for 2016}
    \item{pmMG}{EPA: Magnesium level in Particulate Matter (PM25_MG.ugm3) for 2016}
    \item{pmK}{EPA: Potassium level in Particulate Matter (PM25_K.ugm3) for 2016}
    \item{pmCA}{EPA: Calcium level in Particulate Matter (PM25_CA.ugm3) for 2016}
    \item{pmNH4}{EPA: Ammonium level in Particulate Matter (PM25_NH4.ugm3) for 2016}
    \item{pmNO3}{EPA: Nitrate level in Particulate Matter (PM25_NO3.ugm3) for 2016}
    \item{pmOC}{EPA: Organic Compounds in pmTOT [fine particulate matter] (PM25_OC.ugm3) for 2016}
    \item{pmOM}{EPA: OM compounds in pmTOT (PM25_OM.ugm3) for 2016}
    \item{pmOTHR}{EPA: Other Compounds in pmTOT (PM25_OTHR.ugm3) for 2016}
    \item{pmSO4}{EPA: Sulfate Compounds in pmTOT (PM25_SO4.ugm3) for 2016}
    \item{pmFE}{EPA: Ferrous Compounds in pmTOT (PM25_FE.ugm3) for 2016}
    \item{pmSI}{EPA: Silicon Compounds in pmTOT (PM25_SI.ugm3) for 2016}
    \item{pmTI}{EPA: Titanium Compounds in pmTOT (PM25_TI.ugm3) for 2016}
    \item{pmMN}{EPA: Manganese Compounds in pmTOT (PM25_MN.ugm3) for 2016}
    \item{pmAL}{EPA: Aluminum Compounds in pmTOT (PM25_AL.ugm3) for 2016}
	\item{pmUNSPCRS}{EPA: UNSPCRS Compounds in pmTOT (PM25_UNSPCRS.ugm3) for 2016}
    \item{pmPOA}{EPA: Primary Organic Aerosols in pmTOT (PM25_POA.ugm3) for 2016}
    \item{pmSOA}{EPA: Secondary Organic Aerosols in pmTOT (PM25_SOA.ugm3) for 2016; pmSOA = Avoc+Bvoc}
    \item{pmGLY}{EPA: Glycemic Secondary Organic Aerosols in pmTOT (PM25_GLYSOA.ugm3) for 2016}
    \item{pmOLGB}{EPA: OLGB compounds in pmTOT (PM25_OLGB.ugm3) for 2016}
    \item{pmISOP}{EPA: ISOP compounds in pmTOT (PM25_ISOP.ugm3) for 2016}
    \item{pmEPOX}{EPA: EPOX compounds in pmTOT (PM25_EPOX.ugm3) for 2016}
    \item{pmSQT}{EPA: SQT compounds in pmTOT (PM25_SQT.ugm3) for 2016}
    \item{pmMTN}{EPA: MTN compounds in pmTOT (PM25_MTN.ugm3) for 2016}
    \item{pmMT}{EPA: MT compounds in pmTOT (PM25_MT.ugm3) for 2016}
    \item{pmTOT}{EPA: Total (fine) Particulate Matter (PM25_TOT.ugm3) for 2016}
    \item{SFCtmpK}{EPA: Surface Temperature in Degrees Kelvin for 2016}
    \item{CardioRes}{CDC: Cardio Respiratory Rate (rate I00J98.per100000.cdc) for 2016}
    \item{POPcdc}{CDC: County Population (population.cdc) for 2016}
    \item{POP5yracs}{CDC: 5yracs Population (population.people.5yracs) for 2016}  
    \item{PREMdeath}{Premature Deaths per 100K Residents ...UWPHI for 2018}
    \item{POFHealth}{Poor or Fair Health rate (Poor.or.fair.health) ...UWPHI for 2018}
    \item{PPHdays}{Poor Physical Health days (Poor.physical.health.days) ...UWPHI for 2018} 		
    \item{PMHdays}{Poor Mental Health days (Poor.mental.health.days) ...UWPHI for 2018} 		   		
    \item{LBW}{Low Birth Weight rate (Low.birthweight) ...UWPHI for 2018} 	
    \item{ASmoke}{Adult Smoking Percentage (Adult.smoking) ...UWPHI for 2018} 	    			
    \item{AObes}{Adult Obesity Percentage (Adult.obesity) ...UWPHI for 2018} 	      			
    \item{FEnv}{Food Environment Index (Food.environment.index) ...UWPHI for 2018} 
    \item{PhysInAct}{Physical Inactivity (Physical.inactivity) ...UWPHI for 2018} 
    \item{ExercOPS}{Access to Exercise Opportunities ...UWPHI for 2018} 
    \item{ExsDrink}{Excessive Drinking Rate (Excessive.drinking) ...UWPHI for 2018} 
    \item{AIDrivD}{Alcohol Impaired Driving Deaths ...UWPHI for 2018} 
    \item{STInfect}{Sexually Transmitted Infections ...UWPHI for 2018} 
    \item{TBirths}{Teenage Births (Teen.births) ...UWPHI for 2018} 
    \item{Uninsur}{Uninsured Residences (Uninsured) ...UWPHI for 2018}
    \item{PCDocs}{Primary Care Physicians (Primary.care.physicians) ...UWPHI for 2018}
    \item{Dentists}{Dentists (Dentists) ...UWPHI for 2018}
    \item{PrevntHS}{Preventable Hospital Stays (Preventable.hospital.stays) ...UWPHI for 2018}
    \item{DiabMNT}{Diabetes Monitoring (Diabetes.monitoring) ...UWPHI for 2018}
    \item{MammoSC}{Mammography Screening (Mammography.screening) ...UWPHI for 2018}
    \item{SomCOL}{Some College Education (Some.college) ...UWPHI for 2018}   		
    \item{UnEMP}{Unemployment Rate (Unemployment) ...UWPHI for 2018}    		
    \item{ChildPOV}{Children Living in Poverty (Children.in.poverty) ...UWPHI for 2018} 
    \item{IncomIEQ}{Income Inequality (Income.inequality) ...UWPHI for 2018} 
    \item{ChildSPH}{Children In Single-Parent Households ...UWPHI for 2018} 
    \item{SocASOC}{Social Associations (Social.associations) ...UWPHI for 2018} 
    \item{VioCRM}{Violent Crime Rate (Violent.crime) ...UWPHI for 2018} 
    \item{InjyDths}{Injury Death Rate (Injury.deaths) ...UWPHI for 2018} 
    \item{AirPolpm}{Air Pollution Particulate Matter ...UWPHI for 2018} 
    \item{DrnkWtVi}{Drinking Water Violations (Drinking.water.violations) ...UWPHI for 2018} 
    \item{SevrHOUS}{Severe Housing Problems (Severe.housing.problems) ...UWPHI for 2018} 
    \item{DrivATW}{Driving Alone to Work (Driving.alone.to.work) ...UWPHI for 2018} 
    \item{LComutA}{Long Commute - Driving Alone to Work ...UWPHI for 2018} 
    \item{PAAM}{Premature Age Adjusted Mortality ...UWPHI for 2018} 
    \item{FrqPhysD}{Frequent Physical Distress ...UWPHI for 2018} 
    \item{FrqMentD}{Frequent Mental Distress ...UWPHI for 2018} 
    \item{DiabPrev}{Diabetes Prevalence (Diabetes.prevalence) ...UWPHI for 2018}         
    \item{FoodInSec}{Food Insecurity (Food.insecurity) ...UWPHI for 2018}    
    \item{LimAHFood}{Limited Access to Healthy Foods ...UWPHI for 2018}  
    \item{DrugOdDM}{Drug Overdose Deaths Model predictions ...UWPHI for 2018}  
    \item{InsufSlp}{Insufficient Sleep (Insufficient.sleep) ...UWPHI for 2018} 
    \item{UnInsAds}{Uninsured Adults (Uninsured.adults) ...UWPHI for 2018} 
    \item{UnInsCls}{Uninsured Children (Uninsured.children) ...UWPHI for 2018} 
    \item{HCareCost}{Health Care Costs (Health.care.costs) ...UWPHI for 2018} 
    \item{OthPrimCP}{Other Primary Care Providers ...UWPHI for 2018} 
    \item{MHHIncome}{Median Household Income ...UWPHI for 2018} 
    \item{ChildFRPL}{Children Eligible for Free or Reduced-Price Lunch ...UWPHI for 2018} 
    \item{Population}{County Population (Population) ...UWPHI for 2018} 
    \item{AGELess18}{Residents below 18 Years of Age ...UWPHI for 2018} 
    \item{A65oOVR}{Residents 65 or Older (X..65.and.older) ...UWPHI for 2018} 
    \item{NHispAfA}{Non-Hispanic African-American Residents ...UWPHI for 2018}
    \item{AmINalsN}{American Indian or Alaskan Natives ...UWPHI for 2018}
    \item{Asian}{Asian Residents (X..Asian) ...UWPHI for 2018}
    \item{NHawOPI}{Native Hawaiian and Other Pacific Islanders ...UWPHI for 2018}
    \item{Hispanic}{Hispanic Residents (X..Hispanic) ...UWPHI for 2018}
    \item{NHispWht}{Non-Hispanic White Residents ...UWPHI for 2018}
    \item{LoProEngl}{Low Proficiency in English (not.proficient.in.English) ...UWPHI for 2018}
    \item{Females}{Female Residents ...UWPHI for 2018}
    \item{Rural}{Rural Residents ...UWPHI for 2018}           
    \item{pmOA}{EPA: Organic Aerosols in pmTOT (PM25_OA.ugm3) for 2016}  
    \item{Avoc}{EPA: Anthroprogenic [man-made] Volatile Organic Compounds in pmTOT for 2016}
    \item{pmSEAspry}{EPA: Sea Spray components in pmTOT (PM25_SOAAVOC.ugm3) for 2016}
    \item{pmDUST}{EPA: Dust components in pmTOT (PM25_DUST.ugm3) for 2016}
    \item{pmNH4NO3}{EPA: Ammonium Nitrate components in pmTOT (PM25_NH4NO3.ugm3) for 2016}
    \item{pmSOOT}{EPA: Soot components in pmTOT (PM25_SOOT.ugm3) for 2016}        
    \item{isop}{EPA: SOA Isoprenes (PM25_SOAISOPRENE.ugm3) for 2016}         
    \item{terp}{EPA: SOA Terpenes (PM25_SOATERPENE.ugm3) for 2016}            
    \item{Bvoc}{EPA: Biogenic (natural) Volatile Organic Compounds for 2016; Bvoc = isop + terp}
  }
}
\references{
  Obenchain RL. and Young SS. (2022), EPA Particulate Matter Data - Analyses using Local Control
  Strategy. (24 pages, 22 figures) https://doi.org/10.48550/arXiv.2209.05461

  Pye, H., Ward-Caviness, C., Murphy, B., Appel, K., and Seltzer, K. (2021). Secondary organic
  aerosol association with cardiorespiratory disease mortality in the united states. Nature
  Communications, 12.7215  https://doi.org/10.1038/s41467-021-27484-1
  
  Pye, H. [EPA] (2021), Data For Secondary Organic Aerosol and Cardiorespiratory Disease
  Mortality. https://doi.org/10.5281/zenodo.5713903 
  
  University of Wisconsin, Population Health Institute. https://uwphi.pophealth.wisc.edu
  [UWPHI] UWPHI@med.wisc.edu
  
  Young SS. and Obenchain RL. (2022), "EPA particulate matter data...Analyses
  using Local Control Strategy"  https://doi.org/10.5061/dryad.63xsj3v58 
}
\examples{
    data(pmdata)
    str(pmdata, list.len=122)
}
\keyword{datasets}
