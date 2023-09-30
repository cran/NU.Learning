library(NU.Learning)

# Input plasmode data on Blood Thinner use in Percutaneous Coronary Interventions (PCIs).
# This simulation is based upon data collected at the Lindner Center: Kereiakes et al. (2000).
data(pci15k)
# outcomes:   surv6mo (binary) or cardcost (continuous)
# treatment:  thin  ...is binary

# The following Calculations take about 8 seconds...
# Define Cluster Hierarchy for UNSUPERVISED, nonparametric analyses...
xvars  = c("stent","height","female","diabetic","acutemi","ejfract","ves1proc")
system.time( hclobj <- NUcluster(pci15k, xvars) )     # default clustering method = "ward.D"       
hclobj                                           
plot(hclobj)

# Save NU.Learning basic parameter settings to an environment that will be updated...
NUe = NUsetup(hclobj, pci15k, thin, surv6mo)
ls.str(NUe)

# Compute and Save LTD distributions for several values of K = number of Clusters...
surv0050 = ltdagg(50, NUe)
surv0100 = ltdagg(100, NUe)
surv0200 = ltdagg(200, NUe)
surv0500 = ltdagg(500, NUe)	           # average cluster size ~31 patients 
plot(surv0500, show="ecdf", NUe)
surv0750 = ltdagg(750, NUe)
surv1000 = ltdagg(1000, NUe)	       # average cluster size ~16 patients  
plot(surv1000, show="ecdf", NUe) 

# "Sensitivity Analysis" Summary...
NUcompare(NUe)
# NOTE: LTD Distribution for 500 Clusters appears to Optimize Variance-Bias Trade-Offs... 

# Save and plot Instrumental Variable LAOs for 2 large values of K = number of Clusters...
iv0500 = ivadj(surv0500)
plot(iv0500) 
iv1K = ivadj(surv1000)
plot(iv1K)

# The following Calculations take about 32 seconds...
# Confirm: Does the Observed LRC distribution for 500 clusters truly differ from
# the Random Permutation NULL distribution assuming x-Covariates are Ignorable?
system.time( conf5H <- confirm(surv0500) ) 
conf5H
plot(conf5H)

# The following Calculations take about 52 seconds...
# Simulate crude pmax.value for Kolmogorov-Smirnov D-statistic...
system.time( ksd5H <- KSperm(conf5H) )     
ksd5H
plot(ksd5H)

# Example: "Most-Like-Me" Visualizations for pci15k Patient No. 11870...
# xvars:       "stent","height","female","diabetic","acutemi","ejfract","ves1proc"
xvec11870 = c( 0,      162,     1,       1,         0,        57,       1) 	 
mlme11870 = mlme(NUe, hclobj, surv0500, xvec11870 )
plot(mlme11870, NN = 250) 
# Summary Statistics for 3 values of NN...
mlme.stats(mlme11870, NN = c( 250, 500, 1000 ))

### End of demo(pci15k) ##############################################
