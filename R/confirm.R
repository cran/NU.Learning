"plot.confirm" <- 
function(x, ...)
{
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    nx <- x$nclus
    par(mfrow = c(1,1))
    plot(ecdf(x$dfconf$lstat), verticals=TRUE, do.points=FALSE, ann=FALSE,
        col="gray70", lwd=2, lty=1)  
    lines(ecdf(x$NUdist$lstat), verticals=TRUE, do.points=FALSE, ann=FALSE,
        col="blue2", lwd=2)
    abline(v=0, lty="dashed", col="black")
    if (x$Type == 1)
        title(main = paste("NU.Learning Confirm eCDF Comparison for", nx, "Clusters"), 
            xlab = "Within-Cluster Local Treatment Difference (LTD)", 
            ylab = "Cumulative Probability")
    else
        title(main = paste("NU.Learning Confirm eCDF Comparison for", nx, "Clusters"), 
            xlab = "Within-Cluster Local Rank Correlation (LRC)", 
            ylab = "Cumulative Probability")
}

"print.confirm" <-
function(x, ...)
{
    cat("\nconfirm Object: Compare Observed and NULL Distributions of Local Effect-Sizes...\n")
    cat("   Simulated NULL Distribution uses Random Clusterings of Experimental Units.\n")
    cat("\nData Frame:", x$dframe, "\n")
    cat("Outcome Variable:", x$yvar, "\n")
    cat("Treatment Factor:", x$trtm, "\n")
    cat("Number of Replications:", x$reps, "\n")
    cat("Number of Clusters per Replication:", x$nclus, "\n")
    cat("Number of Random NULL Local Effect-Sizes:", length(x$dfconf[,1]), "\n" )
    cat("\n    Mean Observed Local Effect-Size =", x$NUmean)
    cat("\n    Std. Dev. of Observed Effect-Sizes =", x$NUstde)
    cat("\n    Mean Random NULL Effect-Size =", x$RPmean)
    cat("\n    Std. Dev. of Random Effect-Sizes =", x$RPstde, "\n\n")
    cat("Nonstandard Kolmogorov-Smirnov comparison of Discrete Distributions:\n")
    cat("Observed two-sample KS D-statistic =", x$KSobsD$statistic, "\n\n")
    }

"confirm" <-
function(x, reps=100, seed=12345)
{
    # Compute Random Permutation Distribution of Within-Cluster LTD/LRC estimates. Include
    # less or least relevant comparisons as well as those neutral, more or most relevant. 
    if (missing(x) || (!inherits(x, "ltdagg") && !inherits(x, "lrcagg")))  
        stop("First argument to confirm() must be a ltdagg() or lrcagg() object.")
    if (inherits(x, "lrcagg")) {
        type = 2
	    NUdist <- x$LRCdist
        NUmean <- x$LRCmean
        NUstde <- x$LRCstde
    } else {
        type = 1
        NUdist <- x$LTDdist
        NUmean <- x$LTDmean
        NUstde <- x$LTDstde
    }
    nclus = x$actclust
    names(NUdist) <- c("c", "ID", "y", "t", "lstat")
    units = length(NUdist[,5])
    seed = as.integer(seed)
    olist <- list(hiclus = x$hiclus, dframe = x$dframe, trtm = x$trtm, yvar = x$yvar,
        reps=reps, seed=seed, nclus=nclus, units=units)
    runif(1)		
    set.seed(seed)   # Set seed for Monte Carlo pseudo random sequence...
    for(i in 1:reps) { 
        cperm <- NUdist$c[order(as.vector(rnorm(units)))]  # Permuted Cluster IDs
        pdf <- as.data.frame(cbind(cperm, NUdist[,3:4]))   # y & t columns NOT Permuted
	    names(pdf) <- c("c","y","t")
        if (type == 2) {
            dfLRC <- do.call( rbind, lapply( split(pdf, pdf$c),
    	        function(x) {
    	            x <- na.omit(x)
                    if (length(x[,1]) < 3)
                        LRC = NA
                    else
                        LRC = round( cor(x$y, x$t, method = "spearman"), 8)
                    data.frame(c=x$c[1], LRC=LRC, w=length(x$t))
                } )
            )
            pdf = merge(pdf, dfLRC[,1:2], by.x="c", by.y="c")
        }
        else {
	        dfLTD <- do.call( rbind, lapply( split(pdf, pdf$c),
		        function(x) {
			        n1 = sum(x$t)
			        n0 = length(x$t) - n1
			        if(n1 == 0 || n0 == 0)
                        LTD = NA
                    else
                        LTD = round( sum(as.numeric(x$y * x$t))/n1 - 
                                sum(as.numeric(x$y * (1-x$t)))/n0, 8)
			        data.frame(c=x$c[1], LTD=LTD, w=length(x$t))
		        } )
            )
            pdf = merge(pdf, dfLTD[,1:2], by.x="c", by.y="c")
        }
        names(pdf) <- c("c", "y", "t", "lstat")
        pdf <- as.data.frame(pdf$lstat)
        names(pdf) <- "lstat"		
        if( i == 1 )
            dfconf = pdf
        else
            dfconf = rbind(dfconf, pdf)   
    }
    RPmean <- round(mean(dfconf$lstat, na.rm = TRUE), 8)        # Weighted Average
    RPstde <- round(sqrt(var(dfconf$lstat, na.rm = TRUE)), 8)   # Weighted Std. Dev.
    suppressWarnings( KSobsD <- ks.test(dfconf$lstat, NUdist$lstat) )
    olist <- c(olist, list(Type=type, NUmean=NUmean, NUstde=NUstde, RPmean=RPmean,
        RPstde=RPstde, KSobsD=KSobsD, NUdist=NUdist, dfconf=dfconf))	
    class(olist) <- "confirm"
    olist
}
