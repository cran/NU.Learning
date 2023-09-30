"plot.ivadj" <-
function (x, maxsiz = 0.15, ...) 
{
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(lty = 1, col = "black")
    plot(x$NUtable$PS, x$NUtable$LAO, ann = FALSE, type = "n")
    symbols(x$NUtable$PS, x$NUtable$LAO, circles = sqrt(x$NUtable$w), inches = maxsiz, add = TRUE)
    abline(x$ivfit, lty = 2, lwd=2, col = "red")	
    lines(x$smfit, lty = 2, lwd=2, col = "blue2")
	if (x$Type == 1)
        title(main = paste("IV Adjustment for", length(x$NUtable$LAO), "Clusters"), 
            xlab = "Within-Cluster Treatment Fraction (Propensity)", 
            ylab = "Observed LAO", sub = "Symbol Area proportional to Cluster Size")
    else
        title(main = paste("IV Adjustment for", length(x$NUtable$LAO), "Clusters"), 
            xlab = "Within-Cluster Relative Exposure (Propensity)", 
            ylab = "Observed LAO", sub = "Symbol Area proportional to Cluster Size")			
}

"print.ivadj" <-
function (x, ...) 
{
    cat("\nivadj: Instrumental Variable (IV) Adjustment via Clustering\n")
    cat("\nCluster Tabulation:\n===================\n\n") 
    print(x$NUtable)
    cat("\nSummary of smooth.spline() fit:\n===============================\n\n") 
    print(x$smfit)
    cat("\nSummary of lm() fit:\n====================\n") 	
    print(x$ivsum)
    cat("\nListing of lm() predictions:\n============================\n\n") 	
    print(x$ivpred)
}

"ivadj" <-
function (x) 
{
    if (missing(x) || (!inherits(x, "ltdagg") && !inherits(x, "lrcagg")))  
        stop("First argument to ivadj() must be a ltdagg() or lrcagg() object.")
    if (inherits(x, "lrcagg")) {
        type = 2
	    NUtable <- x$LRCtabl
        expmin <- min(NUtable$PS)
        expmax <- max(NUtable$PS)
        NUtable$PS <- (NUtable$PS - expmin)/(expmax - expmin) # Relative Exposure PS
    } else {
        type = 1
        NUtable = x$LTDtabl
	}
    if (length(NUtable$LAO) < 5) {
        cat("\nIV Adjustment should not be attempted with fewer than 5 Clusters.\n\n")
        return(NULL)
    }
    ivfit <- lm(NUtable$LAO ~ NUtable$PS, weights = NUtable$w)
    ivsum <- summary(ivfit)
    ivpred <- predict(ivfit, data.frame(NUtable$PS), se.fit = TRUE, type="response")
    smfit <- smooth.spline(NUtable$PS, y=NUtable$LAO, w=NUtable$w)
    olist <- list(NUtable = NUtable, Type=type, ivfit=ivfit, ivsum=ivsum,
        ivpred=ivpred, smfit=smfit)
    class(olist) <- "ivadj"
    olist
}
