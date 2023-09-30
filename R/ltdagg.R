"plot.ltdagg" <-
function (x, envir, show = "all", breaks="Sturges", ...) 
{
	if (!is.environment(envir)) 
        stop("The envir argument to plot.ltdagg() must be an existing environment ...from NUsetup().")	
	xmin <- get("LTDmin", envir=envir)
    xmax <- get("LTDmax", envir=envir)
    if (show != "hist" && show != "boxp" && show != "ecdf" && show != "seq")
        show <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (show == "all")
        par(mfrow=c(2,1))
    else
        par(mfrow=c(1,1))
    LTDdist <- na.omit(x$LTDdist$LTD)
    LTDmean <- mean(LTDdist)	
    if (show != "boxp" && show != "ecdf") {   # show "hist" first...
        hist(LTDdist, breaks = breaks, freq = FALSE, xlim = c(xmin, xmax),
            main = paste("LTD Distribution for", x$K, "Clusters"),
            xlab = "Local Treatment Difference (LTD)", ylab = "Density",
            sub = paste(x$infoclus, "Clusters are Informative"))
        abline(v = LTDmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")	  
        if (show == "seq") {
            cat("\nPress Enter to view the Box-Whisker diagram...")
            scan()
        }
	}
    if (show == "seq" || show == "boxp") {    # Vertical "boxp" display...
        boxplot(LTDdist, ann=FALSE, type="l", lty=1, ylim = c(xmin, xmax))
        abline(h = LTDmean, lty = 2, col = "red")
        abline(h = 0, lty = 1, col = "black")           
        title( main = "Box-Whisker diagram of LTD Distribution",
            xlab = paste("K =", x$K, "Clusters Requested"),
            sub = paste(x$infoclus, "Clusters are Informative"),
            ylab = "Local Treatment Difference (LTD)")
        if (show == "seq") {
            cat("\nPress Enter to view the eCDF plot...")
            scan()
        }
    }
    if (show == "all") {   # Horizontal "boxp" display... 
        boxplot(LTDdist, ann=FALSE, type="l", lty=1, ylim = c(xmin, xmax),
            horizontal = TRUE)
        abline(v = LTDmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")	  
        title( main = "Box-Whisker diagram of LTD Distribution",
            xlab = paste("K =", x$K, "Clusters Requested"),
            sub = paste(x$infoclus, "Clusters are Informative"))
        cat("\nPress Enter to view the eCDF plot...")
        scan()
    }
    if (show != "boxp" && show != "hist") {   # show "ecdf" last...
        if (show == "all")
            par(mfrow=c(1,1))
        plot(ecdf(LTDdist), verticals=TRUE, do.points=FALSE, ann=FALSE, col="blue2", pch=46,
            xlim=c(xmin, xmax))
        abline(v = LTDmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")	
        title(main = paste("LTD eCDF for", x$K, "Clusters"),
            ylab = "Probability", xlab = "Local Treatment Difference (LTD)",
            sub = paste(x$infoclus, "Clusters are Informative"))
    }
}

"print.ltdagg" <-
function (x, ...) 
{
    cat("\nltdagg Object: Local Treatment Differences\n")
    cat("Hierarchical Clustering object:", x$hiclus, "\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("Binary Treatment Indicator:", x$trtm, "\n")
    cat("Number of Experimental Units (Patients) =", length(x$LTDdist$y), "\n")
    cat("Number of Clusters Requested:", x$K, "\n")
    if (x$actclust != x$K) 
        cat("Number of Clusters Delivered =", x$actclust, "\n")
    cat("Number of Informative Clusters =", x$infoclus, "\n")
    cat("\n    Number of Experimental Units with LTD estimates =", x$infounits)
    cat("\n    Average of these Local Treatment Differences =", x$LTDmean)
    cat("\n    Standard Deviation of these LTD estimates =", x$LTDstde, "\n\n") 
}

"ltdagg" <-
function (K, envir) 
{
    if (missing(K)) 
        stop("K, the Number of Clusters Requested, must be the first argument to ltdagg().")
    if (!is.environment(envir)) 
        stop("The envir argument to ltdagg() must be an existing environment ...from NUsetup().")
    if (K < 2) { 
        cat("\nCalculations for a single Cluster, K=1, were performed in NUsetup().")
        cat("\nThe K argument to ltdagg() must request at least 2 Clusters.\n\n")
        return(NULL)
    }
    K = floor(K)
    Kmax <- get("Kmax", envir=envir)
    if (K > Kmax) {
        cat("\nThis choice of K=", K, "exceeds our Guideline for the Maximum Number of")
        cat("\nClusters to Request, Kmax=", Kmax, ", ...even in Most Favorable Cases!\n\n")
        return(NULL)    # Too Many Clusters Requested...
    }
    pars <- get("pars", envir=envir)
    hiclus <- get(pars[1])    # get from .GlobalEnv
    dframe <- get(pars[2])    # get from .GlobalEnv
    trtm <- pars[3]
    yvar <- pars[4]
    NumLevels <- get("NumLevels", envir = envir)
    if (NumLevels != 2) {
        cat(paste("\nThe Treatment variable", trtm, "is an Exposure with", NumLevels, "levels."))
        cat("\nLocal Treatment Difference (LTD) analyses are not applicable here.\n\n")
        return(NULL)
    }	
    aggdf <- get("aggdf", envir=envir)
    if ( K %in% aggdf$Blocks ) {
        cat(paste("\nltdagg() calculations for K =", K, "Clusters have already been performed.\n\n"))
        stop("Please specify only a few, NEW values of K within the current NUsetup().")
    }	
    cnam <- paste0("c.", K)
    lcblks <- as.data.frame(as.factor(cutree(hiclus$hclus, k = K)))
    Blocks <- length(table(lcblks))
    N <- length(dframe[,trtm])
    lcmerge <- cbind(1:N, dframe[,yvar], dframe[,trtm], lcblks)
    names(lcmerge) <- c("ID", "y", "t", "c")
    olist <- list(hiclus = pars[1], dframe = pars[2],
        trtm = trtm, yvar = yvar, K = K, actclust = Blocks)
    dfLTD <- do.call( rbind, lapply( split(lcmerge, lcmerge[,4]),
	    function(x) {
	        x <- na.omit(x)       
	        n1 = sum(x$t)
	        n0 = length(x$t) - n1
	        if( n1 == 0 || n0 == 0 ) 
                LTD = NA
            else
                LTD = round( sum(as.numeric(x$y * x$t))/n1 - sum(as.numeric(x$y * (1-x$t)))/n0, 8)
            LAO = round( mean(x$y), 8)
            PS = round( mean(x$t), 8)	
            data.frame(c=x$c[1], LTD=LTD, w=length(x$t), LAO=LAO, PS=PS)
            }
        ) 
    )
    lcmerge = merge(lcmerge, dfLTD[,1:2], by.x="c", by.y="c")
    names(lcmerge) <- c(cnam, "ID", "y", "t", "LTD")
    olist <- c(olist, list(LTDtabl = dfLTD, LTDdist = lcmerge))
    dfLTD = na.omit(dfLTD)			# Exclude Uninformative Clusters in following...
    infoclus = length(dfLTD[,2])    # Number of LTDs != NA
    infounits = sum(dfLTD[,3])      # Sum of corresponding LTD weights
    olist <- c(olist, list(infoclus = infoclus, infounits = infounits ))
    awLTD <- round(mean(lcmerge$LTD, na.rm = TRUE), 8) # Weighted Average LTD
    awstde <- round(sqrt(var(lcmerge$LTD, na.rm = TRUE)), 8) # Weighted Std. Dev.
    olist <- c(olist, list(LTDmean = awLTD, LTDstde = awstde))
    aggnew <- data.frame(Label = "LTD", Blocks = Blocks, LTDmean = awLTD, LTDstde = awstde)
    if ( length(aggdf[,1]) == 1 && aggdf[1,1] == "TEMP" )
        assign("aggdf", aggnew, envir=envir)
    else {
        aggdf <- as.data.frame(rbind(aggdf, aggnew))
        assign("aggdf", aggdf, envir=envir)	
    }
    # get, update & replace boxdf in envir...
    boxdf <- get("boxdf", envir=envir)
    boxnew <- as.data.frame(lcmerge$LTD)
    kdf <- as.data.frame(rep(K, length(boxnew[,1])))
    boxnew <- cbind(boxnew, kdf)
    names(boxnew) <- c("NUstat", "K")
    boxdf <- as.data.frame(rbind(boxdf, boxnew))
    assign("boxdf", boxdf, envir=envir)		
    # finalize olist and envir info...
    localLTDmin <- min(dfLTD$LTD)
    LTDmin <- get("LTDmin", envir=envir)
    if (LTDmin > localLTDmin)
        assign("LTDmin", localLTDmin, envir=envir)
    localLTDmax <- max(dfLTD$LTD)
    LTDmax <- get("LTDmax", envir=envir)
    if (LTDmax < localLTDmax)
        assign("LTDmax", localLTDmax, envir=envir)
    class(olist) <- "ltdagg"
    olist
}	
