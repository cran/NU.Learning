"plot.lrcagg" <-
function (x, envir, show="all", breaks="Sturges", ...) 
{
	if (!is.environment(envir)) 
        stop("The envir argument to plot.lrcagg() must be an existing environment ...from NUsetup().")	
	xmin <- get("LRCmin", envir=envir)
    xmax <- get("LRCmax", envir=envir)
    if (show != "hist" && show != "boxp" && show != "ecdf" && show != "seq")
        show <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (show == "all")
        par(mfrow=c(2,1))
    else
        par(mfrow=c(1,1))
    LRCdist <- na.omit(x$LRCdist$LRC)
    LRCmean <- mean(LRCdist)	
    if (show != "boxp" && show != "ecdf") {   # show "hist" first...
        hist(LRCdist, breaks = breaks, freq = FALSE, xlim = c(xmin, xmax),
            main = paste("LRC Distribution for", x$K, "Clusters"),
            xlab = "Local Rank Correlation (LRC)", ylab = "Density",
            sub = paste(x$infoclus, "Clusters are Informative"))
        abline(v = LRCmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")	
        if (show == "seq") {
            cat("\nPress Enter to view the Box-Whisker diagram...")
            scan()
        }
	}
    if (show == "seq" || show == "boxp") {    # Vertical "boxp" display...
        boxplot(LRCdist, ann=FALSE, type="l", lty=1, ylim = c(xmin, xmax))
        abline(h = LRCmean, lty = 2, col = "red")
        abline(h = 0, lty = 1, col = "black")	          
        title( main = "Box-Whisker diagram of LRC Distribution",
            xlab = paste("K =", x$K, "Clusters Requested"),
            sub = paste(x$infoclus, "Clusters are Informative"),
            ylab = "Local Rank Correlation (LRC)")
        if (show == "seq") {
            cat("\nPress Enter to view the eCDF plot...")
            scan()
        }
    }
    if (show == "all") {   # Horizontal "boxp" display... 
        boxplot(LRCdist, ann=FALSE, type="l", lty=1, ylim = c(xmin, xmax),
            horizontal = TRUE)
        abline(v = LRCmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")  
        title( main = "Box-Whisker diagram of LRC Distribution",
            xlab = paste("K =", x$K, "Clusters Requested"),
            sub = paste(x$infoclus, "Clusters are Informative"))
        cat("\nPress Enter to view the eCDF plot...")
        scan()
    }
    if (show != "boxp" && show != "hist") {   # show "ecdf" last...
        if (show == "all")
            par(mfrow=c(1,1))
        plot(ecdf(LRCdist), verticals=TRUE, do.points=FALSE, ann=FALSE, col="blue2", pch=46,
            xlim=c(xmin, xmax))
        abline(v = LRCmean, lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")			
        title(main = paste("LRC eCDF for", x$K, "Clusters"),
            ylab = "Probability", xlab = "Local Rank Correlation (LRC)",
            sub = paste(x$infoclus, "Clusters are Informative"))
    }
}

"print.lrcagg" <-
function (x, ...) 
{
    cat("\nlrcagg Object: Local Rank Correlations\n")
    cat("Hierarchical Clustering object:", x$hiclus, "\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("Exposure Level Indicator:", x$trtm, "\n")
    cat("Number of Experimental Units (Patients) =", length(x$LRCdist$y), "\n")
    cat("Number of Clusters Requested:", x$K, "\n")
    if (x$actclust != x$K) 
        cat("Number of Clusters Delivered =", x$actclust, "\n")
    cat("Number of Informative Clusters =", x$infoclus, "\n")
    cat("\n    Number of Experimental Units with LRC estimates =", x$infounits)
    cat("\n    Average of these Local Rank Correlations =", x$LRCmean)
    cat("\n    Standard Deviation of LRC estimates =", x$LRCstde, "\n\n") 
}

"lrcagg" <-
function (K, envir) 
{
    if (missing(K)) 
        stop("K, the Number of Clusters Requested, must be the first argument to lrcagg().")
    if (!is.environment(envir)) 
        stop("The envir argument to lrcagg() must be an existing environment ...from NUsetup().")
    if (K < 2) { 
        cat("\nCalculations for a Single Cluster, K=1, were performed in NUsetup().")
        cat("\nThe K argument to lrcagg() must request at least 2 Clusters.\n\n")
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
    NumLevels <- get("NumLevels", envir=envir)
    if (NumLevels == 2) {
        cat("\nThe Treatment variable", trtm, "is an Exposure with only 2 levels.\n")
        cat("Local Rank Correlation (LRC) analyses are not meaningful here.\n\n")
        return(NULL)
    }	
    aggdf <- get("aggdf", envir=envir)
    if ( K %in% aggdf$Blocks ) {
        cat(paste("\nlrcagg() calculations for K =", K, "Clusters have already been performed.\n\n"))
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
    dfLRC <- do.call( rbind, lapply( split(lcmerge, lcmerge[,4]),
	    function(x) {
	        x <- na.omit(x)
            LRC = round( cor(x$y, x$t, method = "spearman"), 8)
            LAO = round( mean(x$y), 8)
            PS = round( mean(x$t), 8)	
            data.frame(c=x$c[1], LRC=LRC, w=length(x$t), LAO=LAO, PS=PS)
            }
        ) 
    )
    lcmerge = merge(lcmerge, dfLRC[,1:2], by.x="c", by.y="c")
    names(lcmerge) <- c(cnam, "ID", "y", "t", "LRC")
    olist <- c(olist, list(LRCtabl = dfLRC, LRCdist = lcmerge))
    dfLRC = na.omit(dfLRC)			# Exclude Uninformative Clusters in the following...
    infoclus = length(dfLRC[,2])    # Number of LRCs != NA
    infounits = sum(dfLRC[,3])      # Sum of corresponding LRC weights
    olist <- c(olist, list(infoclus = infoclus, infounits = infounits ))
    awLRC <- round(mean(lcmerge$LRC, na.rm = TRUE), 8) # Weighted Average LRC
    awstde <- round(sqrt(var(lcmerge$LRC, na.rm = TRUE)), 8) # Weighted Std. Dev.
    olist <- c(olist, list(LRCmean = awLRC, LRCstde = awstde))
    aggnew <- data.frame(Label = "LRC", Blocks = Blocks, LRCmean = awLRC, LRCstde = awstde)
    if ( length(aggdf[,1]) == 1 && aggdf[1,1] == "TEMP" )
        assign("aggdf", aggnew, envir=envir)
    else {
        aggdf <- as.data.frame(rbind(aggdf, aggnew))
        assign("aggdf", aggdf, envir=envir)	
    }
    # get, update & replace boxdf in envir...
    boxdf <- get("boxdf", envir=envir)
    boxnew <- as.data.frame(lcmerge$LRC)
    kdf <- as.data.frame(rep(K, length(boxnew[,1])))
    boxnew <- cbind(boxnew, kdf)
    names(boxnew) <- c("NUstat", "K")
    boxdf <- as.data.frame(rbind(boxdf, boxnew))
    assign("boxdf", boxdf, envir=envir)		
    # finalize olist and envir info...
    localLRCmin <- min(dfLRC$LRC)
    LRCmin <- get("LRCmin", envir=envir)
    if (LRCmin > localLRCmin)
        assign("LRCmin", localLRCmin, envir=envir)
    localLRCmax <- max(dfLRC$LRC)
    LRCmax <- get("LRCmax", envir=envir)
    if (LRCmax < localLRCmax)
        assign("LRCmax", localLRCmax, envir=envir)
    class(olist) <- "lrcagg"
    olist
}	
