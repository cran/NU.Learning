"plot.mlme" <-
function (x, NN = 50, breaks = 50, ...) 
{
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow=c(2,1))
    Type <- x$Type
    outdf <- x$outdf
    xmin = min( outdf$effSiz, na.rm = TRUE )
    if (xmin < 0) xmin = 1.2 * xmin
    xmax = max( outdf$effSiz, na.rm = TRUE )
    if (xmax > 0) xmax = 1.2 * xmax
    totl = length( outdf$effSiz )
    for( i in 1:length(NN) ) {
        if( NN[i] > 0.9 * totl || NN[i] < 10 ) next
        hist(outdf$effSiz, breaks = breaks, freq = TRUE, xlim = c(xmin, xmax),
            main = paste("Local Effect-Sizes for all", totl, "eUnits Available"),
            xlab = paste("Local Effect-Size (", Type, ")"), ylab = "Counts")
        abline(v = mean(outdf$effSiz), lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")
        odfsub <- outdf[c( 1:NN[i] ), ]
        hist(odfsub$effSiz, breaks = "Sturges", freq = TRUE, xlim = c(xmin, xmax),
            main = paste("Distribution for ", NN[i], "eUnits Most-Like You"),
            xlab = paste("Local Effect-Size (", Type, ")"), ylab = "Counts")		
        abline(v = mean(odfsub$effSiz), lty = 2, col = "red")
        abline(v = 0, lty = 1, col = "black")
        if ( i < length(NN) ) {
            cat("\nPress Enter to view next MLMe histogram...")
            scan()
        }
    }
}

"print.mlme" <-
function (x, ...) 
{
    cat("\nmlme Object: Most-Like Me Comparisons...\n\n")
    Type <- x$Type
    outdf <- x$outdf
    xmat <- data.frame( t(x$xvec) )   # Version 1 of xmat...
    names( xmat ) <- x$xvars
    cat("xvec - My TARGET X-Vector:\n")
    print( xmat )
    xmat <- data.frame( t(x$varx) )   # Version 2 of xmat...
    names( xmat ) <- x$xvars
    totl = length( outdf$effSiz )
    cat("X-data Variances for all", totl, "eUnits:\n")
    print( xmat )
    cat("Effect-Size Type:", Type, "\n\n")
    cat(Type, "Summary Statistics...\n")
    print( summary( outdf$effSiz ) )
    cat(Type, "Standard Deviation =", sqrt(var(outdf$effSiz, na.rm = TRUE )), "\n")
    cat("\nFirst 10 Nearest Neighbors ...\n")
    xmat <- outdf[1:10, ]             # Version 3 of xmat...
    print( xmat )
}

"mlme" <-
function (envir, hcl, NUagg, xvec) 
{
    if (missing(envir) || !is.environment(envir)) 
        stop("The envir argument to mlme() must be an existing environment created by NUsetup().")
    if (missing(hcl) || !inherits(hcl, "NUcluster"))  
        stop("The hcl argument to mlme() must be a NUcluster() output object.")
    df = hcl$dframe
    xvars = hcl$xvars
    nox = length(xvars)  # number of xvar names in hcl argument.
    if (missing(NUagg) || (!inherits(NUagg, "ltdagg") && !inherits(NUagg, "lrcagg")))  
        stop("The NUagg argument to mlme() must be a ltdagg() or lrcagg() output object.")
    if (inherits(NUagg, "lrcagg")) {	
        Type = "LRC"
	    NUdist <- NUagg$LRCdist
        NUdist$effSiz <- NUdist$LRC
    }
    else {
        Type = "LTD"
	    NUdist <- NUagg$LTDdist
        NUdist$effSiz <- NUdist$LTD
    }
    nxin = length(xvec)  # number of xvec values input.		
    if ( nxin != nox )
        stop("The length of the xvec argument is", nxin, "...when it should be", nox )
    olist <- list( xvec = xvec, Type = Type )	
    # Create "outdf" data.frame 
    NUdist <- NUdist[order(NUdist$ID),]
    inpdf <- get(NUagg$dframe)
    onams <- c( "OD", "D2", "effSiz", names(inpdf)) 
    outdf <- as.data.frame(cbind(NUdist[,2], NUdist[,5], NUdist[,5], inpdf))
    names(outdf) <- onams
    varx = replicate( nox, 0.0 )  # meaningless initial values
    for( i in 1:nox ) {
        if( is.na(xvec[i]) ) stop("Specified xvec contains at least one NA value.")
        varx[i] = var(outdf[ , xvars[i]], na.rm = TRUE)
    }
    olist <- c(olist, list( xvars = xvars, varx = varx ))
    outdf$D2 <- 0.0
    for( i in 1:nox ) {
        outdf$D2 <- outdf$D2 + (outdf[ , xvars[i]] - xvec[i] )^2 / varx[i]  
    }
    outdf <- outdf[ order(outdf$D2), ]
    outdf$OD <- c( 1:length(outdf$OD) )	
    olist <- c(olist, list( outdf = outdf ))
    class(olist) <- "mlme"
    olist
}
