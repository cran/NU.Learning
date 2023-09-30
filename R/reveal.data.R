"reveal.data" <-
function(x, clus.var="Clus", effe.var="eSiz")
{
    # Prepare for NU.Learning analyses: Form a data.frame by appending the LTD or LRC treatment
    # effect-size measure from ltdagg() or lrcagg() as well as a Cluster membership-number
    # variable to a copy of data.frame specified in NUsetup(). 
    if (missing(x) || (!inherits(x, "ltdagg") && !inherits(x, "lrcagg")))  
        stop("First argument to reveal.data() must be a ltdagg() or lrcagg() output object.")
    if (inherits(x, "lrcagg")) {
        type = 2
	    NUdist <- x$LRCdist
    }
    else {
        type = 1
	    NUdist <- x$LTDdist
    }
    NUdist <- NUdist[order(NUdist$ID),]
    inpdf <- get(x$dframe)
    onams <- c( clus.var, effe.var, names(inpdf)) 
    outdf <- as.data.frame(cbind(NUdist[,1], NUdist[,5], inpdf))
    names(outdf) <- onams
    # class(outdf) <- "data.frame"
    outdf
}
