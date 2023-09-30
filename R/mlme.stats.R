"mlme.stats" <-
function (x, NN = 50, ...) 
{
    if (missing(x) || !inherits(x, "mlme"))  
        stop("First argument to mlme.stats() must be a mlme() object.")
    cat("\nStats for a << Most-Like Me >> Object...\n\n")
    Type <- x$Type
    outdf <- x$outdf
    totl = length( outdf$effSiz )
    xmat <- data.frame( t(x$xvec) )   # Version 1 of xmat...
    names( xmat ) <- x$xvars
    cat("Reference X-Vector:\n")
    print( xmat )
    xmat <- data.frame( t(x$varx) )   # Version 2 of xmat...
    names( xmat ) <- x$xvars
    cat("X-data Variances for all", totl, "eUnits:\n")
    print( xmat )
    cat("Effect-Size Type:", Type, "\n\n")
    cat("Overall", Type, "Summary Statistics...\n")
    print( summary( outdf$effSiz ) )
    cat(Type, "Standard Deviation =", sqrt(var(outdf$effSiz, na.rm = TRUE )), "\n")	
    for( i in 1:length(NN) ) {
        if( NN[i] > 0.9 * totl || NN[i] < 10 ) next
        cat("\nMost-Like-Me Sub-Group", i, "contains", NN[i], "eUnits.\n")
        odfsub <- outdf[c( 1:NN[i] ), ]
        print( summary( odfsub$effSiz ) )
        cat(Type, "Standard Deviation =", sqrt(var(odfsub$effSiz, na.rm = TRUE )), "\n")
    }
}
