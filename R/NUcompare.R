"NUcompare" <-
function (envir) 
{
    if (!is.environment(envir)) 
        stop("The envir argument to NUcompare() must be an existing environment created by NUsetup().")	
    aggdf <- get("aggdf", envir=envir)
    Type <- aggdf[1,1]
    boxdf <- get("boxdf", envir=envir)
    # Sort boxdf into ascending values of its 2nd column, K = # Clusters
    boxdf = boxdf[order(boxdf$K),]
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    boxplot(NUstat ~ K, data=boxdf, ann=FALSE, col="lightgray")
    abline(h = boxdf[1,1], lty = 2, col = "red")
    abline(h = 0, lty = 1, col = "blue2")
    if (Type == "LRC") {    
        title( main = "Box-Whisker comparison of LRC Distributions",
            ylab = "Local Rank Correlation (LRC)", xlab = "K = Number of Clusters")
    } else {
        title( main = "Box-Whisker comparison of LTD Distributions",
            ylab = "Local Treatment Difference (LTD)", xlab = "K = Number of Clusters")
    }
    assign("boxdf", boxdf, envir=envir)  # Save Sorted boxdf
}
