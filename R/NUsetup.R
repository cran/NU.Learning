"NUsetup" <-
function (hclobj, dframe, trex, yvar) 
{
    if (missing(hclobj) || (!inherits(hclobj$hclus, "diana") && 
        !inherits(hclobj$hclus, "hclust"))) 
        stop("First argument to NUsetup must be a diana or hclust object.")	
    hclobj <- deparse(substitute(hclobj))
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("Second argument to NUsetup must be an existing Data Frame.")
    if (missing(trex)) 
        stop("Third argument to NUsetup must name the Treatment or Exposure variable.")
    trex <- deparse(substitute(trex))
    if (!is.element(as.character(trex), dimnames(dframe)[[2]])) 
        stop("Treatment or Exposure must be an existing Data Frame variable.")
    tvec <- dframe[,trex]
    clt <- class(tvec)       # class(clt) is "character"
    if ((clt != "integer") && (clt != "numeric"))
        stop("The levels of the Treatment/Exposure variable must be integer or numeric.")
    Kmax <- floor(length(tvec)/12)    # Guideline: Maximum Number of Clusters = Nobs / 12
    NumLevels <- length(table(tvec))  # Max "levels" in Treatment / Exposure Level indicator
    if (NumLevels < 2)
        stop("Treatment or Exposure Level is identical for all Experimental Units.")
    if (missing(yvar)) 
        stop("Fourth argument to NUsetup must name the y-Outcome variable.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(as.character(yvar), dimnames(dframe)[[2]])) 
        stop("Specified y-Outcome must be an existing Data Frame variable.")
    z <- na.omit(data.frame(cbind(dframe[,yvar], tvec)))
    names(z) <- c("y", "t")
    e <- new.env()
    e$NumLevels <- NumLevels
    e$Kmax <- Kmax
    if (NumLevels > 2) { 
        cat("\nThe Treatment variable is an Exposure with", NumLevels, "levels.")
        cat("\nLocal Treatment Difference (LTD) analyses are not applicable here.")
        cat("\nOnly Local Rank Correlations (LRCs) can be formed Within Clusters.\n\n")
        NUmean = round( cor(z$y, z$t, method = "spearman"), 8)
        e$LRCmin <- min(0, NUmean)
        e$LRCmax <- max(0, NUmean)
        e$aggdf <- data.frame(Label = "TEMP", Blocks = 1, LRCmean = 0, LRCstde = 0)
        }
    else { 
        cat("\nThe Treatment variable has two levels.")
        cat("\nLocal Rank Correlation (LRC) analyses are not applicable here.")
        cat("\nOnly Local Treatment Differences (LTDs) can be formed Within Clusters.\n\n")
        obs <- length(z$t)
        mx <- max(z$t)
        mn <- min(z$t)
        if (mn != 0 || mx != 1) {
          for (i in 1:obs) {
            ifelse(z$t[i] == mx, z$t[i] <- 1, z$t[i] <- 0)
          }
        }
        n1 <- sum(z$t)
        n0 <- obs - n1
        NUmean = round( sum(z$y * z$t)/n1 - sum(z$y * (1-z$t))/n0, 8)
        e$LTDmin <- min(0, NUmean)
        e$LTDmax <- max(0, NUmean)
        e$aggdf <- data.frame(Label = "TEMP", Blocks = 1, LTDmean = 0, LTDstde = 0)		
        }
    e$boxdf <- data.frame(NUmean, 1)
    names(e$boxdf) <- c("NUstat", "K")
    dframe <- deparse(substitute(dframe))
    e$pars <- cbind(hclobj, dframe, trex, yvar)
    e
}
