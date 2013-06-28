plotfit <- function(spot, p=pheno[spot,], showfit=TRUE, legend="topleft", kk=NULL, cex=1, kcol='#ffaa77', lwd=1, xlab="Time after cell seeding (h)", ylab="Number of cells", ...) {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
  }
  
  if (missing(p)) {
    if (exists('pheno')) p <- pheno[spot,]
    else {
      p <- getp0()
      showfit <- FALSE
    }
  }
  
  y <- try(readspot(spot), silent=TRUE)
  if (class(y)=='try-error') cat('warning: cannot read spot=', spot, '\n')
  else {
    nt <- nrow(y)
    fs <- compute.fitstat(y, p)
    yf <- odevaluate(p, nt=nt)
    vt <- yf$t

    cola <- c("#00000077", "#FF000077", "#00CD0077", "#0000FF77")
    colb <- c("#000000", "#FF0000", "#00CD00", "#0000FF")
    if (showfit) {
      matplot(vt, y, type="l", xlab=xlab, ylab=ylab, lty=1, lwd=lwd, cex=cex, col=cola, ...)
      matlines(vt, yf$y, type="l", lty=2, lwd=lwd, cex=cex, col=colb, ...)
    } else {
      matplot(vt, y, type="l", xlab=xlab, ylab=ylab, lty=1, lwd=lwd, cex=cex, col=colb, ...)
    }
    if (!is.null(kk)) {
      for (k in kk) {
        t <- switch(k, "kim"=p["tim"], "kmi"=p["tmi"], "kmp"=p["tmp"], "ka"=p["ta"])
        abline(v=t, col=kcol[k], lwd=lwd)
      }
    }
    if (!is.null(legend)) legend(legend,
                                 legend=c('interphase', 'mitotic', 'polynucleated', 'dead'),
                                 col=c(1,2,3,4), lty=1, lwd=lwd, cex=cex, seg.len=1,
                                 bg="#ffffff")
  }
}

plotk <- function(spot, p=pheno[spot,], kk=c("kim"), height, kcol, lwd=1, xlab="Time after seeding (h)", ylab=kk, type='l', ...) {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    g.kim <- NULL
    g.kmi <- NULL
  }
  
  y <- readspot(spot)
  nt <- nrow(y)
  yf <- odevaluate(p, nt=nt)
  vt <- yf$t
  a <- switch(kk, "kim"=g.kim, "kmi"=g.kmi, 0)
  h <- switch(kk, "kim"=p["him"], "kmi"=p["hmi"], "kmp"=p["hmp"], "ka"=p["ha"])
  t <- switch(kk, "kim"=p["tim"], "kmi"=p["tmi"], "kmp"=p["tmp"], "ka"=p["ta"])
  plot(vt, sig(vt, a, h, t), xlab=xlab, type=type, lwd=lwd, ylab=ylab, ...)
  abline(v=t, col=kcol[kk], lwd=lwd)  
}


