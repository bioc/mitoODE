fitmodel <- function(y, p0, pconst, nfits=1, sd=0, mc.cores=1, best=TRUE) {
  ## use pconst
  assign("g.kim", pconst["g.kim"], 1)
  assign("g.kmi", pconst["g.kmi"], 1)
  assign("g.mit0",  pconst["g.mit0"], 1)
  assign("p.lambda", pconst["p.lambda"], 1)
  
  ## one fit
  if (nfits==1) {
    ## fit model
    par <- try({
      p0 <- getp0(y, p0, sd=sd)
      nls.lm(par=p2par(p0), fn=lscriterion, jac=NULL, y=y, p=p0,
             control=nls.lm.control(maxiter=200))$par
    })
    
    ## return value
    if (class(par)=='try-error' || any(is.na(par))) return(NULL)
    else {
      p <- par2p(par, p0)
      return(c(p, compute.fitstat(y, p)))
    }
  } else {
    ## multiple fits
    pp <- mclapply(1:nfits, function(i) fitmodel(y=y, p0=p0, pconst=pconst, nfits=1, sd=sd), mc.cores=mc.cores)
    zna <- sapply(pp, function(z) any(is.null(pp)))
    pp <- pp[!zna]
    pp <- do.call(rbind, pp)
    if (!is.null(pp)) {
      pp <- pp[order(pp[,'score']),]
      if (best) pp <- pp[1,]
    }
    return(pp)
  }
}

fitspot <- function(spotids, usepheno=TRUE, nfits=4, sd=0.05, mc.cores=4, trace=FALSE) {
   ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
  }
  
  pp <- mclapply(spotids, function(spotid) {
    p <- try({
      y <- readspot(spotid)
      if (usepheno) p0 <- pheno[sample(which(!is.na(pheno[,"score"])), 1),]
      else p0 <- getp0(y)
      p <- fitmodel(y, p0=p0, sd=sd, nfits=nfits)
      if (trace) cat("spotid=", spotid, "score=", p["score"], "\n")
      p
    })
    if (class(p)=="try-error" || is.null(p)) NA
    else p
  }, mc.cores=mc.cores)
  pp <- do.call(rbind, pp)
  rownames(pp) <- spotids
  pp
}

compute.fitstat <- function(y, p) {
  yf <- odevaluate(p, nt=nrow(y))$y
  yfmax <- apply(yf, 2, max)
  names(yfmax) <- c("imax", "mmax", "smax", "amax")
  pen <- sum(penalty(p)^2)
  score <- sum(lscriterion(p2par(p), y, p)^2)
  rss <- score - pen
  c(score=score, rss=rss, pen=pen, yfmax)
}

fresiduals <- function(y, p) {
  odevaluate(p, nt=nrow(y))$y - y
}

lscriterion <- function(par, y, p) {
  p <- par2p(par, p)
  c(as.numeric(fresiduals(y, p))/sqrt(nrow(y)), as.numeric(penalty(p)))
}
