getp0 <- function(y, p0=NULL, sd=0) { 
  ## if not possible, use a predefined p0
  if (is.null(p0)) {
    tmed <- g.tstart+(g.tstop-g.tstart)/2
    p0 <- c(him=0, hmi=0, hmp=0, ha=0,
            tim=tmed, tmi=tmed, tmp=tmed, ta=tmed,
            mu=0.0, i0=40)
  }

  ## cut p0 to par size
  p0 <- p2par(p0)

  ## add noise
  if (sd>0) {
    p0 <- p0 + rnorm(length(p0), sd=sd) ## all
    p0[5:8] <- p0[5:8] + rnorm(4, sd=100*sd) ## extra noise
    p0[10] <- p0[10] + rnorm(1, sd=80*sd) ## extra noise
  }

  ## apply constraints
  constrain(p0)
}

par2p <- function(par, p) {
   p[1:10] <- par
   constrain(p)
}

p2par <- function(p) {
  p[1:10]
}

constrain <- function(p) {
  ## all terms except him and hmi are positive
  p[3:10] <- abs(p[3:10])
  
  ## tim, tmi, ta, tmp
  p[5:8] <- ccirc(p[5:8], g.tstop)

  ## mu
  p[9] <- ccirc(p[9], 1)

  ## return p
  p
}

## constrain z to [0, t]
ccirc <- function(z, t) {
  abs(t-abs(z%%(2*t)-t))
}

sig <- function(t, a, h, tau, gamma=1) {
  a + h / (1+exp(-gamma*(t-tau)))
}

penalty <- function(p) {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    p.lambda <- NULL
  }
  
  z <- c(p[1:4], p[9])
  sqrt(p.lambda * abs(z))
}

odevaluate <- function(p, nt=NULL, oversampling=10) {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    g.mit0 <- NULL
    g.kim <- NULL
    g.kmi <- NULL
  }
  
  ## build time points
  vt <- seq(0, 70, by=0.5/oversampling)

  ## main component
  mu <- p[9]
  vy0 <- c((1-g.mit0)*p[10]*(1-mu), g.mit0*p[10]*(1-mu), 0, 0)
  vp <- c(p[1:10], g.kim, g.kmi)
  y1 <- .Call('rksolve_wrap', vt, vy0=vy0, vp=vp, 4L)

  ## normal contamination
  vy0 <- c((1-g.mit0)*p[10]*(mu), g.mit0*p[10]*(mu), 0, 0)
  vp <- c(rep(0, 10), g.kim, g.kmi)
  y2 <- .Call('rksolve_wrap', vt, vy0=vy0, vp=vp, 4L)
  y <- y1 + y2

  comb <- seq(1, length(vt), by=oversampling)
  vt <- vt[comb]
  y <- y[comb,]
  if (!is.null(nt)) {
     it0 = which(vt>=g.tstart)[1]
     it1 = it0 + nt - 1
     re = it0:it1
     vt = vt[re]
     y = y[re,]
  }
  list(t=vt, y=y)
}
