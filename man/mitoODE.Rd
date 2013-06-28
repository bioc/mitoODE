\name{mitoODE}

\alias{mitoODE}
\alias{getp0}
\alias{fitmodel}

\title{Dynamical modelling of phenotypes in a genome-wide RNAi
live-cell imaging assay}

\description{
The \code{mitoODE} package implements the differential model fitting used
in the ``Dynamical modelling of phenotypes in a genome-wide RNAi
live-cell imaging assay'' paper and includes the code to reproduce the presented figures and results.
See the vignette for details about the model.
}

\usage{
getp0(y, p0=NULL, sd=0)
fitmodel(y, p0, pconst, nfits=1, sd=0, mc.cores=1, best=TRUE)
}

\arguments{
  \item{y}{A matrix with 4 columns containing the cell count time series, such
  as the one returned by \code{readspot}.}
  \item{p0}{A vector of 10 numeric parameters, such as returned by \code{getp0}, containing the initial conditions. See the vignette for details.}
  \item{pconst}{A vector of 10 numeric parameters, containing \code{g.kim}, \code{g.kmi}, \code{g.mit0} and \code{p.lambda}. See the vignette for details.}
  \item{nfits}{The number of fits to be done. Default is 1.}
  \item{sd}{A numeric indicating the standard deviation of some Gaussian noise added to the initial conditions before fitting. Default is 0.}
  \item{mc.cores}{A numeric indicating how many parallel cores should be used to fit the cell count time series. Default is 1.}
  \item{best}{A logical. See return value.}
}

\value{
 \code{getp0} returns a vector of 10 parameters that can be used as initial conditions.
 If \code{TRUE}, \code{fitmodel} returns the best fitted parameters, otherwise returns the all \code{n} fitted parameters.
}

\details{
  The function \code{fitmodel} fits the cell count data to the differential equation model described in the vignette, using
  the penalized Levenberg-Marquardt \code{nls.lm} algorithm.
}

\seealso{
  \code{\link{loadFittedData}}, \code{\link{figure1}}, \code{\link{plotfit}}
}

\author{
  Gregoire Pau, \email{pau.gregoire@gene.com}, 2012
}

\examples{
 ## read spot
 spotid <- 1000
 y <- readspot(spotid)

 ## fit data
 pconst <- c(g.kim=0.025, g.kmi=0.57, g.mit0=0.05, p.lambda=4)
 p0 <- getp0()
 pp <- fitmodel(y, p0, pconst)
 round(pp, 2)

 ## plot fitted data
 plotfit(spotid, pp)
}
