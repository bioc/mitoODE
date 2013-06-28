\name{plotfit}

\alias{plotfit}
\alias{plotk}

\title{Plot fitted values}

\description{
The functions \code{plotfit} and \code{plotk} plot the fitted values and the transition functions across time.
}

\usage{
plotfit(spot, p=pheno[spot,], showfit=TRUE, legend="topleft", kk=NULL, cex=1, kcol='#ffaa77', lwd=1, xlab="Time after cell seeding (h)", ylab="Number of cells", ...)
plotk(spot, p=pheno[spot,], kk=c("kim"), height, kcol, lwd=1, xlab="Time after seeding (h)", ylab=kk, type='l', ...)
}

\arguments{
  \item{spot}{A numeric indicating a spot ID, synchronised with the global
      Mitocheck \code{tab} object.}
  \item{p}{A vector of 10 numeric parameters, such as returned by \code{getp0} or \code{fitmodel}, 
  containing the fitted parameters. Default are the previously fitted parameters in the \code{pheno} matrix.}
  \item{showfit}{A logical. If \code{TRUE}, the fitted data is plotted over the cell count time series. Default is \code{TRUE}.}
  \item{kk}{A character vector containing the transition rates to plot. Values can be \code{kim}, \code{kmi},  \code{kmp} or \code{kd}.}
  \item{kcol}{A named vector of colors for inflection points. Names take the same values as \code{kk}.}
  \item{legend, cex, lwd, xlab, ylab, height, type, ...}{Graphical parameters, see \code{par}.}
}

\value{
  No values.
}

\details{
The function \code{plotfit} plots the cell count times series and the fitted data. The function \code{plotk} plots the transition rate functions. 
}

\seealso{
  \code{\link{fitmodel}}, \code{\link{loadFittedData}}
}

\author{
  Gregoire Pau, \email{pau.gregoire@gene.com}, 2012
}

\examples{
 ## load previously fitted data
 loadFittedData()

 ## plot fitted data
 plotfit(156205)
}
