\name{figures}

\alias{loadFittedData}
\alias{figure1}
\alias{figure2}
\alias{figure3a}
\alias{figure3b}
\alias{figure4}
\alias{pheno}

\title{Reproduce the results of the mitoODE paper}

\description{
The functions \code{loadFittedData} loads all the fitted parameters of the Mitocheck screen. The figure functions contain the code to produce the figures of the paper.
}

\usage{
loadFittedData()
figure1()
figure2()
figure3a()
figure3b()
figure4()
}

\value{
  No values.
}

\details{
The function \code{loadFittedData} loads the data.frame \code{pheno} in the global environment.
The figure functions produce PDF plots in the working directory.
}

\seealso{
  \code{\link{fitmodel}}, \code{\link{plotfit}}
}

\author{
  Gregoire Pau, \email{pau.gregoire@gene.com}, 2012
}

\examples{
 \dontrun{
loadFittedData()
figure1()
figure2()
figure3a()
figure3b()
figure4()
 }
}
