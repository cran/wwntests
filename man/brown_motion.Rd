% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_simulations.R
\name{brown_motion}
\alias{brown_motion}
\title{`brown_motion` Creates at J x N matrix, containing N independent Brownian motion sample paths in
each of the columns.}
\usage{
brown_motion(N, J)
}
\arguments{
\item{N}{the number of independent Brownian motion sample paths to compute.}

\item{J}{the number of steps observed for each sample path (the resolution of the data).}
}
\value{
A J x N matrix containing Brownian motion functional data in the columns.
}
\description{
`brown_motion` Creates at J x N matrix, containing N independent Brownian motion sample paths in
each of the columns.
}
\examples{
b <- brown_motion(250, 50)

}
