% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block_bootstrap.R
\name{block_bootsrap}
\alias{block_bootsrap}
\title{`block_bootstrap` Performs a block bootstrap on the functional data f_data with block size b.}
\usage{
block_bootsrap(f_data, b, B = 300, moving = FALSE)
}
\arguments{
\item{f_data}{the functional data matrix with observed functions in the columns}

\item{b}{the block size (of each block in each bootstrap sample)}

\item{B}{the number of bootstraps samples}

\item{moving}{boolean value specifying whether the block bootstrap should be moving or not. A moving black
bootstrap samples individual functional observations and adds on the consequent block, rather than sampling
blocks of the data.}
}
\value{
Returns a list of B elements, each element being a block bootstrap sample in the same format
as the original functional data f_data.
}
\description{
`block_bootstrap` Performs a block bootstrap on the functional data f_data with block size b.
}
