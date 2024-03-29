% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance_functions.R
\name{covariance_i_j_vec}
\alias{covariance_i_j_vec}
\title{Compute the approximate covariance tensor for lag windows defined by i,j}
\usage{
covariance_i_j_vec(f_data, i, j)
}
\arguments{
\item{f_data}{the functional data matrix with observed functions in the columns}

\item{i, j}{the indices i,j in 1:T that we are computing the covariance for}
}
\value{
A 4-dimensional array, encoding the covariance tensor of the functional data for lag
windows defined by i,j.
}
\description{
`covariance_i_j_vec` computes the approximate covariance tensor of the functional data for lag
windows defined by i,j; a vectorized version of covariance_i_j.
}
