% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gauss_cov}
\alias{gauss_cov}
\title{Computing the Gaussian covariance function for a single distance
measure.}
\usage{
gauss_cov(dists, sigsq, phi)
}
\arguments{
\item{dists}{a numeric matrix representing the distance between spatial
entities.}

\item{sigsq}{the \eqn{\sigma^2} parameter from the Gaussian covariance
function.}

\item{phi}{the \eqn{\phi} parameter from the Gaussian covariance function,
controls the range of the spatial dependence.}
}
\value{
The Gaussian covariance function (for a stationary and
isotropic process) associated with the provided distances (\code{dists})
and the given set of parameters.
}
\description{
Computing the Gaussian covariance function for a single distance
measure.
}
\seealso{
\code{\link{single_exp}}, \code{\link{single_matern}}
}
\keyword{internal}
