\name{myxcorr}
\alias{myxcorr}
\title{ Cross-correlation of two vectors}
\description{
  Compute Matlab-like cross-correlation for two signal vectors:
  signals are zero-padded to the next (upper) close power of 2.
}
\usage{
myxcorr(x, y, maxlag)
}
\arguments{
  \item{x}{ -- first observed sampled discrete signal }
  \item{y}{ -- second observed sampled discrete signal }
  \item{maxlag}{ -- maximum time-lag for correlation }
}

\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }

