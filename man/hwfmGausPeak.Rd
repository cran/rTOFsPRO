\name{hwfmGausPeak}
\alias{hwfmGausPeak}
\title{ Gaussian peak-HWFM estimate }
\description{
  Finds left-half-width for a peak
}
\usage{
hwfmGausPeak(nf, n_mu)
}
\arguments{
  \item{nf}{ -- signal }
  \item{n_mu}{ -- peak position in the signal}
}
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Aux function called by peakWidth function to find peak-width dependence on TOF }
