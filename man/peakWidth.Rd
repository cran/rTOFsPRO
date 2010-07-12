\name{peakWidth}
\alias{peakWidth}
\title{ Peak half-width estimate }
\description{
  Finds peaks and their left half-width accross a TOF record
}
\usage{
peakWidth(s, startsd, nsubdiv, sfr)
}
\arguments{
  \item{s}{ -- TOF spectrum equally sampled in time }
  \item{startsd}{ -- starting point for spectrum subdivision }
  \item{nsubdiv}{ -- number of subdivisions for peak width estimate }
  \item{sfr}{ -- signal fraction (1/sfr) of maximum to include for peak detection}
}
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Auxhiliary function to find peak-width dependenth on TOF and fit parameters for down-sampling }
