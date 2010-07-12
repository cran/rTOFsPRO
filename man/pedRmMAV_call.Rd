\name{pedRmMAV_call}
\alias{pedRmMAV_call}
\title{ Removes peak pedestal }
\description{
  	Estimates pedestal baseline using local minima
	Smoothes baseline estimate using moving average
	Subtracts smooth baseline estimate from signal
}
\usage{
pedRmMAV_call(tofList, pw, nf)
}
\arguments{
  \item{tofList}{-- list of TOF spectra }
  \item{pw}{-- peak half width }
  \item{nf}{-- number fraction of width in window }
}

\author{ Dariya Malyarenko and Maureen Tracy, William and Mary Research Institute, dimaly@wm.edu, mbtrac@wm.edu }
\note{Assumes constant peak width accross a signal record }
