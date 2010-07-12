\name{peakAlignAv}
\alias{peakAlignAv}
\title{Peak alignment to average spectrum}
\description{
  Aligns TOF peaks to the average spectrum peaks
}
\usage{
peakAlignAv(peaks, tofList, noiseL, HWFM)
}
\arguments{
  \item{peaks}{ -- peak positions found in the averaged spectrum }
  \item{tofList}{ -- list of TOF spectra in the same domain as peaks }
  \item{noiseL}{ -- noise amplitude list for the TOF spectra }
  \item{HWFM}{ -- peak half-width (window for alignment) }
}

\author{ Dariya Malyarenko and Maureen Tracy, William and Mary Research Institute, dimaly@wm.edu }
\note{ Algorithm is similar to Cromwell (Matlab) by MD Anderson }

