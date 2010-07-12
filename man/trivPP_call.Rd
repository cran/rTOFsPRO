\name{trivPP_call}
\alias{trivPP_call}
\title{ Trivial peak detection }
\description{
  First difference peak detection
}
\usage{
trivPP_call(tofList, snr, naml, HWFM, fwl, stPk, dataStart)
}
\arguments{
  \item{tofList}{ -- list of TOF signals }
  \item{snr}{ -- signal-to-noise ration for peak detection }
  \item{naml}{ -- noise amplitude list (for tofList)}
  \item{HWFM}{ -- half-width at half-maximum for the model signal wavelet }
  \item{fwl}{ -- filter wavelet length applied for signal processing }
  \item{stPk}{ -- starting point (tof-spectrum index) for peak detection }
  \item{dataStart}{ -- defines what is the absolute start location for the input time points }
}
\author{ Dariya Malyarenko and Maureen Tracy, William and Mary Research Institute, dimaly@wm.edu, mbtrac@wm.edu }
\note{ Algorithm is similar to Cromwell (Matlab) by MD Anderson, but with threshold in repsepct to local minima }

