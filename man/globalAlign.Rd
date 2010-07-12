\name{globalAlign}
\alias{globalAlign}
\title{ Global alignment}
\description{
  Global alignment of TOF-MS spectra by scaling and offsetting
}
\usage{
globalAlign(tof, tofList, tofListMetaData, HWFM, noiseThres, startS, endS, Appli, Detect)
}
\arguments{
  \item{tof}{ -- baseline subtracted and smoothed tofList (to detect global shifts) }
  \item{tofList}{ -- unprocessed data (to apply detected shifts) }
  \item{tofListMetaData}{ -- meta-data structure, from which "offset" and "scale" are used and /or updated }
  \item{HWFM}{ -- peak width in the range of pivot detection }
  \item{noiseThres}{ -- fraction of maximum signal threshold for pivot detection }
  \item{startS}{ -- TOF start (index) for pivot detection range }
  \item{endS}{ -- TOF end (index) for pivot detection range}
  \item{Appli}{ -- apply correction to the tofList (if "true") }
  \item{Detect}{ -- auto-detect correction for the data }
}

\references{ Tracy M.B., et.al.,Proteomics, 2008, 8 (8), p.1530  }
\author{ Maureen Tracy, William and Mary Research Institute, mbtrac@wm.edu }
