\name{offsetBkgr_call}
\alias{offsetBkgr_call}
\title{Constant ADC offset background }
\description{
  Finds constant baseline offset from a specified portion of a spectrum
}
\usage{
offsetBkgr_call(tofList, offsetSamples)
}
\arguments{
  \item{tofList}{ -- list of TOF spectra}
  \item{offsetSamples}{ -- TOF interval for the estimate of constant offset (outised of decay range) }
}
