\name{ExpBkgr_call}
\alias{ExpBkgr_call}
\title{ Analytical background correction}
\description{
  TOF-MS baseline calculation and correction according to analytical model
}
\usage{
ExpBkgr_call(tofList, bsSamples, offsetSamples, modelChoice, dataStart)
}
\arguments{
  \item{tofList}{ -- list of TOF spectra with metadata }
  \item{bsSamples}{ -- TOF indices of data to use for model baseline fit }
  \item{offsetSamples}{ -- TOF interval for the estimate of constant offset }
  \item{modelChoice}{ -- choice of analytical baseline model -- 1 - linear, 2 - exponential 3 - Gaussian }
  \item{dataStart}{ -- starting index for TOF record (from metadata)}
}
\references{ Gatlin-Bunai, C. L., et al., J Proteome Res (2007) 6, p.4517 }
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Type function name and reader the header for details on autofit }

