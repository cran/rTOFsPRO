\name{locjitdet}
\alias{locjitdet}
\title{ Local shift detection }
\description{
  Detects local shifts near peaks between TOF records by correlation
}
\usage{
locjitdet(fs, ps, sw, startS, endS, spectraNames)
}
\arguments{
  \item{fs}{ -- matrix of sampled discrete signals (columns - intensities)}
  \item{ps}{ -- peak set for input signal matrix }
  \item{sw}{ -- maximum lag for cross-correlation }
  \item{startS}{ -- starting index for shift detection range }
  \item{endS}{ -- ending index fro shift detection range}
  \item{spectraNames}{-- names of spectra }
}
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Called by globalAilgn for auto-detection }
