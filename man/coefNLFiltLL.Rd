\name{coefNLFiltLL}
\alias{coefNLFiltLL}
\title{ Nonlinear target filter coefficients for Lorentzian signal}
\description{
  Calculate nonlineartarget filter coefficients for Lorentzian signal wavelet
}
\usage{
coefNLFiltLL(inhw, las, trct)
}
\arguments{
  \item{inhw}{ -- left half-width of an input signal }
  \item{las}{ -- right half-width asymmetry (can be negative)}
  \item{trct}{ -- target wavelet amplitude truncation threshold}
}

\references{ Malyarenko D., et.al., Rapid Commun. Mass Spec., 2006, 20 (11), p.1661  }
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
