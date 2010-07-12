\name{chargeBkgr_call}
\alias{chargeBkgr_call}
\title{ charge accumulation baseline calculation}
\description{
  Calculate MS baseline from charge accumulation.
  (Amplifier capacitor coupling model)
}
\usage{
chargeBkgr_call(tofList, offsetSamples, amp, decay)
}
\arguments{
  \item{tofList}{ -- List of TOF-MS spectra }
  \item{offsetSamples}{ -- TOF interval for the estimate of constant offset }
  \item{amp}{ -- amplitude fraction of accumulated signal }
  \item{decay}{ -- decay constant for the accumulated signal -- exp(-1/tau) }
}

\references{ Malyarenko, D. I., et al., Clin. Chem., 2005, 51(1), p.65 }
\author{ William E Cooke, Professor, Physics, College of William and Mary, wecook@wm.edu}
