\name{msResample}
\alias{msResample}
\title{ Integrative down-sampling }
\description{
 Function resamples mass axis and integrates signal with the step
 determined by the ratio of predicted peak-width for the peak outside
 mass focusing range (given by quadratic polynomial fit) to the
 intitial STD (left half) of the peak in the mass focusing range
}
\usage{
msResample(sigt, t0, masst, sigma0, asym, pfit, mHW, Gnoise)
}
\arguments{
  \item{sigt}{ -- raw equisampled TOF signal }
  \item{t0}{ -- input time ( may start above 1) }
  \item{masst}{ -- mass axis }
  \item{sigma0}{ -- peak-width in the mass-focusing range }
  \item{asym}{ -- right peak asymmetry in mass focusing range (can be negative) }
  \item{pfit}{ -- array of quadrtic fit coefficients [p1,p2,p3] for peak-width dependence on TOF}
  \item{mHW}{ -- sufficient constant point density (peak HW)}
  \item{Gnoise}{ -- Gaussian noise model (if 1) }
}
\references{ Malyarenko, D. I., et al., Clin. Chem., 2005, 51(1), p.65}
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
