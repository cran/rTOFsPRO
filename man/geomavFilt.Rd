\name{geomavFilt}
\alias{geomavFilt}
\title{ Non-linear target filter }
\description{
  Applies non-linear target filter as geometric average of three filtered signals
}
\usage{
geomavFilt(sigIN, trct, fc3list)
}
\arguments{
  \item{sigIN}{ -- input TOF spectrum }
  \item{trct}{ -- amplitude truncatin (precision) threshold }
  \item{fc3list}{  -- list of three non-linear filter coefficient wavelets }
}

\references{ Malyarenko D., et.al., Rapid Commun. Mass Spec., 2006, 20 (11), p.1661  }
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
