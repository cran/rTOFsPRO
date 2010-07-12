\name{findBsEnv}
\alias{findBsEnv}
\title{ Baseline envelope calculation}
\description{
  Finds baseline envelope for analytical fit assuming decaying baseline
}
\usage{
findBsEnv(specDat, idsWindow)
}
\arguments{
  \item{specDat}{ -- TOF-MS spectrum record}
  \item{idsWindow}{ -- down-sampling window for compression }
}
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Called by ExpBkgr function for auto-fit }
