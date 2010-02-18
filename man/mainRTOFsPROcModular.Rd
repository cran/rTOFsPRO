\name{mainRTOFsPROcModular}
\alias{mainRTOFsPROcModular}
\title{ Signal processing for TOF MS }
\description{
	This function links and executes the scripts    
        for signal processing of the linear TOF spectra 	
	to produce aligned peak intensity matrix (and peak positions)   
	that can be used for further statistical analysis
}
\usage{
mainRTOFsPROcModular(OptionsAndParametersFileName)
}
\arguments{
  \item{OptionsAndParametersFileName}{ -- Name of the file containing processing 
	parameters and options }
}
\details{
======================================================================\cr
	PROCESSING STEPS INCLUDE (in a loop over spectra)		\cr
									\cr
	Set up:								\cr
		Load "OptionsAndParameters.txt" file			\cr
									\cr
	Load Data:							\cr
		Load Data and MetaDataSet parameters			\cr
		"multiple files" Loop begins (See note 1)		\cr
									\cr
	Background Correction. Options include:				\cr
		RC (charge accumulation model)				\cr
		Analytical (Gaussian, Exponential, linear) model	\cr
		None (default)						\cr
									\cr
	Preliminary alignment. Options include:				\cr
		Global alignment (using correlation between pivot peaks)\cr
		TOF offset adjustment (from metaData)			\cr
		None							\cr
									\cr
	Downsampling and Filtering. Options include:			\cr
		Integrative downsampling (IDS) only			\cr
		IDS and Optimal Linear Filter (default)			\cr
		IDS and Nonlinear Filter				\cr
		IDS and Matched Filter					\cr
		IDS and MAV Filter					\cr
									\cr
	Pedastal removal. Options include:				\cr
		MAV of local minima 					\cr
		None (default)						\cr
									\cr
	Peak Picking. Options include:					\cr
		Trivial	(1st difference with local minim SNR, Note 2)	\cr
		W & M (Maximum Likelihood Method) (needs extra lib)	\cr
		(call commented off: lib available from wecook@wm.edu)	\cr
									\cr
	Alignment. Options include:					\cr
		Align to average spectra				\cr
		W & M global binning (call commented off: 		\cr
		seprate lib available from wecook@wm.edu)		\cr
									\cr
========================================================================\cr
}
\value{
  \item{alignedPeakList}{ A Structure list of 2 element: "peaks" and "data", where
	"peaks" is a vector of aligned peak positions, and "data" is a list 
	(addressable by spectrumName) of "Intensities","PositionShifts",
	"PositionUncertainties", and"IntensityUncertainties"vectors.}
  \item{alignedPeakListMetaData}{ A data.frame containing string values 
	for experimental meta-data, with rows and columns 
	addressable by spectrumName and attributeName}
}
\references{[1]Malyarenko, D. I., et al., Rapid Commun. Mass Spectrom (2006) 20, p.1670; \cr
		[2]Gatlin-Bunai, C. L., et al., J Proteome Res (2007) 6, p.4517; \cr
		[3]Tracy, M.B., WMBrukerParser R-package (2009): 
	http://cran.r-project.org/web/packages/WMBrukerParser/index.html   }
\author{ Dariya Malyarenko, William and Mary Research Institute, dimaly@wm.edu }
\note{ Best used in "setup" or "tutorial" mode with small data sets < 50 spectra\cr
IMPORTANT: For new data, edit "OptionsAndParameters" file and specify the path to \cr
new working directory, containing data, by changing "TOFfileDIR" (see "Examples")\cr
========================================================================\cr
 Note 1. Due to the large size of data files, in "compression" mode     \cr
 raw tofList data can be loaded and partially processed by IDS from     \cr
 multiple files. The "multiple TOFs loop" includes blocks from "Load    \cr
 data" through "Downsampling and Filtering".At the end of the loop, the \cr
 downsampled spectra and their meta-data are concatenated and the	 \cr
 original tofList is removed.						 \cr
========================================================================\cr
 Note 2. Basic peak picking algorithm adopted from Cromwell toolbox.	 \cr
 (Copyright 2005.) Originally developed by MD Anderson and distributed  \cr
 under "BSD" Mathworks license: bioinformatics.mdanderson.org--   \cr
 cromwell.html. Local minima SNR threshold and baseline-dependent	 \cr
 noise model for each preliminary filtering option was added here.	 \cr
======================================================================\cr
 Note 3. For more details on options and parameters, as well as,        \cr
 individual signal processing routines, and example description	 \cr
 see "rTOFsPRO_overview" manual document in inst -- Doc  			 \cr
======================================================================\cr
}

\seealso{ \code{\link{}}}
\examples{
directory <- system.file("Examples", package = "rTOFsPRO")
OptionsAndParametersFileName <- paste(directory,"/OptionsAndParametersPCAqc11.txt", sep="")
outList <- mainRTOFsPROcModular(OptionsAndParametersFileName)

}
\keyword{IO}

