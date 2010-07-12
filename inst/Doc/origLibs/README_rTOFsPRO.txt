
The rTOFsPRO package is an enhanced equivalent of TOFsPRO Matlab toolbox
(see matlabcentral/fileexchange/24469),with added ability to perform 
preliminary global alignment or initial delay correction in TOF. The
deployed R-package will be available from CRAN after 07/2010. 
See "rTOFsPRO_overview" for package documentation and details.

The "main" function links libraries to perform signal processing of the input
tofList structure (using tofListMetaData parameters) to produce alignedPeakList
output. "alignedPeakList$peaks" comtains aligned peak positions, while
alignedPeakList$data[[i]] includes intensity and uncreatinty information 
for each input spectrum <- tofList[[i]]. (See WMBrukerParser CRAN R-package 
for detailed description of tofList structure and meta-data). Down-sampled 
(compressed) and filtered spectra are saved in tofResampled structure
(similar to tofList).  Options and parameteres for processing are set in
"OptionsAndParameters.txt".

To run the rTOFsPRO R-package examples and libraries, 
(1) move the data example files ("tofListPCAqc11.Rdat" & "tofListMetaDataPCAqc11.Rdat") 
    to your working directorty (one level above rTOFsPRO); 
(2) execute mainRTOFsPROc or mainTOFsPROc modular functions in R-environment:
  >> source(paste("rTOFsPRO/mainRTOFsPROC.txt",sep="/")) # or
  >> source(paste("rTOFsPRO/mainRTOFsPROCModular.txt",sep="/")) 
(3) read through the comments in "OptionsAndParametersPCAqc11.txt" and follow
    instructions to change options and parameters by editing variables in this file
(4) brief description of the input/output and signal processing libraries can
    be found in the header of "main.." functions, and in comments for individual
    processing blocks.

