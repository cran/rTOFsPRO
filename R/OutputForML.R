"OutputForML" <-
function(object2save, choice) {

###-----------------------------------------------------------------
### Function outputForML generates files with array data for MatLab
###
### ALLOWED INPUT and corresponding "choice" parameter:
### "ISP" -	inputSpecList (may be original tofList or tofResampled)
### "TRS" -	tofResampled
### "PL" - 	peakList
### "APL" -	alignedPeakList
###
### NOTE: when "TRS" is run, assumes existance for "rsout"- 
###	  as a global variable in workspace
### DEPENDENCY:
###	 library(R.matlab)
###-----------------------------------------------------------------
 library(R.matlab)

## inputSpectList (either tofList OR tofResampled)

if (identical(choice, "ISP")){ 
	inputSpecList<- object2save
	spectraNames<-names(inputSpecList)
 	numberSpectra<-length(spectraNames)
 	origTimeCounts<-length(inputSpecList[[1]])
 	origTime <- c(1:origTimeCounts)
 	inputSpecArray <-array(0,c(numberSpectra,origTimeCounts))
 	rownames(inputSpecArray)<- spectraNames
 	colnames(inputSpecArray)<- origTime
    
 	for (ii in 1: numberSpectra) {
 		tofTemp<-as.matrix(inputSpecList[[ii]])
 		inputSpecArray[ii,]<-tofTemp
 	}
 	filenameTofList <- paste("inputSpecListForML", ".mat", sep="")
 	filenameOrigTime<- paste("origTimeForML", ".mat", sep="")
 	writeMat(filenameTofList, inputSpecArray=inputSpecArray)
 	writeMat(filenameOrigTime, origTime=origTime)
	print("saved inputSpecListForML.mat and origTimeForML.mat")
 	rm(inputSpecList, tofTemp, inputSpecArray, origTime); gc();
}

## tofResampledAll

if (identical(choice, "TRS")){
	tofResampled <- object2save
	spectraNames<-names(tofResampled)
 	numberSpectra<-length(spectraNames)
 	rsTimeCounts<- length(tofResampled[[1]])#length(rsout$resTime)
 	tofResampledAll <-array(0,c(numberSpectra,rsTimeCounts))
 	rownames(tofResampledAll)<- spectraNames
 	#colnames(tofResampledAll)<- rsout$resTime

 	for (ii in 1: numberSpectra) {
 		tofTemp<-as.matrix(tofResampled[[ii]])
 		tofResampledAll [ii,]<-tofTemp
 	}
 
 	#resampledTime<-rsout$resTime
 	#resampledRate<-rsout$resRate
 	filenameTofResampled <- paste("tofResampledForML", ".mat", sep="")
 	#filenameResTime<- paste("resTimeForML", ".mat", sep="")
 	writeMat(filenameTofResampled, tofResampledAll=tofResampledAll)
 	#writeMat(filenameResTime, resampledTime=resampledTime)
	rm(tofResampled, tofResampledAll); gc();
	print("saved tofResampledForML.mat")
}

## saving peakList

if (identical(choice, "PL")){
	peakList <- object2save
  	numberOfSpectrum<- length(peakList)
  	spectraName <- names(peakList)
  	peakCt <- array(0, numberOfSpectrum)

	# Find number of peaks in each sample
  	for (k in 1:numberOfSpectrum) {
		peakData <- peakList[[spectraName[k]]]
		peakCt[k] <- length(peakData$Positions)
  	}
  	maxPeaks <- max(peakCt);
  	peakPositionArray <- array(0, dim = c(numberOfSpectrum,maxPeaks))
  	peakIntensityArray <- array(0, dim = c(numberOfSpectrum,maxPeaks))

  	for (kk in 1:numberOfSpectrum){
		peakData <- peakList[[spectraName[kk]]]
		peakPositionArray[kk,1:peakCt[kk]] <- peakData$Positions
		peakIntensityArray[kk,1:peakCt[kk]] <- peakData$Intensities
	}

 	filenamePeakListPositions <- paste("peaksForML", ".mat", sep="")
 	filenamePeakListIntensities<- paste("peakIntForML", ".mat", sep="")
 	writeMat(filenamePeakListPositions, peakPositionArray=peakPositionArray)
 	writeMat(filenamePeakListIntensities, peakIntensityArray=peakIntensityArray)
	print("saved peaksForML.mat and peakIntForML.mat")
	rm(peakData, peakCt, peakPositionArray, peakIntensityArray); gc();	
  }

## saving alignedPeakList

if (identical(choice, "APL")){
	alignedPeakList <- object2save
 	alignedPeaks<-alignedPeakList$peaks
 	numberPeaks<-length(alignedPeaks)
 	spectraNames<-names(alignedPeakList$data)
 	numberSpectra<-length(spectraNames)
 	alignedPeakIntArray <-array(0,c(numberSpectra,numberPeaks))
	 rownames(alignedPeakIntArray)<- spectraNames

	 colnames(alignedPeakIntArray)<- alignedPeaks
 	for (ii in 1: numberSpectra) {
		alignedPeakData<-alignedPeakList$data[[spectraNames[ii]]]
 		intensities<-alignedPeakData[["Intensities"]]
 		aplTemp<-as.matrix(intensities)
 		alignedPeakIntArray[ii,]<-aplTemp
 	}
 
 	filenameAlignedPeaks <- paste("alignedPeaksForML", ".mat", sep="")
 	filenameAlignedPeakInt <- paste("alignedPeakIntForML", ".mat", sep="")
 	writeMat(filenameAlignedPeaks, alignedPeaks = alignedPeaks)
 	writeMat(filenameAlignedPeakInt, alignedPeakIntArray = alignedPeakIntArray)
	print("saved alignedPeaksForML.mat and alignedPeakIntForML.mat")
	rm(alignedPeaks, alignedPeakData, intensities,aplTemp, alignedPeakIntArray); gc();	
}

}

