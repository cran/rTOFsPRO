"peakAlignAv" <-
function(peaks,tofList,noiseL,HWFM){		

###----------------------------------------------------------------- 
### Function peakalignAv aligns peaks from tofList spectra with peaks found
###   in averaged spectrum 
###
### INPUT:
### 	peaks - peak positions found in the averaged spectrum;
###	tofList - list of TOF spectra in the same domain (down-sampled or original) as "peaks"
### 	HWFM - peak half-width (window for alignment) -- constant for down-sampled data;
###
### OUTPUT:
### 	alignedPeakList
###
###------------------------------------------------------------------

  alignedPeakList <- list();
  alignedPeakList$peaks <- peaks;
  alignedPeakList$data <- list();

  specName<-names(tofList)
  numSpectra<-length(tofList)
  numPeaks<- length(peaks) 
  peakInt<-array(0,numPeaks) 
  peakShifts<-peakInt 

  zeroVector<-rep(0,numPeaks)

## align

for (i in 1:numSpectra) {
    spectra <- tofList[[specName[i]]]
    nam <- noiseL[[i]]

    for (j in 1:numPeaks) { 	
        w <- spectra[(peaks[j]-HWFM):(peaks[j]+HWFM)];      
       mw<-max(w) 
      imw <- which.max(w)
      x <- diff(w,1,1);		# first difference
	xEnd<-length(x)
       xl <- x[1:xEnd-1]	# left side
       xr <- x[2:xEnd]		# right side
       neartop <- (xl > 0 & xr <= 0) | (xl == 0 & xr < 0)     
      xx <- 1:xEnd
      nt <- xx[neartop]+1 		# tentative maximum     
      if(length(nt)>0) {       
          if(length(nt)>1) { 
            y<-max(w[nt])
            iy<-which.max(w[nt])
           peakInt[j]<- y  
           peakShifts[j]<-nt[iy]-HWFM-1 
           }
          else {   
              peakInt[j]<-w[nt]
              peakShifts[j]<-nt-HWFM-1                
          }    
     } else {
            if(abs(imw-HWFM-1) < HWFM) {  
             peakInt[j] <- mw 
                peakShifts[j] <- imw-HWFM-1 
          }  else {
                 peakInt[j]<-spectra[peaks[j]] 
                 peakShifts[j]<-0  
            } 
      }       
    }     
 
peakInt<-as.numeric(peakInt)
uncertInt <- nam[peaks+peakShifts] # local noise estimate at peak position

## Construct alignedPeakList

 	alignedPeakData <- list("Intensities" = as.vector( peakInt ),
  		"PositionShifts" = as.vector(peakShifts),
  			"PositionUncertainties" = zeroVector + max(abs(peakShifts)),
  				"IntensityUncertainties" = uncertInt )
  
  	alignedPeakList$data[i] <- list( alignedPeakData ); # allow assignment for repeated "names"
  }
  names(alignedPeakList$data) <- specName;
  return(alignedPeakList)
}

