"globalAlign" <-
function(tof,tofList, tofListMetaData, HWFM,noiseThres,startS,endS, Appli, Detect) {

### INPUT:
### tof - baseline subtracted and smoothed tofList (to detect global shifts)
### tofList - unprocessed data (to apply detected shifts)
### tofListMetaData - "offset" and "scale" are used and /or updated
### HWFM - peak width in the range of pivot detection
### noiseThres - fraction of maximum signal threshold for pivot detection
### startS, endS - TOF range (indices) for pivot detection
### Appli - apply correction to the tofList
### Detect - auto-detect correction for the data

Nr <- length(tofList)
spectraName <- names(tofList)
Nc <- 0
for(i in 1:Nr) {
  	len <- length(tof[[i]])
  	if (len > Nc) {
    	Nc <- len
  	}
}

sigtm <- matrix(0, nrow=Nr, ncol=Nc)
sigt0 <- sigtm

digitizerDwell <- as.double(tofListMetaData[1,"instrumentSpecificSettings.timeDelta"]) ;
for (i in 1: Nr) {
    	sigt <- tof[[i]]
    	sig0 <- tofList[[i]]
    	l1 <- length(sigt)
    	l0 <- length(sig0)
    	#print(c(l1,l0)) 
    	sigtm[i, 1:l1] <- sigt;
    	sigt0[i, 1:l0] <- sig0;
    	sigtm[i, is.na(sigt)]<-0  ## zero-fill for NA values 06/22, Dasha
    	sigt0[i, is.na(sig0)]<-0 
	if ( digitizerDwell != as.double(tofListMetaData[ i, "instrumentSpecificSettings.timeDelta" ]) )
	{
		userError( "Spectra do not all have the same digitizer frequency.  Please choose spectra that all have the same experimental settings." );
	}
}

if(Detect){
	sw <- 2 * HWFM
	if ( sw < 2 )
	{
		userError( "The 'Peak Half-Width' must be at least one time tick." );
	}

	if ( noiseThres <= 0 || noiseThres >= 1 )
	{
		userError( "'Threshold Intensity' must be between 0 and 1." );
	}

	### specify the time range to look for pivot peaks

	### try to transfer the list to a matrix and error out if the length are different

	
	dataStart<-as.integer(tofListMetaData[1, "timeOfFlightData.start" ]);
	if ( startS < dataStart )
	{
		userError( "'Start Time' must be between 'dataStart' and the longest spectrum length." );
	}

	if ( endS > Nc )
	{
		userError( "'End Time' must be less than the longest spectrum length." );
	}

	if ( endS <= ( startS + 3 * sw ) )
	{
		userError( "The specified time range ('End Time' - 'Start Time') must be greater than 18 peak half-widths in order to contain at least 3 pivot peaks in the reference spectrum." );
	}

	pivPeaks<- pivPeak(sigtm[1,startS:endS]/max(sigtm[1,startS:endS]), noiseThres, HWFM);
	pivPeaks <- pivPeaks + startS;
	print("Found pivot peaks for global alignment:")
	print(length(pivPeaks))
	#print(pivPeaks)
	if ( length( pivPeaks ) < 3 )
	{
		userError( "Could not find at least 3 peaks in the reference spectrum." );
	}

	locJit <-locjitdet(t(sigtm), pivPeaks, sw, startS, endS, names(tofList));

	pj <- locJit$pjit;  ### pj[i,j] shows jitter for peak "i" in signal "j"
	#print(pj) #debug
	#pjSlope <- rep(0,dim(pj)[2]); pjIntercept <- pjSlope;
} else { pivPeaks <- c(1,2)  # no auto-detection is performed
	pj <-rbind(rep(0,Nr),rep(0,Nr))
}

	t0 <- 1:Nc ## original time axis
	sapl <- 0;

	for (j in 2:(dim(pj)[2])){ ## calculate linear fit coefficients for jitter in "j"-th record
 		ldf <- length(which(pj[,j]!=0)); lndf <- length(pivPeaks)-ldf; 
 		if ( (ldf > 3) & ((ldf-lndf) > 2)) { ## enough data to fit to a line
			z<-line(pivPeaks, pj[,j]);
			pjIntercept<-coef(z)[1];
			pjSlope<-(1.0+coef(z)[2]);
		} else { pjSlope <- as.numeric(tofListMetaData[ j, "timeOfFlightData.scale" ]); 
			pjIntercept <- as.numeric(tofListMetaData[ j, "timeOfFlightData.offset" ]);
		} ## input meta-data values
		#print(c(pjIntercept, pjSlope))  #debug

		if( Detect & !Appli){
			tofListMetaData[ j, "timeOfFlightData.offset" ] <- pjIntercept
			tofListMetaData[ j, "timeOfFlightData.scale" ] <- pjSlope
		}

		if ( Appli ) { ## apply new time scale
        		if (pjSlope!=1 | pjIntercept!=0) {
				tsj<-pjSlope*t0+pjIntercept
				intrp <- approx(tsj, sigt0[j,], t0)  ## appli to "original" (no BS) spectrum
				outSpectrum <-as.vector(intrp$y);
				inputSpectrum <- tofList[[j]]
				#print(c(length(inputSpectrum),length(outSpectrum))) #debug
				outSpectrum <- outSpectrum[1:length(inputSpectrum)]
				outSpectrum[is.na(outSpectrum)] <- inputSpectrum[is.na(outSpectrum)]
				tofList[spectraName[j]] <- list(outSpectrum)
        			sapl <- c(sapl,j);
        		}
        		tofListMetaData[ j, "timeOfFlightData.offset" ] <- pjIntercept
			tofListMetaData[ j, "timeOfFlightData.scale" ] <- pjSlope
        
		} 
 		#print(tofListMetaData[length(tofList), "timeOfFlightData.offset" ])
	}
	if (length(sapl)>1){
		sapl <- sapl[2:length(sapl)]; # spectra with applied global shift 
		print("####    Pre-aligned spectra:")
		print(sapl)
	} else {
		print("####    None of the spectra needed/required pre-alignment    ####")
	}
aOut <- list(tList=tofList, tMeta=tofListMetaData)
rm(sigtm, sigt0, sigt, sig0)
return(as.list(aOut))
}

