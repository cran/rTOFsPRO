"trivPP_call" <-
function(tofList, snr, naml, HWFM, fwl, stPk, dataStart)
{
###----------------------------------------------------
### Find time position of the pivot peaks above SNR 
### input signal
###
### INPUT:
### tofList -  list of observed sampled discrete signals (rows);
### snr - signal-to-noise ration for peak detection
### naml - noise amplitude list (for tofList);
### HWFM – half-width at half-maximum for the model signal wavelet
### fwl - filter wavelet length applied for signal processing
### stPk - starting point (tof-spectrum index) for peak detection (e.g., after matrix hump)
### dataStart - defines what is the absolute start location for the input time points
###
###
### OUTPUT:
### peakList - time point numbers for pivots;
###


spectraName <- names(tofList)

Ns <- length(tofList)

if ( snr <= 0 )
{
	userError( "Signal-to-Noise ratio must be positive." );
}


if ( HWFM < 1 )
{
	userError( "Invalid Peak Half-Width; it must be at least 1.0" );
}

pivLength=0;

## define the results list
peakList <- list()

for (i in 1:Ns){
 spec <- tofList[[i]]
 ## zero-fill for NA values
 spec[is.na(spec)] <- 0 
 peakTime <- spec * 0;
 peakAmps <- peakTime
 timeUncerts <- peakTime
 ampUncerts <- peakTime
 
 Nt <- length(spec)
 ampN <- as.numeric(naml[[i]]);
 stdNoise <- ampN/3;
	
     #print(c(length(ampN),Nt)) #debug
 if (Nt != length(ampN)) {
	userError("Noise amplitude vector is not the same length as signal")
 }
    	qin <- 1:Nt
	k<-qin; #qin[spec/ampN >= snr] ## use noise amplitude 11/08/07 << REV
	if ( length( k ) >= 3 )
	{     ## Step 1: find all local maxima and minima
		x<-diff(spec[k],1,1)
		xl<-x[1:length(x)-1]
		xr<-x[2:length(x)]
		nrtp <- (xl>0 & xr<=0) | (xl==0 & xr<0);
		nrbt<- (xl<0 & xr >=0) | (xl==0 & xr>0)
		xx <- 1:length(x);
		nb <- c(1, xx[nrbt]+1); ## tentative minimum
		nt <- xx[nrtp]+1; ## tentative maximum
		
		## Step 2: Coalesce rising plateaus to find the right-most one
		curt <- 1; curb <- 1; peak <- c(0); bots <- c(1); 
		numtop <- 1; numbot <- 1;
		while ((curt < length(nt)) & (curb < length(nb))) {
			while ((nb[curb]<nt[curt]) & (curb < length(nb))) {
			curb <- curb+1;
			}
		bots[numbot]<-nb[curb];
		numbot<-numbot+1;
			while ((nt[curt]<nb[curb]) & (curt < length(nt))) {
			curt <- curt+1;
			}
		peak[numtop]<-nt[curt];
		numtop <- numtop+1;
		}
		rm(numtop, numbot, curt, curb, nb, nt, nrbt, nrtp)
		leftMin <- c(1); rightMin <- c(1);
		for (i1 in 2:length(peak)) {
		temp <- peak[i1-1]:peak[i1];
		leftMin[i1] <- max(temp[spec[k[temp]]==min(spec[k[temp]])]);
		rightMin[i1] <- min(temp[spec[k[temp]]==min(spec[k[temp]])]);
		}
		rightMin <- c(rightMin[2:length(rightMin)], length(x));
 
		con1 <- abs(spec[k[leftMin]]-spec[k[peak]])/ampN[k[peak]]; ## << REV
		con2 <- abs(spec[k[rightMin]]-spec[k[peak]])/ampN[k[peak]]; ## << REV for ampN vector
		
		peak <- peak[(con1 >= snr) | (con2 >= snr)];
		peak <- peak[spec[k[peak]]/ampN[k[peak]] >= snr];

		if ( length( peak ) > 0 )
		{
			pivots<-k[peak]
		} else
		{
			pivots <- 0;
		}
	} else
	{
		if ( length( k ) > 0 )
		{
			pivots <- k[spec==max(spec[k])];
		} else
		{
			pivots <- 0;
		}
	}	
	## use only peaks no closer than 2*HWFM to the origin (to estimate uncertainties)
	pivots<-pivots[(pivots>2*HWFM)]   ## bug fix for 6086, Dasha 06/22
	pivots <- pivots[which(pivots>stPk & pivots < (length(spec)- ceiling(fwl/2)))]; # trustworthy range

	numberPeaks <- length(pivots)
#print(length(pivots)) #debug
	if ( numberPeaks > 0 )  ## REV << check that stdNoise is vector or add ampN/3?
	{
		peakTime[1:numberPeaks]<-pivots
#print(peakTime) #debug
		peakAmps[1:numberPeaks]<- spec[peakTime[1:numberPeaks]];
		ampUncerts[1:numberPeaks]<-peakTime[1:numberPeaks]*0+ampN[pivots];
		snrPeaks<-peakAmps[1:numberPeaks]/stdNoise[pivots];
	
		for (j in 1:numberPeaks){
			if ((snrPeaks[j]>2) & (peakTime[j]>2*HWFM)){
				y<-1.0/(1.0-2/snrPeaks[j]);
				tempd <- spec[(peakTime[j]-2*HWFM):peakTime[j]]
				timeUncerts[j]<-length(tempd[tempd>=tempd[length(tempd)]/y])
			}else{
				timeUncerts[j]<-HWFM;
			}
		}
	##print("assign_info_list")
		info <- list( 
			"Positions" = peakTime[1:numberPeaks]+dataStart-1,
			"Intensities" = peakAmps[1:numberPeaks],
			"PositionUncertainties" = timeUncerts[1:numberPeaks],
			"IntensityUncertainties" = ampUncerts[1:numberPeaks] )
#print(info$Positions) #debug
	} else
	{       print("no_peaks")
		info <- list( 
			"Positions" = NULL,
			"Intensities" = NULL,
			"PositionUncertainties" = NULL,
			"IntensityUncertainties" = NULL )
	}

	## append results to the list
  	peakList[i] <- list( info )  # allow assignment for repeated names 02/11/10
#print(peakList[[i]]$Positions) #debug

}
names(peakList) <- spectraName;
return(peakList);
#tofList <- NULL;
}

