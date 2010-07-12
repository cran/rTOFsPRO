"ExpBkgr_call" <-
function(tofList, bsSamples, offsetSamples, modelChoice, dataStart)
{

## Usage: Correct Background
## 
## Input Parameters:
## tofList – list of TOF spectra with meta-data
## bsSamples – TOF indices of data to use for model baseline fit 
##             (index is in original, non-truncated record), e.g.: 
##             c([t1:t2],[t3:t4],[t5:t6])
## When automated fit is preferred, provide two values in this field:
## c(Window, corPercent), where 
## idsWindow –- down-sampling step for finding baseline envelope, should be ODD
##           default for Bruker @ 500 MHz sampling rate is 199
## corPercent -– percent correction of baseline envelope for STD of noise
##               default for PCA’07 settings is 0.05 
## offsetSamples –- TOF interval for the estimate of constant offset, e.g.:
##                 c(120000:122000) 
## modelChoice --  choice of analytical baseline model:
##                 1 – linear, 2 – exponential 3 -- Gaussian
##

## use with new meta-data
# dataStart<-as.integer(tofListMetaData[1, "timeOfFlightData.start" ])

if (length(bsSamples)==2) {
	idsWindow <- bsSamples[1];
      corPercent <- bsSamples[2];
} else {
idsWindow <- 199; ## default for Bruker, 500 MHz sampling
corPercent <- 0.05; 
if ( bsSamples[1] <= dataStart )
{
	userError( "'Baseline Start Time' must be greater than 'dataStart' in dbSelector." );
}
}

# Get the input data info from tofList object
spectraName <- names(tofList)
spectraCount <- length(tofList)
#print("Spectra, modelChoice, coef(z)[2], coef(z)[3], coef(z)[1], cOffset")
for (j in 1:spectraCount) {
    spectrum <- tofList[[j]]
    timeLength <- length(spectrum)
	if ( timeLength <= offsetSamples[length(offsetSamples)]  )
	{
		userError( "Time range for offset calculation must be less than the record length." );
	}
        ## estimate and subtract baseline & offset

	t0 <- dataStart:(timeLength+dataStart-1);
	cOffset<-mean(spectrum[offsetSamples]);
	bsTrendModel <- 0*spectrum;

	#cSTD<-sd(spectrum[(timeLength-2000):(timeLength-100)]);
	if (length(bsSamples)<4) {
		#print(idsWindow) ## debug
		envSignl<-findBsEnv(spectrum-cOffset, idsWindow);
		tSamples<-envSignl[,1]; sigSamples<-envSignl[,2];
		#print(tSamples) ## debug
	} else {
		tSamples<- bsSamples ## locations of signal-free regions for baseline model
		sigSamples <- spectrum[tSamples-dataStart]-cOffset;
	}

	# negLength <- length(which(sigSamples <=0));
	if (modelChoice == 1) {
		z <- line(tSamples, sigSamples);
		bsTrendModel <- coef(z)[1]+t0*coef(z)[2];
	} else if (modelChoice == 2) {
		sigSamples[sigSamples < 1] <- 1;
		z<-line(tSamples, log(sigSamples));
		bsTrendModel<- exp(coef(z)[1])*exp(t0*coef(z)[2]);
	} else { ## Gaussian model
		sigSamples[sigSamples < 1] <- 1;
		a1 <- matrix(tSamples, nrow=length(tSamples), ncol=2);
		a1[,1] <- tSamples^2; a1[,2] <- tSamples;
		z<-lsfit(a1,log(sigSamples), intercept=TRUE);
		bsTrendModel<-exp(coef(z)[2]*t0*t0+coef(z)[3]*t0+coef(z)[1]);
	}

#print(c(spectraName[j], modelChoice, coef(z)[2], coef(z)[3], coef(z)[1], cOffset))
	bsTrendModel<-(1-corPercent)*bsTrendModel;

	spectrum<-spectrum-cOffset-bsTrendModel; ## subtract baseline

    ## update the tofList record
    	tofList[spectraName[j]] <- list(spectrum)
}

return(tofList);

}

