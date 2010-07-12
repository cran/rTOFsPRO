"chargeBkgr_call" <-
function(tofList, offsetSamples, amp, decay) 
{

##----------------------------------------------------
## Calculate baseline from charge accumulation
##
## INPUT:
## tofList - list of the observed mass spectra (raw TOF)
## amp - amplitude fraction of accumulated signal
## decay - decay constant for the accumulated signal -- exp(-1/tau)
## offsetSamples –- TOF interval for the estimate of constant offset, e.g.:
##                 c(120000:122000) 
##
## OUTPUT:
## tofList - baseline subtracted lsi of mass spectra
##
## NOTES:
## --1-- Created by Prof.Cooke (wecook@wm.edu) 10/10/04
## --2-- No self-optimization code is included yet, but defaults work
## Slow with long spectra (> 50K points). Faster version with C-function 
## is available
#
# USAGE:
# tofList_noBS <- chargeBkgr_call(tofList,offsetSamples, amp, decay);
# Dependency: none
#

spectraName <- names(tofList)
spectraCount <- length(tofList)

for (j in 1:spectraCount) {
    spectrum <- tofList[[j]]
    lensp <- length(spectrum)
    basecon <- mean(spectrum[offsetSamples]);
    ytemp <- rep(0,lensp);
#amp=0.0005; decay = 0.9988; basecon = 2114;
#amp=0.0007; decay = 0.9989; basecon = 2140;
    for (i in 2:lensp) {
	ytemp[i] <- (1-amp)*decay*(ytemp[i-1])+amp*(spectrum[i]-basecon);
    }
#print(length(y)) #debug
#print(length(ytemp))
    ybs <- ytemp + basecon;
    spectrum <- spectrum-ybs; ## subtract baseline
    tofList[spectraName[j]] <- list(spectrum); ## update the tofList record
}

return(tofList);
}

