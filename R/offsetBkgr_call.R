"offsetBkgr_call" <-
function(tofList, offsetSamples) 
{

##----------------------------------------------------
## Calculate baseline from ADC offset
##
## INPUT:
## tofList - list of the observed mass spectra (raw TOF)
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
# tofList_noBS <- chargeBkgr_call(tofList,offsetSamples);
# Dependency: none
#

spectraName <- names(tofList)
spectraCount <- length(tofList)

for (j in 1:spectraCount) {
    spectrum <- tofList[[j]]
    lensp <- length(spectrum)
    basecon <- mean(spectrum[offsetSamples]);
    spectrum <- spectrum-basecon; ## subtract baseline
    tofList[spectraName[j]] <- list(spectrum); ## update the tofList record
}

return(tofList);
}

