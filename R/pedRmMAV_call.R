"pedRmMAV_call" <-
function(tofList,pw,nf){

###------------------------------------------------------------------ 
### Function pedRmMAV_call manages pedestal removal 
###  It performs the following:
###	Estimates the baseline using local minima
###	Calls mavSmoothing to smooth baseline estimate
###	Subtracts smooth baseline estimate from signal
###
###  INPUT:
###  tofList
### 	pw - peak half-width (used for window width for MAV ~6*HWFM is a 
###       good start) Assumes that peaks are the same width – valid for 
###       resampled data
###	nf – number/fraction of width in window
###
###  OUTPUT:
### 	tofList with pedestals removed
###
###  DEPENDENCIES:
###	mavSmoothing function (from "LibIDSfilts")
###
###------------------------------------------------------------------
ww<-floor(nf*pw)
#print(c(ww,pw))
numSpec<-length(tofList)
lastPt<-length(tofList[[1]])
if (ww < 2) { ww <- 2 }

for (ispec in 1:numSpec) {

sigt<-tofList[[ispec]]
bs0<-array(0,(lastPt-ww-ww))
 for (i in (ww+1):(lastPt-ww)) {
	bs0[i-ww]<-min(sigt[(i-ww):(i+ww)])
  }

zeroVector<-array(0,ww)
bs0<-c(zeroVector+min(sigt[1:ww]),bs0,zeroVector+min(sigt[(lastPt-ww):lastPt]))

bst <- mavSmoothing(bs0,ww)
sigpr <- sigt-bst

#sigpr[1:(ww)]<-0 # Set first ww points to zero if desired

tofList[[ispec]]<-sigpr
}
return(tofList)
}

