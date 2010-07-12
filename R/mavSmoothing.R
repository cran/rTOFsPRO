"mavSmoothing" <-
function(inputSig,pointWindowLength) {

###----------------------------------------------------------------- 
### Function mavSmoothing uses moving average filter to smooth input (e.g., baseline,
### 	removing small scale variations and preserving larger scale "pedestals" 
###
### INPUT:
### 	inputSig - baseline estimated by local minima in spectra
### 	pointWindowLengtht - filter window
###	
### OUTPUT:
### 	bst - smoothed baseline
###------------------------------------------------------------------

 
   L <- pointWindowLength
   
   # create impulse response of the MAV filter;
   h<-1/L*rep(1,L);
   k1<-floor(L/2)
 
   ## get the input info
     inputSpectrum <- inputSig
     Nc <- length(inputSpectrum)
     outSpectrum<-filter(inputSpectrum,h, method="convolution", sides=2, circular=FALSE);
     if ((L/2) != k1){
       outSpectrum<-c(inputSpectrum[1:k1],outSpectrum[(k1+1):(Nc-k1)],inputSpectrum[(Nc-k1+1):Nc]);
     }
     else{
       outSpectrum<-c(inputSpectrum[1:(k1-1)],outSpectrum[k1:(Nc-k1)],inputSpectrum[(Nc-k1+1):Nc]);
     }
     bst <- outSpectrum
   
   return(bst)
 } ###### END of "mavSmoothing" #####

