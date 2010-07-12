"hwfmGausPeak" <-
function(nf,n_mu){
###---------------------------------------------------------
###   Output is a standard deviation (HWFM) for a peak 
### assuming normal distribution fit.
###   SIGMA <- hwfmGausPeak(NF,N_MU) Returns standard
###   deviation SIGMA for a peak at the position N_MU according 
###   to normal pdf with mean, N_MU, for the values in NF.
###   NF is a column vector of signal values
###   For assymetric signals, SIGMA is the left half width. 
###------------------------------------------------------------


fmax<-nf[n_mu];
fhalf<-0.5*fmax;
lhalf<-1;
sigma <- -1;
#print(n_mu)
#print(length(nf))
if ((length(nf)>0) & (n_mu>lhalf)){
while ((n_mu-lhalf >1) & (nf[n_mu-lhalf] > fhalf)) {
  lhalf<-lhalf+1;
#print(n_mu-lhalf)
}
#lhalf
if ((n_mu-lhalf>=1) & (nf[n_mu-lhalf]<=fhalf)){
 
  sigma<-lhalf #/sqrt(2*log(2)); # without Gaussian N-factor, sigma in time-clicks
}
else {
  #disp('Could not find half maximum: input array is too narrow');
  #disp('at the peak'); n_mu
  sigma <- -1;
} 
}
return(sigma)
} ##### END "hwfmGausPeak" #####

