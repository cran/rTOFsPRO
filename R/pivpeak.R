"pivpeak" <-
function(sigt, thres) {

###----------------------------------------------------
### Find time position of the pivot peaks above threshold 
### input signal
###
### INPUT:
### sigt - observed sampled discrete signal;
### thres - intensity threshold for pivot peaks;
###
### OUTPUT:
### pivots - time point numbers for pivots;
###

qin<-1:length(sigt);
k<-qin[sigt>thres];
x<-diff(sigt[k],1,1)
xl<-x[1:length(x)-1]
xr<-x[2:length(x)]
p<-2:length(k)-1;
p1<-p[(xl>0 & xr<=0) | (xl==0 & xr<0)]
pivots<-k[p1]

return(as.vector(pivots))

} ##### END "pivpeak" #####

