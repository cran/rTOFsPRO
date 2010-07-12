"lorLHalf" <-
function(gt, trnct){

tau <- 1/gt; # (sqrt(2)-1)/gt; ## gt is in time clicks, corrected 01/09/08
gt <- ceiling(gt);
lts  <- ceiling(10*gt+2*sqrt(trnct-1)/tau); ## wavelet length with given precision 12/03/07
ts <- as.vector(1:lts); 
tmax <- length(ts)-6*gt;


lfxl <- tau/(1+((ts-tmax)*tau)^2);
lfxl <- lfxl/max(lfxl);

#print(hwfmGausPeak(lfxl,tmax))

lfxl[which(lfxl<1.0/trnct)]=0;
lfxl <- lfxl[1:tmax];

i1 <- which(lfxl==0);
lfxl <- lfxl[(i1[length(i1)]-2*gt):length(lfxl)];

return(as.vector(lfxl));
} ###### END of "lorLHalf" #####

