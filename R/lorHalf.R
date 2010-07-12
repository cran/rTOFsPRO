"lorHalf" <-
function(gt, trnct){

tau <- 1/gt; # (sqrt(2)-1)/gt;  ## gt is in time clicks, corrected 01/09/08
lts <- ceiling(10*gt+2*sqrt(trnct-1)/tau); ## wavelet length with input precision 12/03/07
gt <- ceiling(gt);
ts <- as.vector(1:lts); 
tmax <- 6*gt;

lfx <- tau/(1+((ts-tmax)*tau)^2);
lfx <- lfx/max(lfx);

#print(hwfmGausPeak(lfx,tmax))

lfx[which(lfx<1.0/trnct)]=0;
lfx <- lfx[tmax+1:length(lfx)];

i1 <- which(lfx>0);
lfx <- lfx[1:(i1[length(i1)]+2*gt)];

return(as.vector(lfx));
} ###### END of "lorHalf" #####

