"gausRHalf" <-
function (gt, trnct){

sigma <- gt/sqrt(2*log(2)); ## gt is in time clicks, normalize Gaussian sigma 11/02/07
gt <- ceiling(gt);
lts <- ceiling(30*gt+2*sqrt(2*log(trnct))*sigma); ## wavelet length with given precision 12/03/07
ts <- as.vector(1:lts); 
tmax <- 20*gt;

gfxr <- dnorm(ts, tmax, gt/sqrt(2*log(2))); 
gfxr <- gfxr/max(gfxr);
#print(hwfmGausPeak(gfx,tmax))
gfxr[which(gfxr < (1.0/trnct))] <- 0;
gfxr <- gfxr[tmax+1:length(gfxr)];
i1 <- which(gfxr>0);
gfxr <- gfxr[1:(i1[length(i1)]+6*gt)];

return(as.vector(gfxr));
} ###### END of "gausRHalf" #####

