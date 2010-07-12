"gausHalf" <-
function (gt, trnct){

sigma <- gt/sqrt(2*log(2)); ## gt is in time clicks, normalize sigma for Gaussian HWFM 11/02/07
gt<- ceiling(gt);
lts <- ceiling(30*gt+2*sqrt(2*log(trnct))*sigma);
ts <- as.vector(1:lts); 
tmax <- length(ts)-20*gt;

gfx <- dnorm(ts, tmax, sigma); 
gfx <- gfx/max(gfx);
#print(hwfmGausPeak(gfx,tmax))
gfx[which(gfx < (1.0/trnct))] <- 0;
gfx <- gfx[1:tmax];
i1 <- which(gfx==0);
gfx <- gfx[(i1[length(i1)]-6*gt):length(gfx)];

return(as.vector(gfx));
}  ###### END of "gausHalf" #####

