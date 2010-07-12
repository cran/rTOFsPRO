"myxcorr" <-
function(x, y, maxlag) {

###----------------------------------------------------
### Compute cross-correlation like Matlab.
### Signals are zero-padded to the next (uuper) close power of 2 for the
### longest of the two
###
### INPUT:
### x - observed sampled discrete signal;
### y - observed sampled discrete signal;
### maxlag - maximum tim-lag for correlation;
###
### OUTPUT:
### xcf - cross-correlation function on (-maxlag:maxlag) interval, 2*maxlag+1 length;
###

Nx<-force(length(x))
Ny<-force(length(y))
pow2<-ceiling(log2(max(c(force(2*Nx-1),force(2*Ny-1)))))
xp<-c(x,rep(0,2^pow2-Nx));
yp<-c(y,rep(0,2^pow2-Ny));

Fx<-fft(xp, inverse=FALSE)
Fy<-fft(yp, inverse=FALSE)
xc<-Re(fft(Fx*Conj(Fy), inverse=TRUE))
xcf<-c(xc[force(length(xc)-maxlag+1):force(length(xc))],xc[1:force(maxlag+1)]);

return(as.vector(xcf))
} ###### END of "myxcorr" #####

