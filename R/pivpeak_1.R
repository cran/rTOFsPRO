"pivPeak" <-
function(refsig, thres, phw) {

###----------------------------------------------------
### Find time position of the pivot peaks above threshold 
### input signal
###
### INPUT:
### refsig - observed sampled discrete reference signal nromalized to its maximum = 1;
### thres - intensity threshold for pivot peaks in respect to maximum (e.g. 0.1);
### phw - peak half-width in the pivot detection range
###
### OUTPUT:
### pivots - time point numbers for pivots;
###

qin<-1:length(refsig);
k<-qin[refsig>thres];
x<-diff(refsig[k],1,1)
xl<-x[1:length(x)-1]
xr<-x[2:length(x)]
p<-2:length(k)-1;
p1<-p[(xl>0 & xr<=0) | (xl==0 & xr<0)]
pivots<-k[p1]
n1 <- which(diff(pivots)<phw);
loopn <- 0;
while (length(n1)>0 & loopn < 10){  # combine pivots closer than HW (to get 1/peak)
	loopn <- loopn + 1;
	#print(loopn)
	pcomp <- cbind(pivots[n1],pivots[n1+1]);
	comp <- c(pivots[n1],pivots[n1+1])
	dpiv <- setdiff(pivots,comp)
	pcomp1 <- pcomp[,1];
	d1 <- dim(pcomp)[1]
	for (jj in 1:d1){  # select pivot of max signal for closer than HW pivots
    	pcomp1[jj] <- pcomp[jj,which(refsig[pcomp[jj,]]==max(refsig[pcomp[jj,]]))]
    	}
	pivots <- sort(unique(c(pcomp1,dpiv)))
	n1 <- which(diff(pivots)<phw);
	rm(pcomp,comp,dpiv, pcomp1,d1,jj)
}
rm(qin,k,x,xl,xr,p,p1,n1,loopn)
return(as.vector(pivots))

}

