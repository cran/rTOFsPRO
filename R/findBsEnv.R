"findBsEnv" <-
function(specDat, idsWindow){

t1 <- 1:length(specDat);

trs<-seq(from=1, to=length(t1), by=idsWindow);
k<-1; sigr<-0*trs;

ks <- (idsWindow-1)/2; ## assuming ODD idsWindow
sigr[1]<-sum(specDat[1:ks]); sigr[length(sigr)]<-sum(specDat[(trs[length(trs)]-ks):trs[length(trs)]]);

for (jk in seq(from=trs[2], to=trs[length(trs)-1], by=idsWindow)){
	k <- k+1; sigr[k] <- sum(specDat[(jk-ks):(jk+ks)]);
}

x<- diff(sigr);
xl<-x[1:(length(x)-1)]; xr<-x[2:length(x)];
nearbot <- (xl<0 & xr>=0) | (xl==0 & xr >0);
xx <- 1:length(x);
nb <- c(1, xx[nearbot]+1); ## tentative minimum
smax <- max(sigr[nb]); imax <- which(sigr[nb]==smax);
ts <- nb[imax];
#print(trs[ts]) ## debug
nxt <- imax+1;

while ((nxt <=length(nb)) & (sigr[nb[nxt]]>0.01*sigr[nb[imax]]) & (sigr[nb[nxt]]>1)) {
	if (sigr[nb[nxt]]<sigr[ts[length(ts)]]) {
		ts<-c(ts, nb[nxt]);
	}
	nxt<-nxt+1;
}
tSampl <- trs[ts]; sigSampl <- sigr[ts]/idsWindow;

envBS <- matrix(tSampl, nrow=length(tSampl), ncol=2);
envBS[,1]<-tSampl; envBS[,2]<-sigSampl;

return(as.matrix(envBS));
}

