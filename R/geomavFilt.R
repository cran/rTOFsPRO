"geomavFilt" <-
function(sigIN, trct, fc3list){

wfc1 <- fc3list$fc1; wfc2 <- fc3list$fc2; wfc3 <- fc3list$fc3;

Nw<-max(c(length(wfc1),length(wfc2),length(wfc3)))

fs1 <- filter(sigIN, wfc1, method = "convolution",sides = 2, circular = TRUE);

# for symmetric (around maximum) filters, no shift is needed

mxlag <-floor(length(fs1)/7);
cf<-ccf(sigIN[(3*mxlag+1):(4*mxlag)],as.vector(fs1[(3*mxlag+1):(4*mxlag)]),lag.max=mxlag, type="correlation", plot=FALSE);
sh1<-cf$lag[cf$acf==max(cf$acf)]-1;

#Nw1<-length(wfc1);
#if (Nw1/2==trunc(Nw1/2)) {sh1 <- -Nw1/2;} # for symmetric (around maximum) filters
#else {sh1 <- -(Nw1-1)/2;}
#print(sh1)
fs1 <- circShift(fs1, sh1);

fs2 <- filter(sigIN, wfc2, method = "convolution",sides = 2, circular = TRUE);

# for symmetric (around maximum) filters, no shift is needed

mxlag <-floor(length(fs2)/7);
cf<-ccf(sigIN[(3*mxlag+1):(4*mxlag)],as.vector(fs2[(3*mxlag+1):(4*mxlag)]),lag.max=mxlag, type="correlation", plot=FALSE);
sh1<-cf$lag[cf$acf==max(cf$acf)]-1;
#print(sh1)
fs2 <- circShift(fs2, sh1);

fs3 <- filter(sigIN, wfc3, method = "convolution",sides = 2, circular = TRUE);
# for symmetric (around maximum) filters, no shift is needed

mxlag <-floor(length(fs3)/7);
cf<-ccf(sigIN[(3*mxlag+1):(4*mxlag)],as.vector(fs3[(3*mxlag+1):(4*mxlag)]),lag.max=mxlag, type="correlation", plot=FALSE);
sh1<-cf$lag[cf$acf==max(cf$acf)]-1;
#print(sh1)
fs3 <- circShift(fs3, sh1);

#print(length(fs1))
#print(length(fs2))
#print(length(fs3))

fs<-sign(fs1*fs2*fs3)*((abs(fs1*fs2*fs3))^(1/3.));
i1<-which((fs1*fs2<0) | (fs1*fs3<0) | (fs2*fs3<0) | (fs*fs1<0) | (fs*fs2<0) | (fs*fs3<0)| (fs <0) | (fs1 < 0) | (fs2 < 0) | (fs3 < 0));
fs[i1] <- fs[i1]/(max(fs)*trct);

l1<-length(fs)
fs[(length(fs)-floor(Nw/2)):length(fs)]<-sigIN[(length(sigIN)-floor(Nw/2)):length(sigIN)]; ## correct for filter shift into early points 11/14/07
fs <- fs[1:l1];
sigOUT <- fs;
rm(fs, fs1, fs2, fs3, sigIN);
return(as.vector(sigOUT));

} ###### END of "geomavFilt" #####

