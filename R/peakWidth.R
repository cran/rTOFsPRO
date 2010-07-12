"peakWidth" <-
function(s,startsd, nsubdiv, sfr){
###------------------------------------------ 
### OUTPUT:
###
### imxs - point of the maximum, and estimated width (two columns);
###
### INPUT:
### s - signal equally sampled in time;
### startsd - starting point for subdivision;
### nsubdiv - number of subdivisions for peak width estimate;
### sfr - signal fraction (1/sfr) of maximum to include
###-----------------------------------------------  

Nr <- length(s);
Nc <- ncol(s);
tolrnc <- 0.001
drange_thr  <- sfr #20; ## 1/fraction of maximum signal to include in fit

#source("pivpeak.R");
#source("hwfmGausPeak.R");

psdiv <- floor((Nr-startsd)/nsubdiv);
#print(psdiv); # debug
mxs <- 0; sgmmxs <- 0;
for (k in 1:(nsubdiv-1)){
#startsd+1+(k-1)*psdiv
#startsd+k*psdiv
ss <- s[(startsd+1+(k-1)*psdiv):(startsd+k*psdiv)];
subpoint <- 1:length(ss)
subpeaks <- pivpeak(ss,tolrnc)
#print(k)
mx <- max(ss[subpeaks])
imx <- intersect(subpoint[ss==mx], subpeaks);
if (length(imx)>1){
imx<-imx[length(imx)]
}
pointmx <-imx+startsd+1+(k-1)*psdiv;
pstart <- pointmx-psdiv; pend <- pointmx+psdiv;

sgm <-0;
if ((pstart >=1) & (pend <= length(s)-1)){ 
sgm <- hwfmGausPeak(s[pstart:pend], psdiv+1);
}
if ((sgm > 0) & (!is.element(pointmx, mxs))){
  mxs <- c(mxs,pointmx);
sgmmxs <- c(sgmmxs,sgm);
}
}

mxs <- mxs[2:length(mxs)]; sgmmxs <- sgmmxs[2:length(sgmmxs)];

imxs<- cbind(matrix(mxs, ncol=1),matrix(sgmmxs, ncol=1));

mmx <- max(s[imxs[,1]]);
immx <- which(s[imxs[,1]]==mmx);

## look only for high signals

r1 <- imxs[,1]-imxs[,2]; r2<-imxs[,1]+imxs[,2];
maxsp<-r1*0;
for (k in 1:length(r1)) {
maxsp[k]<-max(s[r1[k]:r2[k]]);
}

kimx <- which((s[imxs[,1]]*imxs[,2]>mmx*imxs[immx,2]/drange_thr) & (s[imxs[,1]]> 0.5*maxsp[]))


#& (s[imxs[,1]]> 0.5*max(s[(imxs[,1]-imxs[,2]):(imxs[,1]+imxs[,2])])));

mxs <- mxs[kimx]; sgmmxs <- sgmmxs[kimx];
imxs <- imxs[kimx,];


t <- 1:Nr;
## PLot the results ##
#plot(t, s*max(sgmmxs)/max(s),type='l');
##lines(t(mxs), sgmmxs, type='p', col=3);
#points(t(mxs), sgmmxs, col=3)

return(as.matrix(imxs));
} ##### END "peakWidth" #####

