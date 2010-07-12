"matchedFilter" <-
function(s, h) {

###----------------------------------------------------
### Produce matched filter smoothing for the input signal
###
### INPUT:
### s - matrix of observed sampled discrete signals;
### wvlt - pure signal wavelet (best fit);
###
### OUTPUT:
### outs - matched filtered signal; 
###

## filter(x, filter, method = c("convolution", "recursive"),
## sides = 2, circular = FALSE)
Ns<-length(s)

Nw<-length(h)

sf <- filter(s, h, method = "convolution",sides = 2, circular = TRUE);

outs<-s;


outs<-sf*sum(s)/sum(sf);

mxlag <-floor(length(outs)/7);
cf<-ccf(s[(3*mxlag+1):(4*mxlag)],as.vector(outs[(3*mxlag+1):(4*mxlag)]),lag.max=mxlag, type="correlation", plot=FALSE);
sh1<-cf$lag[cf$acf==max(cf$acf)]-1;
#print(sh1)


outs <- circShift(outs, sh1);
l1 <- length(outs);
outs[(length(outs)-floor(Nw/3)):length(outs)] <- s[(length(s)-floor(Nw/3)):length(s)] ## correct for filter shift into early points 11/14/07
outs <- outs[1:l1];

rm(sf, s) 
return(as.vector(outs))
} ###### END of "matchedFilter" #####

