"optFilt" <-
function(inSig, wfc) {

Nw <- length(wfc);

fs <- filter(inSig, wfc, method = "convolution",sides = 2, circular = TRUE);
# for symmetric (around maximum) filters, no shift is needed

mxlag <-floor(length(fs)/7);
cf<-ccf(inSig[(3*mxlag+1):(4*mxlag)],as.vector(fs[(3*mxlag+1):(4*mxlag)]),lag.max=mxlag, type="correlation", plot=FALSE);
sh1<-cf$lag[cf$acf==max(cf$acf)]-1;

fs <- circShift(fs, sh1);

l1 <- length(fs);
fs[(length(fs)-floor(Nw/2)):length(fs)] <- inSig[(length(inSig)-floor(Nw/2)):length(inSig)] ## correct for filter shift into early points 11/14/07
fs <- fs[1:l1];

return(as.vector(fs));
} ###### END of "optFilt" #####

