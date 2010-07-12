"trgFilt" <-
function(sWvlt, nu, dWvlt) {

###
### INPUT:
### sWvlt - sampled discrete signal wavelet;
### nu - noise power ratio to signal power;
### dWvlt - desired signal wavelet shifted in respect to input;
###
### OUTPUT:
### wf1 - Wiener filter coefficients in time domain;
###

M <- length(sWvlt);
xc <- myxcorr(as.double(sWvlt), as.double(sWvlt), M-1);

R <- toeplitz(xc[M:length(xc)]);
Rn <- sum(as.double(xc[M:length(xc)]))*diag(1,nrow(R));
R <- R+nu*Rn;

p <- myxcorr(as.double(dWvlt), as.double(sWvlt), M-1);
p <- matrix(as.double(p[M:(2*M-1)]), ncol=1);

# Normalize equation by r0
p <- p/R[1,1]; 
R <- R/R[1,1];

wf1<-solve(R,p);
return(as.vector(wf1));
} ###### END of "trgFilt" #####

