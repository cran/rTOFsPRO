"matchFiltCoefGG" <-
function(inhw, las, trct) {

wvlt <- c(gausHalf(inhw, trct), gausRHalf((inhw+las), trct));
Nw <- length(wvlt);

h<-wvlt[sort(c(1:Nw),decreasing=TRUE)]; # create impulse response of the MATCHED filter;

maxF <- which(h==max(h));

if (maxF>Nw/2){h<-c(h,rep(0,(2*maxF+1-Nw)))} ## symmetrize the filter
else {h <- c(rep(0,(2*maxF+1-Nw)),h)}


return(as.vector(h));
} ###### END of "matchFiltCoefGG" #####

