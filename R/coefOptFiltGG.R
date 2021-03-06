"coefOptFiltGG" <-
function(inhw, las, trct) {

nu <- 0.01;
inpW <- c(gausHalf(inhw, trct), gausRHalf((inhw+las), trct));
trgW <- c(gausHalf(0.8*inhw, trct*1.3), gausRHalf(0.8*(inhw+las), trct*1.3));
trgW <- trgW*sum(inpW)/sum(trgW); # normalize target wavelet ??
inpP <- which(inpW == max(inpW));
trgP <- which(trgW == max(trgW));
shp <- length(1:trgP)+length(inpP:length(inpW))-ceiling(2*inhw)+inpP-trgP; # target shift
fleng <- inpP+shp+length(trgP:length(trgW));
inpW <- c(inpW,rep(0,fleng-length(inpW)));
trgW <- c(trgW,rep(0,fleng-length(trgW)));

wfc <- trgFilt(inpW,nu,circShift(trgW, shp));

maxW <- which(wfc==max(wfc)); 
hw <- floor(length(wfc)/2);
wfc <- circShift(wfc, (hw-maxW)); ## make a filter symmetric about maximum
## no need to truncate if no secondary maximum is present
#wfc <- coefTrunc(wfc);

return(as.vector(wfc));
} ###### END of "coefOptFiltLL" #####

