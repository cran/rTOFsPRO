"coefOptFiltLL" <-
function(inhw, las, trct) {

nu <- 0.001;
inpW <- c(lorLHalf(inhw, trct), lorHalf((inhw+las), trct));
trgW <- c(lorLHalf(0.8*inhw, trct*1.4), lorHalf(0.8*(inhw+las), trct*1.4));
trgW <- trgW*sum(inpW)/sum(trgW); # normalize target wavelet ??
inpP <- which(inpW == max(inpW));
trgP <- which(trgW == max(trgW));
shp <- length(1:trgP)+length(inpP:length(inpW))-ceiling(2*inhw)+inpP-trgP; # target shift
fleng <- inpP+shp+length(trgP:length(trgW));
inpW <- c(inpW,rep(0,fleng-length(inpW)));
trgW <- c(trgW,rep(0,fleng-length(trgW)));

wfc <- trgFilt(inpW,nu,circShift(trgW, shp));
wfc <- coefTrunc(wfc);

return(as.vector(wfc));
} ###### END of "coefOptFiltLL" #####

