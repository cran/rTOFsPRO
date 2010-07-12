"coefNLFiltGL" <-
function(inhw, las, trct) {

nu1 <- 0.01; nu2 <- 0.001; nu3 <- 0.0001;
inpW<-c(gausHalf(inhw, trct), lorHalf((inhw+las), trct));
trgW1<-c(gausHalf(0.2*inhw, trct*5), lorHalf(0.2*(inhw+las), trct*5));
trgW2<-c(gausHalf(0.5*inhw, trct*2), lorHalf(0.5*(inhw+las), trct*2));
trgW1 <- trgW1*sum(inpW)/sum(trgW1); # normalize target wavelet ??
trgW2 <- trgW2*sum(inpW)/sum(trgW2); # normalize target wavelet ??

inpP <- which(inpW == max(inpW));
trgP <- which(trgW2 == max(trgW2));
shp <- length(1:trgP)+length(inpP:length(inpW))-ceiling(2*inhw)+inpP-trgP; # target shift
fleng <- inpP+shp+length(trgP:length(trgW2));
inpW <- c(inpW,rep(0,fleng-length(inpW)));
trgW1 <- c(trgW1,rep(0,fleng-length(trgW1)));
trgW2 <- c(trgW2,rep(0,fleng-length(trgW2)));

wfc1 <- trgFilt(inpW,nu1,circShift(trgW1, shp));
wfc2 <- trgFilt(inpW,nu2,circShift(trgW1, shp));
wfc3 <- trgFilt(inpW,nu3,circShift(trgW2, shp));
wfc1 <- coefTrunc(wfc1);
wfc2 <- coefTrunc(wfc2);
wfc3 <- coefTrunc(wfc3);

nlfCoef3 <- list(fc1=wfc1, fc2=wfc2, fc3=wfc3)

return(as.list(nlfCoef3))
} ###### END of "coefNLFiltGL" #####

