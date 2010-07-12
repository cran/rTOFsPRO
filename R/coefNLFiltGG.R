"coefNLFiltGG" <-
function(inhw, las, trct) {

nu1 <- 0.01; nu2 <- 0.001; nu3 <- 0.0001;
inpW<-c(gausHalf(inhw, trct), gausRHalf((inhw+las), trct));
trgW1<-c(gausHalf(0.2*inhw, trct*5), gausRHalf(0.2*(inhw+las), trct*5));
trgW2<-c(gausHalf(0.5*inhw, trct*2), gausRHalf(0.5*(inhw+las), trct*2));
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
maxW <- which(wfc1==max(wfc1)); 
hw <- floor(length(wfc1)/2);
wfc1 <- circShift(wfc1, (hw-maxW)); ## make a filter symmetric about maximum

wfc2 <- trgFilt(inpW,nu2,circShift(trgW1, shp));
maxW <- which(wfc2==max(wfc2)); 
hw <- floor(length(wfc2)/2);
wfc2 <- circShift(wfc2, (hw-maxW)); ## make a filter symmetric about maximum

wfc3 <- trgFilt(inpW,nu3,circShift(trgW2, shp));
maxW <- which(wfc3==max(wfc3)); 
hw <- floor(length(wfc3)/2);
wfc3 <- circShift(wfc3, (hw-maxW)); ## make a filter symmetric about maximum

## no need to truncate if no secondary maximum is present
#wfc1 <- coefTrunc(wfc1);
#wfc2 <- coefTrunc(wfc2);
#wfc3 <- coefTrunc(wfc3);

nlfCoef3 <- list(fc1=wfc1, fc2=wfc2, fc3=wfc3)

return(as.list(nlfCoef3))
}  ###### END of "coefNLFiltGG" #####

