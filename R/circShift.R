"circShift" <-
function(dataV, shiftN) {

if (shiftN>0) {

shdat <- c(dataV[sort(c((length(dataV)-shiftN+1):length(dataV)),decreasing=TRUE)],dataV[1:(length(dataV)-shiftN)]);

} else if (shiftN < 0) {

shdat <- c(dataV[(-shiftN+1):length(dataV)],dataV[1:(-shiftN)]);

}  else {shdat <- dataV}

return(as.vector(shdat));
} ###### END of "circShift" #####

