"coefTrunc" <-
function(fcoef){

maxF <- which(fcoef==max(fcoef));
hw <- floor(length(fcoef)/2);
ftrn <- circShift(fcoef, (hw-maxF)); ## symmetrize the coefficients around maximum
flen <- floor(length(fcoef)/3)+20;
ftrn <- ftrn[(hw-floor(flen/2)):(hw+floor(flen/2))];
return(as.vector(ftrn));
} ###### END of "coefTrunc" #####

