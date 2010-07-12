"noiseEstmtr2" <-
function(spect)
{

#----------------------------------------------------
# Estimate noise level (STD) for TOF record (assuming constant peak width)
#
# INPUT:
# spect - spectrum signal (resampled or range of constant peak width);
#
# OUTPUT:
# ns - noise estimate assuming Gaussian (random) noise;
#
# USAGE:
# ns <- noiseEstmtr2 (spect);
# Dependency: none
#

 
## Part 1  Noise Level Estimation

len <-length(spect); 

#A <- spect[1:(len-1)]; B <- spect[2:len]; d1 <- A - B; 
d1 <- -diff(spect,1,1); 

spect <- sort(d1); 
#print(length(spect)) # debug
tstart <- min(which(spect>0))-round(0.05*len);
tend <- tstart+round(0.1*len);
#print(c(tstart,tend)) #debug
spect <- spect[tstart:tend];
tt <- 1:length(spect); #tt <- t(tt);
linef <- lsfit(tt, spect, intercept  = TRUE); # line(tt,spect);  
fcoefs <- coef(linef)
ns <- len*fcoefs[2]/3.55; 

return(ns);
}

