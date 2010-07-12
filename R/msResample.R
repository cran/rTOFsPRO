"msResample" <-
function(sigt, t0, masst, sigma0, asym, pfit, mHW, Gnoise){
		
###--------------------------------------------------------------  
### Function resamples mass axis and integrates signal with the step
### determined by the ratio of predicted STD for the peak outside
### mass focusing range (given by quadratic polynomial fit) to the
### intitial STD (left half) of the peak in the mass focusing range
### as measured by gpdf_hw
###
### INPUT:
### 
### sigt - raw equisampled TOF signal;
### t0 - input time ( may start above 1)
### masst - mass axis, corresponding to sigt; ## hard coded
### sigma0 - value of STD for the instrumental function peak in the
### mass focusing range (in time points, can be fractional);
### asym - right peak asymmetry in mass focusing range (can be negative);
### pfit - is the array of polyfit coefficients [p1,p2,p3] for the
### expected peak STD  outside mass focusing range; ## hard-coded
### mHW - sufficient constant point density (peak HW)
###
### Output:
###
### sigr - resampled, integrated TOF signal;
### massr - resampled mass axis, corresponding to sigr;
### rrr - integer resampling rate
### ii - resampled time
### npikw - resampled peak width
### nasym - resampled asymmetry
###-------------------------------------------------------------------------


#t1<-1:length(sigt);  -- data Start?? Read from input

rr<-(pfit[1]*t0*t0+pfit[2]*t0+pfit[3])/sigma0;
rr0<-rr; rrl0 <- length(rr);
tt<-1:length(rr);
stp<-round(sigma0/mHW);
npikw <- sigma0; nasym <- asym;
#print("const_rstep")
#print(stp)
if (stp >= 3){
  if ((stp/2.0-trunc(stp/2.0))==0) { stp <- stp-1; }
  trs <- seq (from=1, to=length(tt), by=stp);
  k1 <- 1; sigrStp <- rep(0, length(trs));
  ks <- (stp-1)/2;
  sigrStp[1] <- sum(sigt[1:ks]);
  ltrs <- length(trs);
  sigrStp[ltrs] <- sum(sigt[(trs[ltrs]-ks):trs[ltrs]]);
  for (jk in seq(from=trs[2], to=trs[ltrs-1], by=stp)) {
     k1 <- k1+1; sigrStp[k1] <- sum(sigt[(jk-ks):(jk+ks)]);
  }
npikw <- npikw/stp; nasym <- nasym/stp;
masst <- masst[trs];
sigt <- sigrStp;
rr0 <- rr0[trs];
rr <- rr[trs];
tt <- tt[trs];  ## tt or t0 ???
}

t1 <- 1:length(tt);

mintt<-eval(which(rr0==min(rr0)))

rr[(rr-floor(rr))<0.5]<-floor(rr[(rr-floor(rr))<0.5]);

rr[-((rr-floor(rr))<0.5)] <- ceiling(rr[-((rr-floor(rr))<0.5)]);

#print("minimum_RSrate")
#print(min(rr))

if (min(rr)>1) { ### YOU MAY comment the warning off -- it will run, but unreasonable PARS
#warning("Check that your constant HW parameter is NOT less than minimum HW of the quadratic fit")
}
ir <- which(rr>=2); ir <- ir[ir>mintt];

#print("how_many_irS")
#print(length(ir))
#print(ir[1])

if (length(ir)<1){
  sigr <- sigt;
  massr <- masst;
  rrr <- rep(1, length(tt));
  ii <- tt;
} else {
  #print("onset_of_sampling")
  #print(ir[1])  ## onset of sampling
  mr <- max(rr[(ir[1]):length(rr)])  ## maximum sampling rate for this record
  #print("max_rr")
  #print(mr)
  sigr <- sigt[1:(ir[1]-2)]; 
  #print(length(sigr))
  sigr[length(sigr)] <- sigr[length(sigr)]+0.5*sigt[ir[1]-1];
  massr <- masst[1:(ir[1]-2)];
  rrr <- rr[1:(ir[1]-2)];
  ii <- 1:(ir[1]-2);
  rsp <- ir[1]; rst <- ir[1];
  if (mr<3){  # case of mr=2
     rr[(ii[length(ii)]+1):length(rr)]<-2;
     #rst <- c(rst, ii[length(ii)]+1);
     #print("length_t1")
     #print(length(t1))
     if (t1[length(t1)]>(ir[1])){ #+1?
     ii<-c(ii, t1[seq(from=ir[1], to=length(t1), by=2)]);
     } else {
     ii <- c(ii, ir[1])
     #print(ii[length(ii)])
     }
     
  }else{ ## mr > 2
  for (jj in (rr[ir[1]]+1):mr) { ## start index resampling for mr>2
     ich <- which(rr>=jj);
     ich <- ich[ich>ir[1]];
     rsp <- c(rsp,ich[1]);
     rr[(ii[length(ii)]+jj-1):rsp[length(rsp)-1]]<-rr[rsp[length(rsp)-1]];
     rst <- c(rst, ii[length(ii)]+jj-1);
     ii<-c(ii, t1[seq(from=(ii[length(ii)]+jj-1), to=(rsp[length(rsp)]-jj), by=(jj-1))]);
  }  ## finished index resampling
##### HERE: treat tails for mr>2
  rr[(ii[length(ii)]+mr):rsp[length(rsp)]] <- rr[rsp[length(rsp)]];
  rst <- c(rst, (ii[length(ii)]+mr)); rst <- rst[2:length(rst)];
  
  strt <- ii[length(ii)]+mr;
  endd <- t1[length(t1)]-ceiling(mr/2);
  if (endd > (strt+mr)) {
  ii <- c(ii, t1[seq(from=strt, to=endd, by=mr)]);
  } else { ii <- c(ii, t1[strt])}
  
  }
  
  if (rst[1]<length(ii)) {
  irs <- ii[rst[1]:length(ii)];
  } else {
  irs <- rst[1]
  }
  
  massr <- c(massr, masst[irs]);
  rrr <- c(rrr, rr[irs],rr[irs[length(irs)]]);
  sr <- rep(0,length(irs)-1);
  
  if(length(irs)>2){
  for (ri in 1:(length(irs)-1)) {  ## do integrative resampling of signal intensity
     #print("inside_int_resampl")
     #print(ri)
     rri <- rr[irs[ri]];
     rri1 <- rr[irs[ri+1]];
     if ((rri/2.0-trunc(rri/2.0)) == 0) { ## even resampling rate
       k <- rri/2;
       sri <- sum(sigt[(irs[ri]-k+1):(irs[ri]+k-1)])+0.5*sigt[irs[ri]-k];
       if ((rri1/2.0-trunc(rri1/2.0))!=0) { ## next resampling rate is odd 
         sri <- sri + sigt[irs[ri]+k];
       } else { sri <- sri + 0.5*sigt[irs[ri]+k];}

     } else { ## odd resampling rate
       k <- (rri-1)/2;
       sri <- sum(sigt[(irs[ri]-k):(irs[ri]+k)]);
       if ((rri1/2.0-trunc(rri1/2.0))==0) { ## next resampling rate is even
         sri <- sri + 0.5*sigt[irs[ri]+k+1];
       }

     }
     
     sr[ri] <- sri;

  }  ## end loop for integrative resampling
  
  sigr <- c(sigr, sr);
  }  ## end-if length(irs)>2
  
  #print("after_resample_loop")
  #print(length(sigr))
  #print(length(ii))
  
  ## resample last point assuming that next rate is the same
  ri <- length(irs); rri <- rr[irs[ri]];
  if (length(irs)>1){
  if ((rri/2.0-trunc(rri/2.0)) == 0) { ## even resampling rate
    k <- rri/2;
    sri <- sum(sigt[(irs[ri]-k+1):(irs[ri]+k-1)]);
    sri <- sri +0.5*sigt[irs[ri]-k] + 0.5*sigt[irs[ri]+k];
  } else { ## odd resampling rate
    k <- (rri-1)/2;
    sri <- sum(sigt[(irs[ri]-k):(irs[ri]+k)]);
  }
  ## add last resampled points
  sigr <- c(sigr, sri, mean(sigr[(length(sigr)-k):length(sigr)]));
  massr <- c(massr, masst[ii[length(ii)-1]]);
  }else {
  sigr <- c(sigr, sigt[irs])
  massr <- c(massr, massr[irs])
  }
  
  if (ii[length(ii)]< t1[length(t1)]){
  rrr <- c(rrr, rrr[length(rrr)]);
  ii <- c(ii, t1[length(t1)]);
  massr <- c(massr, masst[t1[length(t1)]]);
  sigr <- c(sigr, mean(sigt[ii[length(ii)-1]:ii[length(ii)]])); 
  } 
  
#print("before_ramp")
#print(length(sigr))

sigr[is.na(sigr)]<-sigt[t1[ii[is.na(sigr)]]] ## correct for missing edge values ??

## scale by sqrt(rate) for Gaussian noise
  if (Gnoise==1){    
    sigr[(ir[1]-2):length(sigr)]=sigr[(ir[1]-2):length(sigr)]/sqrt(rrr[(ir[1]-2):length(rrr)]);
#print("IDS scaled for Gnoise: check to back the peak- amps after detection")
  }  

  ## find ramps
  rmpi <- which(diff(rrr)==1); rmpi <- rmpi+1; ## ramp positions
  #print("ramp_positions")
  #print(rmpi)
  rmpa <- sigr[rmpi]-sigr[rmpi-2]; ## ramp amplitudes
  ## estimate noise amplitude (2*std) near signal minimum
  sgnl <- sigr[rmpi[1]:(length(sigr)-2*mr-1)];
  #print("length_signl")
  #print(length(sgnl))
  pst <- 1; ## intialize
  mi <- which(sgnl==min(sgnl)); pst <- rmpi[1]+mi;
  if (length(pst)>1) {pst <- pst[1]} ## correct for mulitple minima ??
  
  if ((pst > (rmpi[1]+3*mr)) & (pst < (length(sigr)-3*rrr[length(sigr)]))){
  nsa <- 2*sd(sigr[(pst-2*rrr[pst]):(pst+2*rrr[pst])]);
  }else{  ## pst is too close to 1st ramp or to signal-end: estimate "pre-ramp" noise
  sgnl <- sigr[(rmpi[1]-7*npikw-2*mr):(rmpi[1]-2*mr-1)];
  mi <- which(sgnl==min(sgnl)); pst <- rmpi[1]-7*npikw-2*mr+mi;
  nsa <- 2*sd(sigr[(pst-2*mr):(pst+2*mr)]);
  }
  #print("noise_amplitude"); print(nsa);
  i1 <- which((rmpa>nsa) & ((rmpa+nsa)>(sigr[rmpi+1]-sigr[rmpi-1])));
  #ii[rmpi[i1]]; ## ramps to correct
  #print(t0[tt[ii[rmpi]]])  ## all ramps in original time 

  if (length(i1)>=1){  ## at least one ramp was detected
  for (cr in 1:length(i1)) { ## correct for ramps
     si <- rmpi[i1[cr]];
     trmp <- ii[rmpi[i1[cr]]];
     trmp1 <- ii[rmpi[i1[cr]]-2];
     rrmp <- rrr[rmpi[i1[cr]]]; ## resampling rate at the ramp
     
     lcr <- trunc(length(which(rr0[trmp:length(rr0)]<rrmp))/rrmp);
     wcr <- 1/lcr;
     # length of points to correct to the right of ramp 
     lcr <- abs(floor(lcr - 2*nsa/(wcr*rmpa[i1[cr]])));
     # weights for points to correct to the right of ramp
     wcr <- 1/lcr;
     
     if(rrmp > 1) {  # safety check for rate = 1 at the ramp added 03/20/10
     	lcl <- trunc(length(which(rr0[mintt:trmp1]>(rrmp-1)))/(rrmp-1));
     } else { lcl <- 1 } 
     wcl <- 1/lcl;
     # length of points to correct to the left of ramp
     lcl <- abs(floor(lcl-2*nsa/(wcl*rmpa[i1[cr]])));
     # weights for points to correct to the right of ramp
     wcl <- 1/lcl; 
     rrl <- length(which(rrr==rrmp)); 
     rll <- length(which(rrr==(rrmp-1)));
 
  
     	if (lcr>0.5*rrl) {lcr <- floor(0.5*rrl); wcr <- 1/lcr}
     	if (lcr >1){
       		dcr <- 0.5*rmpa[i1[cr]]*(1-wcr*(1:lcr));
       		sts <- si:(si+lcr-1);
       		sigr[sts] <- sigr[sts]-dcr;
     	}

     	if (lcl>0.5*rll) {lcl <- floor(0.5*rll); wcl <- 1/lcl}
    	 if (lcl > 1) {
       		dcr <- 0.5*rmpa[i1[cr]]*(1-wcl*(1:lcl));
       		sts <- seq(from=(si-2), to=(si-2-lcl+1), by=-1);
       		sigr[sts] <- sigr[sts]+dcr;
     	}

  }  ## finish loop for ramp correction
 } ## safe check that at least one ramp exists

#print("after_ramp")
#print(c(length(sigr), length(ii)))
#print(ir[1])
 ii <- tt[ii];
#print(ii[ir[1]])
ii[ir[1]:(length(ii)-1)]<-ii[(ir[1]+1):length(ii)] ## correct 1-point shift in resampling onset 11/19/07
#print(ii[ir[1]])
ii[length(ii)]<-rrl0; ## correct copy of end value 11/19/07
#print(ii[(length(ii)-2):length(ii)])  

} ### safety check: if resampling is at all needed


if (length(sigr)<length(ii)) {ii <- ii[1:length(sigr)]} ## check for correct length to interpolate

rst0<-t0[ii]; 
if (!exists("rst")) {rst <- ii} # if no down-sampling was performed
#print('INFO: *resSignal*, *resMass*, *resRate* and *resTime* fields of the list are returned')
resampleOut <- list(resSignal=sigr, resMass=massr, resRate=rrr, resTime=rst0, resPikw=npikw, resAsym=nasym)
rm(ii, rst, sigr, rrr)
return(as.list(resampleOut))
}  ###### END of "msResample" #####

