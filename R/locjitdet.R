"locjitdet" <-
function(fs, ps, sw, startS, endS, spectraNames){

###-----------------------------------------------------------
### Compute jitter  and triggering error in time points for the set of signals
### at the pivit peak positions.
### Jitter is measured by the shifts of the maxima for the local
### cross-correlation function in the vicinity of a peak.
### Jitter is detected in resspect to the FIRST (reference) signal in the set.
###
### INPUT:
### fs - matrix of sampled discrete signals (columns - intensities)
### ps - peak set for fs precomputed eg. with "pivpiks"
### sw - maximum lag for cross-correlation (shift window ~ of twice
### peak sigma, estimated eg. with "peakwidth_pf");
###
### OUTPUT:
### out$pjit = fsj (i,j) - jitter in point numbers for the corresponding peak
### "i" in the signal "j" as measured by the local cross-correlation to the reference signal;
### out$trerr = avshift(j) - average triggering error for the j-th record
###-----------------------------------------------------------

### "sw" is a window for local correlation function (maximum lag)
### sw=10 is recommended for SELDI data (focused from 2 to 15 kDa)

	cw<-3*sw; # is a half-window in a spectrum for the local X-corr calculation
	ps<-as.vector(ps)
	
	l0<-length(ps)

	### remove peaks for which correlation cannot be calculated with the required window
	maxpoint<-eval(nrow(fs))
	ps[ps<(1+cw)]<-0               ### remove peaks closer than correlation window to the record's zero
	ps[ps>(maxpoint-eval(cw)-1)]<-0 ### remove peaks closer than correlation window to the record's end
	t<-ps[ps!=0]
	ps<-t
	
	Nr<-length(ps);
	if ( Nr < 3 )
	{
		userError( "Could not find at least 3 peaks in the reference spectrum." );
	}

	Nc<-ncol(fs); ### if we have multiple spectra

	if (Nr<l0) { print('WARNING: less peaks used -- some peaks outside allowed correlation window')}

	fsj<-matrix(rep(0,Nr*Nc),nrow=Nr,ncol=Nc);  ### initialize jitter table for the peaks
	msj<-fsj; ### initialize maximum correlation coefficient table for the peaks

	enhstart<-eval(cw-floor(sw/2)) ### peak enhancement range for the reference signal
	enhend<-eval(cw+floor(sw/2))
	refsig <- fs[,1]/max(fs[startS:endS,1]);
	lmax <- trunc(length(startS:endS)/2);
	#lmax <- length(startS:endS)-1;
	if (Nr>0)
	{
		for (j in 2:Nc)
		{   
			### calculate jitter for all peaks
			nsigJ <- fs[,j]/max(fs[startS:endS,j]);
			cfg <- ccf(refsig[startS:endS], nsigJ[startS:endS], lag.max=lmax, type="correlation", plot=FALSE);
			cw <- 3*abs(cfg$lag[cfg$acf==max(cfg$acf)]) + 3*sw;
			
			cwStart <- cw + 1;
			cwEnd <- length( refsig ) - cw;
			
			peaksWindowStart <- ps[ ps < cwStart ];
			peaksWindowEnd <- ps[ ps > cwEnd ];
			
			windowStartError <- length( peaksWindowStart ) > 0;
			windowEndError <- length( peaksWindowEnd ) > 0;
			
			if ( windowStartError || windowEndError )
			{
				msg <- paste( "\nWhile aligning spectrum '", spectraNames[j], "' to reference spectrum '", spectraNames[1], "':\n\n", sep = "" );
				checkRefSpecMsg <- "\nPlease check the reference spectrum.";
				if ( cwStart >= cwEnd )
				{
					userError( paste( msg, "For all pivot peaks, the correlation window goes beyond the range of the spectrum data.", checkRefSpecMsg, sep = "" ) );
				}
				
				msg <- paste( msg, "For some pivot peaks, the correlation window goes beyond the range of the spectrum data.", sep = "" );
				msg <- paste( msg, checkRefSpecMsg, sep = "" );
				msg <- paste( msg, "\n\nThe peaks that caused this problem are ", sep = "" );
				
				correctionStart <- "\nIf the reference spectrum is good, then adjusting ";
				correctionEnd <- " will exclude these peaks and should resolve this error.";
				
				if ( windowStartError && windowEndError )
				{
					msg <- paste( msg, "<= ", peaksWindowStart[ length( peaksWindowStart ) ], " and >= ", peaksWindowEnd[ 1 ], ".", sep = "" );
					msg <- paste( msg, correctionStart, "'Start Time' to ", cwStart, " and 'End Time' to ", cwEnd, correctionEnd, sep = "" );
				} else if ( windowStartError )
				{
					msg <- paste( msg, "<= ", peaksWindowStart[ length( peaksWindowStart ) ], ".", sep = "" );
					msg <- paste( msg, correctionStart, "'Start Time' to ", cwEnd, correctionEnd, sep = "" );
				} else
				{
					msg <- paste( msg, ">= ", peaksWindowEnd[ 1 ], ".", sep = "" );
					msg <- paste( msg, correctionStart, "'End Time' to ", cwEnd, correctionEnd, sep = "" );
				}
				
				userError( msg );
			}
			
			for (i in 1:Nr)
			{
				startp<-eval(ps[i]-cw)
				endp<-eval(ps[i]+cw)
				
				refc<-refsig[startp:endp];

				fsc<-nsigJ[startp:endp];

				cf<-ccf(refc,fsc,lag.max=trunc(cw/2), type="correlation", plot=FALSE);
				mm<-max(cf$acf);
			    fsj[i,j]<-cf$lag[cf$acf==max(cf$acf)];

			    #if ((fsj[i,j]==trunc(cw/2)) | (fsj[i,j]==-trunc(cw/2))) {fsj[i,j]=0;} 
		    	msj[i,j]<-mm;
			}
 		}
	} else {print('WARNING: EMPTY JITTER MATRIX -- all peaks outside allowed correlation window')}

	avshift<-as.vector(rep(0,Nc-1))

	for (jj in 1:(Nc-1))
	{
		signjj<-sign(fsj[,jj+1])
		lneg<-force(length(signjj[signjj<0]))
		lpos<-force(length(signjj[signjj>0]))
		lzer<-force(length(signjj[signjj==0]))
		votesign<-max(c(lneg, lpos, lzer))

		avneg<-sum(fsj[signjj<0,jj+1])
		if (lneg>0){avneg<-avneg/lneg}

		avpos<-sum(fsj[signjj>0,jj+1])
		if (lpos>0){avpos<-avpos/lpos}

		avzer=sum(fsj[signjj==0,jj+1])
		if (lzer>0){avzer<-avzer/lzer}

		avshift[jj]<-avneg*(lneg==votesign)+avpos*(lpos==votesign)+avzer*(lzer==votesign)
	}

	jitout <- list(pjit=fsj[,], trerr=avshift)

	return(as.list(jitout))
	#return(as.vector(fsj));
}

