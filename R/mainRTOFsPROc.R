"mainRTOFsPROc" <-
function(OptionsAndParametersFileName) {
##======================================================================##
##======================================================================##
##									##
##		MAIN RTOFsPROc CODE  021210				##
##									##
## 	This script links and executes the functions of rTOFsPRO	## 
##	package for signal processing of the linear TOF spectra 	##
##	to produce aligned peak intensity matrix (and peak positions)   ##
##	that can be used for further statistical analysis		##
##									##
##   		>>> best used in PRODUCTION mode <<<			## 
##  	 >>> for larger data sets (as memory permits) <<<		##
##   		    >>> no "debug" option <<<				##
##									##
## OUTPUT:								##
##									##
## alignedPeakList and alignedPeakListMetaData structures containing	## 
## on peak intensities and positions and their uncertainties, as well	##
## as clinical and instrumental meta-data.				##
## "..._Test.RDat" files save auxhiliary info and data structures,	## 
## if desired (e.g. tofResampled and peakList) -- see documentation	##
##									##
## INPUT:								##
## tofList and tofListMetaData structures (can be generated from raw	##
## e.g., using WMBrukerParser R-package -- see documentation for	## 
## detailed description)						##
## Signal processing parameters are read from "OptionsAndParameters"	## 
## file, and can be edited there. 					##
##									##
##	PROCESSING STEPS INCLUDE (in a loop over spectra)		##
##									##
##	Set up:								##
##		Load "OptionsAndParameters.txt" file			##
##		Load processing libraries				##
##									##
##	Load Data:							##
##		Load Data and MetaDataSet parameters			##
##		"multiple files" Loop begins (See note 1)		##
##									##
##	Background Correction. Options include:				##
##		RC (charge accumulation model)				##
##		Analytical (Gaussian, Exponential, linear) model	##
##		None (default)						##
##									##
##	Preliminary alignment. Options include:				##
##		Global alignment (using correlation between pivot peaks)##
##		TOF offset adjustment (from metaData)			##
##		None							##
##									##
##	Downsampling and Filtering. Options include:			##
##		Integrative downsampling (IDS) only			##
##		IDS and Optimal Linear Filter (default)			##
##		IDS and Nonlinear Filter				##
##		IDS and Matched Filter					##
##		IDS and MAV Filter					##
##									##
##	Pedastal removal. Options include:				##
##		MAV of local minima 					##
##		None (default)						##
##									##
##	Peak Picking. Options include:					##
##		Trivial	(1st difference with local minim SNR, Note 2)	##
##		W & M (Maximum Likelihood Method) (needs extra lib)	##
##		(call commented off: lib available from wecook@wm.edu)	##
##									##
##	Alignment. Options include:					##
##		Align to average spectra				##
##		W & M global+binning (call commented off: 		##
##		seprate lib available from wecook@wm.edu)		##	
##									##
#========================================================================#
# Note 1. Due to the large size of data files, in "compression" mode     #
# raw tofList data can be loaded and partially processed by IDS from     #
# multiple files. The "multiple TOFs loop" includes blocks from "Load    # 
# data" through "Downsampling and Filtering".At the end of the loop, the #
# downsampled spectra and their meta-data are concatenated and the	 # 
# original tofList is removed.						 #
#========================================================================#
# Note 2. Basic algorithm adopted from Matlab-based Cromwell toolbox.	 #	
# (Copyright 2005.) Originally developed by MD Anderson and distributed  #
# under "BSD" Mathworks license: http://bioinformatics.mdanderson.org/   #
# cromwell.html. Local minima SNR threshold and baseline-dependent	 #
# noise model for each preliminary filtering option was added here.	 #
##======================================================================##

############################  Set Up Block  ##############################

### Set memory limit and clean up space ###

	memory.limit(1800)
	rm(list=ls())
	gc()

### Load Options and Parameters File ###

	# source(paste("rTOFsPRO/OptionsAndParametersPCAqc33.txt",sep="/")) 
	source(OptionsAndParametersFileName)
#TOFfileDIR = system.file("Examples", package ="rTOFsPRO")
	setwd(TOFfileDIR) # working directory for IO

	if (monitorMemory){
		startMemoryUsed<-memory.size(max=FALSE)
		print("Memory Used Before loading data")
		print(startMemoryUsed)
	}
	if (averageOnly) {
		print("Will produce average spectrum *only*: *no processing* requested")
	}

####### Load TOF spectra and metaData from multiple files (LOOP START) #######

#MBT#	setwd(TOFfileDIR)  ## set working directory to load data files from, 
			   ## and save the results (if desired, otherwise default is used)

Nfiles <-length(rTofFileNames); 

if (Nfiles != length(rTofMetaFileNames)) { 
	userError("Spectra files must have corresponding meta-data files: check <OptionsAndParameters>")
}


for (rdl in 1:Nfiles){  ## For multiple TOF data files (e.g., multiple biprocessor RobotRuns)
			## from the same experiment (same instrumental settings)
        rTofFileName <- rTofFileNames[rdl]  # reads from "OptionsAndParameters" file
	rTofMetaFileName <-rTofMetaFileNames[rdl]
	load(rTofFileName)
	load(rTofMetaFileName)
print("---------------------------"); 
print("Processing input data file:"); print(rTofFileName); 
print("---------------------------");
 
### Initialization ###

	dataStart <- 1 + as.integer(tofListMetaData[1, "timeOfFlightData.start"]);# same for all spectra
	sampleRate <- as.numeric(tofListMetaData[1, "instrumentSpecificSettings.timeDelta" ]);
	Tav1 <-  min(as.integer(tofListMetaData[,"timeOfFlightData.end"]))
	if (rdl > 1) { Tav <- min(c(Tav, Tav1))} else { Tav <- Tav1} # length of average spectrum
	if (Tav < osS[length(osS)]) {
		print("Baseline offset calculation range is outside of a spectrum length:"); 
		userError("Choose osS < min(timeOfFlightData.end) in <OptionsAndParameters>")
	}

### Get information from tofList ###

	spectraName <- names(tofList)
	spectraCount <- length(spectraName)
	Ns <- length(tofList[[1]])
	t0 <- dataStart:(Ns + dataStart - 1) 
 
### perform only for raw tof data!

########################### End Set Up Block  ############################
	if (monitorMemory){
		memoryUsed<-memory.size(max=FALSE)
		print("Memory Used After loading data")
		print(memoryUsed)
	}

########################## End Load Data Block  ##########################


#################### Preliminary Alignment Block  #########################

# Call appropriate preliminary alignment procedure

	if (identical(paMethod, "globalAlign")) {

### Apply global alignment ###
		
		print("Preliminary Alignment Method: Global Align")
		if (detect) {
			tof <- tofList; # initialize list of bs-spectra
			##need to subtract "rough" baseline and smooth *beforehand* to detect shifts, 
			##but apply to "unsubtracted" spectra
			if (rdl > 1 & exists("savl")) { ## in case of mutiple input data files
				tof <- c(savl, tof);  ## use previous average spectrum as ref.
				tofList <- c(savl, tofList); 
				tofListMetaData <- rbind(tofListMetaData[1,], tofListMetaData);
				rownames(tofListMetaData)[1] <- names(savl)
				
			} 
			tof <- ExpBkgr_call(tof, bsS,osS,mCH, dataStart) # subtract rough baseline
			for (ia in 1:spectraCount){
				insp <- tof[[ia]]
				smsp <- mavSmoothing(insp,2*pikw0) # smooth with MAV
				tof[[ia]] <- as.vector(smsp)
			} 
		}

	aList <- globalAlign(tof, tofList, tofListMetaData, pikw0,noiseThres,startS,endS, appli, detect)
        tofList <- aList$tList; tofListMetaData <- aList$tMeta;
	rm(aList); gc();  ## clean up memory  
	if (rdl > 1 & detect & exists("savl")) { ## in case of mutiple data files
		lentl <- length(tofList); lenmd <- dim(tofListMetaData)[1]
		tofList <- tofList[2:lentl];  ## remove previous average spectrum as ref.
		tof <- tof[2:lentl];
		tofListMetaData <- tofListMetaData[2:lenmd,] 
	} 
	rm(tof); gc(); ## clean memory

####### End global alignment #######

		} else if (identical(paMethod, "delay")) {

####### Begin Offset Adjust #######

		print("Preliminary Alignment Method: Time Delay Adjustment")
	## assumes the same sampling rate for all records, but different time-zero for some

	if(!exists("zero0") ) {
		v1 <- as.numeric(tofListMetaData[,"instrumentSpecificSettings.timeZero" ])
		zero0<- max(v1)
	} 
	#print(v1)
	#print(sampleRate)
        adsp <- 0;
	spectraName <- names(tofList)
	spectraCount <- length(tofList)

	for (i in 1:spectraCount){
		zero<- as.numeric(tofListMetaData[i,"instrumentSpecificSettings.timeZero" ])
		if(zero!=zero0){
			delay<-round((zero0-zero)/sampleRate);
			adsp<-c(adsp,i);
	   		spectrum <- tofList[[spectraName[i]]]
			shsp <- circShift(spectrum, delay)
			tofList[[spectraName[i]]] <- as.vector(shsp) 
		}
	} 

	if (length(adsp)>1){
	adsp <- adsp[2:length(adsp)]; # spectra with applied offset adjustment 
	print("####   Pre-adjusted spectra:")
	print(adsp)
	} else {
	print("####   None of the spectra needed pre-adjustment    ####")
	}
	inputSpectrum <- tofList[[1]] 
	Ns <- length(inputSpectrum)

### End Offset Adjust ###

		} else {
		print("Preliminary Alignment Method: none")
	} # end preliminary alignment options


################## End Preliminary Alignment Block  #######################

###################  Construct Average Spectrum  ########################

	if (rdl == 1){
		sav <- tofList[[1]][1:Tav];
		Navs <- 0;
	} else {
		sav <- sav[1:Tav]*Navs+tofList[[1]][1:Tav];
	}
	if (spectraCount > 1){
		for (i in 2:spectraCount) {
    			sav <- sav+tofList[[i]][1:Tav]
		}
	} else { i <- 1 }
	Navs <- spectraCount + Navs;
	sav <- sav/Navs;
	cofsav <- mean(sav[osS]);  # constant offset for average-spectrum
	savl <- list();
	savl[[1]] <- sav; # make into tofList-like strcuture
	names(savl) <- "average"
	savbs <- savl

########################## End Construct Mean Spectrum Block  ###################

if (!averageOnly){
###################  Mean Spectrum Background Correction Block  ########################

# Call appropriate baseline subtraction procedure

	if (identical(bsMethod, "RC")) {
	  
		### Start RC-Background Correction ###
		print("Background Subtraction Method: RC")
		savbs <- chargeBkgr_call(savl,osS,amp,decay); # baseline-subtracted average spectrum

		### End RC-Background Correction ###

		} else if (identical(bsMethod, "Analytical")) {

		### Start Analytical Background Correction ###

		print("Background Subtraction Method: Analytical function fit");
		if (mCH == 1) { print("# Linear baseline model #")}
		else if (mCH == 2) {print("# Exponential baseline model #")}
		else {print("# Gaussian baseline model #")}
		savbs <- ExpBkgr_call(savl, bsS,osS,mCH, dataStart)

		### End Analytical-Background Correction ###

		} else if (identical(bsMethod, "ADCoffset")) {
		print("Background Subtraction Method: ADC offset");
		savbs <- offsetBkgr_call(savl, osS)
		} else { # no BS-subtraction

		print("Background Subtraction Method/Default: none")
		
	} # end baseline subtraction options

	bsav <- savl[[1]]-savbs[[1]]; #baseline of average spectrum

################## End Mean Spectrum Background Correction Block  ######################

#################### Mean Spectrum Downsampling/Filtering Block  #########################

	### MZ axis calculation ###
	if (useMeta) { # use instrumental meta-data for calibration
		tms <- (as.numeric(tofListMetaData[1,"param.T0"])+(t0-1)*as.numeric(tofListMetaData[1,"param.TDelta"]))/as.numeric(tofListMetaData[1,"param.U"])
		masst <- as.numeric(tofListMetaData[1,"param.c1"])*tms*tms+ as.numeric(tofListMetaData[1,"param.c0"])*tms+as.numeric(tofListMetaData[1,"param.c2"]);
	} else {
		tms <- (nDelay+(t0-1)*as.numeric(tofListMetaData[1,"param.TDelta"]))/as.numeric(tofListMetaData[1,"param.U"])
		masst <- A*tms*tms+ B*tms + C;
	}
	trnct <- 2.0^(instrPrc+1) # signal wavelet truncation number

	## pikw0,lasym0 - original input peak shape parameters from file

	pfit <- as.vector(c(p1,p2,p3))
	minHW <- mpikw; # for now, hard code minimum sufficient peakHW (point density),same scale as "pikw"

	## down-sample average spectrum
		msrec <- savbs[[1]]
		cOffset<-mean(msrec[(length(msrec)-2000):(length(msrec)-100)]);# PAR: Need to make parameters
		msrec<-msrec-cOffset;
		rsavout<-msResample(msrec,t0, masst, pikw0, lasym0, pfit, minHW, Gnoise);
		pikw<-rsavout$resPikw; fHW0 <- ceiling(pikw)+1; ## IDS/MAV signal HW
		lasym<-rsavout$resAsym; nsk <- ceiling(lasym); fHW <- fHW0; iir <- rsavout$resTime;
		rsav <- rsavout$resSignal; 
		nsav <- noiseEstmtr2(rsav)
		fwl <- 0; # initial filter length
		namavl <- savbs; namavl[[1]]<- (bsav[iir]-cofsav)*nsav;
		fsavl <- savbs; 

	### Calculate appropriate filter coefficients,length and filtered HW #####
	### filter average spectrum and record its filter noise amplitude "namavl"

	if ( identical(algo, "none") ){  ## added processing for "rsav"
		print("Filter: none")
		fsavl[[1]] <- rsav
		##noise amplitude threshold for down-sampled signal
    		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]] <- (sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*3; # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]]<-(rsav*0+nsav)*3;  # noise amplitude  
    		}
		fwl <- 0; # filter length
		fHW <- fHW;
	} else if (identical(algo, "optimalLinear")) {
		print("Filter: Optimal Linear")
		if (identical(wavelet, "LL")) {
			print("Wavelet Shape: Lorentzian")
			wof <- coefOptFiltLL(fHW, nsk, trnct);
			filtSig <- optFilt(rsav, wof);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra LL-filter phase shift
		} else if (identical(wavelet, "GG")) {
			print("Wavelet Shape: Gaussian")
			wof <- coefOptFiltGG(fHW, nsk, trnct);
			filtSig <- optFilt(rsav, wof);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		} else if (identical(wavelet, "GL")) {
			print("Wavelet Shape: lhalf-Gaussian-rhalf-Lorentzian")
			wof <- coefOptFiltGL(fHW, nsk, trnct);
			filtSig <- optFilt(rsav, wof);
			fsavl[[1]] <- circShift(filtSig, 1); # correction
		} else { print("No wavelet function was selected: using <GG> as default")
			wavelet <- "GG" # default
			wof <- coefOptFiltGG(fHW, nsk, trnct);
			filtSig <- optFilt(rsav, wof);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		}

    		fsavl[[1]][length(fsavl[[1]])] <- 0; # end-point correction 
    		## noise amplitude threshold for OLF
		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]]<-(sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*3*(max(fsavl[[1]])/max(rsav))/3.5; # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]]<-(rsav*0+nsav)*3*(max(fsavl[[1]])/max(rsav))/3.5;    # noise amplitude  
    		} 
		fwl <- length(wof); # filter length
		fHW <- 0.9*fHW; # updated filtered HW
	} else if (identical(algo, "matched")) {
		print("Filter: Matched")
		if (identical(wavelet, "LL")) {
			print("Wavelet Shape: Lorentzian")
			mfc <- matchFiltCoefLL(fHW, nsk, trnct)
			filtSig <- matchedFilter(rsav, mfc);
			fsavl[[1]] <- filtSig
		} else if (identical(wavelet, "GG")) {
			print("Wavelet Shape: Gaussian")
			mfc <- matchFiltCoefGG(fHW, nsk, trnct)
			filtSig <- matchedFilter(rsav, mfc);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		} else if (identical(wavelet, "GL")) {
			print("Wavelet Shape: lhalf-Gaussian-rhalf-Lorentzian")
			mfc <- matchFiltCoefGL(fHW, nsk, trnct)
			fsavl[[1]] <- matchedFilter(rsav, mfc);
		} else { print("No wavelet function was selected: using <GG> as default")
			wavelet <- "GG" # default
			mfc <- matchFiltCoefGG(fHW, nsk, trnct)
			filtSig <- matchedFilter(rsav, mfc);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		}
    		## noise amplitude threshold for matched
		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]]<-(sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*0.5*(max(fsavl[[1]])/max(rsav)); # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]]<-(rsav*0+nsav)*0.5*(max(fsavl[[1]])/max(rsav));    # noise amplitude  
    		} 
		fwl <- length(mfc); # filter length
		fHW <- 1.4*fHW; # updated filtered HW
	} else if (identical(algo, "nonLinear")) {
		print("Filter: Nonlinear")
		if (identical(wavelet, "LL")) {
			print("Wavelet Shape: Lorentzian")
			coef3list <- coefNLFiltLL(fHW, nsk, trnct);
			filtSig <- geomavFilt(rsav, trnct, coef3list);
			fsavl[[1]] <- filtSig
		} else if (identical(wavelet, "GG")) {
			print("Wavelet Shape: Gaussian")
			coef3list <- coefNLFiltGG(fHW, nsk, trnct);
			filtSig <- geomavFilt(rsav, trnct, coef3list);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		} else if (identical(wavelet, "GL")) {
			print("Wavelet Shape: lhalf-Gaussian-rhalf-Lorentzian")
			coef3list <- coefNLFiltGL(fHW, nsk, trnct);
			fsavl[[1]] <- geomavFilt(rsav, trnct, coef3list);
		} else { print("No wavelet function was selected: using <GG> as default")
			wavelet <- "GG" # default
			coef3list <- coefNLFiltGG(fHW, nsk, trnct);
			filtSig <- geomavFilt(rsav, trnct, coef3list);
			fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		}
    		## noise amplitude threshold for nonlinear
		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]]<-(sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*1.5*(max(fsavl[[1]])/max(rsav)); # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]]<-(rsav*0+nsav)*1.5*(max(fsavl[[1]])/max(rsav));    # noise amplitude  
    		} 
		wfc1 <- coef3list$fc1; wfc2 <- coef3list$fc2; wfc3 <- coef3list$fc3;
		fwl <- max(c(length(wfc1),length(wfc2),length(wfc3))); # filter length
		fHW <- 0.6*fHW; # updated filtered HW
	} else if (identical(algo, "MAV")) {
		print("Filter: MAV")
		fsavl[[1]]<- mavSmoothing(rsav, fHW0);
		## noise amplitude threshold for MAV
		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]] <- (sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*3*(max(fsavl[[1]])/max(rsav))/sqrt(fHW0); # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]] <- (rsav*0+nsav)*3*(max(fsavl[[1]])/max(rsav))/sqrt(fHW0);  # noise amplitude   
    		}
		fwl <- fHW0; # filter length
		#fHW <- 1.1*fHW; # ??updated filtered HW
	} else {
		print("NO Filter specified: using OLF with Gaussian wavelet as Default")
		wof <- coefOptFiltGG(fHW, nsk, trnct);
		filtSig <- optFilt(rsav, wof);
		fsavl[[1]] <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
		fsavl[[1]][length(fsavl[[1]])] <- 0; # end-point correction 
    		## noise amplitude threshold for OLF
		if (nsdepBS) { #noise dependent on baseline: Poisson distribution
    			namavl[[1]]<-(sqrt(abs(bsav[iir]-cofsav)*nsav)+nsav)*3*(max(fsavl[[1]])/max(rsav))/3.5; # noise amplitude
                } else { # stationary noise (independent of baseline)
    			namavl[[1]]<-(rsav*0+nsav)*3*(max(fsavl[[1]])/max(rsav))/3.5;    # noise amplitude  
    		} 
		fwl <- length(wof); # filter length
		fHW <- 0.9*fHW; # updated filtered HW
	} ### end filter and wavelet options and average spectrum processing

#################### End Mean Spectrum Downsampling/Filtering Block  #########################

#####################  Mean Spectrum Pedestal Removal Block  ###########################

	# Call appropriate pedestal removal procedure

	if (identical(prMethod, "lmMAV")) {

	### Start local minima MAV for Pedestal Removal ###

		print("Pedestal Removal Method: local minima MAV")
		fsavl <- pedRmMAV_call(fsavl,fHW,nf) # remove pedestal in average spectrum

	### End Moving Average Filter for Pedestal Removal ###

		}  else  {
		print("Pedestal Removal Method: none")
		
	} # end pedestal removal options

	if (saveOut) {
		save(savl, savbs, file = "sav_bs_Test.Rdat")
		rm(savl, savbs); gc();
	}

#################### End Mean Spectrum Pedastal Removal Block  #########################
      
	tofResampled <- tofList; naml <- tofList;  ## intitialize resampled and noise amplitude lists
	tof <- list(); length(tof) <- 1; # initialze "processed" 1-spectrum list

	for (k in 1:spectraCount){  ## loop through tofList spectra
		##################  Background Correction Block  ########################

		# Call appropriate baseline subtraction procedure for each spectrum
		tof[1] <- tofList[k]
		names(tof) <- names(tofList)[k]
		if (identical(bsMethod, "RC")) {
	  
		### Start RC-Background Correction ###
	  	
                tof <- chargeBkgr_call(tof,osS,amp,decay);

		### End RC-Background Correction ###

		} else if (identical(bsMethod, "Analytical")) {

		### Start Analytical Background Correction ###

	 	tof <- ExpBkgr_call(tof, bsS,osS,mCH, dataStart)

		### End Analytical-Background Correction ###

		} else if (identical(bsMethod, "ADCoffset")) {

 		tof <- offsetBkgr_call(tof, osS)

		} else { # no BS-subtraction
		tof <- tof;
		
	} # end baseline subtraction options
#tofbs<-tof;

		################## End Background Correction Block  ######################
		
		#################### Downsampling/Filtering Block  #########################

		inputSpectrum <- tof[[1]] # baseline-subtracted
		msrec <- inputSpectrum

		## subtract extra constant baseline
		cOffset<-mean(msrec[(length(msrec)-2000):(length(msrec)-100)]);# PAR: Need to make parameters
		msrec<-msrec-cOffset;

	### Down-sample each spectrum
		rsout<-msResample(msrec,t0, masst, pikw0, lasym0, pfit, minHW, Gnoise);
		#pikw<-rsout$resPikw; fHW <- ceiling(pikw)+1; ## IDS signal HW
		#lasym<-rsout$resAsym; iir <- rsout$resTime; ## already assigned from average above
		srs <- rsout$resSignal; 
		
		ns <- noiseEstmtr2(rsout$resSignal) # STD of noise for this spectrum
		nam <- (srs*0+ns)*3;  # initialize noise amplitude 
		mbs <- tofList[[k]]-tof[[1]]; # spectrum baseline

	## apply appropriate filter
		if ( identical(algo, "none") ){
			   filtSig <- srs;
			   ##noise amplitude threshold for down-sampled signal
    			   if (nsdepBS) { #noise dependent on baseline: Poisson distribution
			   	cof <- mean(tofList[[k]][osS]) 
    			   	nam <- (sqrt(abs(mbs[iir]-cof)*ns)+ns)*3; # noise amplitude
                           } else { # stationary noise (independent of baseline)

    			     	nam<-(filtSig*0+ns)*3;  # noise amplitude  
    			   }
 		}  else if (identical(algo, "matched")) {
			   filtSig <- matchedFilter(srs, mfc);
			   if (identical(wavelet, "GG")) {
        			filtSig <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
    			   } ## noise amplitude threshold for matched filter
			   if (nsdepBS) { #noise dependent on baseline for "matched"
			   	cof <- mean(tofList[[k]][osS]) 
    			   	nam <- (sqrt(abs(mbs[iir]-cof)*ns)+ns)*0.5*max(filtSig)/max(srs); # noise amplitude
                           } else { # stationary noise (independent of baseline)
    			     	nam <- (srs*0+ns)*0.5*max(filtSig)/max(srs);  # noise amplitude   
    			   }

		} else if (identical(algo, "nonLinear")) {
			   filtSig <- geomavFilt(srs, trnct, coef3list);
			   if (identical(wavelet, "GG")) {
        			filtSig <- circShift(filtSig, 1); # correction of extra GG-filter phase shift
    			   }
    			   ## noise amplitude threshold for non-linear filter
			   if (nsdepBS) { #noise dependent on baseline for "nonlinear": Poisson distribution
			   	cof <- mean(tofList[[k]][osS]) 
    			   	nam<-(sqrt(abs(mbs[iir]-cof)*ns)+ns)*1.5*max(filtSig)/max(srs); # noise amplitude
                           } else { # stationary noise (independent of baseline)
    			     	nam<-(srs*0+ns)*1.5*(max(filtSig)/max(srs));    # noise amplitude  
    			   }
  
		} else if (identical(algo, "MAV")) {
			   filtSig <- mavSmoothing(rsout$resSignal, fHW0);
			   ## noise amplitude threshold for MAV
			   if (nsdepBS) { #noise dependent on baseline: Poisson distribution
			   	cof <- mean(tofList[[k]][osS]) 
    			   	nam <- (sqrt(abs(mbs[iir]-cof)*ns)+ns)*3*(max(filtSig)/max(srs))/sqrt(fHW); # noise amplitude
                           } else { # stationary noise (independent of baseline)
    			     	nam <- (srs*0+ns)*3*(max(filtSig)/max(srs))/sqrt(fHW);  # noise amplitude   
    			   }

		} else {  ## We use OLF as default (if no filter option specified)
			   filtSig <- optFilt(srs, wof);
        		   filtSig<- circShift(filtSig, 1); # correction of phase shift	
    			   filtSig[length(filtSig)] <- 0; # end-point correction 
    			   ## noise amplitude threshold for OLF
			   if (nsdepBS) { #noise dependent on baseline: Poisson distribution
			   	cof <- mean(tofList[[k]][osS]) 
    			   	nam<-(sqrt(abs(mbs[iir]-cof)*ns)+ns)*3*(max(filtSig)/max(srs))/3.5; # noise amplitude
                           } else { # stationary noise (independent of baseline)
    			     	nam<-(srs*0+ns)*3*(max(filtSig)/max(srs))/3.5;    # noise amplitude  
    			   } 
		}
	#####################  Pedestal Removal Block  ###########################
		tof[[1]] <- filtSig
		# Call appropriate pedestal removal procedure
#tofbpr<-tof;
		if (identical(prMethod, "lmMAV")) {

		### Start local minima MAV for Pedestal Removal ###

			tof <-pedRmMAV_call(tof,fHW,nf)
			filtSig <- tof[[1]]

		### End Moving Average Filter for Pedestal Removal ###

		}  else  {
			filtSig <- tof[[1]]
		
		} # end pedestal removal options
#tofpr<-tof;

	#################### End Pedastal Removal Block  #########################


		tofResampled[[k]] <- filtSig; naml[[k]] <- nam;

	} ## End of spectrum loop
	fHW <- ceiling(fHW); # make integer
	rm(filtSig, mbs, nam,srs, ns,iir, msrec, inputSpectrum, rsav, nsav); gc();

############## End all tofList processing (bs, ids, filter, pr) Block  #################

	## update master lists for data compression (multiple input files) mode
	if (rdl == 1){
		masterRTOFl <- tofResampled
		masterMeta <- tofListMetaData
		masterNaml <- naml
	} else {
		masterRTOFl <- c(masterRTOFl, tofResampled);
		masterMeta <- rbind(masterMeta, tofListMetaData)
		masterNaml <- c(masterNaml, naml);
	}

	if (monitorMemory){
		memoryUsedBeforeCleanUp<-memory.size(max=FALSE)
		print("Memory Used Before cleaning up data")
		print(memoryUsedBeforeCleanUp)
	}

	## remove large tofList data from memory
	rm(tofList, tofResampled, tof, naml) ## 
	gc()
	if (monitorMemory){
		memoryUsedAfterCleanUp<-memory.size(max=FALSE)
		print("Memory Used After cleaning up data")
		print(memoryUsedAfterCleanUp)
	}
   } # end (not for "average only")
} ## end data reduction for raw TOF data

####### END LOOP: Load TOF spectra and metaData from multiple files and compress #######

## Final (compressed) data structures for all input files (concatinated): 

if (!averageOnly){
	tofResampled <- masterRTOFl; 
	tofListMetaData <- masterMeta;
	naml <- masterNaml;
	spectraName <- names(tofResampled)  # spectra names for all data
	spectraCount <- length(spectraName)  # all spectra from multiple files

	rm(masterNaml, masterRTOFl, masterMeta, sav); gc(); # clean-up memory from "master"-lists

	if (saveOut){
		save(tofResampled,file="tofResampledAllProc_Test.Rdat")
		save(tofListMetaData,file="tofListMeta_Test.Rdat")
	}

#######################  Peak Detection Block  ###########################

	if (nsdepBS) { print("# Baseline-dependent noise threshold used for peak detection #")
	} else { print("# Constant noise threshold used for peak detection #")}


	# Call appropriate peak picking procedure

	if  (identical(ppMethod, "Trivial")) {


	### Start Trivial Peak Picker Block ###

		print("Peak Picking Method: Trivial")
		peakLav <- trivPP_call(fsavl, snr, namavl, fHW, fwl, startPP, dataStart)
		if (identical(alMethod, "AlignToAverage") & noPL) {
		peakL <- peakLav
		} else {
		peakL <- trivPP_call(tofResampled, snr, naml, fHW, fwl, startPP, dataStart)
		}

	### END Trivial Peak Picker block ###

	} else if (identical(ppMethod, "WM")) {

	### Start W&M Peak Picker Block ###

		print("Peak Picking Method: WM")

	## This assumes data is in TOF, explanation: 
	## tBegin <- tStartIn*sampleRate
	## tStop <- tEndIn*sampleRate
	## rsout$resTime[5800] = time index that was mapped to 5800
	## newStart <- rsout$resTime[min(which(rsout$resTime>tStartIn))]

		tBegin <- min(which(rsout$resTime>tStart))
		tStop <- max(which(rsout$resTime<tEnd))	
#		tBegin <- min(which(rsout$resTime>tStartIn))
#	 	tStop <- length(tofResampled[[1]])
#		tStop <- 3000
		threshold <- snr # use the same as for other detection methods
#		Control <- list()
#		kGroup<- 1
#		Control[[1]] <- list(labGroup="labGroup",include = TRUE, width=width,sNR=threshold,sampleRate=1,offset=0,slope=1,tof=tofList)

	 	width<-fHW # where pikw<-rsout$resPikw
#	  	sampleRate <- Control[[kGroup]]$sampleRate
#	  	width <- sampleRate*Control[[kGroup]]$width

#		peakL<-PP_call(tofResampled,base,width,threshold,tBegin,tStop,Noise,NoiseModel,bGetArrays,dataStart)

	### END W&M Peak Picker Block ###

	} else  {
		print("Peak Picking Method: none")
	} # end peak picking options

	### If peak picking is not performed, processing is terminated

	if (exists("peakL") & length(peakL)==spectraCount) {

		##### Update peakList in original time and peakListMetaData ###

		peakList<-list()
		length(peakList) <- spectraCount
		names(peakList) <- names(peakL)
	
		## Generated abbreviated peakListMetaData if it doesn't exist
		if(!exists("peakListMetaData")){peakListMetaData<-plMeta(tofListMetaData)}
	
			## re-Set default offset and scale parameters in metaData, if needed
			#peakListMetaData[ , "peakData.offset" ] <- 0
			#peakListMetaData[ , "peakData.scale" ] <- 1
	
		if (identical(ppMethod, "WM")) { ## NOTE: peak list will be in original time
			peakList <- peakL$peakList
			### Correct peakTimes 	  
			## To interpolate the peakList "time" values fro W&M peak picker 
			##(which allows non-integer position)
			## rsout$resTime is an vector where rsout$resTime[5000] = time index 
			## that was mapped to 5000
			## so peakValue=5000, true peakTime is rsout$resTime[peakValue]
			## so peakValue=5000.4, true peakTime is:
			# rsout$resTime[as.integer(peakValue)] + 0.4*rsout$resTime[as.integer(peakValue+1)]
			for (k in 1:spectraCount) {  ## correct original-time positions
				peakData <- peakList[[k]];
				peakIndex <- as.integer(peakData$Positions)
				peakDiff <- peakData$Positions - peakIndex
				peakDelta <- rsout$resTime[1+peakIndex]-rsout$resTime[peakIndex]
				peakData$Positions <- rsout$resTime[peakIndex] + peakDiff*peakDelta
				peakList[[k]] <- peakData;
				rm(peakData, peakIndex, peakDiff, peakDelta)
			}
		} else {
			peakList <- peakL
		}	

		rm(peakL); gc(); # clean-up memory
		### End Update peakList in original and PeakListMetaData ###

		if (saveOut) {
			save(peakList,file="peakList_Test.Rdat")
			save(peakListMetaData,file="peakListMeta_Test.Rdat")
		}

	}  ## end if peakList detection was performed
###################### End Peak Picking Block  ###########################

#######################  alignment Block  #############################

	alignedPeakListMetaData<-aplMeta(tofListMetaData)

	# Call appropriate alignment procedure

	if (identical(alMethod, "AlignToAverage")) {

	### Start Align to Average Spectra ###

		print("Alignment Method: Align to Average")

		if (align2Training) {   # align to the peak list determined from training
		# provide DOWN-SAMPLED aligned peak list (aplTrain) file for training data
		# NOTE: choose "DS" domain to generate and save when training
			load(aplTrainFileName) # .Rdat file name containing "aplTrain" object
			#print(aplTrain$peaks)
			peaksav<-aplTrain$peaks # peak positions from training	
			#print(peaksav)
                } else {
			peaksav<-peakLav[[1]]$Positions # peak positions from average spectrum
		}
		
		alignedPeakList<- peakAlignAv(peaksav,tofResampled,naml,fHW) 
#		alignedPeakList<-peakAlignAv(peaksav,tofList,naml,fHW)


	### END Align to Average Spectra ###

	}  else if (identical(alMethod, "WM")) {

	### START W & M Aligner ###
	### Note: requires peakList

		print("Alignment Method: WM")
		print("Requires preliminary peakList detection for all spectra")
		width <- fHW
		if (!exists("windmin")) windmin <- min(3*width,maxShift) ## Smallest window to use for binning
		if (!exists("windmax")) windmax <- min(8*width,maxShift) ## Starting window size
		if (exists("peakList") & length(peakList)==spectraCount) {
#		results3 <- AlignWM_call(tofResampled,t0,peakList,peakListMetaData,alignFraction,windmin,windmax,detectFraction,ShiftYes,SlopeYes)
#		alignedPeakList <- results3$alignedPeakList
#		alignedpeakListMetaData <-results3$alignedPeakListMetaData
#		rm(results3); gc();
		} else { print("ERROR: peakList is missing: no alignment could be performed")}

### END W & M Aligner ###

	} else  {
		print("Alignment Method: none")
		
	} # end alignment options

	## scale peak signals back to "integrated" intensities for Gaussian noise 07/01/09
	if (exists("alignedPeakList")){
		if (Gnoise == 1) { # Gaussian noise (sqrt-scaled IDS spectra)
			rrr <- rsavout$resRate
			pks <- alignedPeakList$peaks # in IDS domain
    			for (jd in 1:spectraCount) {
				dmp <- alignedPeakList$data[[jd]]$Intensities;
        			dmp <- dmp*sqrt(rrr[pks]); # square root of down-sampling window
				alignedPeakList$data[[jd]]$Intensities <- as.vector(dmp);
    			}
		}
		rm(rrr,pks,dmp); gc();
	## Convert aligned peak positions to original time

	#### Check domain (revise -- interpolate if using WM-alignment)
		peaksRS <- as.numeric(alignedPeakList$peaks)
		peaksOT <- rsavout$resTime[peaksRS]
		peaksMZ <- masst[peaksOT]
		pks3d <- cbind(peaksRS, peaksOT, peaksMZ);
		avspr1 <- fsavl[[1]];

		if (identical(alPeakDomain, "DS")) {
			dmn <- rep("down-sampled-TOF", spectraCount)
			alignedPeakListMetaData[,"alignedPeakData.domain"] <- dmn

		} else if (identical(alPeakDomain, "MZ")) {
			alignedPeakList$peaks<-peaksMZ
			dmn <- rep("mass-to-charge", spectraCount)
			alignedPeakListMetaData[,"alignedPeakData.domain"] <- dmn
		} else {
			alignedPeakList$peaks<-peaksOT
			dmn <- rep("original-TOF", spectraCount)
			alignedPeakListMetaData[,"alignedPeakData.domain"] <- dmn
		}
		print("---------------------------"); 
		print("---------------------------"); 

		if (saveOut) {
			save(alignedPeakList,file="alignedPeakList_Test.Rdat")
			save(alignedPeakListMetaData,file="alignedPeakListMeta_Test.Rdat")
			save(avspr1, rsavout, pks3d, tms, masst, t0, file = "avSpec_peaks3d_tms_mz_Test.Rdat")
			print("saved ..._Test.Rdat files")
			rm(fsavl, rsavout, pks3d, tms, masst, t0); gc();
		}
		if (save4Mlb) { #output for MatLab if desired 
			if (exists("tofResampled")){ 
			OutputForML(tofResampled, "TRS") 
		}
		if (exists("peakList")){ 
			OutputForML(peakList, "PL") 
		}
		if (exists("alignedPeakList")){ 
			OutputForML(alignedPeakList, "APL") 
		}
		if (exists("inputSpecList")){ # placeholder for "debug"-spectra
			OutputForML(peakList, "ISP") 
		}
		
	}
	rm(peaksRS, peaksOT, peaksMZ, avspr1, i, jd, k, rsout, nsk); gc();
	}
}

#######################  End Alignment Block  #############################
if (exists("alignedPeakList") & exists("alignedPeakListMetaData") & !averageOnly){
return(list(alignedPeakList=alignedPeakList, alignedPeakListMetaData= alignedPeakListMetaData ))  
} else {return(savl)}


} # end of main function 

