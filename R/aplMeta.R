"aplMeta" <-
function(tofListMetaData) {

###-----------------------------------------------------------------
###
### Function alpMeta uses tofListMetaData to generate alignedPeakListMetaData
###
### INPUT
###	toflistMetaData
###
### OUTPUT
###	alignedPeakListMetaData 
###----------------------------------------------------------------- 
## Create MetaData
  
APLmetaDataFields = c("acquisitionInfo.refId", 
"acquisitionInfo.spotName", "alignedPeakData.domain", "alignedPeakData.offset", 
"alignedPeakData.scale", "replicateNumber", "sampleInfo.groupName", 
"sampleInfo.sampleDescription", "sampleInfo.sampleName", "sampleInfo.sampleSource")

commonFieldNames <- intersect(colnames(tofListMetaData), APLmetaDataFields)

numberOfMetaDataFields <- length(APLmetaDataFields)

spectraNameList<-dimnames(tofListMetaData)[1]
numSpectra<-dim(tofListMetaData)[1]
spectraName<-rownames(tofListMetaData)

alignedPeakListMetaData<-array(0,c(numSpectra,numberOfMetaDataFields))
dimnames(alignedPeakListMetaData)[2]<-list(APLmetaDataFields)
dimnames(alignedPeakListMetaData)[1]<-spectraNameList

## assign new (not-common) fields here:
alignedPeakListMetaData[,"alignedPeakData.scale"]<-tofListMetaData[,"timeOfFlightData.scale"]
alignedPeakListMetaData[,"alignedPeakData.offset"]<-tofListMetaData[,"timeOfFlightData.offset"]
alignedPeakListMetaData[,"alignedPeakData.domain"]<-tofListMetaData[,"timeOfFlightData.domain"]

for (ij in 1:length(commonFieldNames)) { ## assign common fields
	cfi <- commonFieldNames[ij]
	alignedPeakListMetaData[,cfi]<-tofListMetaData[,cfi]
}


return(alignedPeakListMetaData)
}

