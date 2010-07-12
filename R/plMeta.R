"plMeta" <-
function(tofListMetaData) {

###-----------------------------------------------------------------
###
### Function plMeta uses tofListMetaData to generate peakListMetaData
###
### INPUT
###	toflistMetaData
###
### OUTPUT
###	peakListMetaData 
###----------------------------------------------------------------- 
## Create MetaData
  
PLmetaDataFields <- c("acquisitionInfo.refId", 
"acquisitionInfo.spotName", "peakData.domain", "peakData.offset", 
"peakData.scale", "replicateNumber", "sampleInfo.groupName", 
"sampleInfo.sampleDescription", "sampleInfo.sampleName", "sampleInfo.sampleSource"
)

commonFieldNames <- intersect(colnames(tofListMetaData), PLmetaDataFields)
numberOfMetaDataFields <- length(PLmetaDataFields)
spectraNameList<-dimnames(tofListMetaData)[1]
numSpectra<-dim(tofListMetaData)[1]
spectraName<-rownames(tofListMetaData)

peakListMetaData<-array(0,c(numSpectra,numberOfMetaDataFields))
dimnames(peakListMetaData)[2]<-list(PLmetaDataFields)
dimnames(peakListMetaData)[1]<-spectraNameList

## assign new (not-common) fields here:
peakListMetaData[,"peakData.scale"]<-tofListMetaData[,"timeOfFlightData.scale"]
peakListMetaData[,"peakData.offset"]<-tofListMetaData[,"timeOfFlightData.offset"]
peakListMetaData[,"peakData.domain"]<-tofListMetaData[,"timeOfFlightData.domain"]

for (i in 1:length(commonFieldNames)) { ## assign common fields
	cfi <- commonFieldNames[i]
	peakListMetaData[,cfi]<-tofListMetaData[,cfi]
}

return(peakListMetaData)
}

