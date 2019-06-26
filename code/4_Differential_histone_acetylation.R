# Script to run differential histone acetylation analyses
# Gokul Ramaswami

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(nlme)

# define directories
dataDir <- "../data/histone_acetylation/"
processedDataDir <- "../processed_data/differential_acetylation/"
metaDataDir <- "../metadata/"
resultsDir <- "../results/differential_acetylation/"

############### FUNCTIONS #################

runlme <- function(thisdat,expression) {
  lme1 <- eval(parse(text=expression));
  ##Get the summary of the model
  smodel = summary(lme1);
  return(smodel)
}

# Function to run DAR analysis
runDAR <- function(datAce,datMeta,processedDataDir,analysisName) {

	## Set up covariates for linear mixed effects model
	biolrep <- as.numeric(as.factor(datMeta[,"BrainID"]))
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region_Fixed"]))-1
	CET <- as.numeric(datMeta[,"CET_Filled"])
	Bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	peakNum <- as.numeric(datMeta[,"PeakNum"])
	FRIP <- as.numeric(datMeta[,"FRIPFract"])
	Dup <- as.numeric(datMeta[,"DupFract"])
	Align <- as.numeric(datMeta[,"AlignFract"])

	varmat <- cbind(condition, age, sex, region, CET, Bank, peakNum, FRIP, Dup, Align)
	rownames(varmat) <- rownames(datMeta)
		
	## LME
	if (!file.exists(paste0(processedDataDir,analysisName,".Rdata"))) {
		Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datAce),ncol=ncol(varmat))
		rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datAce)
		colnames(Bmat) <- paste("beta",colnames(varmat),sep=".")
		colnames(SEmat) <- paste("SE",colnames(varmat),sep=".")
		colnames(Pmat) <- paste("p",colnames(varmat),sep=".")

		for (i in 1:nrow(datAce)) {
    		if (i %% 1000 == 0) {cat(paste("On peak ",i,"\n",sep=""))}
    		thisAce <- as.numeric(datAce[i,])
    		expression_model <- paste("lme(thisAce ~ ",paste(colnames(varmat),collapse=" + "),", rand = ~1|biolrep, data = thisdat)",sep="")
    		designmatrix <- data.frame(thisAce, varmat, biolrep)

    		lme1.out <- try(runlme(designmatrix,expression_model),silent=F);

    		if (substr(lme1.out[1],1,5)!="Error") {
        		tabOut <- lme1.out$tTable
        		Bmat[i,] <- tabOut[-c(1),"Value"]
        		SEmat[i,] <- tabOut[-c(1),"Std.Error"]
        		Pmat[i,] <- tabOut[-c(1),"p-value"]
    		} else {
        		cat('Error in LME of Peak',rownames(datAce)[i],"id",'\n')
        		cat('Setting P-value=NA,Beta value=NA, and SE=NA\n');
        		Bmat[i,] <- SEmat[i,] <- Pmat[i,] <- NA;
    		}
		}

		DAR <- cbind(Bmat, SEmat, Pmat)
		DAR <- as.data.frame(DAR)
		DAR[,"p.condition.fdr"] <- p.adjust(DAR[,"p.condition"],method="BH")
		save(DAR, file=paste0(processedDataDir,analysisName,".Rdata"))
	} else {
		load(paste0(processedDataDir,analysisName,".Rdata"))
	}

	DAR = DAR[complete.cases(DAR),]

	DAR
}

###########################################

# Load in histone acetylation data
load(paste0(dataDir,"Histone_acetylation_data_GH.RData")) # datAce_GH

# Load in histone acetylation metadata
datMeta_acetylation <- read.table(paste0(metaDataDir,"SuppTable1_acetylation.txt"),sep="\t",row.names=1,header=TRUE)


## All ASD vs Control
sampleSubset <- datMeta_acetylation$Outlier=="N"
DAR_All <- runDAR(datAce_GH[,sampleSubset], datMeta_acetylation[sampleSubset,], processedDataDir, analysisName="All")
## Disparate ASD subtype vs Control
sampleSubset <- datMeta_acetylation$Outlier=="N" & datMeta_acetylation$Subtype!="Convergent_ASD"
DAR_Disparate <- runDAR(datAce_GH[,sampleSubset], datMeta_acetylation[sampleSubset,], processedDataDir, analysisName="Disparate")
## Convergent ASD subtype vs Control
sampleSubset <- datMeta_acetylation$Outlier=="N" & datMeta_acetylation$Subtype!="Disparate_ASD"
DAR_Convergent <- runDAR(datAce_GH[,sampleSubset], datMeta_acetylation[sampleSubset,], processedDataDir, analysisName="Convergent")



# Write out DAR summary statistics
write.table(DAR_All,file=paste0(resultsDir,"Differential_histone_acetylation_statistics_AllASD.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(DAR_Convergent,file=paste0(resultsDir,"Differential_histone_acetylation_statistics_Convergent.txt"),sep="\t",quote=F,row.names=T,col.names=T)





#####################  SESSION INFO ###########################

# R Under development (unstable) (2016-03-02 r70268)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.6 (El Capitan)
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] nlme_3.1-125
# 
# loaded via a namespace (and not attached):
# [1] grid_3.3.0      lattice_0.20-33





