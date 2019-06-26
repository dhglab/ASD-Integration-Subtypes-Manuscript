# Script to run differential DNA methylation analyses
# Gokul Ramaswami

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(nlme)

# define directories
dataDir <- "../data/DNA_methylation/"
processedDataDir <- "../processed_data/differential_methylation/"
metaDataDir <- "../metadata/"
resultsDir <- "../results/differential_methylation/"

############### FUNCTIONS #################

runlme <- function(thisdat,expression) {
  lme1 <- eval(parse(text=expression));
  ##Get the summary of the model
  smodel = summary(lme1);
  return(smodel)
}

# Function to run DMG analysis
runDMG <- function(datMeth,datMeta,processedDataDir,analysisName) {

	## Set up covariates for linear mixed effects model
	biolrep <- as.numeric(as.factor(datMeta[,"newBrainBankID"]))
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"BrainRegion_update"]))-1
	batch <- as.numeric(datMeta[,"Batch"])-1
	bank <- as.numeric(as.factor(datMeta[,"BrainCentre.M"]))-1
	age <- as.numeric(datMeta[,"Age"])
	CET <- as.numeric(datMeta[,"CET"])

	varmat <- cbind(condition, sex, region, batch, bank, age, CET)
	rownames(varmat) <- rownames(datMeta)
		
	## LME
	if (!file.exists(paste0(processedDataDir,analysisName,".Rdata"))) {
		Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datMeth),ncol=ncol(varmat))
		rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datMeth)
		colnames(Bmat) <- paste("beta",colnames(varmat),sep=".")
		colnames(SEmat) <- paste("SE",colnames(varmat),sep=".")
		colnames(Pmat) <- paste("p",colnames(varmat),sep=".")

		for (i in 1:nrow(datMeth)) {
    		if (i %% 1000 == 0) {cat(paste("On gene ",i,"\n",sep=""))}
    		thisMeth <- as.numeric(datMeth[i,])
    		expression_model <- paste("lme(thisMeth ~ ",paste(colnames(varmat),collapse=" + "),", rand = ~1|biolrep, data = thisdat)",sep="")
    		designmatrix <- data.frame(thisMeth, varmat, biolrep)

    		lme1.out <- try(runlme(designmatrix,expression_model),silent=F);

    		if (substr(lme1.out[1],1,5)!="Error") {
        		tabOut <- lme1.out$tTable
        		Bmat[i,] <- tabOut[-c(1),"Value"]
        		SEmat[i,] <- tabOut[-c(1),"Std.Error"]
        		Pmat[i,] <- tabOut[-c(1),"p-value"]
    		} else {
        		cat('Error in LME of Gene',rownames(datMeth)[i],"id",'\n')
        		cat('Setting P-value=NA,Beta value=NA, and SE=NA\n');
        		Bmat[i,] <- SEmat[i,] <- Pmat[i,] <- NA;
    		}
		}

		DMG <- cbind(Bmat, SEmat, Pmat)
		DMG <- as.data.frame(DMG)
		DMG[,"p.condition.fdr"] <- p.adjust(DMG[,"p.condition"],method="BH")
		save(DMG, file=paste0(processedDataDir,analysisName,".Rdata"))
	} else {
		load(paste0(processedDataDir,analysisName,".Rdata"))
	}

	DMG = DMG[complete.cases(DMG),]

	DMG
}

###########################################

# Load in DNA methylation data
load(paste0(dataDir,"DNA_methylation_data_GH.RData")) # datMeth_Prom_GH, datMeth_GB_GH

# Load in DNA methylation metadata
datMeta_methylation <- read.table(paste0(metaDataDir,"SuppTable1_methylation.txt"),sep="\t",row.names=1,header=TRUE)


### Promoters
## All ASD vs Control
sampleSubset <- datMeta_methylation$Outlier=="N"
DMG_Prom_All <- runDMG(datMeth_Prom_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="Prom_All")
## Disparate ASD subtype vs Control
sampleSubset <- datMeta_methylation$Outlier=="N" & datMeta_methylation$Subtype!="Convergent_ASD"
DMG_Prom_Disparate <- runDMG(datMeth_Prom_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="Prom_Disparate")
## Convergent ASD subtype vs Control
sampleSubset <- datMeta_methylation$Outlier=="N" & datMeta_methylation$Subtype!="Disparate_ASD"
DMG_Prom_Convergent <- runDMG(datMeth_Prom_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="Prom_Convergent")

### Gene Bodies
## All ASD vs Control
sampleSubset <- datMeta_methylation$Outlier=="N"
DMG_GB_All <- runDMG(datMeth_GB_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="GB_All")
## Disparate ASD subtype vs Control
sampleSubset <- datMeta_methylation$Outlier=="N" & datMeta_methylation$Subtype!="Convergent_ASD"
DMG_GB_Disparate <- runDMG(datMeth_GB_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="GB_Disparate")
## Convergent ASD subtype vs Control
sampleSubset <- datMeta_methylation$Outlier=="N" & datMeta_methylation$Subtype!="Disparate_ASD"
DMG_GB_Convergent <- runDMG(datMeth_GB_GH[,sampleSubset], datMeta_methylation[sampleSubset,], processedDataDir, analysisName="GB_Convergent")

# Write out DMG summary statistics
write.table(DMG_Prom_All,file=paste0(resultsDir,"Differential_DNA_methylation_statistics_Promoter_AllASD.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(DMG_Prom_Convergent,file=paste0(resultsDir,"Differential_DNA_methylation_statistics_Promoter_Convergent.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(DMG_GB_All,file=paste0(resultsDir,"Differential_DNA_methylation_statistics_GB_AllASD.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(DMG_GB_Convergent,file=paste0(resultsDir,"Differential_DNA_methylation_statistics_GB_Convergent.txt"),sep="\t",quote=F,row.names=T,col.names=T)





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







