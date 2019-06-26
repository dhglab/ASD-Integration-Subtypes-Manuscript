# Script to run differential miRNA expression analyses
# Gokul Ramaswami

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(nlme)

# define directories
dataDir <- "../data/miRNA_expression/"
processedDataDir <- "../processed_data/differential_miRNA_expression/"
metaDataDir <- "../metadata/"
resultsDir <- "../results/differential_miRNA_expression/"

############### FUNCTIONS #################

runlme <- function(thisdat,expression) {
  lme1 <- eval(parse(text=expression));
  ##Get the summary of the model
  smodel = summary(lme1);
  return(smodel)
}

# Function to run DEG analysis
runDEG <- function(datExpr,datMeta,processedDataDir,analysisName) {

	## Set up covariates for linear mixed effects model
	biolrep <- as.numeric(as.factor(datMeta[,"Brain.ID"]))
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age..yrs."])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	RIN <- as.numeric(datMeta[,"RIN"])
	bank <- as.numeric(as.factor(datMeta[,"Brain.bank"]))-1
	ExonProp <- as.numeric(datMeta[,"Proportion.of.exonic.mRNA.reads"])
	Depth <- log10(as.numeric(datMeta[,"Sequencing.depth"]))
	PMI <- as.numeric(datMeta[,"PMI..hrs."])

	varmat <- cbind(condition, age, sex, region, RIN, bank, ExonProp, Depth, PMI)
	rownames(varmat) <- rownames(datMeta)
		
	## LME
	if (!file.exists(paste0(processedDataDir,analysisName,".Rdata"))) {
		Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(varmat))
		rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datExpr)
		colnames(Bmat) <- paste("beta",colnames(varmat),sep=".")
		colnames(SEmat) <- paste("SE",colnames(varmat),sep=".")
		colnames(Pmat) <- paste("p",colnames(varmat),sep=".")

		for (i in 1:nrow(datExpr)) {
    		if (i %% 1000 == 0) {cat(paste("On gene ",i,"\n",sep=""))}
    		thisExpr <- as.numeric(datExpr[i,])
    		expression_model <- paste("lme(thisExpr ~ ",paste(colnames(varmat),collapse=" + "),", rand = ~1|biolrep, data = thisdat)",sep="")
    		designmatrix <- data.frame(thisExpr, varmat, biolrep)

    		lme1.out <- try(runlme(designmatrix,expression_model),silent=F);

    		if (substr(lme1.out[1],1,5)!="Error") {
        		tabOut <- lme1.out$tTable
        		Bmat[i,] <- tabOut[-c(1),"Value"]
        		SEmat[i,] <- tabOut[-c(1),"Std.Error"]
        		Pmat[i,] <- tabOut[-c(1),"p-value"]
    		} else {
        		cat('Error in LME of Gene',rownames(datExpr)[i],"id",'\n')
        		cat('Setting P-value=NA,Beta value=NA, and SE=NA\n');
        		Bmat[i,] <- SEmat[i,] <- Pmat[i,] <- NA;
    		}
		}

		DEG <- cbind(Bmat, SEmat, Pmat)
		DEG <- as.data.frame(DEG)
		DEG[,"p.condition.fdr"] <- p.adjust(DEG[,"p.condition"],method="BH")
		save(DEG, file=paste0(processedDataDir,analysisName,".Rdata"))
	} else {
		load(paste0(processedDataDir,analysisName,".Rdata"))
	}

	DEG = DEG[complete.cases(DEG),]

	DEG
}

###########################################

# Load in miRNA data
load(paste0(dataDir,"miRNA_expression_data_GH.RData")) # datExpr_miRNA_GH

# Load in mRNA metadata
datMeta_miRNA <- read.table(paste0(metaDataDir,"SuppTable1_miRNA.txt"),sep="\t",row.names=1,header=TRUE)


## Disparate ASD subtype vs Control
sampleSubset <- datMeta_miRNA$Outlier=="N" & datMeta_miRNA$Subtype!="Convergent_ASD"
DEG_Disparate <- runDEG(datExpr_miRNA_GH[,sampleSubset], datMeta_miRNA[sampleSubset,], processedDataDir, analysisName="Disparate")

## Convergent ASD subtype vs Control
sampleSubset <- datMeta_miRNA$Outlier=="N" & datMeta_miRNA$Subtype!="Disparate_ASD"
DEG_Convergent <- runDEG(datExpr_miRNA_GH[,sampleSubset], datMeta_miRNA[sampleSubset,], processedDataDir, analysisName="Convergent")
# Write out DEG summary statistics
write.table(DEG_Convergent,file=paste0(resultsDir,"Differential_miRNA_expression_statistics.txt"),sep="\t",quote=F,row.names=T,col.names=T)






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





