# Script to classify samples missing in at least one dataset into the two SNF groups
# Gokul Ramaswami

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(SNFtool)

# define directories
dataDir <- "../data/"
metaDataDir <- "../metadata/"
resultsDir <- "../results/"

############### FUNCTIONS #################

regressCovariates_mRNA <- function(datExpr,datMeta,numSeqStatPCs) {

	# regress out technical covariates (RIN, Brain Bank, Batch, Sequencing Stats) and Biological (Age, Sex, Region) from data 
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	RIN <- as.numeric(datMeta[,"RIN"])
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	
	regvars <- as.data.frame(cbind(condition,age,sex,region,RIN,bank))
	
	# Check if Batch should be included
	if (length(unique(datMeta[,"SeqBatch"])) > 1) {
		if (length(datMeta$SeqBatch[datMeta$SeqBatch=="batch2"])) { # Check if batch2 samples are present
			regvars = cbind(regvars,as.numeric(datMeta[,"SeqBatch"] == "batch2"))
			colnames(regvars)[ncol(regvars)] <- "batch2"
		}
		## don't include batch3 as covariate - did not include in SNF
	}
	
	if (numSeqStatPCs > 0) {
		for (i in 1:numSeqStatPCs) {
			regvars = cbind(regvars,as.numeric(datMeta[,paste0("seqStatPC",eval(i))]))
			colnames(regvars)[ncol(regvars)] <- paste0("seqStatPC",eval(i))
		}
	}	
	
	## Run the regression and make the adjusted FPKM matrix via matrix multiplication                                                                                                                                       

	X <- model.matrix(as.formula(paste0("~",paste(colnames(regvars),collapse="+"))), data = regvars)
	rownames(X) <- rownames(datMeta)
	Y <- datExpr
	
	beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
	
	# regress out beta values from data
	to_regress <- (as.matrix(X[,3:(ncol(regvars)+1)]) %*% (as.matrix(beta[3:(ncol(regvars)+1),])))
	datExpr.reg <- datExpr - t(to_regress)

	## This datExpr.reg is now a technical variable corrected matrix.
	rownames(datExpr.reg) <- rownames(datExpr)
	colnames(datExpr.reg) <- rownames(datMeta)
	
	datExpr.reg

}

regressCovariates_miRNA <- function(datExpr,datMeta) {

	# regress out technical covariates (RIN, Brain Bank, Exon proportion) and Biological (Age, Sex, Region) from data                                                                                                                                  
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	RIN <- as.numeric(datMeta[,"RIN"])
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	ExonProp <- as.numeric(datMeta[,"Proportion.of.exonic.mRNA.reads"])
	SeqDepth <- log10(datMeta[,"Sequencing.depth"])
	PMI <- as.numeric(datMeta[,"PMI"])

	regvars <- as.data.frame(cbind(condition,age,sex,region,RIN,bank,ExonProp,SeqDepth,PMI))
	
	## Run the regression and make the adjusted FPKM matrix via matrix multiplication                                                                                                                                       

	X <- model.matrix(~condition+age+sex+region+RIN+bank+ExonProp+SeqDepth+PMI, data = regvars)
	rownames(X) <- rownames(datMeta)
	Y <- datExpr
	
	beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
	
	# regress out beta values from data
	to_regress <- (as.matrix(X[,3:(ncol(regvars)+1)]) %*% (as.matrix(beta[3:(ncol(regvars)+1),])))
	datExpr.reg <- datExpr - t(to_regress)

	## This datExpr.reg is now a technical variable corrected matrix.
	rownames(datExpr.reg) <- rownames(datExpr)
	colnames(datExpr.reg) <- rownames(datMeta)
	
	datExpr.reg
}

regressCovariates_Ace <- function(datAce,datMeta) {

	# regress out covariates (Age, Region, Brain Bank, Number of Peaks, FRIP) from data                                                                                                                                  
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	FRIP <- as.numeric(datMeta[,"FRIPFract"])
	Dup <- as.numeric(datMeta[,"DupFract"])
	Align <- as.numeric(datMeta[,"AlignFract"])
	CET <- datMeta[,"CET_Fill_Missing"]
	

	regvars <- as.data.frame(cbind(condition, sex, region, age, CET, bank, FRIP, Dup, Align))
	
	## Run the regression and make the adjusted FPKM matrix via matrix multiplication                                                                                                                                       

	X <- model.matrix(as.formula(paste0("~",paste(colnames(regvars),collapse="+"))), data = regvars)
	rownames(X) <- rownames(datMeta)
	Y <- datAce
	
	beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
	
	# regress out beta values from data
	to_regress <- (as.matrix(X[,3:(ncol(regvars)+1)]) %*% (as.matrix(beta[3:(ncol(regvars)+1),])))
	datAce.reg <- datAce - t(to_regress)

	## This datAce.reg is now a technical variable corrected matrix.
	rownames(datAce.reg) <- rownames(datAce)
	colnames(datAce.reg) <- rownames(datMeta)
	
	datAce.reg
}

regressCovariates_Meth <- function(datMeth,datMeta) {

	# regress out technical covariates (Brain Bank, Batch, CET) and Biological (Age, Sex, Region) from data                                                                                                                                  
	condition <- 2-as.numeric(as.factor(datMeta[,"Diagnosis"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	batch <- as.numeric(as.factor(datMeta[,"Batch"]))-1
	CET <- as.numeric(datMeta[,"CET"])

	regvars <- as.data.frame(cbind(condition,age,sex,region,bank,batch,CET))
	
	## Run the regression and make the adjusted FPKM matrix via matrix multiplication                                                                                                                                       

	X <- model.matrix(as.formula(paste0("~",paste(colnames(regvars),collapse="+"))), data = regvars)
	rownames(X) <- rownames(datMeta)
	Y <- datMeth
	
	beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
	
	# regress out beta values from data
	to_regress <- (as.matrix(X[,3:(ncol(regvars)+1)]) %*% (as.matrix(beta[3:(ncol(regvars)+1),])))
	datMeth.reg <- datMeth - t(to_regress)

	## This datMeth.reg is now a technical variable corrected matrix.
	rownames(datMeth.reg) <- rownames(datMeth)
	colnames(datMeth.reg) <- rownames(datMeta)
	
	datMeth.reg
}

Run_PCA <- function (data, nPC = 5) {
  pCdat <- prcomp(t(data), center=FALSE,scale=FALSE);
  topPCs <- pCdat$x[,1:nPC];
  # Calculate variance explained by each PC
  varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
  topVar <- varExp[1:nPC]
  colnames(topPCs) <- paste("Data\n", colnames(topPCs), " (", signif(100 * topVar[1:nPC], 2), "%)", sep = "")
  return(topPCs)
}

### Train logistic regression model on 3 data types
TrainModel_3DataTypes <- function(datMat1,datMat2,datMat3,datName1,datName2,datName3,samples_Group1,samples_Group2,sample_Zscores,training_response) {

	sampleSubset <- colnames(datMat1)[colnames(datMat1) %in% colnames(datMat2)] # datMat1/datMat2 shared
	sampleSubset = sampleSubset[sampleSubset %in% colnames(datMat3)] # datMat1/datMat2/datMat3 shared
	sampleSubset = sampleSubset[!(sampleSubset %in% c(samples_Group1,samples_Group2))] # datMat1/datMat2/datMat3 shared, not in training set
	
	if (length(sampleSubset) > 0) {
		# train classifier
		train <- cbind(sample_Zscores[c(samples_Group1,samples_Group2),c(datName1,datName2,datName3)], training_response) # training set
		model <- glm(training_response ~ get(datName1) + get(datName2) + get(datName3), data = train, family = binomial) # Logistic regression model

		# Try different cutoffs to maximize cross-validation
		cutoffs <- seq(0,1,0.05)
		## exhaustive leave 1 out cross-validation
		accurate_pred <- rep(0,length(cutoffs))
		for (i in 1:length(c(samples_Group1,samples_Group2))) {
			model_CV <- glm(training_response ~ get(datName1) + get(datName2) + get(datName3), data = train[-c(i),], family = binomial)
			model_CV_predict <- predict(model_CV,newdata=sample_Zscores[c(samples_Group1,samples_Group2)[i],c(datName1,datName2,datName3)], type="response")
			for ( j in 1:length(cutoffs) ) {
				cutoff_curr <- cutoffs[j]
				if ( (model_CV_predict < cutoff_curr & training_response[i] == 0) | (model_CV_predict >= cutoff_curr & training_response[i] == 1) ) {
					accurate_pred[j] = accurate_pred[j] + 1
				}
			} 
		}
		names(accurate_pred) <- cutoffs
		cutoff_best <- names(sort(accurate_pred))[length(accurate_pred)]
		
		# Predict missing subtypes
		prediction <- predict(model,newdata=sample_Zscores[sampleSubset,c(datName1,datName2,datName3)],type="response") # test set
		Group1 <- names(prediction[prediction < as.numeric(cutoff_best)])
		Group2 <- names(prediction[prediction >= as.numeric(cutoff_best)])
	
		list(Group1,Group2)
		
	} else {
		list(c(),c())
	}

}

### Train logistic regression model on 2 data types
TrainModel_2DataTypes <- function(datMat1,datMat2,datName1,datName2,samples_Group1,samples_Group2, sample_Zscores,training_response,sample_union_unAssigned) {

	sampleSubset <- colnames(datMat1)[colnames(datMat1) %in% colnames(datMat2)] # datMat1/datMat2 shared
	sampleSubset = sampleSubset[sampleSubset %in% sample_union_unAssigned] # Only look at samples without an assigned group

	if (length(sampleSubset) > 0) {
		# train classifier
		train <- cbind(sample_Zscores[c(samples_Group1,samples_Group2),c(datName1,datName2)], training_response) # training set
		model <- glm(training_response ~ get(datName1) + get(datName2), data = train, family = binomial) # Logistic regression model

		# Try different cutoffs to maximize cross-validation
		cutoffs <- seq(0,1,0.05)
		## exhaustive leave 1 out cross-validation
		accurate_pred <- rep(0,length(cutoffs))
		for (i in 1:length(c(samples_Group1,samples_Group2))) {
			model_CV <- glm(training_response ~ get(datName1) + get(datName2), data = train[-c(i),], family = binomial)
			model_CV_predict <- predict(model_CV,newdata=sample_Zscores[c(samples_Group1,samples_Group2)[i],c(datName1,datName2)], type="response")
			for ( j in 1:length(cutoffs) ) {
				cutoff_curr <- cutoffs[j]
				if ( (model_CV_predict < cutoff_curr & training_response[i] == 0) | (model_CV_predict >= cutoff_curr & training_response[i] == 1) ) {
					accurate_pred[j] = accurate_pred[j] + 1
				}
			} 
		}
		names(accurate_pred) <- cutoffs
		cutoff_best <- names(sort(accurate_pred))[length(accurate_pred)]

		# Predict missing subtypes
		prediction <- predict(model,newdata=sample_Zscores[sampleSubset,c(datName1,datName2)],type="response") # test set
		Group1 <- names(prediction[prediction < as.numeric(cutoff_best)])
		Group2 <- names(prediction[prediction >= as.numeric(cutoff_best)])
	
		list(Group1,Group2)
		
	} else {
		list(c(),c())
	}

}

### Train logistic regression model on 1 data type
TrainModel_1DataType <- function(datMat1,datName1,samples_Group1,samples_Group2,sample_Zscores, training_response,sample_union_unAssigned) {

	sampleSubset <- colnames(datMat1)
	sampleSubset = sampleSubset[sampleSubset %in% sample_union_unAssigned] # Only look at samples without an assigned group

	if (length(sampleSubset) > 0) {
		# train classifier
		train <- as.data.frame(cbind(sample_Zscores[c(samples_Group1,samples_Group2),c(datName1)], training_response)) # training set
		colnames(train)[1] <- datName1
		model <- glm(training_response ~ get(datName1), data = train, family = binomial) # Logistic regression model

		# Try different cutoffs to maximize cross-validation
		cutoffs <- seq(0,1,0.05)
		## exhaustive leave 1 out cross-validation
		accurate_pred <- rep(0,length(cutoffs))
		for (i in 1:length(c(samples_Group1,samples_Group2))) {
			model_CV <- glm(training_response ~ get(datName1), data = train[-c(i),], family = binomial)
			model_CV_predict <- predict(model_CV,newdata=sample_Zscores[c(samples_Group1,samples_Group2)[i],], type="response")
			for ( j in 1:length(cutoffs) ) {
				cutoff_curr <- cutoffs[j]
				if ( (model_CV_predict < cutoff_curr & training_response[i] == 0) | (model_CV_predict >= cutoff_curr & training_response[i] == 1) ) {
					accurate_pred[j] = accurate_pred[j] + 1
				}
			} 
		}
		names(accurate_pred) <- cutoffs
		cutoff_best <- names(sort(accurate_pred))[length(accurate_pred)]

		# Predict missing subtypes
		prediction <- predict(model,newdata=sample_Zscores[sampleSubset,],type="response") # test set
		Group1 <- names(prediction[prediction < as.numeric(cutoff_best)])
		Group2 <- names(prediction[prediction >= as.numeric(cutoff_best)])
	
		list(Group1,Group2)
		
	} else {
		list(c(),c())
	}

}

##################################################################


## mRNA data
# Load in mRNA data
load(paste0(dataDir,"mRNA_expression/","mRNA_expression_data_GH.RData")) # datExpr_mRNA_GH
# Load in mRNA metadata
datMeta_mRNA <- read.table(paste0(metaDataDir,"SuppTable1_mRNA.txt"),sep="\t",row.names=1,header=TRUE)
# Load in iASD vs Control differential expression statistics from Parikshak et al. 2016
datDiff_mRNA_Parik <- read.table(paste0(metaDataDir,"Parikshak_2016/","Parikshak_iASD_vs_CTL_DEG.txt"),sep="\t", header=TRUE,row.names=1)


## miRNA data
# Load in miRNA data
load(paste0(dataDir,"miRNA_expression/","miRNA_expression_data_GH.RData")) # datExpr_miRNA_GH
# Load in mRNA metadata
datMeta_miRNA <- read.table(paste0(metaDataDir,"SuppTable1_miRNA.txt"),sep="\t",row.names=1,header=TRUE)
# Load in ASD vs Control differential expression statistics from Wu et al. 2016
datDiff_miRNA_Wu <- read.table(paste0(metaDataDir,"Wu_2016/","Wu_ASD_vs_CTL_DEG.txt"),sep="\t", header=TRUE,row.names=1)


## DNA methylation data - Promoters
# Load in DNA methylation data
load(paste0(dataDir,"DNA_methylation/","DNA_methylation_data_GH.RData")) # datMeth_Prom_GH
# Load in DNA methylation metadata
datMeta_methylation <- read.table(paste0(metaDataDir,"SuppTable1_methylation.txt"),sep="\t",row.names=1,header=TRUE)
# Load in initial differential methylation statistics: all ASD vs control
datDiff_Meth_Prom <- read.table(paste0(resultsDir,"differential_methylation/","Differential_DNA_methylation_statistics_Promoter_AllASD.txt"), sep="\t",header=TRUE,row.names=1)


## Histone acetylation data
# Load in histone acetylation data
load(paste0(dataDir,"histone_acetylation/","Histone_acetylation_data_GH.RData")) # datAce_GH
# Fix the Acetylation sample names (lowercase the "BA")
colnames(datAce_GH) <- gsub("_BA","_ba",colnames(datAce_GH))
# Load in histone acetylation metadata
datMeta_acetylation <- read.table(paste0(metaDataDir,"SuppTable1_acetylation.txt"),sep="\t",row.names=1,header=TRUE)
# Load in initial differential acetylation statistics: all ASD vs control
datDiff_Ace <- read.table(paste0(resultsDir,"differential_acetylation/","Differential_histone_acetylation_statistics_AllASD.txt"), sep="\t",header=TRUE,row.names=1)


### Adjustment of datasets for technical and biological covariates
## mRNA data
# add seqPCs to Metadata for mRNA
datSeq <- data.matrix(datMeta_mRNA[,c(18:32)])
datSeq[,c(1:8)] = log10(datSeq[,c(1:8)]) # Log transform the large read count numbers
datSeqNorm <- t(scale(datSeq))
PC.datSeq <- prcomp(datSeqNorm);
num_topPCs <- 5 # Number of top seq stat PCs
topPC.datSeq <- PC.datSeq$rotation[,1:num_topPCs];
colnames(topPC.datSeq) <- paste0("SeqSV",c(1:num_topPCs)) ## Compute sequencing stats PCA
for (i in 1:num_topPCs) {
	datMeta_mRNA[,paste0("seqStatPC",eval(i))] <- as.numeric(topPC.datSeq[,i])
}
# Regress out technical and biological covariates
datExpr_mRNA_GH.reg <- regressCovariates_mRNA(datExpr_mRNA_GH,datMeta_mRNA,numSeqStatPCs=5)
# Standard normalize the data
datExpr_mRNA_GH.reg.stdNorm <- t(standardNormalization(t(datExpr_mRNA_GH.reg)))

## miRNA data
# Regress out technical and biological covariates
datExpr_miRNA_GH.reg <- regressCovariates_miRNA(datExpr_miRNA_GH,datMeta_miRNA)
# Standard normalize the data
datExpr_miRNA_GH.reg.stdNorm <- t(standardNormalization(t(datExpr_miRNA_GH.reg)))

## DNA methylation data
# Regress out technical and biological covariates
datMeth_Prom_GH.reg <- regressCovariates_Meth(datMeth_Prom_GH,datMeta_methylation)
# Standard normalize the data
datMeth_Prom_GH.reg.stdNorm <- t(standardNormalization(t(datMeth_Prom_GH.reg)))

## Histone acetylation data
# Regress out technical and biological covariates
datAce_GH.reg <- regressCovariates_Ace(datAce_GH,datMeta_acetylation)
# Standard normalize the data
datAce_GH.reg.stdNorm <- t(standardNormalization(t(datAce_GH.reg)))


## Subset to Differential features
# mRNA DEG (FDR of 10%)
DX_signif_mRNA <- rownames(datDiff_mRNA_Parik[datDiff_mRNA_Parik[,"p.condition.FDR"] < 0.1,])
datExpr_mRNA_GH.reg.stdNorm.DEG <- datExpr_mRNA_GH.reg.stdNorm[DX_signif_mRNA,]
# miRNA DEG (FDR of 10%)
DX_signif_miRNA <- rownames(datDiff_miRNA_Wu[datDiff_miRNA_Wu[,"FDR.adjusted.Pval.diagnosis"] < 0.1,])
datExpr_miRNA_GH.reg.stdNorm.DEG <- datExpr_miRNA_GH.reg.stdNorm[DX_signif_miRNA,]
# DNA methylation promoter DMG (FDR of 10%)
DX_signif_Meth <- rownames(datDiff_Meth_Prom[datDiff_Meth_Prom[,"p.condition.fdr"] < 0.1,])
datMeth_Prom_GH.reg.stdNorm.DMG <- datMeth_Prom_GH.reg.stdNorm[DX_signif_Meth,]
# Histone acetylation DAR (FDR of 20%)
DX_signif_Ace <- rownames(datDiff_Ace[datDiff_Ace[,"p.condition.fdr"] < 0.2,])
datAce_GH.reg.stdNorm.DAR <- datAce_GH.reg.stdNorm[DX_signif_Ace,]


### Calculate scaled PC1 (Z scores) on the differential features
## mRNA
mRNA_DEG_PC1 <- Run_PCA(datExpr_mRNA_GH.reg.stdNorm.DEG)[,1]
# Transform PC1 loadings into Z scores
DEG_Loadings_mRNA_PC1_Z <- as.numeric(scale(as.numeric(mRNA_DEG_PC1)))
names(DEG_Loadings_mRNA_PC1_Z) <- colnames(datExpr_mRNA_GH.reg.stdNorm.DEG)
## miRNA
miRNA_DEG_PC1 <- Run_PCA(datExpr_miRNA_GH.reg.stdNorm.DEG)[,1]
# Transform PC1 loadings into Z scores
DEG_Loadings_miRNA_PC1_Z <- as.numeric(scale(as.numeric(miRNA_DEG_PC1)))
names(DEG_Loadings_miRNA_PC1_Z) <- colnames(datExpr_miRNA_GH.reg.stdNorm.DEG)
## DNA methylation
Meth_DMG_PC1 <- Run_PCA(datMeth_Prom_GH.reg.stdNorm.DMG)[,1]
# Transform PC1 loadings into Z scores
DMG_Loadings_Meth_PC1_Z <- as.numeric(scale(as.numeric(Meth_DMG_PC1)))
names(DMG_Loadings_Meth_PC1_Z) <- colnames(datMeth_Prom_GH.reg.stdNorm.DMG)
## Histone acetylation
Ace_DAR_PC1 <- Run_PCA(datAce_GH.reg.stdNorm.DAR)[,1]
# Transform PC1 loadings into Z scores
DAR_Loadings_Ace_PC1_Z <- as.numeric(scale(as.numeric(Ace_DAR_PC1)))
names(DAR_Loadings_Ace_PC1_Z) <- colnames(datAce_GH.reg.stdNorm.DAR)


## Orient the PC loadings so that control is negative and ASD is positive
# mRNA
mRNA_beta <- summary(lm(DEG_Loadings_mRNA_PC1_Z ~ datMeta_mRNA[,"Diagnosis"]))$coefficients[2,1]
if (mRNA_beta > 0) {
	DEG_Loadings_mRNA_PC1_Z = -1 * DEG_Loadings_mRNA_PC1_Z
}
# miRNA
miRNA_beta <- summary(lm(DEG_Loadings_miRNA_PC1_Z ~ datMeta_miRNA[,"Diagnosis"]))$coefficients[2,1]
if (miRNA_beta > 0) {
	DEG_Loadings_miRNA_PC1_Z = -1 * DEG_Loadings_miRNA_PC1_Z
}
# DNA methylation
Meth_beta <- summary(lm(DMG_Loadings_Meth_PC1_Z ~ datMeta_methylation[,"Diagnosis"]))$coefficients[2,1]
if (Meth_beta > 0) {
	DMG_Loadings_Meth_PC1_Z = -1 * DMG_Loadings_Meth_PC1_Z
}
# Histone acetylation
Ace_beta <- summary(lm(DAR_Loadings_Ace_PC1_Z ~ datMeta_acetylation[,"Diagnosis"]))$coefficients[2,1]
if (Ace_beta > 0) {
	DAR_Loadings_Ace_PC1_Z = -1 * DAR_Loadings_Ace_PC1_Z
}


## Load in SNF group assignments
SNF_groups <- read.table(paste0(resultsDir,"SNF/","SNF_group_assignments.txt"), row.names=1,header=FALSE)
Samps_group_1 <- names(as.matrix(SNF_groups)[SNF_groups[,1]=="SNF_group_1",])
Samps_group_2 <- names(as.matrix(SNF_groups)[SNF_groups[,1]=="SNF_group_2",])


## Combine Z scores together into single matrix
samples_union <- union(names(DAR_Loadings_Ace_PC1_Z),union(names(DMG_Loadings_Meth_PC1_Z), union(names(DEG_Loadings_mRNA_PC1_Z),names(DEG_Loadings_mRNA_PC1_Z))))
## Fill in missing Z-scores with NAs (samples are missing from that dataset)
# mRNA
missing_mRNA <- samples_union[!(samples_union %in% names(DEG_Loadings_mRNA_PC1_Z))]
missing_mRNA_values <- rep(NA,length(missing_mRNA))
DEG_Loadings_mRNA_PC1_Z_filled <- c(DEG_Loadings_mRNA_PC1_Z,missing_mRNA_values)
names(DEG_Loadings_mRNA_PC1_Z_filled) <- c(names(DEG_Loadings_mRNA_PC1_Z),missing_mRNA)
DEG_Loadings_mRNA_PC1_Z_filled = DEG_Loadings_mRNA_PC1_Z_filled[samples_union]
# mRNA
missing_miRNA <- samples_union[!(samples_union %in% names(DEG_Loadings_miRNA_PC1_Z))]
missing_miRNA_values <- rep(NA,length(missing_miRNA))
DEG_Loadings_miRNA_PC1_Z_filled <- c(DEG_Loadings_miRNA_PC1_Z,missing_miRNA_values)
names(DEG_Loadings_miRNA_PC1_Z_filled) <- c(names(DEG_Loadings_miRNA_PC1_Z),missing_miRNA)
DEG_Loadings_miRNA_PC1_Z_filled = DEG_Loadings_miRNA_PC1_Z_filled[samples_union]
# DNA methylation
missing_Meth <- samples_union[!(samples_union %in% names(DMG_Loadings_Meth_PC1_Z))]
missing_Meth_values <- rep(NA,length(missing_Meth))
DMG_Loadings_Meth_PC1_Z_filled <- c(DMG_Loadings_Meth_PC1_Z,missing_Meth_values)
names(DMG_Loadings_Meth_PC1_Z_filled) <- c(names(DMG_Loadings_Meth_PC1_Z),missing_Meth)
DMG_Loadings_Meth_PC1_Z_filled = DMG_Loadings_Meth_PC1_Z_filled[samples_union]
# Histone acetylation
missing_Ace <- samples_union[!(samples_union %in% names(DAR_Loadings_Ace_PC1_Z))]
missing_Ace_values <- rep(NA,length(missing_Ace))
DAR_Loadings_Ace_PC1_Z_filled <- c(DAR_Loadings_Ace_PC1_Z,missing_Ace_values)
names(DAR_Loadings_Ace_PC1_Z_filled) <- c(names(DAR_Loadings_Ace_PC1_Z),missing_Ace)
DAR_Loadings_Ace_PC1_Z_filled = DAR_Loadings_Ace_PC1_Z_filled[samples_union]
## Join into a data frame
sample_Zscores <- as.data.frame(cbind(DEG_Loadings_mRNA_PC1_Z_filled,DEG_Loadings_miRNA_PC1_Z_filled,DMG_Loadings_Meth_PC1_Z_filled, DAR_Loadings_Ace_PC1_Z_filled))
colnames(sample_Zscores) <- c("mRNA","miRNA","Meth","Ace")


### Fit logistic regression classifiers to extrapolate groups of samples with missing data types
### Samples with all 4 data types: Samps_group_1, Samps_group_2
## Set up training factor response
training_response <- c(rep(0,length(Samps_group_1)),rep(1,length(Samps_group_2))) 

### Samples with 3 data types: 1. mRNA, miRNA, Meth   2. mRNA, miRNA, Ace   3. mRNA, Meth, Ace   4. miRNA, Meth, Ace

######### 1. mRNA, miRNA, Meth
group3_type1 <- TrainModel_3DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datExpr_miRNA_GH.reg.stdNorm.DEG,datMeth_Prom_GH.reg.stdNorm.DMG,datName1="mRNA",datName2="miRNA",datName3="Meth",Samps_group_1,Samps_group_2,sample_Zscores,training_response)

######### 2. mRNA, miRNA, Ace
group3_type2 <- TrainModel_3DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datExpr_miRNA_GH.reg.stdNorm.DEG,datAce_GH.reg.stdNorm.DAR,datName1="mRNA",datName2="miRNA",datName3="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response)

######### 3. mRNA, Meth, Ace
group3_type3 <- TrainModel_3DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datMeth_Prom_GH.reg.stdNorm.DMG,datAce_GH.reg.stdNorm.DAR,datName1="mRNA",datName2="Meth",datName3="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response)

######### 4. miRNA, Meth, Ace
group3_type4 <- TrainModel_3DataTypes(datExpr_miRNA_GH.reg.stdNorm.DEG,datMeth_Prom_GH.reg.stdNorm.DMG,datAce_GH.reg.stdNorm.DAR,datName1="miRNA",datName2="Meth",datName3="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response)

### Samples with 2 data types: 1. mRNA, miRNA   2. mRNA, Meth   3. mRNA, Ace   4. miRNA, Meth   5. miRNA, Ace   6. Meth, Ace
## Filter samples that have already been assigned
sample_union <- rownames(sample_Zscores)
sample_union_unAssigned <- sample_union[!(sample_union %in% c(Samps_group_1,Samps_group_2,group3_type1[[1]],group3_type1[[2]],group3_type2[[1]],group3_type2[[2]],group3_type3[[1]],group3_type3[[2]],group3_type4[[1]],group3_type4[[2]]))]

######### 1. mRNA, miRNA
group2_type1 <- TrainModel_2DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datExpr_miRNA_GH.reg.stdNorm.DEG,datName1="mRNA",datName2="miRNA",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 2. mRNA, Meth
group2_type2 <- TrainModel_2DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datMeth_Prom_GH.reg.stdNorm.DMG,datName1="mRNA",datName2="Meth",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 3. mRNA, Ace
group2_type3 <- TrainModel_2DataTypes(datExpr_mRNA_GH.reg.stdNorm.DEG,datAce_GH.reg.stdNorm.DAR,datName1="mRNA",datName2="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 4. miRNA, Meth
group2_type4 <- TrainModel_2DataTypes(datExpr_miRNA_GH.reg.stdNorm.DEG,datMeth_Prom_GH.reg.stdNorm.DMG,datName1="miRNA",datName2="Meth",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 5. miRNA, Ace
group2_type5 <- TrainModel_2DataTypes(datExpr_miRNA_GH.reg.stdNorm.DEG,datAce_GH.reg.stdNorm.DAR,datName1="miRNA",datName2="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 6. Meth, Ace
group2_type6 <- TrainModel_2DataTypes(datMeth_Prom_GH.reg.stdNorm.DMG,datAce_GH.reg.stdNorm.DAR,datName1="Meth",datName2="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)


### Samples with 1 data type: 1. mRNA   2. miRNA   3. Meth   4. Ace
## Filter samples that have already been assigned
sample_union_unAssigned = sample_union_unAssigned[!(sample_union_unAssigned %in% c(group2_type1[[1]],group2_type1[[2]],group2_type2[[1]],group2_type2[[2]],group2_type3[[1]],group2_type3[[2]],group2_type4[[1]],group2_type4[[2]],group2_type5[[1]],group2_type5[[2]],group2_type6[[1]],group2_type6[[2]]))]

######### 1. mRNA
group1_type1 <- TrainModel_1DataType(datExpr_mRNA_GH.reg.stdNorm.DEG,datName1="mRNA",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 2. miRNA
group1_type2 <- TrainModel_1DataType(datExpr_miRNA_GH.reg.stdNorm.DEG,datName1="miRNA",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 3. Meth
group1_type3 <- TrainModel_1DataType(datMeth_Prom_GH.reg.stdNorm.DMG,datName1="Meth",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

######### 4. Ace
group1_type4 <- TrainModel_1DataType(datAce_GH.reg.stdNorm.DAR,datName1="Ace",Samps_group_1,Samps_group_2,sample_Zscores,training_response,sample_union_unAssigned)

#### Combine together all the group assignments
Group1_combined <- c(Samps_group_1,group3_type1[[1]],group3_type2[[1]],group3_type3[[1]],group3_type4[[1]],group2_type1[[1]],group2_type2[[1]],group2_type3[[1]],group2_type4[[1]],group2_type5[[1]],group2_type6[[1]],group1_type1[[1]],group1_type2[[1]],group1_type3[[1]],group1_type4[[1]])
Group2_combined <-c(Samps_group_2,group3_type1[[2]],group3_type2[[2]],group3_type3[[2]],group3_type4[[2]],group2_type1[[2]],group2_type2[[2]],group2_type3[[2]],group2_type4[[2]],group2_type5[[2]],group2_type6[[2]],group1_type1[[2]],group1_type2[[2]],group1_type3[[2]],group1_type4[[2]])


## Translate the group assignments into subtype assignments
# Control = Control
# Group 1 ASD = Disparate Subtype
# Group 2 ASD = Convergent Subtype

ASD_samps <- union(rownames(datMeta_mRNA[datMeta_mRNA[,"Diagnosis"]=="ASD",]),union(rownames(datMeta_miRNA[datMeta_miRNA[,"Diagnosis"]=="ASD",]),union(rownames(datMeta_methylation[datMeta_methylation[,"Diagnosis"]=="ASD",]),rownames(datMeta_acetylation[datMeta_acetylation[,"Diagnosis"]=="ASD",]))))

subtype_assign <- cbind(c(Group1_combined,Group2_combined),c(rep("Group1",length(Group1_combined)),rep("Group2",length(Group2_combined))))
subtype_assign[,2][!(subtype_assign[,1] %in% ASD_samps)] <- "Control"
subtype_assign[,2][subtype_assign[,2] =="Group1"] <- "Disparate_ASD"
subtype_assign[,2][subtype_assign[,2] =="Group2"] <- "Convergent_ASD"

## Write out subtype assignments
write.table(subtype_assign,file=paste0(resultsDir,"SNF/","Classified_sample_assignments.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")





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
# [1] SNFtool_2.2.1
# 
# loaded via a namespace (and not attached):
# [1] heatmap.plus_1.3







