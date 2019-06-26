# Script to run SNF
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
	condition <- 2-as.numeric(as.factor(datMeta[,"ASD.CTL"]))
	age <- as.numeric(datMeta[,"Age"])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"RegionID"]))-1
	RIN <- as.numeric(datMeta[,"RIN"])
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	
	regvars <- as.data.frame(cbind(condition,age,sex,region,RIN,bank))
	
	# Check if Batch should be included
	if (length(unique(datMeta[,"SeqBatch"])) > 1) {
		if (length(datMeta$SeqBatch[datMeta$SeqBatch=="batch2"])) { # Check if batch2 samples are present
			regvars = cbind(regvars,as.numeric(datMeta[,"SeqBatch"] == "batch2"))
			colnames(regvars)[ncol(regvars)] <- "batch2"
		}
		## don't include batch3 as covariate - only 1 sample
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
	age <- as.numeric(datMeta[,"Age..yrs."])
	sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
	region <- as.numeric(as.factor(datMeta[,"Region"]))-1
	RIN <- as.numeric(datMeta[,"RIN"])
	bank <- as.numeric(as.factor(datMeta[,"Brain.bank"]))-1
	ExonProp <- as.numeric(datMeta[,"Proportion.of.exonic.mRNA.reads"])
	SeqDepth <- log10(datMeta[,"Sequencing.depth"])
	PMI <- as.numeric(datMeta[,"PMI..hrs."])

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
	region <- as.numeric(as.factor(datMeta[,"Region_Fixed"]))-1
	bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
	peakNum <- as.numeric(datMeta[,"PeakNum"])
	FRIP <- as.numeric(datMeta[,"FRIPFract"])
	Dup <- as.numeric(datMeta[,"DupFract"])
	Align <- as.numeric(datMeta[,"AlignFract"])
	CET <- datMeta[,"CET_Filled"]
	

	regvars <- as.data.frame(cbind(condition, sex, region, age, CET, bank, peakNum, FRIP, Dup, Align))
	
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
	region <- as.numeric(as.factor(datMeta[,"BrainRegion_update"]))-1
	bank <- as.numeric(as.factor(datMeta[,"BrainCentre.M"]))-1
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

Run_SNF <- function(datExpr_mRNA, datExpr_miRNA, datMeth, datAce, K, alpha, T, C, datMeta) {

	# SNF - clustering
	Dist_mRNA <- SNFtool::dist2(as.matrix(t(datExpr_mRNA)),as.matrix(t(datExpr_mRNA)));
	Dist_miRNA <- SNFtool::dist2(as.matrix(t(datExpr_miRNA)),as.matrix(t(datExpr_miRNA)));
	Dist_Meth <- SNFtool::dist2(as.matrix(t(datMeth)),as.matrix(t(datMeth)));
	Dist_Ace <- SNFtool::dist2(as.matrix(t(datAce)),as.matrix(t(datAce)));
	
	# Affinity matrices
	W_mRNA <- affinityMatrix(Dist_mRNA, K, alpha)
	W_miRNA <- affinityMatrix(Dist_miRNA, K, alpha)
	W_Meth <- affinityMatrix(Dist_Meth, K, alpha)
	W_Ace <- affinityMatrix(Dist_Ace, K, alpha)
	# Fuse graphs together
	W_Fused = SNF(list(W_mRNA,W_miRNA,W_Meth,W_Ace), K, T)
	# Clustering
	group_Fused <- spectralClustering(W_Fused, C); 
	
	# Orient SNF group 2 to be the "convergent ASD" group to match manuscript
	group_Fused_1_CTL <- length(group_Fused[group_Fused == 1 & datMeta[,"ASD.CTL"] == "CTL"])
	group_Fused_2_CTL <- length(group_Fused[group_Fused == 2 & datMeta[,"ASD.CTL"] == "CTL"])
	if (group_Fused_2_CTL > group_Fused_1_CTL) {
		group_Fused = (C + 1) - group_Fused
	}
	
	group_Fused
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


## Subset mRNA, miRNA, Meth, Ace datasets to common samples
shared_mRNA_Meth_samps <- intersect(colnames(datExpr_mRNA_GH),colnames(datMeth_Prom_GH))
shared_mRNA_miRNA_Meth_samps <- intersect(shared_mRNA_Meth_samps,colnames(datExpr_miRNA_GH))
shared_mRNA_miRNA_Meth_Ace_samps <- intersect(shared_mRNA_miRNA_Meth_samps, colnames(datAce_GH))

## Filter datasets to common samples
# mRNA
datExpr_mRNA_filtered <- datExpr_mRNA_GH[,shared_mRNA_miRNA_Meth_Ace_samps]
datMeta_mRNA_filtered <- datMeta_mRNA[shared_mRNA_miRNA_Meth_Ace_samps,]
# miRNA
datExpr_miRNA_filtered <- datExpr_miRNA_GH[,shared_mRNA_miRNA_Meth_Ace_samps]
datMeta_miRNA_filtered <- datMeta_miRNA[shared_mRNA_miRNA_Meth_Ace_samps,]
# Meth
datMeth_Prom_filtered <- datMeth_Prom_GH[,shared_mRNA_miRNA_Meth_Ace_samps]
datMeta_methylation_filtered <- datMeta_methylation[shared_mRNA_miRNA_Meth_Ace_samps,]
# Ace
datAce_filtered <- datAce_GH[,shared_mRNA_miRNA_Meth_Ace_samps]
datMeta_acetylation_filtered <- datMeta_acetylation[shared_mRNA_miRNA_Meth_Ace_samps,]


### Adjustment of datasets for technical and biological covariates
## mRNA data
# add seqPCs to Metadata for mRNA
datSeq <- data.matrix(datMeta_mRNA_filtered[,c(18:32)])
datSeq[,c(1:8)] = log10(datSeq[,c(1:8)]) # Log transform the large read count numbers
datSeqNorm <- t(scale(datSeq))
PC.datSeq <- prcomp(datSeqNorm);
num_topPCs <- 5 # Number of top seq stat PCs
topPC.datSeq <- PC.datSeq$rotation[,1:num_topPCs];
colnames(topPC.datSeq) <- paste0("SeqSV",c(1:num_topPCs)) ## Compute sequencing stats PCA
for (i in 1:num_topPCs) {
	datMeta_mRNA_filtered[,paste0("seqStatPC",eval(i))] <- as.numeric(topPC.datSeq[,i])
}
# Regress out technical and biological covariates
datExpr_mRNA_filtered.reg <- regressCovariates_mRNA(datExpr_mRNA_filtered,datMeta_mRNA_filtered,numSeqStatPCs=5)
# Standard normalize the data
datExpr_mRNA_filtered.reg.stdNorm <- t(standardNormalization(t(datExpr_mRNA_filtered.reg)))

## miRNA data
# Regress out technical and biological covariates
datExpr_miRNA_filtered.reg <- regressCovariates_miRNA(datExpr_miRNA_filtered,datMeta_miRNA_filtered)
# Standard normalize the data
datExpr_miRNA_filtered.reg.stdNorm <- t(standardNormalization(t(datExpr_miRNA_filtered.reg)))

## DNA methylation data
# Regress out technical and biological covariates
datMeth_Prom_filtered.reg <- regressCovariates_Meth(datMeth_Prom_filtered,datMeta_methylation_filtered)
# Standard normalize the data
datMeth_Prom_filtered.reg.stdNorm <- t(standardNormalization(t(datMeth_Prom_filtered.reg)))

## Histone acetylation data
# Regress out technical and biological covariates
datAce_filtered.reg <- regressCovariates_Ace(datAce_filtered,datMeta_acetylation_filtered)
# Standard normalize the data
datAce_filtered.reg.stdNorm <- t(standardNormalization(t(datAce_filtered.reg)))


### Similarity Network Fusion 
## Set all the SNF parameters:
K = 20;		# number of neighbors
alpha = 0.5;  	# hyperparameter
T = 15; 	# Number of Iterations)
C = 2; # Number of clusters

## Subset to Differential features
# mRNA DEG (FDR of 10%)
DX_signif_mRNA <- rownames(datDiff_mRNA_Parik[datDiff_mRNA_Parik[,"p.condition.FDR"] < 0.1,])
datExpr_mRNA_filtered.reg.stdNorm.DEG <- datExpr_mRNA_filtered.reg.stdNorm[DX_signif_mRNA,]
# miRNA DEG (FDR of 10%)
DX_signif_miRNA <- rownames(datDiff_miRNA_Wu[datDiff_miRNA_Wu[,"FDR.adjusted.Pval.diagnosis"] < 0.1,])
datExpr_miRNA_filtered.reg.stdNorm.DEG <- datExpr_miRNA_filtered.reg.stdNorm[DX_signif_miRNA,]
# DNA methylation promoter DMG (FDR of 10%)
DX_signif_Meth <- rownames(datDiff_Meth_Prom[datDiff_Meth_Prom[,"p.condition.fdr"] < 0.1,])
datMeth_Prom_filtered.reg.stdNorm.DMG <- datMeth_Prom_filtered.reg.stdNorm[DX_signif_Meth,]
# Histone acetylation DAR (FDR of 20%)
DX_signif_Ace <- rownames(datDiff_Ace[datDiff_Ace[,"p.condition.fdr"] < 0.2,])
datAce_filtered.reg.stdNorm.DAR <- datAce_filtered.reg.stdNorm[DX_signif_Ace,]

# Run SNF
SNF_assignments <- Run_SNF(datExpr_mRNA_filtered.reg.stdNorm.DEG,datExpr_miRNA_filtered.reg.stdNorm.DEG, datMeth_Prom_filtered.reg.stdNorm.DMG,datAce_filtered.reg.stdNorm.DAR,K, alpha, T, C,datMeta_mRNA_filtered)
SNF_assignments = paste0("SNF_group_",SNF_assignments)
names(SNF_assignments) <- colnames(datExpr_mRNA_filtered.reg.stdNorm.DEG)

# Write out SNF group assignments
write.table(as.data.frame(cbind(names(SNF_assignments),SNF_assignments)),file=paste0(resultsDir,"SNF/","SNF_group_assignments.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=F)




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





