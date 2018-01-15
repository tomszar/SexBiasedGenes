#This code will remove genes with low counts from the database
#Because it will do that for the whole dataset, it need a bug amount of memory
#Run this in cluster

####LIBRARIES------------------------------------------------------------------------------------------------####
library(dplyr)
library(DESeq2)
library(edgeR)

####LOAD DATA BASES------------------------------------------------------------------------------------------####
samples <- read.csv("storage/home/tug156/code/SubjSample.csv")
counts  <- read.csv("storage/home/tug156/code/GeneReads.csv", row.names = 1)

#Remove low counts
cpm1   <- cpm(counts[1:100,])
#keep those rows (genes) with cpm values greater than 145
keep   <- rowSums(cpm1>1)>=median(summary(samples$SMTSD))
sum(keep)
counts <- counts[keep,]
remove(cpm1)

#Save file
write.csv(counts, file = "HighGeneReads.csv", quote = FALSE, col.names = FALSE)

####DESEQ DATABASE-------------------------------------------------------------------------------------------####
#conditions           <- samples[,c(2:4,6)]
#rownames(conditions) <- samples[,1]
#conditions$GENDER    <- as.factor(conditions$GENDER)
#DesMat <- DESeqDataSetFromMatrix(countData = counts, 
#                                 colData   = conditions,
#                                 design    = ~ GENDER + AGE)
#DesMat

#musclecond           <- musclesamp[,c(3:4)]
#rownames(musclecond) <- musclesamp[,1]
#musclecond$GENDER    <- as.factor(musclecond$GENDER)
#MatMuscle <- DESeqDataSetFromMatrix(countData = musclecount, 
#                                    colData   = musclecond,
#                                    design    = ~ GENDER + AGE)
#MatMuscle

####PCA------------------------------------------------------------------------------------------------------####
#rld <- rlogTransformation(MatAdip, blind=TRUE)
#plotPCA(rld, intgroup=c("GENDER", "AGE"))

####NORMALIZATION AND DE GENES-------------------------------------------------------------------------------####
#MatAdip <- estimateSizeFactors( MatAdip ) 
#sizeFactors(MatAdip)
#MatAdip <- estimateDispersions( MatAdip )
#plotDispEsts( MatAdip , cex.lab=1.5 , ymin=0.01)
#MatAdip <- nbinomWaldTest(MatAdip)
#resadip <- results(MatAdip, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
#table(resadip$padj < 0.05)
#resadipSig <- subset(resadip, padj < 0.05)
#plotCounts(MatAdip, gene=which.min(resadip$padj), intgroup="GENDER")
#save(MatAdip, file="MatAdip.RData")

#MatMuscle <- estimateSizeFactors( MatMuscle ) 
#sizeFactors(MatMuscle)
#MatMuscle <- estimateDispersions( MatMuscle )
#plotDispEsts( MatMuscle , cex.lab=1.5 , ymin=0.01)
#MatMuscle <- nbinomWaldTest(MatMuscle)
#resmuscle <- results(MatMuscle, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
#table(resmuscle$padj < 0.05)
#resmuscleSig <- subset(resmuscle, padj < 0.05)
#plotCounts(MatMuscle, gene=which.min(resmuscle$padj), intgroup="GENDER")
#save(MatMuscle, file="MatMuscle.RData")

####LIMMA ANALYSIS-------------------------------------------------------------------------------------------####
#adipdesign <- model.matrix(~ as.factor(adiposesamp[,3]))
#colnames(adipdesign) <- c("inter", "sexf")
#contr.matrix <- makeContrasts(
#  FvsM = inter-sexf, 
#  levels = colnames(adipdesign))
#contr.matrix

#v    <- voom(adiposecount, adipdesign, plot=TRUE)
#vfit <- lmFit(v, adipdesign)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#vfit <- eBayes(vfit)
#plotSA(vfit, main="Final model: Mean???variance trend")
#vfit <- treat(vfit, lfc=0)
#summary(decideTests(vfit))
