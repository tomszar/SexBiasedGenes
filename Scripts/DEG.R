####LIBRARIES------------------------------------------------------------------------------------------------####
library(ggplot2)
library(dplyr)
library(DESeq2)
library(edgeR)
library(limma)

####LOAD DATA BASES------------------------------------------------------------------------------------------####
setwd('..')
path <- getwd()
setwd(paste(path, "/DataBases", sep = ""))
samples <- read.csv("SubjSampleSimple.csv")
counts  <- read.table("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct", skip=2, header=T, row.names=1)

####CLEAN DATA-----------------------------------------------------------------------------------------------####
samples <- droplevels(filter(samples, SMTSD != '' & SMTSD != "Cells - Leukemia cell line (CML)"))
#Remove tissue with one sex
for(l in 1:length(levels(samples$SMTSD))){
  dat <- filter(samples, SMTSD == levels(samples$SMTSD)[l])
  if(mean(dat$GENDER) == 1 |  mean(dat$GENDER) == 2){
    print(levels(samples$SMTSD)[l])
    samples <- filter(samples, SMTSD != levels(samples$SMTSD)[l])
  }
}
remove(dat)

samples <- droplevels(samples)
counts  <- counts[ , colnames(counts) %in% as.character(samples$SAMPID)]
samples <- samples[as.character(samples$SAMPID) %in% colnames(counts), ]

#Samples count
#Order table
samples <- within(samples, SMTSD <- factor(SMTSD, levels=names(sort(table(SMTSD), 
                                                                    decreasing=TRUE))))

median(summary(samples$SMTSD))
ggplot(samples, aes(SMTSD)) + geom_bar(aes(fill=as.factor(GENDER)), position = position_stack(reverse = TRUE)) + 
  coord_flip() + geom_hline(yintercept = median(summary(samples$SMTSD))) + 
  theme(legend.position = "top")

#Remove low counts
cpm1   <- cpm(counts)
#keep those rows (genes) with cpm values greater than 145
keep   <- rowSums(cpm1>1)>=median(summary(samples$SMTSD))
sum(keep)
counts <- counts[keep,]
remove(cpm1)

####DESEQ DATABASE-------------------------------------------------------------------------------------------####
adiposecond           <- adiposesamp[,c(3:4)]
rownames(adiposecond) <- adiposesamp[,1]
adiposecond$GENDER    <- as.factor(adiposecond$GENDER)
MatAdip <- DESeqDataSetFromMatrix(countData = adiposecount, 
                                 colData   = adiposecond,
                                 design    = ~ GENDER + AGE)
MatAdip

musclecond           <- musclesamp[,c(3:4)]
rownames(musclecond) <- musclesamp[,1]
musclecond$GENDER    <- as.factor(musclecond$GENDER)
MatMuscle <- DESeqDataSetFromMatrix(countData = musclecount, 
                                    colData   = musclecond,
                                    design    = ~ GENDER + AGE)
MatMuscle

####PCA------------------------------------------------------------------------------------------------------####
rld <- rlogTransformation(MatAdip, blind=TRUE)
plotPCA(rld, intgroup=c("GENDER", "AGE"))

####NORMALIZATION AND DE GENES-------------------------------------------------------------------------------####
MatAdip <- estimateSizeFactors( MatAdip ) 
sizeFactors(MatAdip)
MatAdip <- estimateDispersions( MatAdip )
plotDispEsts( MatAdip , cex.lab=1.5 , ymin=0.01)
MatAdip <- nbinomWaldTest(MatAdip)
resadip <- results(MatAdip, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
table(resadip$padj < 0.05)
resadipSig <- subset(resadip, padj < 0.05)
plotCounts(MatAdip, gene=which.min(resadip$padj), intgroup="GENDER")
save(MatAdip, file="MatAdip.RData")

MatMuscle <- estimateSizeFactors( MatMuscle ) 
sizeFactors(MatMuscle)
MatMuscle <- estimateDispersions( MatMuscle )
plotDispEsts( MatMuscle , cex.lab=1.5 , ymin=0.01)
MatMuscle <- nbinomWaldTest(MatMuscle)
resmuscle <- results(MatMuscle, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
table(resmuscle$padj < 0.05)
resmuscleSig <- subset(resmuscle, padj < 0.05)
plotCounts(MatMuscle, gene=which.min(resmuscle$padj), intgroup="GENDER")
save(MatMuscle, file="MatMuscle.RData")

####LIMMA ANALYSIS-------------------------------------------------------------------------------------------####
adipdesign <- model.matrix(~ as.factor(adiposesamp[,3]))
colnames(adipdesign) <- c("inter", "sexf")
contr.matrix <- makeContrasts(
  FvsM = inter-sexf, 
  levels = colnames(adipdesign))
contr.matrix

v    <- voom(adiposecount, adipdesign, plot=TRUE)
vfit <- lmFit(v, adipdesign)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
vfit <- eBayes(vfit)
plotSA(vfit, main="Final model: Mean???variance trend")
vfit <- treat(vfit, lfc=0)
summary(decideTests(vfit))
