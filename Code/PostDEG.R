####LIBRARIES------------------------------------------------------------------------------------------------####
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

####LOAD DATA BASES------------------------------------------------------------------------------------------####
setwd('..')
path <- getwd()
setwd(paste(path, "/Data Bases", sep = ""))
load("MatAdip.RData")
load("MatMuscle.RData")

####SIGNIFICANT GENES----------------------------------------------------------------------------------------####
resadip <- results(MatAdip, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
table(resadip$padj < 0.05)
resadipSig <- subset(resadip, padj < 0.05)

test <- results(MatAdip, lfcThreshold = 1, format = "GRanges")



resmuscle <- results(MatMuscle, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2"))
table(resmuscle$padj < 0.05)
resmuscleSig <- subset(resmuscle, padj < 0.05)

####ANNOTATIONS----------------------------------------------------------------------------------------------####
#rownames(resadipSig) <- gsub("\\..*","", DEGsHSCCMP$gene_id)

geneids <- rownames(resadipSig)



test  <- mapIds(org.Hs.eg.db,
                     keys=geneids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


ensembl            <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
DEGsHSCCMP[,18]    <- getBM(filters="ensembl_gene_id" , attributes= c("external_gene_name"), 
                            values = DEGsHSCCMP$gene_id, mart = ensembl)
