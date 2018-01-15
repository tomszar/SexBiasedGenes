####DIFFERENTIAL GENE EXPRESSION SCRIPT--------------------------------------------------------------------------------------------------

#####Libraries needed----------------------------------------------------------------------------------------------------------------####
library(ggplot2)
library(dplyr)
library(DESeq2) 
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db) #The last FOUR are Bioconductor packages. To install them, you'll first need to run the following commands:
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")
#biocLite("edgeR")
#biocLite("AnnotationDbi")
#biocLite("org.Hs.eg.db")

####LOAD DATA BASES------------------------------------------------------------------------------------------------------------------#### 
samples <- read.csv("SubjSampleSimple.csv") #Database with sample information
counts  <- read.table("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct", skip=2, header=T, row.names=1) #Database with read counts
#This will take a while to load

####SUBSET SAMPLES-------------------------------------------------------------------------------------------------------------------####
brainsamp    <- droplevels(filter(samples, SMTSD == 'ENTER TISSUE HERE')) #Selecting the appropriate tissue
braincount   <- select(counts, which(colnames(counts) %in% as.character(brainsamp$SAMPID))) #Filtering count data by tissue
sampincount  <- as.data.frame(cbind(colnames(braincount)))
colnames(sampincount) <- "SAMPID"
brainsamp  <- left_join(sampincount, brainsamp) 

####CLEAN DATA-----------------------------------------------------------------------------------------------------------------------####
#Remove low counts
cpm        <- cpm(braincount) #Transforming raw count to count per mil (cpm) values
keep.exprs <- rowSums(cpm>1)>=ncol(cpm) #Selecting rows with cpm values higher than number of columns (at least one cpm per sample)
braincount <- braincount[keep.exprs, ] #Subsetting count database

####DESEQ DATABASE-------------------------------------------------------------------------------------------------------------------####
braincond           <- brainsamp[,c(3:4)] #Create covariate conditions (Gender and age)
rownames(braincond) <- brainsamp[,1] 
braincond$GENDER    <- as.factor(braincond$GENDER) 
MatBrain <- DESeqDataSetFromMatrix(countData = braincount,  #Creating deseq database with counts and covariates (braincond)
                                  colData   = braincond,
                                  design    = ~ GENDER + AGE)
MatBrain 

####NORMALIZATION AND DE GENES-------------------------------------------------------------------------------------------------------####
#This is the core of differential gene expression analysis. Parameters can be changed according to different criteria (whether we want
#to be more or less stringent in differential gene expression estimation)
MatBrain <- estimateSizeFactors( MatBrain ) #The following three lines of code are the proper estimation
MatBrain <- estimateDispersions( MatBrain ) #This takes a while
MatBrain <- nbinomWaldTest(MatBrain)
resbrain <- results(MatBrain, pAdjustMethod = "BH", contrast = c("GENDER", "1", "2")) #Export the results, with contrast between males and 
#females. P adjustment can be changed, this case is BH method
table(resbrain$padj < 0.05) #Number of differentially expressed genes
resbrainpSig <- subset(resbrain, padj < 0.05) #Subset results to those differentially expressed only
plotCounts(MatBrain, gene=which.min(resbrain$padj), intgroup="GENDER")

####ANNOTATIONS----------------------------------------------------------------------------------------------------------------------####
geneids <- rownames(resbrainpSig) #Getting gene IDs, and removing the number after the dor
geneids <- gsub("\\..*","",geneids)

test  <- mapIds(org.Hs.eg.db,
                keys=geneids,
                column="SYMBOL",
                keytype="ENSEMBL",
                multiVals="first")


ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
test    <- getBM(filters="ensembl_gene_id" , attributes= c("external_gene_name"), 
                            values = geneids, mart = ensembl)

