##LIBRARIES
library(WGCNA)
library(DESeq2)
library(edgeR)

options(stringsAsFactors = FALSE);
# Caution: skip this line if you run RStudio or other third-party R environments.
enableWGCNAThreads()

##DATABASES
setwd('..')
path <- getwd()
setwd(paste(path, "/Results/GTEx", sep = ""))
samples <- read.csv("SubjSample.csv")
counts  <- read.csv(gzfile("GeneReads.csv.gz"), row.names=1)

#TESTING ONE TISSUE
SamplesByTissue <- split(samples, samples$SMTSD)
GeneCounts <- counts[,colnames(counts) %in% SamplesByTissue[[1]]$SAMPID] #Select reads of the selected tissue
keep       <- rowSums(cpm(GeneCounts)>1) >= min(summary(as.factor(SamplesByTissue[[1]]$GENDER)))
symbols    <- as.character(counts[keep,1])
GeneCounts <- GeneCounts[keep,] #Removing low counts from GeneCounts file
insample   <- SamplesByTissue[[1]]$SAMPID %in% colnames(GeneCounts)
cond       <- SamplesByTissue[[1]][insample, 3:4] #Create covariate conditions (Gender and age)
rownames(cond) <- SamplesByTissue[[1]][insample, 1] #Adding sample ID as rowname 
cond$GENDER    <- as.factor(cond$GENDER) #making gender a factor
#Order cond and GeneCounts to match each other
cond       <- cond[order(rownames(cond)), ]
GeneCounts <- GeneCounts[, order(colnames(GeneCounts))]
GeneCounts <- varianceStabilizingTransformation(as.matrix(GeneCounts))
GeneCounts <- t(GeneCounts)

#Choosing soft threshold beta
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft    <- pickSoftThreshold(GeneCounts, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#One-step network construction for large datasets
setwd(paste(path, "/Results/WGCNA", sep = ""))
bwnet = blockwiseModules(GeneCounts,
                         power = 7, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "AdiposeTOM-blockwise",
                         verbose = 3)
table(bwnet$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
