#!/bin/bash
#First move to the database folder where imputed genotypes are
cd ../DataBases/ImputedGenotypes
mkdir croppedfiles

#Then, for each plink file, cropp the chromosome file using the GeneRanges
generanges="../../Results/WGCNA/Adipose/GeneRanges.txt"
for file in *.bed
do
	inputname=$(echo $file | cut -f 1 -d '.' | tr '/' '\\')
	echo $inputname
	counter="${inputname##*_}"
	echo $counter
	outname="cropped_chr${counter}"
	./plink --bfile $inputname --extract range $(echo $generanges) --maf --make-bed --out $outname
done

mv cropped*.* croppedfiles/
