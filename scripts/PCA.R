library(WGCNA) ##This package needs to be downloaded from bioconductor
library(data.table)
library(RColorBrewer)

####Read in the following files:
## meta: metadata - make sure sample names match sample names from barcode file used in demultiplexing
## data: SNP table output from pipeline (vcftools -012 format)
## snps: SNP positions (this is also output from vcftools -012)
## inds: sample names (also vcftools -012 output)
setwd("~/Documents/DIRECTORY")

meta <- read.delim("meta.txt",header=T)
data <- data.frame(fread("SNP_MATRIX.012",na.strings="-1",header=F),row.names=1)
inds <- read.delim("SNP_INDIVIDUALS.012.indv",header=F)
rownames(data) <- inds[,1]
snps <- read.delim("SNP_POSITIONS.012.pos",header=F)
ordermeta <- meta[match(rownames(data),meta$Sample),]

###Filter uninformative SNPs
##Filters SNPs with only one genotype present
genNum <- apply(data,2,function(x) length(table(x)))
informative <- data[,genNum>1]
dim(informative)

###OPTIONAL: Filter out "bad" individuals: These may have been outliers or samples you think there was a problem with. You can also use the last two lines to filter out potential migrants as long as that is marked in the metadata
# badones <- c("SAMPLE1",
			# "SAMPLE2",
			# "SAMPLE2")
# informative <- informative[-which(rownames(informative)%in%badones),]
# ordermeta <- meta[match(rownames(informative),rownames(meta)),]
# informative <- informative[ordermeta$Group!="migrant",]
# ordermeta <- ordermeta[ordermeta$Group!="migrant",]

##Fill in missing values with median genotype for that SNP
for (i in 1:ncol(informative)) {
	print (i)
	snp <- informative[,i]
	snp[is.na(snp)] <- median(snp,na.rm=T)
	informative[,i] <- snp
}


####PCA
pca <- prcomp(informative)
colors <- c(brewer.pal(12,"Set3"),brewer.pal(8,"Set1"),brewer.pal(8,"Set2"))
par(mfrow=c(1,2),tck=0.025,mar=c(4,4,1,1),mgp=c(1.5,0.1,0))
plot(pca$x[,1],pca$x[,2],pch=19,col=labels2colors(ordermeta$State,colorSeq=colors),
     cex=1.5,cex.lab=1.7,xlab="PC1",ylab="PC2")
plot(1:10,1:10,type="n",axes=F,xlab="",ylab="")
legend("left",bty="n",cex=1.2,legend=unique(ordermeta$State),fill=unique(labels2colors(ordermeta$State,colorSeq=colors)))

