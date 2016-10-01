library(WGCNA)
library(data.table)

####Setup metadata
setwd("~/Documents/MigratoryBirds/SampleMetadata/")

barcodes <- read.delim("~/Documents/MigratoryBirds/RAD/readcounts_05.16/YWAR_readcounts.txt")
names(barcodes)
meta <- read.delim("YWAR_RADplates_04.12.16_final.txt")
names(meta)
meta <- meta[!is.na(meta$Plate.),]
meta$barcodes <- barcodes$Barcode[match(meta$Field_Number,barcodes$Sample)]
meta$readcount <- barcodes$NoDups[match(meta$Field_Number,barcodes$Sample)]
meta$totreads <- barcodes$Reads[match(meta$Field_Number,barcodes$Sample)]
rownames(meta) <- paste("YWAR",meta$Plate.,"_sample_",meta$barcodes,sep="")
meta <- meta[,c(1,2,3,4,6,7,8,9,10,12,19,20,21)]


###Genotype data
setwd("~/Documents/MigratoryBirds/RAD/")

data <- data.frame(fread("YWAR_premature.012",na.strings="-1",header=F),row.names=1)
inds <- read.delim("YWAR_premature.012.indv",header=F)
rownames(data) <- inds[,1]
snps <- read.delim("YWAR_premature.012.pos",header=F)

##Filter out individuals with low coverage
goodGT <- apply(data,1,function(x) length(na.omit(x)))
goodGT <- goodGT/ncol(data)
#goodGT_ordered <- goodGT[match(rownames(meta),rownames(data))]
#plot(meta$totreads,goodGT_ordered,ylab="% SNPs genotyped (n=252,404)",
#     xlab="# reads sequenced")

filtered <- data[goodGT>0.5,]
dim(filtered)
rownames(filtered)
badones <- c("YWAR1_sample_ACAGATTC",
			"YWAR2_sample_AGTACAAG",
			"YWAR1_sample_AGTGGTCA",
			"YWAR1_sample_CACCTTAC",
			"YWAR1_sample_CCATCCTC")
filtered <- filtered[-which(rownames(filtered)%in%badones),]
ordermeta <- meta[match(rownames(filtered),rownames(meta)),]
filtered <- filtered[ordermeta$Group!="migrant",]
ordermeta <- ordermeta[ordermeta$Group!="migrant",]

###Filter uninformative SNPs
genNum <- apply(filtered,2,function(x) length(table(x)))
#indCov <- apply(filtered,2,function(x) length(na.omit(x)))
informative <- filtered[,genNum>1]
dim(informative)
MAF <- colMeans(informative,na.rm=T)
pos <- snps[genNum>1,]

###PCA
# for (i in 231000:ncol(informative)) {
	# print (i)
	# snp <- informative[,i]
	# snp[is.na(snp)] <- median(snp,na.rm=T)
	# informative[,i] <- snp
# }


noNA <- na.omit(t(informative[,MAF>0.05]))
dim(noNA)
pca <- prcomp(t(noNA))

par(mar=c(5,4,1,1),cex.lab=1.3,mfrow=c(1,4),mgp=c(1.5,0,0),tck=0.025,fg="white",col.axis="white")
colors <- brewer.pal(12,"Set3")[c(1:8,10:12)]
par(mfrow=c(1,1))
plot(-pca$x[,1],pca$x[,3],pch=19,col=labels2colors(ordermeta$State,colorSeq=colors),
     xlab="",ylab="",cex=1.5,xlim=range(-pca$x[,1]),axes=F,cex.lab=2)
axis(1,pos=0,at=seq(-100,100,by=5),lwd=3,cex.axis=1)
axis(2,pos=0,at=seq(-100,100,by=5),lwd=3,cex.axis=1)
text(1.5,14,"PC3",cex=1.3)
text(13,-1,"PC1",cex=1.3)
box(lwd=2)

#write.table(as.matrix(informative),file="YWAR_median_genotypes.txt",quote=F,row.name=T,col.name=T)
#text(pca$x[,1],pca$x[,3],rownames(ordermeta),cex=0.5)
plot(ordermeta$Long~pca$x[,1],pch=19,col=labels2colors(ordermeta$State,colorSeq=colors),xlab="PC1",ylab="Longitude")
plot(ordermeta$Lat~pca$x[,3],pch=19,col=labels2colors(ordermeta$State,colorSeq=colors),xlab="PC2",ylab="Latitude")
par(mar=c(5,1,1,1))
plot(1:10,1:10,type="n",axes=F,xlab="",ylab="")
legend("left",bty="n",cex=1.5,legend=unique(ordermeta$State),fill=unique(labels2colors(ordermeta$State,colorSeq=colors)))
summary(lm(ordermeta$Lat~pca$x[,2]))
#summary(lm(ordermeta$Lat[ordermeta$State!="AZ"]~pca$x[ordermeta$State!="AZ",3]))

##Training set: Half locations and top n SNPs (based on Near_Town-level frequencies)
aggGen <- aggregate(informative,list(Loc=ordermeta$Near_Town),mean)
aggMeta <- aggregate(ordermeta[,8:9],list(Loc=ordermeta$Near_Town),mean)


##Training set: Half individuals and top n SNPs (48 from PC1 and PC3)
par(mar=c(5,4,1,4),cex.lab=1.3,mfrow=c(2,2),mgp=c(1.5,0.2,0),tck=0.025)
trainers <- sample(1:ncol(noNA),73,replace=F)
n=96
pca <- prcomp(t(noNA[,trainers]))
summary(lm(ordermeta$Long[trainers]~pca$x[,1]))
summary(lm(ordermeta$Lat[trainers]~pca$x[,3]))
PC1snps <- rownames(pca$rotation)[abs(pca$rotation[,1])>min(tail(sort(abs(pca$rotation[,1])),(n/2)))]
PC3snps <- rownames(pca$rotation)[abs(pca$rotation[,3])>min(tail(sort(abs(pca$rotation[,3])),(n/2)))]
goodSNPs <- unique(c(PC1snps,PC3snps))
panel1 <- noNA[goodSNPs,]
pcaP1 <- prcomp(t(panel1))
pcframe <- data.frame(cbind(ordermeta,pcaP1$x))

##predict longitude
longlm <- lm(Long~PC1,data=pcframe[trainers,])
predictors <- pcframe[-trainers,]
predlong <- predict.lm(longlm,newdata=predictors)
plot(predlong~pcframe$Long[-trainers],col=rgb(0.5,0,0.5,alpha=0.5),pch=19,
     ylab="Predicted Longitude",xlab="Observed Longitude")
abline(0,1,lwd=2)
hist(predlong-pcframe$Long[-trainers],xlab="Predicted-Observed Longitude",main="")

##predict latitude
latlm <- lm(Lat~PC2,data=pcframe[trainers,])
predlat <- predict.lm(latlm,newdata=predictors)
plot(predlat~pcframe$Lat[-trainers],col=rgb(0.5,0,0.5,alpha=0.5),pch=19,
     ylab="Predicted Latitude",xlab="Observed Latitude")
abline(0,1,lwd=2)
hist(predlat-pcframe$Lat[-trainers],xlab="Predicted-Observed Latitude",main="")

##calculate distance between predicted & observed
library(geosphere)
distframe <- data.frame(cbind(eLong=predlong,eLat=predlat,
                              oLong=pcframe$Long[-trainers],oLat=pcframe$Lat[-trainers]))
distance <- apply(distframe,1,function(x) distm(x[1:2],x[3:4],fun=distHaversine))
distance <- (distance/1000)*0.621 ##This converts to miles
hist(distance)
summary(distance)


###Create structure file
structurenames <- rep(rownames(informative),each=1)
structurepops <- rep(as.factor(ordermeta$Population),each=1)
structuretypes <- as.numeric(rep(as.factor(ordermeta$Group),each=1))
structurelocs <- as.numeric(rep(as.factor(ordermeta$Near_Town),each=1))

structuregens <- c()
for (i in 1:nrow(informative)) {
  print(i)
	row <- informative[i,]
	row[is.na(row)] <- -9
	hap1 <- row
	hap2 <- row
	hap2[hap2==1] <- 0
	hap2[hap2==2] <- 1
	hap1[hap1==2] <- 1
	structuregens <- rbind(structuregens,hap1,hap2)
}
fakestring <- rep(1,length(structurepops))
structureframe <- data.frame(cbind(structurenames,structurepops,structuretypes,structurelocs,fakestring,fakestring,structuregens))
write.table(structureframe,file="YWAR_premature_structure.str",sep="\t",quote=F,col.names=F,row.names=F)
