library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

#Find working directory, and set to bam directory if not there
getwd()	
#setwd("target dir")
list.files()

#Set "search rate" for bam files (how many reads at a time it will look for)
files<-list.files(pattern=".bam$");
file.exists(files)
bamfiles <- BamFileList(files, yieldSize=200000)


#Prepare GTF file for counting, change file name depending on species
seqinfo(bamfiles[1])

#Create transcript db and exon coordinates for each gene
#setwd("GTF directory")
(txdb <- makeTxDbFromGFF("/external/mgmt3/genome/scratch/Sibille/dnewton/PITT_tetrad_SCT_RNAseq/Reference/genes/Homo_sapiens.GRCh38.98.gtf", format="gtf", circ_seqs=character())) 

(ebg <- exonsBy(txdb, by="gene"))
#Setwd(Back to bam directory)

#this is the counting step - run overnight 
register(MulticoreParam(workers=10))

#IMPORTANT FOR SINGLE V.S. PAIRED END#####################################################################################
se <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="Union",singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE)

se <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="Union",singleEnd=TRUE,ignore.strand=TRUE,fragments=FALSE)
#####################################################################################################################

#Check the proper # of columns are in the se
se

dim(se) 
assayNames(se)
head(assay(se), 3) 

HT <- colSums(assay(se)) # this provides the total number of reads mapped to exons for each sample

# outputs the object RR to a .csv, this also provides the order of the samples - important for integrating with the pTable, make sure your pTable has the names and order in its rows as in the assay file. Save as pTable.csv and upload to .bam directory of cluster.
write.csv(HT,file="HTDepth_trim.csv")

#reordered pTable based on order of HTDepth
rowRanges(se)
pTable<-read.csv("pTable.csv")
colnames(se)
rownames(pTable)
rownames(pTable) <- pTable[,5]
rownames(pTable) <- gsub("_", "-", rownames(pTable))

#subset pTable just for current flow cell - if performing a single run, not batches, this is unnecessary
pTable <- pTable[rownames(pTable) %in% colnames(se),]
pTable <- pTable[order(match(rownames(pTable), colnames(se))),]


#Check for proper subsetting
rownames(pTable) == colnames(se) #returned 79 TRUEs, good to go
pTable<-DataFrame(pTable)
(colData(se) <- DataFrame(pTable))
colData(se) #merge was successful, all trait data is in se


#######Get feature data and download from Ensembl BioMart, currently done manually on https://www.ensembl.org/biomart/martview, future implementation with biomaRt package would be better
featureData <- rownames(se)
sum(duplicated(featureData)) #if this is not 0, there are duplicates and something went weird above, check reference genome steps
write.csv(featureData,file="featureData.csv")
#############

#no Tx-specific data right now, just merging
mData <- read.table("martquery_1104155809_616.txt", sep = ",", header=TRUE)
mData <- mData[!duplicated(mData$Gene.stable.ID),]
noAnnotations <- setdiff(featureData, mData$Gene.stable.ID)
withAnnotations <- intersect(featureData, mData$Gene.stable.ID)
mTable1 <- mData[which(withAnnotations %in% mData$Gene.stable.ID),]
length(noAnnotations)

#Skip if length of NoAnnotations is 0, this steps reconsiles non annotated Gene_Stable_IDs, more of an issue for Mouse than Human reference genomes
###
mTable2 <- as.data.frame(matrix(nrow=0, ncol=8)) #nrow is number of rows in noAnnotations, and ncol is number of cols in other metadata table
colnames(mTable2) <- colnames(mTable1)
mTable2$Gene.stable.ID <- noAnnotations
colnames(mTable1)==colnames(mTable2) # check that these are the same
x <- rbind(mTable1,mTable2)
dim(x) #check that the number of genes (rows) and metadata (columns) is correct
sum(duplicated(x$Gene.stable.ID)) #if this is not 0, there are duplicates
y <- x[match(as.factor(featureData), x$Gene.stable.ID),]
sum(as.factor(featureData)== y$Gene.stable.ID) # returns 60623 - good to go
rownames(y) <- y$Gene.stable.ID
###

#Check dimensions line up
sum(as.character(rownames(y)) == as.character(rownames(mcols(se)))) #good to go, they match up

#Add metadata for genes
mcols(se) <- cbind(mcols(se), y)

se #metadata is present


se$Cell.Type <- relevel(as.factor(se$Cell.Type), "Pyr-L2n3")

#Useful group variable for some clustering analyses, not necessarily for actual differential expression
se$group <- as.factor(paste0(se$Subject.Group, se$Cell.Type)) #make group variable for design
se$group <- gsub("-", "_", se$group)
se$group <- as.factor(se$group)


dds0 <- DESeqDataSet(se, design = ~ group)

#### Simple cutoff for inital screening purposes of 30 raw reads across all samples, removes spurious "genes"
keep <- rowSums(counts(dds0)) > 30
sum(keep)
dds0 <- dds0[keep,]

#dds0 <- DESeqDataSet(se, design = ~ group) #Optional - may need to re-design after filtering

#use variane stablizing transformation for normalization, much faster than rlog transformation and controls heteroscedastity very well (plots below)
register(MulticoreParam(workers=10))
vsd <- vst(dds0,blind=FALSE)


par(mfrow = c( 1, 2 ) )
library(ggplot2)
write.csv(plotPCA(vsd, intgroup = c("Cell.Type"), returnData=TRUE), file="Novaseq_PCAdata(min30count)_vst.csv")

pdf("plot_PCA_Celltype_Novaseq5(min10count)_vst.pdf")
plotPCA(vsd, intgroup = c("Cell.Type"))
dev.off()

#Get normalized data

register(MulticoreParam(workers=10))
write.csv(assay(vsd), file="Normalized_vst_Counts_Novaseq5.csv")

dds <- DESeq(dds0, parallel=TRUE)

ncounts <- counts(dds0, normalized=TRUE)
write.csv(ncounts,"Normalized_Counts_Novaseq5.csv")

######checking variance plots of each normalization method

library(vsn)
ntd <- normTransform(dds)


untransplot <- meanSdPlot(ncounts)

ntdplot <- meanSdPlot(assay(ntd)) 
ntdplot2 <- ntdplot$gg + scale_y_continuous(limits=c(0,8))

vstplot <- meanSdPlot(assay(vsd))
vstplot2 <- vstplot$gg + scale_y_continuous(limits=c(0,8))


#Export differentia=
pdf("Novaseq_untransformed_variance.pdf")
untransplot$gg
dev.off()

pdf("Novaseq_log2_transformed_variance.pdf")
ntdplot2
dev.off()

pdf("Novaseq_vst_variance.pdf")
vstplot2
dev.off()

pdf("Novaseq_vst_noscale_variance.pdf")
vstplot$gg
dev.off()


##############Removing outliers and merging runs 1-5, BAM files were the result of trimming: 5'-3bp(11 including 8bp adapter), 3'-15bp
###Specific to PITT Tetrad dataset - 5 NovaSeq flowcells
load("se_nova1.rData")
load("se_nova2.rData")
load("se_nova3.rData")
load("se_nova4.rData")
load("se_nova5.rData")

outliers <- read.csv("MiSeq Outlier List current.csv")
pTable <- read.csv("pTable.csv")

se_p1_cl <- se_nova1[,!(colnames(se_nova1) %in% outliers$Subject[outliers$Pool==1])]
se_p2_cl <- se_nova2[,!(colnames(se_nova2) %in% outliers$Subject[outliers$Pool==2])]

seFullData <- cbind(se_p1_cl, se_p2_cl, se_nova3, se_nova4, se_nova5)
###


seFullData$group <- as.factor(paste0(seFullData$Subject.Group, seFullData$Cell.Type)) #make group variable for design
seFullData$group <- gsub("-", "_", seFullData$group)
seFullData$group <- as.factor(seFullData$group)


dds0 <- DESeqDataSet(seFullData, design = ~ group)

keep <- rowSums(counts(dds0)) > 30
sum(keep)
dds0 <- dds0[keep,]



#dds0 <- DESeqDataSet(se, design = ~ group) #Optional, may need to re-design after filtering
register(MulticoreParam(workers=10))

vsd <- vst(dds0,blind=FALSE)

par(mfrow = c( 1, 2 ) )
library("ggplot2")
write.csv(plotPCA(vsd, intgroup = c("Cell.Type"), returnData=TRUE), file="MiSeq_poolFullDatadata(min30count)_vst.csv")

pdf("plot_PCA_CelltypeMiseqPoolFullData(min30count)_vst.pdf")
plotPCA(vsd, intgroup = c("Cell.Type"))
dev.off()

#Get normalized data

register(MulticoreParam(workers=10))
write.csv(assay(vsd), file="Normalized_vst_Counts_NOVA_poolFullData.csv")

write.csv(assay(rld), file="Normalized_rld_Counts_NOVA_poolFullData.csv")


dds <- DESeq(dds0, parallel=TRUE)

ncounts <- counts(dds, normalized=TRUE)
write.csv(ncounts,"Normalized_Counts_NOVA_poolFullData-no5outliers.csv")

##Do variance plots for merged data as well

library(vsn)
ntd <- normTransform(dds)


untransplot <- meanSdPlot(ncounts)

ntdplot <- meanSdPlot(assay(ntd)) 
ntdplot2 <- ntdplot$gg + scale_y_continuous(limits=c(0,8))

vstplot <- meanSdPlot(assay(vsd))
vstplot2 <- vstplot$gg + scale_y_continuous(limits=c(0,8))

# rldplot <- meanSdPlot(assay(rld)) 
# rldplot2 <- rldplot$gg + scale_y_continuous(limits=c(0,8))


pdf("FullData_untransformed_variance.pdf")
untransplot$gg
dev.off()

pdf("FullData_log2_transformed_variance.pdf")
ntdplot2
dev.off()

pdf("FullData_vst_variance.pdf")
vstplot2
dev.off()

pdf("FullData_vst_noscale_variance.pdf")
vstplot$gg
dev.off()




#Testing batch effect
#1. Overlap of genes "identified" in each cell-type across pools (>30 ncounts in each CT)
#2. Overlap of top 500/1000 variable genes
#3. Correlations of mean gene expression (by CT)

#1. Overlap of genes "identified" in each cell-type across pools (>30 ncounts in each CT)
#2. Overlap of top 500/1000 variable genes
library(DESeq2)
library(Rsamtools);
library("GenomicFeatures");
library("GenomicAlignments");
library("BiocParallel");
library("pheatmap");
library("RColorBrewer");

options(stringsAsFactors = FALSE)
options(scipen = 999)
load("se_nova1.rData")
load("se_nova2.rData")
load("se_nova3.rData")
load("se_nova4.rData")
load("se_nova5.rData")

outliers <- read.csv("MiSeq Outlier List current.csv")
pTable <- read.csv("pTable.csv")

se_nova1 <- se_nova1[,!(colnames(se_nova1) %in% outliers$Subject[outliers$Pool==1])]
se_nova2 <- se_nova2[,!(colnames(se_nova2) %in% outliers$Subject[outliers$Pool==2])]

##For Multi-Exec function in MobaxTerm
dds0 <- DESeqDataSet(se_nova1, design = ~ group)
dds0 <- DESeqDataSet(se_nova2, design = ~ group)
dds0 <- DESeqDataSet(se_nova3, design = ~ group)
dds0 <- DESeqDataSet(se_nova4, design = ~ group)
dds0 <- DESeqDataSet(se_nova5, design = ~ group)

keep <- rowSums(counts(dds0)) > 30
sum(keep)
dds0 <- dds0[keep,]

#Initially did >30 criteria within each CT, further tested >100 for bias by low-expressing genes
dds0ident <- dds0
dds0ident <- estimateSizeFactors(dds0ident)

ncounts <- counts(dds0ident, normalized=TRUE)

keepPV <- rowSums(ncounts[,grep("PV", colnames(ncounts))]) > 100
keepSST <- rowSums(ncounts[,grep("SST", colnames(ncounts))]) > 100
keepVIP <- rowSums(ncounts[,grep("VIP", colnames(ncounts))]) > 100
keepPyrU <- rowSums(ncounts[,grep("L2n3", colnames(ncounts))]) > 100
keepPyrL <- rowSums(ncounts[,grep("L5n6", colnames(ncounts))]) > 100

keep <- list(keepPV, keepSST, keepVIP, keepPyrU, keepPyrL)
keepV <- Reduce("|", keep)

ncounts[!keepPV,grep("PV", colnames(ncounts))] <- NA
ncounts[!keepSST,grep("SST", colnames(ncounts))] <- NA
ncounts[!keepVIP,grep("VIP", colnames(ncounts))] <- NA
ncounts[!keepPyrU,grep("L2n3", colnames(ncounts))] <- NA
ncounts[!keepPyrL,grep("L5n6", colnames(ncounts))] <- NA

ncounts_rm<- ncounts[keepV,]



#Seperately for each pool now
library(Vennerable)
library(gplots)
library(venn)

#Pool 1
PVgenes1 <- (as.character(rownames(ncounts)[keepPV]))
PYRUgenes1 <- (as.character(rownames(ncounts)[keepPyrU]))
PYRLgenes1 <- (as.character(rownames(ncounts)[keepPyrL]))
SSTgenes1 <- (as.character(rownames(ncounts)[keepSST]))
VIPgenes1 <- (as.character(rownames(ncounts)[keepVIP]))

identList1 <- list(PVgenes1, PYRUgenes1, PYRLgenes1, SSTgenes1, VIPgenes1)
# names(identList1) <- c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP")
# vtable1 <- venn::venn(list("PV"=PVgenes1, "PYR2/3"=PYRUgenes1, "PYR5/6"=PYRLgenes1, "SST"=SSTgenes1, "VIP"=VIPgenes1))
# vtable1 <- venn(identList1, ggplot=TRUE)

pdf("Pool1 identified genes Venn - 100.pdf")
par(cex=1.5)
venn(identList1, snames=c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP"), zcolor=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE)
dev.off()

save(PVgenes1, PYRUgenes1, PYRLgenes1, SSTgenes1, VIPgenes1, file="pool_1_identified_genes - 100.rData")

#top variable genes
p1rowvars <- as.data.frame(ncounts)
p1rowvars$rowvar <- apply(ncounts, 1, var, na.rm=TRUE)
p1top500 <- rownames(p1rowvars[order(p1rowvars$rowvar, decreasing = TRUE),])[1:500]
p1top1000 <- rownames(p1rowvars[order(p1rowvars$rowvar, decreasing = TRUE),])[1:1000]

save(p1top500, p1top1000, file="p1_top_variable_genes.rData")

#Pool 2
PVgenes2 <- (as.character(rownames(ncounts)[keepPV]))
PYRUgenes2 <- (as.character(rownames(ncounts)[keepPyrU]))
PYRLgenes2 <- (as.character(rownames(ncounts)[keepPyrL]))
SSTgenes2 <- (as.character(rownames(ncounts)[keepSST]))
VIPgenes2 <- (as.character(rownames(ncounts)[keepVIP]))

identList1 <- list(PVgenes2, PYRUgenes2, PYRLgenes2, SSTgenes2, VIPgenes2)
# names(identList1) <- c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP")
# vtable1 <- venn::venn(list("PV"=PVgenes2, "PYR2/3"=PYRUgenes2, "PYR5/6"=PYRLgenes2, "SST"=SSTgenes2, "VIP"=VIPgenes2))
# vtable1 <- venn(identList1, ggplot=TRUE)

pdf("Pool2 identified genes Venn - 100.pdf")
par(cex=1.5)
venn(identList1, snames=c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP"), zcolor=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE)
dev.off()
save(PVgenes2, PYRUgenes2, PYRLgenes2, SSTgenes2, VIPgenes2, file="pool_2_identified_genes - 100.rData")

#top variable genes
p2rowvars <- as.data.frame(ncounts)
p2rowvars$rowvar <- apply(ncounts, 1, var, na.rm=TRUE)
p2top500 <- rownames(p2rowvars[order(p2rowvars$rowvar, decreasing = TRUE),])[1:500]
p2top2000 <- rownames(p2rowvars[order(p2rowvars$rowvar, decreasing = TRUE),])[1:1000]

save(p2top500, p2top2000, file="p2_top_variable_genes.rData")



#Pool 3
PVgenes3 <- (as.character(rownames(ncounts)[keepPV]))
PYRUgenes3 <- (as.character(rownames(ncounts)[keepPyrU]))
PYRLgenes3 <- (as.character(rownames(ncounts)[keepPyrL]))
SSTgenes3 <- (as.character(rownames(ncounts)[keepSST]))
VIPgenes3 <- (as.character(rownames(ncounts)[keepVIP]))

identList1 <- list(PVgenes3, PYRUgenes3, PYRLgenes3, SSTgenes3, VIPgenes3)
# names(identList1) <- c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP")
# vtable1 <- venn::venn(list("PV"=PVgenes3, "PYR2/3"=PYRUgenes3, "PYR5/6"=PYRLgenes3, "SST"=SSTgenes3, "VIP"=VIPgenes3))
# vtable1 <- venn(identList1, ggplot=TRUE)

pdf("Pool3 identified genes Venn - 100.pdf")
par(cex=1.5)
venn(identList1, snames=c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP"), zcolor=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE)
dev.off()

save(PVgenes3, PYRUgenes3, PYRLgenes3, SSTgenes3, VIPgenes3, file="pool_3_identified_genes - 100.rData")

#top variable genes
p3rowvars <- as.data.frame(ncounts)
p3rowvars$rowvar <- apply(ncounts, 1, var, na.rm=TRUE)
p3top500 <- rownames(p3rowvars[order(p3rowvars$rowvar, decreasing = TRUE),])[1:500]
p3top3000 <- rownames(p3rowvars[order(p3rowvars$rowvar, decreasing = TRUE),])[1:1000]

save(p3top500, p3top3000, file="p3_top_variable_genes.rData")




#Pool 4
PVgenes4 <- (as.character(rownames(ncounts)[keepPV]))
PYRUgenes4 <- (as.character(rownames(ncounts)[keepPyrU]))
PYRLgenes4 <- (as.character(rownames(ncounts)[keepPyrL]))
SSTgenes4 <- (as.character(rownames(ncounts)[keepSST]))
VIPgenes4 <- (as.character(rownames(ncounts)[keepVIP]))

identList1 <- list(PVgenes4, PYRUgenes4, PYRLgenes4, SSTgenes4, VIPgenes4)
# names(identList1) <- c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP")
# vtable1 <- venn::venn(list("PV"=PVgenes4, "PYR2/3"=PYRUgenes4, "PYR5/6"=PYRLgenes4, "SST"=SSTgenes4, "VIP"=VIPgenes4))
# vtable1 <- venn(identList1, ggplot=TRUE)

pdf("Pool4 identified genes Venn - 100.pdf")
par(cex=1.5)
venn(identList1, snames=c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP"), zcolor=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE)
dev.off()

save(PVgenes4, PYRUgenes4, PYRLgenes4, SSTgenes4, VIPgenes4, file="pool_4_identified_genes - 100.rData")

#top variable genes
p4rowvars <- as.data.frame(ncounts)
p4rowvars$rowvar <- apply(ncounts, 1, var, na.rm=TRUE)
p4top500 <- rownames(p4rowvars[order(p4rowvars$rowvar, decreasing = TRUE),])[1:500]
p4top4000 <- rownames(p4rowvars[order(p4rowvars$rowvar, decreasing = TRUE),])[1:1000]

save(p4top500, p4top4000, file="p4_top_variable_genes.rData")




#Pool 5
PVgenes5 <- (as.character(rownames(ncounts)[keepPV]))
PYRUgenes5 <- (as.character(rownames(ncounts)[keepPyrU]))
PYRLgenes5 <- (as.character(rownames(ncounts)[keepPyrL]))
SSTgenes5 <- (as.character(rownames(ncounts)[keepSST]))
VIPgenes5 <- (as.character(rownames(ncounts)[keepVIP]))

identList1 <- list(PVgenes5, PYRUgenes5, PYRLgenes5, SSTgenes5, VIPgenes5)
# names(identList1) <- c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP")
# vtable1 <- venn::venn(list("PV"=PVgenes5, "PYR2/3"=PYRUgenes5, "PYR5/6"=PYRLgenes5, "SST"=SSTgenes5, "VIP"=VIPgenes5))
# vtable1 <- venn(identList1, ggplot=TRUE)

pdf("Pool5 identified genes Venn - 100.pdf")
par(cex=1.5)
venn(identList1, snames=c("PV", "PYR-2/3", "PYR-5/6", "SST", "VIP"), zcolor=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE)
dev.off()

save(PVgenes5, PYRUgenes5, PYRLgenes5, SSTgenes5, VIPgenes5, file="pool_5_identified_genes - 100.rData")

#top variable genes
p5rowvars <- as.data.frame(ncounts)
p5rowvars$rowvar <- apply(ncounts, 1, var, na.rm=TRUE)
p5top500 <- rownames(p5rowvars[order(p5rowvars$rowvar, decreasing = TRUE),])[1:500]
p5top5000 <- rownames(p5rowvars[order(p5rowvars$rowvar, decreasing = TRUE),])[1:1000]

save(p5top500, p5top5000, file="p5_top_variable_genes.rData")


#Comparing across runs, within each cell-type
library(Vennerable)
library(gplots)
library(venn)
options(stringsAsFactors = FALSE)
options(scipen = 999)

load("p1_top_variable_genes.rData")
load("p2_top_variable_genes.rData")
load("p3_top_variable_genes.rData")
load("p4_top_variable_genes.rData")
load("p5_top_variable_genes.rData")
load("pool_1_identified_genes.rData")
load("pool_2_identified_genes.rData")
load("pool_3_identified_genes.rData")
load("pool_4_identified_genes.rData")
load("pool_5_identified_genes.rData")

PVList <- list(PVgenes1, PVgenes2, PVgenes3, PVgenes4, PVgenes5)

pdf("PV identified genes Venn (across pools).pdf")
par(cex=1.5)
venn(PVList, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


PYRUList <- list(PYRUgenes1, PYRUgenes2, PYRUgenes3, PYRUgenes4, PYRUgenes5)

pdf("PYR2-3 identified genes Venn (across pools).pdf")
par(cex=1.5)
venn(PYRUList, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


PYRLList <- list(PYRLgenes1, PYRLgenes2, PYRLgenes3, PYRLgenes4, PYRLgenes5)

pdf("PYR5-6 identified genes Venn (across pools).pdf")
par(cex=1.5)
venn(PYRLList, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


SSTList <- list(SSTgenes1, SSTgenes2, SSTgenes3, SSTgenes4, SSTgenes5)

pdf("SST identified genes Venn (across pools).pdf")
par(cex=1.5)
venn(SSTList, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


VIPList <- list(VIPgenes1, VIPgenes2, VIPgenes3, VIPgenes4, VIPgenes5)

pdf("VIP identified genes Venn (across pools).pdf")
par(cex=1.5)
venn(VIPList, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()

top500list <- list(p1top500, p2top500, p3top500, p4top500, p5top500)
top1000list <- list(p1top1000, p2top2000, p3top3000, p4top4000, p5top5000)

pdf("Top 500 most variable genes Venn (across pools).pdf")
par(cex=1.5)
venn(top500list, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


pdf("Top 1000 most variable genes Venn (across pools).pdf")
par(cex=1.5)
venn(top1000list, snames=c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5"), par=TRUE)
dev.off()


#3. Correlations of mean gene expression (by CT)
library(DESeq2)
library(Rsamtools);
library("GenomicFeatures");
library("GenomicAlignments");
library("BiocParallel");
library("pheatmap");
library("RColorBrewer");

load("seFullData_5p-3bp_3p-15bp.rData")

coldata <- colData(seFullData)
coldata$Pool <- c(rep(1, 75), rep(2, 75), rep(3, 81), rep(4, 80), rep(5, 68))
colData(seFullData) <- coldata


coldataTR <- coldata[,c("X", "Pool")]
ncounts <- read.csv("Normalized_vst_Counts_NOVA_poolFullData.csv")

ncountst <- t(ncounts)
colnames(ncountst) <- ncountst[1,]
ncountst <- ncountst[-1,]
ncountst <- as.data.frame(ncountst)
ncountst <- sapply(ncountst, as.numeric)

ncountsTT <- t(ncountst)
ncountsTT <- as.data.frame(ncountsTT)
names(ncountsTT) <- names(ncounts)[-1]

p1 <- ncountsTT[,c(1:75)]
p2 <- ncountsTT[,c(76:150)]
p3 <- ncountsTT[,c(151:231)]
p4 <- ncountsTT[,c(232:311)]
p5 <- ncountsTT[,c(312:379)]

p1$mean <- apply(p1, 1, mean)
p2$mean <- apply(p2, 1, mean)
p3$mean <- apply(p3, 1, mean)
p4$mean <- apply(p4, 1, mean)
p5$mean <- apply(p5, 1, mean)

CrossPoolncounts <- cbind(p1$mean, p2$mean, p3$mean, p4$mean, p5$mean)
CrossPoolncounts <- as.data.frame(CrossPoolncounts)
names(CrossPoolncounts) <- c("Pool1", "Pool2", "Pool3", "Pool4", "Pool5")

library(ggplot2)
library(reshape2)
library(dplyr)
library(GGally)

CrossPoolMelt <- melt(CrossPoolncounts)

scattermatrix <- ggpairs(CrossPoolncounts)

rownames(CrossPoolncounts) <- colnames(ncountst)

#Save counts, and which tetrads in which pools
write.csv(CrossPoolncounts, file="Cross-pool count comparisons.csv")

pdf("Cross-pool correlation matrix.pdf")
scattermatrix
dev.off()





