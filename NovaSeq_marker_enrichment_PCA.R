library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(Rcmdr)
library(reshape2)
library(funtimes)
library(ClusterR)
library(mclust)
library(plotly)
library(car)
library(stringr)
options(stringsAsFactors = FALSE)
SEfun <- function(x){
  x <- sd(x)/(sqrt(length(x)))
}

#Load and clean data of duplicates/unwanted genes
ncounts <- read.csv("Normalized_Counts_NOVA_poolFullData-no5outliers.csv")
#metadata from BioMart
mData <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/MiSeq/Flow cell 1/martquery_1104155809_616.txt")
pTable <- read.csv("pTable.csv")
#Make sample ID name consistent
pTable$X <- gsub("-", "\\.", pTable$X)
pTable$X <- gsub("_", "\\.", pTable$X)
pTable$Cell.Type <- gsub("-", "\\.", pTable$Cell.Type)
pTable <- pTable[,c(3:5)]

mData <- mData[!duplicated(mData$Gene.name),]
mData <- as.data.frame(mData[mData$Gene.stable.ID %in% ncounts$X,])

ncounts <- ncounts[!duplicated(ncounts$X),]
ncounts <- ncounts[which(mData$Gene.stable.ID %in% ncounts$X),]

ncounts2 <- merge(mData, ncounts, by.x="Gene.stable.ID", by.y="X")

rownames(ncounts2) <- ncounts2$Gene.name
ncounts2 <- ncounts2[,c(9:ncol(ncounts2))]

#User-defined cell-type markers of interest - these were derived from the Allen Brain human cell-types database
markerlist <- c("Slc17a7", "Sst", "Pvalb", "Vip","Gad1", "Gad2","Stmn2", "Lhx6", "Cux2", "Themis")
markerlist <- toupper(markerlist)

markers <- ncounts2[rownames(ncounts2) %in% markerlist,]

markers <- as.data.frame(t(markers))
markers$Gene.Name <- rownames(markers)
markers <- merge(pTable, markers, by.x="X", by.y="Gene.Name", all.y=TRUE)

#Make key variables factors for plotting
markers$Cell.Type <- as.factor(markers$Cell.Type)
markers$Subject.Group <- as.factor(markers$Subject.Group)

#Make individual plots since the count scales will be different across genes, not ideal for facet plots
for(i in 1:length(markerlist)){
  workingGene <- toupper(markerlist[i])
  assign(paste0(markerlist[i], "_plot"), 
         ggplot(markers, aes_string(x="Cell.Type", y=workingGene, fill="Cell.Type")) + 
           scale_fill_manual(values=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"))  + 
           geom_boxplot(outlier.alpha = 0,alpha=0.6, width=0.6) + 
           guides(fill=guide_legend(title=NULL)) + 
           ggtitle(str_to_title(markerlist[i])) + 
           theme_classic() +
           theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "none") + 
           expand_limits(y=0))
}

#Arrange in grid nicely
grid.arrange(SLC17A7_plot, SST_plot, PVALB_plot, VIP_plot, STMN2_plot, GAD1_plot, GAD2_plot, LHX6_plot, CUX2_plot, THEMIS_plot, nrow=2)

#PCA plot

PCA_data <- read.csv("Novaseq5_PCAdata(min10count)_vst.csv")

#Make names nicer for legend
PCA_data$Cell.Type <- gsub("Pyr-L2n3", "L2/3 PYR", PCA_data$Cell.Type)
PCA_data$Cell.Type <- gsub("Pyr-L5n6", "L5/6 PYR", PCA_data$Cell.Type)
unique(PCA_data$Cell.Type)

#refactor
PCA_data$Cell.Type <- factor(PCA_data$Cell.Type, levels=c("PVALB", "L2/3 PYR", "L5/6 PYR", "SST", "VIP"))
PCA_data$Pool
PCA_data$Cell.Type

#plot with nicer output than default DEseq2 PCA plots
ggplot(PCA_data, aes(x=PC1, y=PC2, colour=Cell.Type)) +
  geom_point(size=4, alpha=.5) +
  scale_colour_manual(values=c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF")) +
  theme_classic() +
  xlab("PC1 (45% variance)") +
  ylab("PC2 (19% variance)")


#quick k-means to check clusters are actually seperated 
PCA_KM <- PCA_data[,2:3]
rownames(PCA_KM) <- PCA_data$X
PCA_kmeans <- kmeans(PCA_KM, centers=4, iter.max=10000, nstart=10000)
kmeans_res <- as.data.frame(PCA_kmeans$cluster, row.names = names(PCA_kmeans$cluster)) # 100% accuracy 