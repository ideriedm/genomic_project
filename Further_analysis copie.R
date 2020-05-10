# tutorial : https://f1000research.com/articles/5-1408/v3

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)   
library(tidyr)
library(stringi)

library(limma) 
library(edgeR) 

# RNAseq data:
rnaseq <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/data_RNA_Seq_v2_expression_median.txt")
# Only use rnaseq$Entrez_Gene_Id as gene ID in this work => delete Hugo_Symbol
rnaseq$Hugo_Symbol <- NULL 

rnaseq_genes <- rnaseq$Entrez_Gene_Id  

head(rnaseq)

# Gene lengths for all 20440 IDs
gene_lenght <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/dtg_gene_lengths_final.txt")
# Order by rnaseq_genes
gene_lenght <- gene_lenght[match(rnaseq_genes, gene_lenght$entrez),] 

# Get metadata about samples in data_clinical_sample.txt
metadf <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/data_clinical_sample.txt",skip=4)
colnames(metadf)

unique(metadf$ONCOTREE_CODE)

# Separate ILC from IDC samples
ILC_samples <-  metadf[metadf$ONCOTREE_CODE == "ILC"]$SAMPLE_ID
IDC_samples <-  metadf[metadf$ONCOTREE_CODE == "IDC"]$SAMPLE_ID

# select ILC and IDC samples
rnaseq <-rnaseq[,(names(rnaseq) %in% ILC_samples )|(names(rnaseq) %in% IDC_samples )| (names(rnaseq)=="Entrez_Gene_Id"),with=FALSE] 

rnaseq_mat <- as.matrix(rnaseq[,grepl( "TCGA" , names( rnaseq ) ),with=FALSE]) 

row.names(rnaseq_mat) <- rnaseq$Entrez_Gene_Id

colnames(rnaseq_mat) 

sample_types <- data.table(colnames(rnaseq_mat))

# col V1 : TCGA-A1-A0SD-01
# col 2 = type : ILC/IDC
sample_types[sample_types$V1 %in% ILC_samples ,"type"] <- "ILC"
sample_types[sample_types$V1 %in% IDC_samples ,"type"] <- "IDC"

# Counts from all samples were stored in a single file
x <- DGEList(counts=rnaseq_mat,genes=gene_lenght,samples=sample_types, group=sample_types$type )
class(x)

# 20440 genes, 616 TCGA biopsy
dim(x) 

# vector of TCGA-..-....-01
# change the names to shorter names : substring(colnames(x), 1, nchar(colnames(x))-3) 
samplenames <- colnames(x) 

# vector of ILC/IDC
group <- x$samples$group 



###################### gene annotations ##############################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Homo.sapiens")
library(Homo.sapiens) 

# contains the Entrez_gene_Id : e.g. 83440
geneid <- rownames(x) 
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
dim(genes) # 21514 3

head(genes) # 3 columns : ENTREZID, SYMBOL, TXCHROM

# keeping only the first occurrence of each gene ID 
genes <- genes[!duplicated(genes$ENTREZID),]

# Gene order is the same in both the annotation and the data object. 


x


#############################    Data pre-processing  ##################################

# cpm function in edgeR
cpm <- cpm(x) # counts per million
lcpm <- cpm(x, log=TRUE) # log2-counts per million

########################### Mean and Median, library size #######################################


L <- mean(x$samples$lib.size) * 1e-6 
M <- median(x$samples$lib.size) * 1e-6 
c(L, M) # 19.05765 18.94148 for Ciriello et. Al, 2015 

summary(lcpm)



#########################   Removing lowly expressed genes ###################################

table(rowSums(x$counts==0)==616) # 616 samples 

# filterByExpr function in the edgeR
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) # we removed lowly expressed genes

# Fig1 from Tutorial
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



#############################  Normalising gene expression distributions #########################

# calcNormFactors function in edgeR
x <- calcNormFactors(x, method = "TMM")
# With DGEList-objects, normalisation factors are stored in x$samples$norm.factors
x$samples$norm.factors

################################### fig 2 #########################################
# Subset of samples selected
# diff are amplified for better visual demonstration
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(oma=c(6,2,2,2),mfrow=c(1,2),pch=16) # fit labels in figures 

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm[,1:10], las=2, main="") 
title(main="A. Unnormalised data", ylab="Log-cpm")

x2 <- calcNormFactors(x2)
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm[,1:10], las=2, main="") 
title(main="B. Normalised data", ylab="Log-cpm")


################################# Unsupervised clustering of samples ###############################

# Fig 3
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- group
View(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
# multi-dimensional scaling (MDS) plot => DE detected before formal tests
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")


library(Glimma)
# Glimma-plots with MDS-Plot.html 
glMDSPlot(lcpm, labels=group, groups=x$samples[,c(2,5)], launch=FALSE)



############################## Differential expression analysis ###############################

# Linear models are fitted to the data, assumption : data normally distributed. 

# col1 = groupIDC, col2 = groupILC
# design matrix
design <- model.matrix(~0+group)
# col1 = IDC, col2 = ILC
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  ILCvsIDC = ILC-IDC,   # or ILCvsIDC = IDC-ILC, 
  levels = colnames(design))
contr.matrix

############################ Removing heteroscedascity from count data ##########################
# Fig 4a 
# voom mean-variance trend
v <- voom(x, design, plot=TRUE)  
v

######################### Fitting linear models for comparisons of interest ##################

# Fig 4b 
# resultant mean-variance trend
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit) 

######################### Examining the number of DE genes #############################

summary(decideTests(efit))


# treat : to calculate p-values from empirical 
# Bayes moderated t-statistics with a minimum log-FC requirement
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
# ILC - IDC : up : 234

de.common <- which(dt[,1]!=0)
ILC_vs_IDC_up <- which(dt[,1]==1)
length(de.common)

head(tfit$genes$entrez[de.common], n=30)

# for entrez - gene name correspondence
dftmp <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/data_RNA_Seq_v2_expression_median.txt")
# gene named CDH1 after HUGO denomination = 999 under entrez_gene_id
dftmp[dftmp$Hugo_Symbol=="CDH1"]$Entrez_Gene_Id 

############################# Examining individual DE genes from top to bottom ##########################

idc.vs.ilc <- topTreat(tfit, coef=1, n=Inf)
head(idc.vs.ilc)

############################# Graphical representations of DE results ###############

# fig6
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13), hl.col = c("chocolate1","green") ) 

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
        side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)


############################ Heatmap of DE genes #############################
library(gplots)
idc.vs.ilc.topgenes <- idc.vs.ilc$entrez[1:100]
i <- which(v$genes$entrez %in% idc.vs.ilc.topgenes) # indices

# hugo name of 100 DE genes
rnaseq <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/data_RNA_Seq_v2_expression_median.txt")
gene_hugo <- rnaseq[ILC_vs_IDC_up]$Hugo_Symbol
write.table(data.frame(gene_hugo), file = "hugo.csv")

mycol <- colorpanel(1000,"darkgreen","white","deeppink")

heatmap.3(lcpm[i,], scale="row", main= "ILC vs IDC DE genes", 
          labRow=v$genes$entrez[i], labCol=group, key = TRUE,
          col=mycol, trace="none", density.info="none", dendrogram="column")


######################### PCA ################################

library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)


# i from DE_genes = indices from top DE ILCvsIDC genes
# lcpm[i,] contains RNAseq data
RNA_data.pca <- prcomp(t(lcpm[i,]), center = TRUE,scale. = TRUE)

summary(RNA_data.pca)
str(RNA_data.pca)
ggbiplot(RNA_data.pca,group=x$samples$group,ellipse = TRUE, circle = TRUE,var.axes=FALSE) +
  ggtitle("PCA of TCGA dataset comparing 100 top DE genes") +
  theme_minimal()+
  theme(legend.position = "bottom")

########################## tSNE ####################################

# install.packages("Rtsne")
# install.packages("irlba")
library(Rtsne)
library(ggplot2)
library(cowplot)
library(irlba)

tsne <- Rtsne(t(lcpm[i,]), dims = 2, perplexity=3, partial_pca=TRUE, verbose=TRUE, max_iter = 500)
colors = rainbow(length(unique(group)),start = .9, end = .3)
names(colors) = unique(group)
plot(tsne$Y,col="white", main = "tSNE analysis")
text(tsne$Y, labels=group, col=colors[group])

