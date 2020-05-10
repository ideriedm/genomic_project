library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringi)
library(BBmisc)

######################### Reproduction of Fig. 5B of Ciriello et al., 2015 #######################################


# RPPA data 20,440 proteins x 819 samples
rppa <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/data_rppa.txt")
# Creation of a correspondance table of 673 samples :
# col V1 : 01/06
# col V2 : TCGA-EW-A2FR
# col sample : TCGA-EW-A2FR-01
sample_ids <- as.data.table(lapply(stri_split_regex(stri_reverse(names(rppa[ , grepl( "TCGA" , names( rppa ) ), with = FALSE ])), 
                                                    pattern = '[-\\s]+', n = 2), stri_reverse))
sample_ids <- sample_ids[, data.table(t(.SD), keep.rownames=TRUE)]
sample_ids$sample <- paste0(sample_ids$V2,"-",sample_ids$V1) 

# Column 1 (= Composite.Element.REF), contains the name of the proteins 
# to upper case => all in capital letters
rppa$Composite.Element.REF <- toupper(rppa$Composite.Element.REF) 

# Two lines of column 1 have 3 values of IDs separated by "|", other lines have only 2. We will clean these 2 lines:
rppa[rppa$Composite.Element.REF == toupper("CHEK1|Chk1_pS296;CHEK1|CHK1_pS296"),"Composite.Element.REF"] <- "CHEK1|CHK1_PS296"
rppa[rppa$Composite.Element.REF == toupper("CDKN2A|P16INK4A;CDKN2A|p16_INK4a"),"Composite.Element.REF"] <- "CDKN2A|P16_INK4A"

# First column with IDs separated by "|"  are split into 2 columns with IDs
rppa <- separate(data = rppa, col = Composite.Element.REF, into = c("ID1", "ID2"), sep = "\\|")

# ILC class : from mmc9 = supplementary table 8
# Classes = reactive, immune, proliferative
ilc_class <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/mmc9_selected.txt",sep="\t")
# Copy V2 to the same col named as in ilc_class (= Sample ID)
sample_ids$`Sample ID` <- sample_ids$V2  
# Merge the sample_ids to the ilc_class according to the col Sample ID
ilc_matched <- merge(x=ilc_class, y=sample_ids, by="Sample ID") 
# Keep only the col sample & 60 Gene-classifier Class Assignment
ilc_matched <- ilc_matched[,c("sample", "60 Gene-classifier Class Assignment")] 
# transpose as a list
ilc_matched_t <- transpose(ilc_matched[,c("60 Gene-classifier Class Assignment")]) 
# col names according to the sample col, to match the RPPA names
names(ilc_matched_t) <- ilc_matched$sample 

# Take the col in rppa that are in ilc
rppa_ilc <- rppa[, (names(rppa) %in% names(ilc_matched_t))|(names(rppa)=="ID1")|(names(rppa)=="ID2"), with=FALSE]
# bind the ILC_subtypes to the rppa_ilc, as a new first row
rppa_ilc <- rbind(ilc_matched_t, rppa_ilc, fill=TRUE) 

# col names = TCGA-..-....-01, n=93
# Contains only 1 line with the 3 ILC subtypes
annotations <- rppa_ilc[1,]
annotations$ID2 <- NULL
annotations$ID1 <- NULL

# col V1 : "1:Reactive-like","2:Immune-related" or "3:Proliferative" 
# col V2 : TCGA-..-....-01, n=93
annotations2 <-data.table(t(annotations),names(annotations))
unique(annotations2$V1)
annotations2[annotations2$V1 == "Immune-related"]$V1  <- "2:Immune-related" 
annotations2[annotations2$V1 == "Reactive-like"]$V1  <- "1:Reactive-like" 
annotations2[annotations2$V1 == "Proliferative"]$V1  <- "3:Proliferative" 
# Ordered according to V1 : 
# first : "1:Reactive-like" , second : "2:Immune-related" and third : "3:Proliferative"
annotations2 <- annotations2[order(V1)]

# To reproduce Figure 5B, we need to select the 37 proteins we are interested in:
fig5bRows <- c("C-KIT","PKC-ALPHA","PKC-ALPHA_PS657","BETA-CATENIN","E-CADHERIN","MYH11","14-3-3_EPSILON","P70S6K","RAPTOR","EIF4G","HSP70","P62-LCK-LIGAND","ASNS","STAT5-ALPHA","PRAS40_PT246","MTOR_PS2448","MEK1","MEK1_PS217_S221","PKC-PAN_BETAII_PS660","EGFR_PY1173","MIG-6","SMAD4","14-3-3_EPSILON","MYH11","RAD51","CD49B","C-KIT","CYCLIN_E1","FOXM1","	PCNA","CHK1_PS345","RAD50","RAD51","XRCC1","BRCA2","SF2","SRC","P27","CD49B","FIBRONECTIN","MAPK_PT202_Y204","MEK1","MEK1_PS217_S221","PKC-ALPHA_PS657", "PKC-PAN_BETAII_PS660","SRC_PY527","	YB-1_PS102")
rppa_ilc_fig5b <- rppa_ilc[(rppa_ilc$ID1 %in% fig5bRows)|(rppa_ilc$ID2 %in% fig5bRows),names(rppa_ilc),with=FALSE]
rppa_ilc_fig5b <- rppa_ilc_fig5b[order(match(rppa_ilc_fig5b$ID2, fig5bRows))]


#################### ILC subtypes #################################################################################

rppa_ilc_fig5b$ID1 <- NULL
row.names(rppa_ilc_fig5b) <- rppa_ilc_fig5b$ID2 # put ID2 as character

rppa_ilc_fig5b <- setcolorder(rppa_ilc_fig5b, annotations2$V2)

rppa_ilc_fig5b_1 <- data.matrix(rppa_ilc_fig5b)
rppa_ilc_fig5b_1 <- rppa_ilc_fig5b_1[,-94] # Remove column with NA


#################### RPPA subtypes #################################################################################

#RPPA subtype: from TCGA_supp_info_samples_patients.txt (reactive, non-reactive, undefined)
rppa_class <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/TCGA_supp_info_samples_patients.txt",sep="\t")
unique(rppa_class$`RPPA Clusters`) # write the name of the different RPPA clusters


rppa_class$RPPA_subtype <- rppa_class$`RPPA Clusters`
# ReacI, ReacII = reactive
rppa_class[(rppa_class$`RPPA Clusters`== "ReacI")| (rppa_class$`RPPA Clusters`== "ReacII")]$RPPA_subtype <- "Reactive"
# The rest, except NA = non-reactive
rppa_class[(rppa_class$`RPPA Clusters`== "Basal") | (rppa_class$`RPPA Clusters`== "X") | (rppa_class$`RPPA Clusters`== "LumA/B") | (rppa_class$`RPPA Clusters`== "Her2")| (rppa_class$`RPPA Clusters`== "LumA")]$RPPA_subtype <- "non-Reactive"
# NA = undefined
rppa_class[is.na(rppa_class$`RPPA Clusters`)]$RPPA_subtype <- "Undefined"

unique(rppa_class$RPPA_subtype)

rppa_class$`Sample ID` <- rppa_class$`Complete TCGA ID`
# Merge the sample_ids to the ilc_class according to the col Sample ID
rppa_matched <- merge(x=rppa_class, y=sample_ids, by="Sample ID") 
# Keep only the 2 col, sample = ID with 01 
rppa_matched <- rppa_matched[,c("sample", "RPPA_subtype")] 

rppa_ilc_matched_t <- transpose(rppa_matched[,c("RPPA_subtype")]) # transpose as a list
# col names according to the sample col, to match the RPPA names
colnames(rppa_ilc_matched_t) <- rppa_matched$sample 
rppa_ilc_matched_t <- rppa_ilc_matched_t[, names(rppa_ilc_matched_t) %in% names(ilc_matched_t), with = FALSE] 

# Contains 49 samples
# col V1 = 1:reactive-like,...
# col V2 = TCGA-A2-A0EX-01
# col color = black, red, green, according to their ILC subtype
annotations3 <- annotations2[annotations2$V2 %in% names(rppa_ilc_matched_t)]
annotations3[annotations3$V1 == "1:Reactive-like","color"] <- "black"
annotations3[annotations3$V1 == "2:Immune-related","color"] <- "red"
annotations3[annotations3$V1 == "3:Proliferative","color"] <- "green"


rppa_ilc_matched_t <- setcolorder(rppa_ilc_matched_t, annotations3$V2)
rppa_ilc_matched_t2 <- data.table(t(rppa_ilc_matched_t),names(rppa_ilc_matched_t) )

rppa_ilc_matched_t2[rppa_ilc_matched_t2$V1 == "Reactive","color"] <- "black"
rppa_ilc_matched_t2[rppa_ilc_matched_t2$V1 == "Undefined","color"] <- "grey"
rppa_ilc_matched_t2[rppa_ilc_matched_t2$V1 == "non-Reactive","color"] <- "white"

# The heatmap protein expression levels will go from cyan to red, through white
mycol <- colorpanel(1000,"cyan","white","red")

#################### Bottom heatmap #################################################################################
#RPPA score for cell cycle, DNA damage response, RAS MAPK: from mmc2.txt 
rppa_scores <- fread("/Users/ines/Library/Mobile Documents/com~apple~CloudDocs/SV/SV-MA2/genomics & bioinfo/project/brca_tcga_pub2015/important/mmc2.txt",skip = 2)
rppa_scores <- rppa_scores[, c("Case.ID", "Cell cycle score", "DNA damage response score", "Ras/MAPK score")] # keep 4 col

rppa_scores$`Sample ID` <- rppa_scores$Case.ID
# Merge the sample_ids to the ilc_class according to the col Sample ID
rppa_scores <- merge(x=rppa_scores, y=sample_ids, by="Sample ID") 
# Keep only the 2 col, sample = ID with 01
rppa_scores <- rppa_scores[,c("sample", "Cell cycle score", "DNA damage response score", "Ras/MAPK score")] 

rppa_scores_matched_t <- transpose(rppa_scores[,c("Cell cycle score", "DNA damage response score", "Ras/MAPK score")]) # transpose as a list
# col names according to the sample col, to match the RPPA names
colnames(rppa_scores_matched_t) <- rppa_scores$sample 

rppa_scores_matched_t <- rppa_scores_matched_t[, names(rppa_scores_matched_t) %in% names(rppa_ilc_matched_t), with = FALSE]
rppa_scores_matched_t$ID2 <- c("Cell cycle","DNA damage response", "Ras/MAPK")
rppa_scores_matched_t[is.na(rppa_scores_matched_t)] <- 0

rppa_scores_matched_t <- setcolorder(rppa_scores_matched_t, annotations3$V2)
rppa_ilc_fig5b_2 <- rppa_ilc_fig5b[,(names(rppa_ilc_fig5b) %in% colnames(rppa_ilc_matched_t))|(names(rppa_ilc_fig5b)=="ID2"), with=FALSE]



################################## plotting with Heatmap.3 ############################################

library("gplots")
library("devtools")
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

rppa_ilc_fig5b_2 <- rppa_ilc_fig5b[,(names(rppa_ilc_fig5b) %in% colnames(rppa_ilc_matched_t))|(names(rppa_ilc_fig5b)=="ID2"), with=FALSE]
row.names(rppa_ilc_fig5b_2) <- row.names(rppa_ilc_fig5b)


prob_matrix=data.matrix(rppa_ilc_fig5b_2[,1:49,with=FALSE]) 
prot_names=row.names(rppa_ilc_fig5b_2)
patient_ids=names(rppa_ilc_fig5b_2[,1:49,with=FALSE])
rownames(prob_matrix)= prot_names
colnames(prob_matrix)= patient_ids

# cell cycle & co
prob_matrix2 = data.matrix(rppa_scores_matched_t[,1:49,with=FALSE])
rownames(prob_matrix2)= rppa_scores_matched_t$ID2
prot_names2 = rppa_scores_matched_t$ID2
colnames(prob_matrix2)= patient_ids



#Create color bars
protclass_colors=c(rep("yellow",7),rep("chocolate1",6),rep("pink",6),rep("blueviolet",5),rep("darkgreen",10),rep("darkblue",3))

subtype_colors=annotations3$color 
rppa_colors= rppa_ilc_matched_t2$color 
rlab=t(protclass_colors)
clab=cbind(rppa_colors,subtype_colors)

rownames(rlab)=c("Protein Class")
colnames(clab)=c("RPPA subtypes","ILC subtypes")

mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
main_title="Heatmap 5B"
par(cex.main=1)
heatmap.3(prob_matrix, hclustfun=myclust, distfun=mydist, na.rm = TRUE,scale="row",dendrogram="none", margins=c(6,12),
          Rowv=FALSE, Colv=FALSE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE,symkey=FALSE,
          density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=prot_names, cexRow=1, col=mycol,
          ColSideColorsSize=3, RowSideColorsSize=2, KeyValueName="Protein expression", colsep = c(23,42), rowsep = c(7,13,19,24,34,37))
legend(x=-0.1,y=0.7,legend=c("1:Reactive-like","2:Immune-related","3:Proliferative","RPPA Reactive","RPPA non-Reactive","RPPA Undefined", "RL high", "RL low", "IR high", "IR low", "Pro high", "Pro low"),
       fill=c("black","red","green","black","white","grey","yellow","chocolate1","pink","blueviolet","darkgreen","darkblue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

mycol <- colorpanel(1000,"blue","white","orange")
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
main_title="Heatmap 5B 2"
par(cex.main=1)

heatmap.3(prob_matrix2, hclustfun=myclust, distfun=mydist, na.rm = TRUE,scale="row",dendrogram="none", 
          margins=c(17,8),Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE,symkey=FALSE,density.info="none", 
          trace="none", main=main_title, labCol=FALSE, labRow=prot_names2, cexRow=1, col=mycol, 
          KeyValueName="Signature Score", colsep = c(23,42))

