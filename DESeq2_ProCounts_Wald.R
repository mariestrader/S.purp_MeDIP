# This script uses DESeq2 to analyze methylation counts of promoters
# This script also outputs files for GOMWU analysis and DAPC
# Strader et al 2019, Journal of Experimental Marine Biology and Ecology
# Parental environments alter DNA methylation in offspring of the purple sea urchin, Strongylocentrotus purpuratus

# note on naming scheme. Samples were initially labeled using (B) Blue: non-upwelling conditions (650 uatm, 17C); (R) Red: upwelling conditions (1330 uatm, 13C), with the first letter denoting the maternal environment and the second letter denoting the developmental environment (B: 450 uatm, R: 1050 uatm). The naming scheme was updated for the manuscript for clarity. BB=NL, BR=NH, RB=UL, RR=UH

#---------------- load libraries ------------------
library('DESeq2') #DESeq2_1.18.1
library('arrayQualityMetrics') #arrayQualityMetrics_3.34.0
library('vegan') #vegan_2.5-3
library('rgl') #rgl_0.99.16 
library('ape') #ape_5.2
library('pheatmap') #pheatmap_1.0.10
library('dplyr') #dplyr_0.7.7
library('adegenet') #adegenet_2.1.1
library('Rmisc') #Rmisc_1.5
library('rlang') #rlang_0.3.0.1
library('ggplot2') #ggplot2_3.1.0

#---------------- set directory ------------------
#setwd()

#---------------- read in counts data ------------------
pcountsALL=read.table("MEDIP_proCounts_10July2018.txt", header=T) #version mapped only to CDS regions

pcounts=as.data.frame(pcountsALL[-c(2:6) ]) #removing top rows with general count numbers
row.names(pcounts)<-pcounts$Geneid
pcounts$Geneid=NULL

length(pcounts[,1]) #29923

colnames(pcounts)=c("BB1_BL","BB1_G","BB1_PR","BB2_BL","BB2_G","BB2_PR","BB3_BL","BB3_G","BB3_PR","BR1_BL","BR1_G","BR1_PR","BR2_BL","BR2_PR","BR3_BL","BR3_G","BR3_PR","RB1_BL","RB1_G","RB1_PR","RB2_BL","RB2_G","RB2_PR","RB3_BL","RB3_G","RB3_PR","RR1_BL","RR1_G","RR1_PR","RR2_BL","RR2_G","RR2_PR","RR3_BL","RR3_G","RR3_PR")

summary(colSums(pcounts))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 126003  379926  460429  473915  519338  898705

totalCounts=colSums(pcounts)
totalCounts
#BB1_BL  BB1_G BB1_PR BB2_BL  BB2_G BB2_PR BB3_BL  BB3_G BB3_PR BR1_BL  BR1_G BR1_PR BR2_BL BR2_PR BR3_BL 
#432234 877839 596985 521128 437977 487299 476244 440651 692467 602793 620167 283123 381500 368240 378351 
# BR3_G BR3_PR RB1_BL  RB1_G RB1_PR RB2_BL  RB2_G RB2_PR RB3_BL  RB3_G RB3_PR RR1_BL  RR1_G RR1_PR RR2_BL 
#543294 359608 898705 545381 482068 460429 126003 466068 437693 497304 362945 461329 363598 517547 515494 
# RR2_G RR2_PR RR3_BL  RR3_G RR3_PR 
#297790 364724 440137 407782 442143 

min(totalCounts) #126003
mean(totalCounts) #473915.4
max(totalCounts)  #898705

#---------------- Get promoter count means and remove low promoters ------------------
mns = apply(pcounts, 1, mean)

pcounts=pcounts[mns>3,] #get rid of promoters that show little or no methylation

table(mns > 3)
#FALSE  TRUE 
# 15370 14553

totalCountsFilt=colSums(pcounts)
totalCountsFilt
#BB1_BL  BB1_G BB1_PR BB2_BL  BB2_G BB2_PR BB3_BL  BB3_G BB3_PR BR1_BL  BR1_G BR1_PR BR2_BL BR2_PR BR3_BL 
#423108 842261 569480 507685 431376 476580 464580 417819 663109 567419 599532 276070 367798 351993 366986 
# BR3_G BR3_PR RB1_BL  RB1_G RB1_PR RB2_BL  RB2_G RB2_PR RB3_BL  RB3_G RB3_PR RR1_BL  RR1_G RR1_PR RR2_BL 
#526279 349806 859717 527721 470280 443084 122030 455322 423668 478303 349115 446730 351146 489546 500716 
# RR2_G RR2_PR RR3_BL  RR3_G RR3_PR 
#290622 353000 424130 394200 424065  

min(totalCountsFilt) #122030
max(totalCountsFilt) #859717
mean(totalCountsFilt) #457293.6

#---------------- set up sample metadata ------------------
colData <- data.frame(row.names= colnames(pcounts), bucket= c("BB1","BB1","BB1","BB2","BB2","BB2","BB3","BB3","BB3","BR1","BR1","BR1","BR2","BR2","BR3","BR3","BR3","RB1","RB1","RB1","RB2","RB2","RB2","RB3","RB3","RB3","RR1","RR1","RR1","RR2","RR2","RR2","RR3","RR3","RR3"), stage=c("BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR"))

colData$treat = substr(colData$bucket, 1,2)
colData$treat_maternal = substr(colData$bucket, 1,1)
colData$treat_dev = substr(colData$bucket, 2,2)

#---------------- identify sample outliers ------------------
#dds<-DESeqDataSetFromMatrix(countData=pcounts, colData=colData, design= ~ treat+stage)
#vsd=varianceStabilizingTransformation(dds, blind=T)
#e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
#arrayQualityMetrics(e, intgroup=c("stage"), force=T, outdir= "report_for_PRO_vsd_MEDIP_07102018")

#Was determined there are 4 ouliers: BB1_BL, BB2_G, RB2_G, RR2_G

#---------------- remove outliers from both gcounts and colData ------------------
#remove outliers from both gcounts and colData
pcounts$BB1_BL<-NULL
pcounts$BB2_G<-NULL
pcounts$RB2_G<-NULL
pcounts$RR2_G<-NULL

colData=as.data.frame(colData[-c(1,5,22,31), ])

#---------------- Wald Tests ------------------

dds<-DESeqDataSetFromMatrix(pcounts,

	colData = colData, 

	design = formula(~ treat_maternal + treat_dev + stage ))


meth.rld=rlog(dds)

meth.coldata = colData

#### Wald test for maternal treatment
ddsM_RB<-DESeq(dds, minReplicatesForReplace=Inf) 
resM_RB<-results(ddsM_RB, contrast=c('treat_maternal', 'R', 'B')) 
mcols(resM_RB,use.names=TRUE)
summary(resM_RB)
#out of 14553 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 185, 1.3% 
#LFC < 0 (down)   : 200, 1.4% 
#outliers [1]     : 5, 0.034% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

#### Wald test for developmental treatment
ddsD_RB<-DESeq(dds, minReplicatesForReplace=Inf) 
resD_RB<-results(ddsD_RB, contrast=c('treat_dev', 'R', 'B')) 
mcols(resD_RB,use.names=TRUE)
summary(resD_RB)
#out of 14553 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2, 0.014% 
#LFC < 0 (down)   : 2, 0.014% 
#outliers [1]     : 5, 0.034% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

####Wald test for stages
dds_BL_G<-DESeq(dds, minReplicatesForReplace=Inf) 
res_BL_G<-results(dds_BL_G, contrast=c('stage', 'BL', 'G')) 
mcols(res_BL_G,use.names=TRUE)
summary(res_BL_G)
#out of 14553 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 5, 0.034% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

dds_BL_PR<-DESeq(dds, minReplicatesForReplace=Inf)
res_BL_PR<-results(dds_BL_PR, contrast=c('stage', 'BL', 'PR')) 
mcols(res_BL_PR,use.names=TRUE)
summary(res_BL_PR)
#out of 14553 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 18, 0.12% 
#LFC < 0 (down)   : 4, 0.027% 
#outliers [1]     : 5, 0.034% 
#low counts [2]   : 1129, 7.8% 
#(mean count < 4)

dds_G_PR<-DESeq(dds, minReplicatesForReplace=Inf) 
res_G_PR<-results(dds_G_PR, contrast=c('stage', 'G', 'PR')) 
mcols(res_G_PR,use.names=TRUE)
summary(res_G_PR)
#out of 14553 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 5, 0.034% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

save(dds,resM_RB,resD_RB,res_BL_G,res_BL_PR,res_G_PR,meth.rld, meth.coldata,pcounts, file="MeDIP_DESeq2_PRO_Wald.Rdata")

#---------------- Make output files for rld and p values; get #DMGs for 0.05 cutoff ------------------
ll=load("MeDIP_DESeq2_PRO_Wald.Rdata")

meth.rld.df = assay(meth.rld)
colnames(meth.rld.df) = colnames(pcounts)

vals=(cbind(resM_RB$stat, resM_RB$pvalue, resM_RB$padj, resD_RB$stat, resD_RB$pvalue, resD_RB$padj, res_G_PR$stat, res_G_PR$pvalue, res_G_PR$padj, res_BL_G$stat, res_BL_G$pvalue, res_BL_G$padj, res_BL_PR$stat, res_BL_PR$pvalue, res_BL_PR$padj)) # collect pvalues for each comparison
colnames(vals)=c("stat_resM_RB","pval_resM_RB", "padj_resM_RB","stat_resD_RB","pval_resD_RB", "padj_resD_RB","stat_res_G_PR","pval_res_G_PR", "padj_res_G_PR","stat_res_BL_G","pval_res_BL_G", "padj_res_BL_G","stat_res_BL_PR","pval_res_BL_PR", "padj_res_BL_PR") # name columns
rldpvals=as.data.frame(cbind(meth.rld.df,vals)) # combine RLD with pvals

#get number of DMGs with padj < 0.01
#maternal condition
length((rldpvals[rldpvals$padj_resM_RB<= 0.01 & !is.na(rldpvals$padj_resM_RB), ])[,1]) #156

#developmental condition
length((rldpvals[rldpvals$padj_resD_RB<= 0.01 & !is.na(rldpvals$padj_resD_RB), ])[,1]) #1

#Gastrula v. Prism 
length((rldpvals[rldpvals$padj_res_G_PR<= 0.01 & !is.na(rldpvals$padj_res_G_PR), ])[,1]) #0

#Blastula v. Prism 
length((rldpvals[rldpvals$padj_res_BL_PR<= 0.01 & !is.na(rldpvals$padj_res_BL_PR), ])[,1]) #4

#Blastula v. Gastrula 
length((rldpvals[rldpvals$padj_res_BL_G<= 0.01 & !is.na(rldpvals$padj_res_BL_G), ])[,1]) #0

write.csv(rldpvals,file="MeDIP_Spurp_PRO_rldpvals.csv",quote=F, row.names=T)

#---------------- principal coordinate analysis ------------------
ll=load("MeDIP_DESeq2_PRO_Wald.Rdata")

meth.rld.df = assay(meth.rld)
colnames(meth.rld.df) = colnames(pcounts)

dd.veg=vegdist(t(meth.rld.df), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

conditions=meth.coldata
conditions$treat_maternal=as.factor(conditions$treat_maternal)
conditions$treat_dev=as.factor(conditions$treat_dev)

dd.pcoa
#$values
#   Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1    50.688293   0.14745346  0.133166238 0.1474535      0.1331662
#2    29.573446   0.08602986  0.099832904 0.2334833      0.2329991
#3    18.300417   0.05323635  0.083166238 0.2867197      0.3161654
#4    17.210215   0.05006493  0.072055127 0.3367846      0.3882205


#---------------- plotting principal coordinate analysis ------------------

### definte maternal conditions for plotting
colorsMat <- c("coral", "turquoise4")
conditions$treat_maternal_c <- colorsMat[as.factor(conditions$treat_maternal)]
conditions$treat_maternal=as.factor(conditions$treat_maternal)
conditions$treat_maternal = factor(conditions$treat_maternal, levels = c('B', 'R'), labels = c('Non-Upwelling', 'Upwelling'))

plot(scores[,1], scores[,2], col=conditions$treat_maternal_c,pch=19, xlab="PCo1", ylab="PCo2", main="PCoA by Maternal Condition")
ordispider(scores,conditions$treat_maternal,label=F)
ordiellipse(scores,conditions$treat_maternal,label=F )
legend("topleft", legend = levels(conditions$treat_maternal), col=colorsMat, pch = 19)


### definte developmental conditions for plotting 
colorsMat <- c("coral", "turquoise4")
conditions$treat_dev_c <- colorsMat[as.factor(conditions$treat_dev)]
conditions$treat_dev=as.factor(conditions$treat_dev)
conditions$treat_dev = factor(conditions$treat_dev, levels = c('B', 'R'), labels = c('Non-Upwelling', 'Upwelling'))

plot(scores[,1], scores[,2], col=conditions$treat_dev_c,pch=19, xlab="PCo1", ylab="PCo2", main="PCoA by Larval Condition")
ordispider(scores,conditions$treat_dev,label=F)
ordiellipse(scores,conditions$treat_dev,label=F )
legend("topleft", legend = levels(conditions$treat_dev), col=colorsMat, pch = 19)


### definte developmental stage for plotting 
colorsStage <- c("coral", "turquoise4", "gold")
conditions$stage_c <- colorsStage[as.factor(conditions$stage)]
conditions$stage=as.factor(conditions$stage)
conditions$stage = factor(conditions$stage, levels = c('BL', 'G', "PR"), labels = c('Blastula', 'Gastrula', 'Prism'))

plot(scores[,1], scores[,2], col=conditions$stage_c,pch=19, xlab="PCo1", ylab="PCo2", main="PCoA by Developmental Stage")
ordispider(scores,conditions$stage,label=F)
ordiellipse(scores,conditions$stage,label=F )
legend("topleft", legend = levels(conditions$stage), col=colorsStage, pch = 19)

#---------------- multivariate statistics ------------------
ll=load("MeDIP_DESeq2_PRO_Wald.Rdata")
rld=assay(meth.rld)
ad=adonis(t(rld)~treat_maternal+treat_dev+stage,data=meth.coldata,method="manhattan") 
ad
#adonis(formula = t(rld) ~ treat_maternal + treat_dev + stage,      data = meth.coldata, method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#treat_maternal  1  19875603 19875603  1.7932 0.05782  0.004 **
#treat_dev       1  11563178 11563178  1.0432 0.03364  0.320   
#stage           2  24137313 12068656  1.0888 0.07022  0.204   
#Residuals      26 288181820 11083916         0.83833          
#Total          30 343757913                  1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

labs=c("Maternal_condition","Developmental_condition","Stage","Residuals")
cols=c("skyblue","green2","coral","grey80")

labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="GBM")

#---------------- volcano plots ------------------
ll=load("MeDIP_DESeq2_PRO_Wald.Rdata")

with(resM_RB, plot(log2FoldChange, -log10(pvalue), pch=20, main="Maternal Condition Padj<0.01", xlim=c(-5,5), ylim=c(1,45)))
with(subset(resM_RB, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))

with(resD_RB, plot(log2FoldChange, -log10(pvalue), pch=20, main="Developmental Condition Padj<0.01", xlim=c(-5,5),ylim=c(1,45)))
with(subset(resD_RB, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))

with(res_BL_PR, plot(log2FoldChange, -log10(pvalue), pch=20, main="BL v. PR Padj<0.01", xlim=c(-5,5),ylim=c(1,45)))
with(subset(res_BL_PR, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))

