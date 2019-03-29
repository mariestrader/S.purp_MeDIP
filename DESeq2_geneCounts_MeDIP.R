# This script uses DESeq2 to analyze methylation counts of genes
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
gcountsALL=read.table("MEDIP_geneCounts_6July2018.txt", header=T) 
#Mapped to Strongylocentrotus_purpuratus.Spur_3.1.dna.toplevel using hisat2, sam files cleaned, sorted, dups removed with picard tools, counted with featurecounts -t gene -g gene_id counted multimapped reads and overlaps

gcounts=as.data.frame(gcountsALL[-c(2:6) ]) #removing top rows with general count numbers
row.names(gcounts)<-gcounts$Geneid
gcounts$Geneid=NULL

length(gcounts[,1]) #30240

colnames(gcounts)=c("BB1_BL","BB1_G","BB1_PR","BB2_BL","BB2_G","BB2_PR","BB3_BL","BB3_G","BB3_PR","BR1_BL","BR1_G","BR1_PR","BR2_BL","BR2_PR","BR3_BL","BR3_G","BR3_PR","RB1_BL","RB1_G","RB1_PR","RB2_BL","RB2_G","RB2_PR","RB3_BL","RB3_G","RB3_PR","RR1_BL","RR1_G","RR1_PR","RR2_BL","RR2_G","RR2_PR","RR3_BL","RR3_G","RR3_PR")

summary(colSums(gcounts))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1173097 3535692 4140541 4274867 4794318 8110571

totalCounts=colSums(gcounts)
totalCounts
# BB1_BL   BB1_G  BB1_PR  BB2_BL   BB2_G  BB2_PR  BB3_BL   BB3_G  BB3_PR  BR1_BL   BR1_G  BR1_PR  BR2_BL 
#3542891 7138590 5373064 4663930 3528492 4256346 4140541 4165290 5966408 5543711 5483689 2753554 3587477 
# BR2_PR  BR3_BL   BR3_G  BR3_PR  RB1_BL   RB1_G  RB1_PR  RB2_BL   RB2_G  RB2_PR  RB3_BL   RB3_G  RB3_PR 
#2899708 3429055 4968640 3616494 8110571 4852457 4736180 4069137 1173097 4591530 4049210 4870062 3741686 
# RR1_BL   RR1_G  RR1_PR  RR2_BL   RR2_G  RR2_PR  RR3_BL   RR3_G  RR3_PR 
#4270730 3400495 4529580 4714959 2466288 3499958 4081959 3414140 3990436 

min(totalCounts) #1173097
mean(totalCounts) #4274867
max(totalCounts)  #8110571

#---------------- Get density of methylation across genes (are there 2 peaks?) yes! ------------------

glength=as.data.frame(gcountsALL[-c(2:5,7:41) ])
colnames(glength)=c("annot", "length")

mns = as.data.frame(apply(gcounts, 1, mean))
row.names(mns)->mns$annot
colnames(mns)=c("means", "annot")

m=merge(mns, glength, by='annot',all.x=T)
m$means<-m$means+0.00001
m$class=-log(m$means/m$length)
dim(m)#30240     4
hist(m$class, breaks=100)
gbm_class=data.frame(cbind("gene"=m$annot,"class"=m$class))
write.csv(gbm_class,file=".../gbm_class.csv",quote=F, row.names=F)

#---------------- Get gene count means and remove low genes ------------------
mns = apply(gcounts, 1, mean)
gcounts=gcounts[mns>3,] #get rid of genes that show little or no methylation
table(mns > 3)
#FALSE  TRUE 
# 5893 24347 
dim(gcounts) #24347    35

totalCountsFilt=colSums(gcounts)
totalCountsFilt
# BB1_BL   BB1_G  BB1_PR  BB2_BL   BB2_G  BB2_PR  BB3_BL   BB3_G  BB3_PR  BR1_BL   BR1_G  BR1_PR  BR2_BL 
#3539981 7124707 5363262 4659516 3526018 4251950 4136327 4157350 5954841 5531505 5476589 2751190 3582640 
# BR2_PR  BR3_BL   BR3_G  BR3_PR  RB1_BL   RB1_G  RB1_PR  RB2_BL   RB2_G  RB2_PR  RB3_BL   RB3_G  RB3_PR 
#2893091 3425375 4962700 3612910 8094731 4846625 4732099 4063090 1171630 4587433 4044546 4863194 3736563 
# RR1_BL   RR1_G  RR1_PR  RR2_BL   RR2_G  RR2_PR  RR3_BL   RR3_G  RR3_PR 
#4265121 3396052 4518217 4710323 2462965 3495980 4076445 3408342 3982961

min(totalCountsFilt) #1171630
max(totalCountsFilt) #8094731
mean(totalCountsFilt) #4268751

#---------------- set up sample metadata ------------------
colData <- data.frame(row.names= colnames(gcounts), bucket= c("BB1","BB1","BB1","BB2","BB2","BB2","BB3","BB3","BB3","BR1","BR1","BR1","BR2","BR2","BR3","BR3","BR3","RB1","RB1","RB1","RB2","RB2","RB2","RB3","RB3","RB3","RR1","RR1","RR1","RR2","RR2","RR2","RR3","RR3","RR3"), stage=c("BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR","BL","G","PR"))

colData$treat = substr(colData$bucket, 1,2)
colData$treat_maternal = substr(colData$bucket, 1,1)
colData$treat_dev = substr(colData$bucket, 2,2)

#---------------- identify sample outliers ------------------
#dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ treat+stage)
#vsd=varianceStabilizingTransformation(dds, blind=T)
#e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
#arrayQualityMetrics(e, intgroup=c("stage"), force=T, outdir= "report_for_genes_vsd_MEDIP_07082018")

#Was determined there are 3 ouliers: BB1_BL, BB2_G, RR2_G

#---------------- remove outliers from both gcounts and colData ------------------
gcounts$BB1_BL<-NULL
gcounts$BB2_G<-NULL
gcounts$RR2_G<-NULL

colData=as.data.frame(colData[-c(1,5,31), ])

#---------------- Wald Tests ------------------

dds<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = formula(~ treat_maternal + treat_dev + stage ))

meth.rld=rlog(dds)

meth.coldata = colData

#### Wald test for maternal treatment
ddsM_RB<-DESeq(dds, minReplicatesForReplace=Inf) 
resM_RB<-results(ddsM_RB, contrast=c('treat_maternal', 'R', 'B')) 
mcols(resM_RB,use.names=TRUE)
summary(resM_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 589, 2.4% 
#LFC < 0 (down)   : 642, 2.6% 
#outliers [1]     : 5, 0.021% 
#low counts [2]   : 0, 0% 
#(mean count < 2)


#### Wald test for developmental treatment
ddsD_RB<-DESeq(dds, minReplicatesForReplace=Inf) 
resD_RB<-results(ddsD_RB, contrast=c('treat_dev', 'R', 'B')) 
mcols(resD_RB,use.names=TRUE)
summary(resD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 5, 0.021% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

#### Wald test for stages
dds_BL_G<-DESeq(dds, minReplicatesForReplace=Inf) 
res_BL_G<-results(dds_BL_G, contrast=c('stage', 'BL', 'G')) 
mcols(res_BL_G,use.names=TRUE)
summary(res_BL_G)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 7, 0.029% 
#LFC < 0 (down)   : 7, 0.029% 
#outliers [1]     : 5, 0.021% 
#low counts [2]   : 13685, 56% 
#(mean count < 68)

dds_BL_PR<-DESeq(dds, minReplicatesForReplace=Inf)
res_BL_PR<-results(dds_BL_PR, contrast=c('stage', 'BL', 'PR'))
mcols(res_BL_PR,use.names=TRUE)
summary(res_BL_PR)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 45, 0.18% 
#LFC < 0 (down)   : 105, 0.43% 
#outliers [1]     : 5, 0.021% 
#low counts [2]   : 9909, 41% 
#(mean count < 35)

dds_G_PR<-DESeq(dds, minReplicatesForReplace=Inf) 
res_G_PR<-results(dds_G_PR, contrast=c('stage', 'G', 'PR'))
mcols(res_G_PR,use.names=TRUE)
summary(res_G_PR)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1, 0.0041% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 5, 0.021% 
#low counts [2]   : 0, 0% 
#(mean count < 2)


save(dds,resM_RB,resD_RB,res_BL_G,res_BL_PR,res_G_PR,meth.rld, meth.coldata,gcounts, file="MeDIP_DESeq2_full_Wald_genes.Rdata")

#---------------- Look at treatment differneces in methylation at each dev stage ------------------

#---------------- Blastula ------------------
gcountsB=as.data.frame(select(gcounts, ends_with("BL")))
colDataB=subset(meth.coldata, stage=="BL")

dds<-DESeqDataSetFromMatrix(gcountsB,

	colData = colDataB, 

	design = formula(~ treat_maternal + treat_dev))

######### test for maternal
dds.mb<-DESeq(dds, minReplicatesForReplace=Inf)
BresM_RB=results(dds.mb, contrast=c('treat_maternal', 'R', 'B')) 
mcols(BresM_RB,use.names=TRUE)
summary(BresM_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 138, 0.57% 
#LFC < 0 (down)   : 146, 0.6% 
#outliers [1]     : 15, 0.062% 
#low counts [2]   : 0, 0% 
#(mean count < 0)


######### test for developmental
dds.db<-DESeq(dds, minReplicatesForReplace=Inf)
BresD_RB=results(dds.db, contrast=c('treat_dev', 'R', 'B')) 
mcols(BresD_RB,use.names=TRUE)
summary(BresD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 4, 0.016% 
#LFC < 0 (down)   : 10, 0.041% 
#outliers [1]     : 15, 0.062% 
#low counts [2]   : 3774, 16% 
#(mean count < 9)


#---------------- Gastrula ------------------
gcountsG=as.data.frame(select(gcounts, ends_with("G")))
colDataG=subset(meth.coldata, stage=="G")

dds<-DESeqDataSetFromMatrix(gcountsG,

	colData = colDataG, 

	design = formula(~ treat_maternal + treat_dev))

######### test for maternal
dds.mg<-DESeq(dds, minReplicatesForReplace=Inf)
GresM_RB=results(dds.mg, contrast=c('treat_maternal', 'R', 'B')) 
mcols(GresM_RB,use.names=TRUE)
summary(GresM_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 128, 0.53% 
#LFC < 0 (down)   : 104, 0.43% 
#outliers [1]     : 8, 0.033% 
#low counts [2]   : 0, 0% 
#(mean count < 1)

######### test for developmental
dds.dg<-DESeq(dds, minReplicatesForReplace=Inf)
GresD_RB=results(dds.dg, contrast = c('treat_dev', 'R', 'B'))
mcols(GresD_RB,use.names=TRUE)
summary(GresD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 10, 0.041% 
#LFC < 0 (down)   : 4, 0.016% 
#outliers [1]     : 8, 0.033% 
#low counts [2]   : 5193, 21% 
#(mean count < 12)


#---------------- Prism ------------------
gcountsP=as.data.frame(select(gcounts, ends_with("PR")))
colDataP=subset(meth.coldata, stage=="PR")

dds<-DESeqDataSetFromMatrix(gcountsP,

	colData = colDataP, 

	design = formula(~ treat_maternal + treat_dev))

######### test for maternal
dds.mp<-DESeq(dds, minReplicatesForReplace=Inf)
PresM_RB=results(dds.dg, contrast = c('treat_maternal', 'R', 'B'))
mcols(PresM_RB,use.names=TRUE)
summary(PresM_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 128, 0.53% 
#LFC < 0 (down)   : 104, 0.43% 
#outliers [1]     : 8, 0.033% 
#low counts [2]   : 0, 0% 
#(mean count < 1)


######### test for developmental
dds.dp<-DESeq(dds, minReplicatesForReplace=Inf)
PresD_RB=results(dds.dp, contrast = c('treat_dev', 'R', 'B'))
summary(PresD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2, 0.0082% 
#LFC < 0 (down)   : 1, 0.0041% 
#outliers [1]     : 8, 0.033% 
#low counts [2]   : 0, 0% 
#(mean count < 2)

#---------------- Look at differences in developmental treatment between each maternal treatment 

# subset samples that come from upwelling maternal condition
gcountsU=as.data.frame(select(gcounts, starts_with("R")))
colDataU=subset(meth.coldata, treat_maternal=="R")

dds<-DESeqDataSetFromMatrix(gcountsU,

	colData = colDataU, 

	design = formula(~ treat_dev))

######### test for developmental
dds.u<-DESeq(dds, minReplicatesForReplace=Inf)
UresD_RB=results(dds.u, contrast=c('treat_dev', 'R', 'B'))
mcols(UresD_RB,use.names=TRUE)
summary(UresD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 9, 0.037% 
#low counts [2]   : 0, 0% 
#(mean count < 0)

# subset samples that come from non-upwelling maternal condition
gcountsN=as.data.frame(select(gcounts, starts_with("B")))
colDataN=subset(meth.coldata, treat_maternal=="B")

dds<-DESeqDataSetFromMatrix(gcountsN,

	colData = colDataN, 

	design = formula(~ treat_dev))

######### test for developmental
dds.n<-DESeq(dds, minReplicatesForReplace=Inf)
NresD_RB=results(dds.n, contrast=c('treat_dev', 'R', 'B')) 
mcols(NresD_RB,use.names=TRUE)
summary(NresD_RB)
#out of 24347 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 1, 0.0041% 
#outliers [1]     : 10, 0.041% 
#low counts [2]   : 0, 0% 
#(mean count < 1)


#---------------- Make output files for rld and p values; get #DMGs for 0.05 cutoff ------------------
ll=load("MeDIP_DESeq2_full_Wald_genes.Rdata")

meth.rld.df = assay(meth.rld)
colnames(meth.rld.df) = colnames(gcounts)

vals=(cbind(resM_RB$stat, resM_RB$pvalue, resM_RB$padj, resD_RB$stat, resD_RB$pvalue, resD_RB$padj, res_G_PR$stat, res_G_PR$pvalue, res_G_PR$padj, res_BL_G$stat, res_BL_G$pvalue, res_BL_G$padj, res_BL_PR$stat, res_BL_PR$pvalue, res_BL_PR$padj)) # collect pvalues for each comparison
colnames(vals)=c("stat_resM_RB","pval_resM_RB", "padj_resM_RB","stat_resD_RB","pval_resD_RB", "padj_resD_RB","stat_res_G_PR","pval_res_G_PR", "padj_res_G_PR","stat_res_BL_G","pval_res_BL_G", "padj_res_BL_G","stat_res_BL_PR","pval_res_BL_PR", "padj_res_BL_PR") # name columns
rldpvals=as.data.frame(cbind(meth.rld.df,vals)) # combine RLD with pvals

#get number of DMGs with padj < 0.01
#maternal condition
length((rldpvals[rldpvals$padj_resM_RB<= 0.01 & !is.na(rldpvals$padj_resM_RB), ])[,1]) #605

#developmental condition
length((rldpvals[rldpvals$padj_resD_RB<= 0.01 & !is.na(rldpvals$padj_resD_RB), ])[,1]) #0

#Gastrula v. Prism 
length((rldpvals[rldpvals$padj_res_G_PR<= 0.01 & !is.na(rldpvals$padj_res_G_PR), ])[,1]) #1

#Blastula v. Prism 
length((rldpvals[rldpvals$padj_res_BL_PR<= 0.01 & !is.na(rldpvals$padj_res_BL_PR), ])[,1]) #22

#Blastula v. Gastrula 
length((rldpvals[rldpvals$padj_res_BL_G<= 0.01 & !is.na(rldpvals$padj_res_BL_G), ])[,1]) #0

write.csv(rldpvals,file="MeDIP_Spurp_rldpvals.csv",quote=F, row.names=T)

ll=read.csv("MeDIP_Spurp_rldpvals.csv")
llsig = ll[order(ll$pval_resM_RB),]
llsig=llsig[llsig$pval_resM_RB<= 0.05 & !is.na(llsig$pval_resM_RB), ]
row.names(llsig)=llsig$X
llsig=as.data.frame(llsig[-c(1,34:48) ],) #collect wanted data columns for output
dim(llsig) #3123   32

write.csv(llsig,file="SigGenes_Maternal0.05.csv",quote=F, row.names=T)

#---------------- principal coordinate analysis ------------------
ll=load("MeDIP_DESeq2_full_Wald_genes.Rdata")

meth.rld.df = assay(meth.rld)
colnames(meth.rld.df) = colnames(gcounts)

dd.veg=vegdist(t(meth.rld.df), "manhattan")
div.dd.veg=dd.veg/1000

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

conditions=meth.coldata
conditions$treat_maternal=as.factor(conditions$treat_maternal)
conditions$treat_dev=as.factor(conditions$treat_dev)

dd.pcoa
#$values
#   Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1    94.636783   0.22620636  0.129911135 0.2262064      0.1299111
#2    54.482743   0.13022783  0.097653071 0.3564342      0.2275642
#3    25.012229   0.05978569  0.081524039 0.4162199      0.3090882
#4    20.507393   0.04901797  0.070771350 0.4652378      0.3798596
#5    17.006008   0.04064875  0.062706834 0.5058866      0.4425664

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
#treat_maternal  1  27427677 27427677 2.09978 0.06556  0.015 *
#treat_dev       1  10961329 10961329 0.83917 0.02620  0.631  
#stage           2  27297537 13648768 1.04491 0.06525  0.348  
#Residuals      27 352678282 13062159         0.84299         
#Total          31 418364825                  1.00000         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

labs=c("Maternal_condition","Developmental_condition","Stage","Residuals")
cols=c("skyblue","green2","coral","grey80")

labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="GBM")

#---------------- OUTPUT FOR GO-MWU TESTS ------------------

### for maternal condition
resM_RB$lp=-log(resM_RB$pvalue,10)
resM_RB$lp[resM_RB$log2FoldChange<0]=-resM_RB$lp[resM_RB$log2FoldChange<0]

maternal_lpv=data.frame(cbind("gene"=row.names(resM_RB),"lp"=resM_RB$lp))
maternal_fc=data.frame(cbind("gene"=row.names(resM_RB),"fc"=resM_RB$log2FoldChange))
maternal_stat=data.frame(cbind("gene"=row.names(resM_RB),"stat"=resM_RB$stat))
write.csv(maternal_lpv,file="ReanalysisJuly2018/GOMWU/maternal_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(maternal_stat,file="ReanalysisJuly2018/GOMWU/maternal_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(maternal_fc,file="ReanalysisJuly2018/GOMWU/maternal_wald_fc_mrna.csv",quote=F, row.names=F)

### output for fishers test for maternal condition
resM_RB$padj_go[resM_RB$padj<=0.01]<-1
resM_RB$padj_go[resM_RB$padj>0.01]<-0 
resM_RB$padj_go=as.factor(resM_RB$padj_go) 
summary(resM_RB$padj_go)
#    0     1 
#23737   605 
resM_RB=resM_RB[complete.cases(resM_RB),]
resM_RB$padj_go=as.character(resM_RB$padj_go)
maternal_fishers=data.frame(cbind("gene"=row.names(resM_RB),"padj_GO"=resM_RB$padj_go))
write.csv(maternal_fishers,file="ReanalysisJuly2018/GOMWU/maternal_wald_lpv_mrna.csv",quote=F, row.names=F)

# for larval condition
resD_RB$lp=-log(resD_RB$pvalue,10)
resD_RB$lp[resD_RB$log2FoldChange<0]=-resD_RB$lp[resD_RB$log2FoldChange<0]

dev_lpv=data.frame(cbind("gene"=row.names(resD_RB),"lp"=resD_RB$lp))
dev_fc=data.frame(cbind("gene"=row.names(resD_RB),"fc"=resD_RB$log2FoldChange))
dev_stat=data.frame(cbind("gene"=row.names(resD_RB),"stat"=resD_RB$stat))
write.csv(dev_lpv,file="ReanalysisJuly2018/GOMWU/larval_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(dev_stat,file="ReanalysisJuly2018/GOMWU/larval_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(dev_fc,file="ReanalysisJuly2018/GOMWU/larval_wald_fc_mrna.csv",quote=F, row.names=F)

#for comparision between BL and G
res_BL_G$lp=-log(res_BL_G$pvalue,10)
res_BL_G$lp[res_BL_G$log2FoldChange<0]=-res_BL_G$lp[res_BL_G$log2FoldChange<0]

BL_G_lpv=data.frame(cbind("gene"=row.names(res_BL_G),"lp"=res_BL_G$lp))
BL_G_fc=data.frame(cbind("gene"=row.names(res_BL_G),"fc"=res_BL_G$log2FoldChange))
BL_G_stat=data.frame(cbind("gene"=row.names(res_BL_G),"stat"=res_BL_G$stat))
write.csv(BL_G_lpv,file="ReanalysisJuly2018/GOMWU/BL_G_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(BL_G_stat,file="ReanalysisJuly2018/GOMWU/BL_G_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(BL_G_fc,file="ReanalysisJuly2018/GOMWU/BL_G_wald_fc_mrna.csv",quote=F, row.names=F)

#for comparision between G and PR
res_G_PR$lp=-log(res_G_PR$pvalue,10)
res_G_PR$lp[res_G_PR$log2FoldChange<0]=-res_G_PR$lp[res_G_PR$log2FoldChange<0]

G_PR_lpv=data.frame(cbind("gene"=row.names(res_G_PR),"lp"=res_G_PR$lp))
G_PR_fc=data.frame(cbind("gene"=row.names(res_G_PR),"fc"=res_G_PR$log2FoldChange))
G_PR_stat=data.frame(cbind("gene"=row.names(res_G_PR),"stat"=res_G_PR$stat))
write.csv(G_PR_lpv,file="ReanalysisJuly2018/GOMWU/G_PR_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(G_PR_stat,file="ReanalysisJuly2018/GOMWU/G_PR_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(G_PR_fc,file="ReanalysisJuly2018/GOMWU/G_PR_wald_fc_mrna.csv",quote=F, row.names=F)

#for comparision between BL and PR
res_BL_PR$lp=-log(res_BL_PR$pvalue,10)
res_BL_PR$lp[res_BL_PR$log2FoldChange<0]=-res_BL_PR$lp[res_BL_PR$log2FoldChange<0]

BL_PR_lpv=data.frame(cbind("gene"=row.names(res_BL_PR),"lp"=res_BL_PR$lp))
BL_PR_fc=data.frame(cbind("gene"=row.names(res_BL_PR),"fc"=res_BL_PR$log2FoldChange))
BL_PR_stat=data.frame(cbind("gene"=row.names(res_BL_PR),"stat"=res_BL_PR$stat))
write.csv(BL_PR_lpv,file="ReanalysisJuly2018/GOMWU/BL_PR_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(BL_PR_stat,file="ReanalysisJuly2018/GOMWU/BL_PR_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(BL_PR_fc,file="ReanalysisJuly2018/GOMWU/BL_PR_wald_fc_mrna.csv",quote=F, row.names=F)


#---------------- volcano plots ------------------
lnames=load("MeDIP_DESeq2_full_Wald_genes.Rdata")

with(resM_RB, plot(log2FoldChange, -log10(pvalue), pch=20, main="Maternal Condition Padj<0.01", xlim=c(-5,5), ylim=c(1,45)))
with(subset(resM_RB, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))

with(resD_RB, plot(log2FoldChange, -log10(pvalue), pch=20, main="Developmental Condition Padj<0.01", xlim=c(-5,5),ylim=c(1,45)))
with(subset(resD_RB, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(res_BL_PR, plot(log2FoldChange, -log10(pvalue), pch=20, main="Blastula v. Prism Condition Padj<0.01", xlim=c(-5,5),ylim=c(1,45)))
with(subset(res_BL_PR, padj<0.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))


#---------------- DAPC ------------------
meth.rld.df=assay(meth.rld)
degs<-rownames(meth.rld.df)

a.rld<-rld[,c(grep("BB",colnames(rld)), grep("RR",colnames(rld)))]
a.rld.supp<-rld[,c(grep("BR",colnames(rld)),grep("RB",colnames(rld)))]

##Calling PCs
pcp=prcomp(t(a.rld[degs,]), retx=TRUE, center=TRUE, scale=TRUE)
scores=pcp$x
screeplot(pcp,bstick=T)

clus=find.clusters(t(a.rld[degs,]),max.n.clus=15, n.clust=2, n.pca=4)

colnames(rld)
clus$grp=c(substr(colnames(a.rld), start=1,stop=2))

#Discriminate between maternal conditions (B or R), keep 1 DF, use PCs to account for 80% of var, sets up axis for upwelling mom-like GBM and nonupwelling mom-like GBM
dp=dapc(t(a.rld[degs,]),clus$grp, n.da=1, perc.pca=80)

scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

#cluster based on dev treatment
pred.sup<-predict.dapc(dp,newdata=(t(a.rld.supp[degs,])))

names(pred.sup)
pred.sup$assign
colnames(a.rld.supp)

test<-dp
test$ind.coord<-pred.sup$ind.scores
test$posterior<-pred.sup$posterior
test$assign<-pred.sup$assign
#test$grp<-as.factor(test$assign)
test$grp<-as.factor(substr(colnames(a.rld.supp), start =1, stop=2))

scatter(test,bg="white",scree.da=FALSE,legend=TRUE,solid=.4,ylim=c(0,0.5),xlim=c(-6,5)) 

dpc=data.frame(rbind(dp$ind.coord,pred.sup$ind.scores))

write.csv(dpc,"MeDIP_DFA.csv",quote=F)

#---------------- prepare morphometric data ------------------

#setwd("/Users/mariestrader/Dropbox/MEDIP/PPhoxMorphometrics/morph") #change to where your files are stored
gas=read.csv("PurplepHox_Gastrula_edited.csv") #load in raw file
hb=read.csv("PurplepHox_HatchedBlastula_corrected_edited.csv") #load in raw file
pr=read.csv("PurplepHox_Prism_corrected_edited.csv") #load in raw file

#### Gastrula
gas$Diameter=as.numeric(gas$Diameter)

all.G=summarySE(gas,measurevar="Diameter",groupvars=c("tubeLabel"), na.rm=T)
all.G

all.G$tubeLabel=sub("G_","", all.G$tubeLabel)
all.G$tubeLabel=paste0(all.G$tubeLabel, "_G")
colnames(all.G)=c("bucket","N","length","sd","se","ci")

#### Hatched blastula
hb$Diameter=as.numeric(hb$length)
all.BL=summarySE(hb,measurevar="length",groupvars=c("tube_label"), na.rm=T)
all.BL
all.BL$tube_label=paste0(all.BL$tube_label, "_BL")
colnames(all.margBL)=c("bucket","N","length","sd","se","ci")

#### Prism
pr$length=as.numeric(pr$length)
all.PR=summarySE(pr,measurevar="length",groupvars=c("Bucket"), na.rm=T)
all.PR
all.PR$Bucket=paste0(all.PR$Bucket, "_PR")
colnames(all.PR)=c("bucket","N","length","sd","se","ci")

#Combine length data for all stages
length=rbind(all.BL,all.G,all.PR)
length$stage<-as.factor(substr(length$bucket, start =5, stop=6))

write.csv(length,"morphData_dapc.csv",quote=F)

#---------------- correlate morphometric data with DAPC ------------------

dpc=read.csv("MeDIP_DFA.csv")
length=read.csv("/Users/mariestrader/Dropbox/MEDIP/morphData_dapc.csv")
colnames(dpc)=c("bucket","LD1")
m=merge(length, dpc, by='bucket', all.x=T)
m$treat<-as.factor(substr(m$bucket, start =1, stop=2))
m=m[complete.cases(m),]
m$treat_maternal = factor(substr(m$treat, 1,1))
m$treat_larval = factor(substr(m$treat, 2,2))

treat_colors <- c("blue", "turquoise4", "coral","red")
m$treat_c <- treat_colors[as.factor(m$treat)]

plot(m$length~m$LD1, pch=20, col=m$treat_c, xlim=c(-6,5),xlab="Discriminant Function", ylab="Embryo Length in mm")
lm1=lm(m$length~m$LD1)
summary(lm1)
#lm(formula = m$length ~ m$LD1)
#
#Residuals:
#      Min        1Q    Median        3Q       Max 
#-0.011570 -0.006462 -0.001613  0.006509  0.014196 
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.1137131  0.0014845  76.603   <2e-16 ***
#m$LD1       0.0018398  0.0007783   2.364    0.025 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.008257 on 29 degrees of freedom
#Multiple R-squared:  0.1616,	Adjusted R-squared:  0.1327 
#F-statistic: 5.588 on 1 and 29 DF,  p-value: 0.02499
abline(lm1)
legend("topleft", legend = levels(m$treat), col=treat_colors, pch = 19)



