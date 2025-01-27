
library(future.apply)
library(dplyr)
library(limma)
library(Rsubread)
library(trqwe)
library(ggplot2)
library(ChIPQC)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(limma)
library(Rsubread)
library(GenomicRanges)
library(DESeq2)
library(ggrepel)
library(data.table)
library(dplyr)
library(limma)
library(Rsubread)
library(GenomicRanges)



***************************************************************************************

hs_TEs <- mcreadRDS("/mnt/data/user_data/yiman/workshop/RNAseq_ref/TEtranscripts/hg38_rmsk_TE.rds",mc.cores=20)
hs_TEs$id <- paste0(hs_TEs[,1],"_",hs_TEs[,2],"_",hs_TEs[,3])
hs_TEs <- hs_TEs[!duplicated(hs_TEs$id),]
hs_TEs <- toGRanges(hs_TEs)

NF391.rmigg <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.peaks.rmIgG.csv")
NF391.rmigg.ol <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.peaks.rmIgG.ol.csv")
nrow(NF391.rmigg)
nrow(NF391.rmigg.ol)
NF391.rmigg <- toGRanges(NF391.rmigg)
NF391.rmigg.ol <- toGRanges(NF391.rmigg.ol)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)

NF391_vs_IgG_up <- read.table(file = "/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391_vs_IgG_up.bed",
    sep="\t",header=FALSE)
NF391_vs_IgG_up <- NF391_vs_IgG_up[,1:3]
colnames(NF391_vs_IgG_up) <- c("chr","start","end")
NF391_vs_IgG_up <- toGRanges(NF391_vs_IgG_up)

NF391_vs_IgG_up.TE <- annoPeaks(NF391_vs_IgG_up, 
                       annoData = hs_TEs,
                       bindingType='fullRange', 
                       bindingRegion=c(-1,1),
                       ignore.peak.strand=TRUE, 
                       select='bestOne')
pie1(table(NF391_vs_IgG_up.TE$insideFeature))
4573/10589
write.csv(as.data.frame(NF391_vs_IgG_up.TE),"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.csv")

NF391_vs_IgG_up.TE <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.csv")
NF391_vs_IgG_up.TE$class_id <- as.character(NF391_vs_IgG_up.TE$class_id)
NF391_vs_IgG_up.TE$TE.gene.class <- "others"
# NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE")),]$TE.gene.class <- paste0(NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE")),]$class_id,"-",NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE")),]$family_id)
NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$TE.gene.class <- NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$class_id
table(NF391_vs_IgG_up.TE$TE.gene.class)
write.table(NF391_vs_IgG_up.TE, file = "/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.bed",
    sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

pie1(table(NF391_vs_IgG_up.TE$insideFeature))
pie1(table(NF391_vs_IgG_up.TE$family_id))
pie1(table(NF391_vs_IgG_up.TE$class_id))
pie1(table(NF391_vs_IgG_up.TE$TE.gene.class))

pdf("/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/TE.class.NF391.peaks.final.pdf")
pie1(table(NF391_vs_IgG_up.TE$TE.gene.class))
dev.off()


cd /mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/bw_files/

computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
-R /mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.bed \
-S IgG.bw combine_all_p391.bw \
--numberOfProcessors 30 --skipZeros \
--missingDataAsZero \
-o /mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/center.NF391.vs.IgG.up.TE.mat.gz 

plotHeatmap -m /mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/center.NF391.vs.IgG.up.TE.mat.gz \
 -out /mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/center.NF391.vs.IgG.up.TE.pdf  \
 --colorList  'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' \
 --missingDataColor "white"

####################################################

NF391_vs_IgG_up.TE <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.csv")
NF391_vs_IgG_up.TE$class_id <- as.character(NF391_vs_IgG_up.TE$class_id)
NF391_vs_IgG_up.TE$TE.gene.class <- "other_TE"
NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$TE.gene.class <- NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$class_id
table(NF391_vs_IgG_up.TE$TE.gene.class)
TE.res <- NF391_vs_IgG_up.TE$TE.gene.class
TE.res <- as.character(TE.res)
NF391.peak.all <- data.frame(ID=c(TE.res,rep("Non-TE_peaks",10589-4573)))
pie1(table(NF391.peak.all$ID))

pdf("/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/all.class.NF391.peaks.final.pdf")
pie1(table(NF391.peak.all$ID))
dev.off()

####################################################

NF391_vs_IgG_up.TE <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391.vs.IgG.up.TE.anno.csv")
NF391_vs_IgG_up.TE$class_id <- as.character(NF391_vs_IgG_up.TE$class_id)
NF391_vs_IgG_up.TE$TE.gene.class <- "other_TE"
NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$TE.gene.class <- NF391_vs_IgG_up.TE[which(NF391_vs_IgG_up.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$class_id
table(NF391_vs_IgG_up.TE$TE.gene.class)
TE.res <- NF391_vs_IgG_up.TE$TE.gene.class
TE.res <- as.character(TE.res)

NF391_vs_IgG_up <- read.table(file = "/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391_vs_IgG_up.bed",
    sep="\t",header=FALSE)
NF391_vs_IgG_up$id <- paste0(NF391_vs_IgG_up$V1,"_",NF391_vs_IgG_up$V2,"_",NF391_vs_IgG_up$V3)
NF391_vs_IgG_up.TE$id <- paste0(NF391_vs_IgG_up.TE$seqnames,"_",NF391_vs_IgG_up.TE$start,"_",NF391_vs_IgG_up.TE$end)
NF391_vs_IgG_up.TE$id <- as.character(NF391_vs_IgG_up.TE$id)
nonTE.id <- setdiff(NF391_vs_IgG_up$id,NF391_vs_IgG_up.TE$id)
nrow(NF391_vs_IgG_up)
nrow(NF391_vs_IgG_up.TE)
length(nonTE.id)

nonTE.peak <- NF391_vs_IgG_up[which(NF391_vs_IgG_up$id %in% nonTE.id),]
nrow(nonTE.peak)
nonTE.peak <- nonTE.peak[,1:3]
colnames(nonTE.peak) <- c("chr","start","end")

nonTE.peak.anno <- annotatePeak(GRanges(nonTE.peak), tssRegion=c(-1, 1),
                         TxDb=txdb, annoDb="org.Hs.eg.db",verbose=FALSE)
nonTE.peak.anno <- as.data.frame(nonTE.peak.anno)
nonTE.peak.anno$annotation_tidy <- nonTE.peak.anno$annotation
nonTE.peak.anno[grep("Intron",nonTE.peak.anno$annotation),]$annotation_tidy <- "Intron"
nonTE.peak.anno[grep("Exon",nonTE.peak.anno$annotation),]$annotation_tidy <- "Exon"
nonTE.peak.anno[grep("Downstream",nonTE.peak.anno$annotation),]$annotation_tidy <- "Downstream"
table(nonTE.peak.anno$annotation_tidy)
nonTE.peak.anno$id <- paste0(nonTE.peak.anno$seqnames,"_",nonTE.peak.anno$start,"_",nonTE.peak.anno$end)

NF391_vs_IgG_up$anno <- "Unknown"
match_indices <- match(NF391_vs_IgG_up$id, NF391_vs_IgG_up.TE$id)
NF391_vs_IgG_up$anno[!is.na(match_indices)] <- NF391_vs_IgG_up.TE$TE.gene.class[match_indices[!is.na(match_indices)]]
match_indices <- match(NF391_vs_IgG_up$id, nonTE.peak.anno$id)
NF391_vs_IgG_up$anno[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(NF391_vs_IgG_up$anno)

NF391_vs_IgG_up$anno11 <- "Unknown"
match_indices <- match(NF391_vs_IgG_up$id, NF391_vs_IgG_up.TE$id)
NF391_vs_IgG_up$anno11[!is.na(match_indices)] <- "TEs"
match_indices <- match(NF391_vs_IgG_up$id, nonTE.peak.anno$id)
NF391_vs_IgG_up$anno11[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(NF391_vs_IgG_up$anno11)
write.csv(NF391_vs_IgG_up,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391_vs_IgG_up_anno_final.csv")

NF391_vs_IgG_up <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/NF391_vs_IgG_up_anno_final.csv")
pie1(table(NF391_vs_IgG_up$anno))
pie1(table(NF391_vs_IgG_up$anno11))
tmp <- as.data.frame(table(NF391_vs_IgG_up$anno11))
tmp$pct <- tmp$Freq / sum(tmp$Freq)

pdf("/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/all.class.NF391.peaks.final.final.pdf")
pie1(table(NF391_vs_IgG_up$anno11))
dev.off()


######################ATAC##############################

library(ChIPpeakAnno)

ATAC_DEG <- read.table(file = "/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.bed")
colnames(ATAC_DEG)[1:3] <- c("chr","start","end")
ATAC_DEG <- ATAC_DEG[,1:3]
ATAC_DEG <- toGRanges(ATAC_DEG) #7314

hs_TEs <- mcreadRDS("/mnt/data/user_data/yiman/workshop/RNAseq_ref/TEtranscripts/hg38_rmsk_TE.rds",mc.cores=20)
hs_TEs$id <- paste0(hs_TEs[,1],"_",hs_TEs[,2],"_",hs_TEs[,3])
hs_TEs <- hs_TEs[!duplicated(hs_TEs$id),]
hs_TEs <- toGRanges(hs_TEs)

ATAC.TE <- annoPeaks(ATAC_DEG, 
                       annoData = hs_TEs,
                       bindingType='fullRange', 
                       bindingRegion=c(-1,1),
                       ignore.peak.strand=TRUE, 
                       select='bestOne')
6783/7314

write.csv(as.data.frame(ATAC.TE),"/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.TE.anno.csv")

ATAC.TE <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.TE.anno.csv")
ATAC.TE$class_id <- as.character(ATAC.TE$class_id)
ATAC.TE$TE.gene.class <- "other_TE"
ATAC.TE[which(ATAC.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$TE.gene.class <- ATAC.TE[which(ATAC.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$class_id
table(ATAC.TE$TE.gene.class)
TE.res <- ATAC.TE$TE.gene.class
TE.res <- as.character(TE.res)

ATAC_DEG <- read.table(file = "/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.bed")
colnames(ATAC_DEG)[1:3] <- c("chr","start","end")
ATAC_DEG <- ATAC_DEG[,1:3]

ATAC_DEG$id <- paste0(ATAC_DEG$chr,"_",ATAC_DEG$start,"_",ATAC_DEG$end)
ATAC.TE$id <- paste0(ATAC.TE$seqnames,"_",ATAC.TE$start,"_",ATAC.TE$end)
ATAC.TE$id <- as.character(ATAC.TE$id)
nonTE.id <- setdiff(ATAC_DEG$id,ATAC.TE$id)
nrow(ATAC_DEG)
nrow(ATAC.TE)
length(nonTE.id)

nonTE.peak <- ATAC_DEG[which(ATAC_DEG$id %in% nonTE.id),]
nrow(nonTE.peak)
nonTE.peak <- nonTE.peak[,1:3]
colnames(nonTE.peak) <- c("chr","start","end")

nonTE.peak.anno <- annotatePeak(GRanges(nonTE.peak), tssRegion=c(-1, 1),
                         TxDb=txdb, annoDb="org.Hs.eg.db",verbose=FALSE)
nonTE.peak.anno <- as.data.frame(nonTE.peak.anno)
nonTE.peak.anno$annotation_tidy <- nonTE.peak.anno$annotation
nonTE.peak.anno[grep("Intron",nonTE.peak.anno$annotation),]$annotation_tidy <- "Intron"
nonTE.peak.anno[grep("Exon",nonTE.peak.anno$annotation),]$annotation_tidy <- "Exon"
nonTE.peak.anno[grep("Downstream",nonTE.peak.anno$annotation),]$annotation_tidy <- "Downstream"
table(nonTE.peak.anno$annotation_tidy)
nonTE.peak.anno$id <- paste0(nonTE.peak.anno$seqnames,"_",nonTE.peak.anno$start,"_",nonTE.peak.anno$end)

ATAC_DEG$anno <- "Unknown"
match_indices <- match(ATAC_DEG$id, ATAC.TE$id)
ATAC_DEG$anno[!is.na(match_indices)] <- ATAC.TE$TE.gene.class[match_indices[!is.na(match_indices)]]
match_indices <- match(ATAC_DEG$id, nonTE.peak.anno$id)
ATAC_DEG$anno[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(ATAC_DEG$anno)

ATAC_DEG$anno11 <- "Unknown"
match_indices <- match(ATAC_DEG$id, ATAC.TE$id)
ATAC_DEG$anno11[!is.na(match_indices)] <- "TEs"
match_indices <- match(ATAC_DEG$id, nonTE.peak.anno$id)
ATAC_DEG$anno11[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(ATAC_DEG$anno11)
write.csv(ATAC_DEG,"/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.TE.anno.final.csv")

ATAC_DEG <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/DEG.peaks.4samples.TE.anno.final.csv")
pie1(table(ATAC_DEG$anno))
pie1(table(ATAC_DEG$anno11))
tmp <- as.data.frame(table(ATAC_DEG$anno11))
tmp$pct <- tmp$Freq / sum(tmp$Freq)

pdf("/mnt/data/user_data/yiman/workshop/ATAC_Seq/IPM/IPM_A549/outs/all.class.ATAC.peaks.final.final.pdf")
pie1(table(ATAC_DEG$anno11))
dev.off()

######################H3K9me3##############################

library(data.table)
library(ChIPpeakAnno)

H3K9_DEG <- fread(file = "/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_broad/DEG/H3K9me3.DEG.peaks.p0.05.FC1.bed")
H3K9_DEG <- as.data.frame(H3K9_DEG)
colnames(H3K9_DEG)[1:3] <- c("chr","start","end")
H3K9_DEG <- H3K9_DEG[,1:3]
H3K9_DEG <- toGRanges(H3K9_DEG) #3903

H3K9.TE <- annoPeaks(H3K9_DEG, 
                       annoData = hs_TEs,
                       bindingType='fullRange', 
                       bindingRegion=c(-1,1),
                       ignore.peak.strand=TRUE, 
                       select='bestOne')
3249/3903

write.csv(as.data.frame(H3K9.TE),"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/H3K9.DEG.peaks.p0.05.FC1.TE.anno.csv")

H3K9.TE <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/H3K9.DEG.peaks.p0.05.FC1.TE.anno.csv")
H3K9.TE$class_id <- as.character(H3K9.TE$class_id)
H3K9.TE$TE.gene.class <- "other_TE"
H3K9.TE[which(H3K9.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$TE.gene.class <- H3K9.TE[which(H3K9.TE$class_id %in% c("DNA","LINE","LTR","SINE","Satellite")),]$class_id
table(H3K9.TE$TE.gene.class)
TE.res <- H3K9.TE$TE.gene.class
TE.res <- as.character(TE.res)

H3K9_DEG <- as.data.frame(H3K9_DEG)
colnames(H3K9_DEG)[1:3] <- c("chr","start","end")
H3K9_DEG <- H3K9_DEG[,1:3]

H3K9_DEG$id <- paste0(H3K9_DEG$chr,"_",H3K9_DEG$start,"_",H3K9_DEG$end)
H3K9.TE$id <- paste0(H3K9.TE$seqnames,"_",H3K9.TE$start,"_",H3K9.TE$end)
H3K9.TE$id <- as.character(H3K9.TE$id)
nonTE.id <- setdiff(H3K9_DEG$id,H3K9.TE$id)
nrow(H3K9_DEG)
nrow(H3K9.TE)
length(nonTE.id)

nonTE.peak <- H3K9_DEG[which(H3K9_DEG$id %in% nonTE.id),]
nrow(nonTE.peak)
nonTE.peak <- nonTE.peak[,1:3]
colnames(nonTE.peak) <- c("chr","start","end")

nonTE.peak.anno <- annotatePeak(GRanges(nonTE.peak), tssRegion=c(-1, 1),
                         TxDb=txdb, annoDb="org.Hs.eg.db",verbose=FALSE)
nonTE.peak.anno <- as.data.frame(nonTE.peak.anno)
nonTE.peak.anno$annotation_tidy <- nonTE.peak.anno$annotation
nonTE.peak.anno[grep("Intron",nonTE.peak.anno$annotation),]$annotation_tidy <- "Intron"
nonTE.peak.anno[grep("Exon",nonTE.peak.anno$annotation),]$annotation_tidy <- "Exon"
nonTE.peak.anno[grep("Downstream",nonTE.peak.anno$annotation),]$annotation_tidy <- "Downstream"
table(nonTE.peak.anno$annotation_tidy)
nonTE.peak.anno$id <- paste0(nonTE.peak.anno$seqnames,"_",nonTE.peak.anno$start,"_",nonTE.peak.anno$end)

H3K9_DEG$anno <- "Unknown"
match_indices <- match(H3K9_DEG$id, H3K9.TE$id)
H3K9_DEG$anno[!is.na(match_indices)] <- H3K9.TE$TE.gene.class[match_indices[!is.na(match_indices)]]
match_indices <- match(H3K9_DEG$id, nonTE.peak.anno$id)
H3K9_DEG$anno[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(H3K9_DEG$anno)

H3K9_DEG$anno11 <- "Unknown"
match_indices <- match(H3K9_DEG$id, H3K9.TE$id)
H3K9_DEG$anno11[!is.na(match_indices)] <- "TEs"
match_indices <- match(H3K9_DEG$id, nonTE.peak.anno$id)
H3K9_DEG$anno11[!is.na(match_indices)] <- nonTE.peak.anno$annotation_tidy[match_indices[!is.na(match_indices)]]
table(H3K9_DEG$anno11)
write.csv(H3K9_DEG,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/H3K9.DEG.peaks.p0.05.FC1.TE.anno.final.csv")

H3K9_DEG <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/H3K9.DEG.peaks.p0.05.FC1.TE.anno.final.csv")
pie1(table(H3K9_DEG$anno))
pie1(table(H3K9_DEG$anno11))
tmp <- as.data.frame(table(H3K9_DEG$anno11))
tmp$pct <- tmp$Freq / sum(tmp$Freq)


pdf("/mnt/data/user_data/yiman/workshop/CUT-TAG/IPM/IPM_A549/MACS3_new/hs_q0.05_narrow/DEG/all.class.H3K9me3.peaks.final.final.pdf")
pie1(table(H3K9_DEG$anno11))
dev.off()


