setwd("/home/data/vip9t07/livercan")
#read data
data <- read.csv("liver_fpkm.csv",header = T )
#colnames(data) <- c("gene_id",	"hu_con_4",	"hu_con_5",	"lung_1",	"lung_2",	"gene_name","gene_chr",	"gene_start",	"gene_end",	"gene_strand",	"gene_length",	"gene_biotype",	"gene_description",	"tf_family")
#data <- filter(data,!duplicated(data$gene_name))
#data <- data[!duplicated(data$gene_name),]##to remove duplicated values
#mean data according to duplicated rows
data1 <- data[,c(2:19)]
library(dplyr)
data1 <- data1 %>%
  group_by(gene_name) %>%
  summarise_all(mean)
data1 <- as.data.frame(data1)

row.names(data1) <- data1[,1]
data1 <- data1[,-1]
datafilt <- data1

#fpkm to tpm
expMatrix <- datafilt
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
#write.csv(tpms,file = "liver_exprset.csv")

#differential analysis
#to assign order
library(limma)
exprSet <- tpms
exprSet <- exprSet[,c(1:6,15:17)]

#group
group_list <- factor(c(rep("ctrl",6),rep("cirrhosis",3)), levels = c("ctrl","cirrhosis"),ordered = F)
group_list
table(group_list)

#data correction
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#to test if data needed transfer
exprSet <- log2(exprSet+1)
dat <- exprSet

################ !!!ATTENTION!!! 2 methods for differential analysis!!!
# 1. without makecontrast
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg1 <- topTable(fit,coef=2,adjust='BH',number = Inf)
# 2. with makecontrast
design2 <- model.matrix(~0+group_list)
colnames(design2) = levels(factor(group_list))
rownames(design2) = colnames(group_list)
design2
cont.matrix <- makeContrasts('cirrhosis-ctrl', levels = design2)
cont.matrix
fit2=lmFit(dat,design2)
fit2=contrasts.fit(fit2,cont.matrix)
fit2=eBayes(fit2)
options(digits = 4)
deg2=topTable(fit2,adjust='BH',number = Inf)
head(deg2) 
# get deg
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}

#differential expression genes
if(T){
  logFC_t=1
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  
  save(DEG,file = 'anno_DEG.Rdata')
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
}
write.csv(DEG,file = "NC-vs-CONTROL.csv")
##GSEA
gene=data.frame(DEG$symbol,DEG$logFC,stringsAsFactors=FALSE)
geneID=select(org.Hs.eg.db,keys=DEG$symbol,columns="ENTREZID",keytype="SYMBOL")
geneID=na.omit(geneID)
colnames(gene) <- c("SYMBOL", "log2FoldChange")
tmp=left_join(geneID,gene,by="SYMBOL")
genelist=tmp$log2FoldChange
names(genelist)=tmp$ENTREZID
head(genelist)
genelist_sort=sort(genelist,decreasing = T)
head(genelist_sort)
#by GO
gse_ALL <- gseGO(genelist_sort,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont="ALL",
                 pvalueCutoff = 0.9)
head(gsemf)
library(enrichplot)
gseaplot2(gse_ALL, 1,title = gsemf@result$Description[1],pvalue_table = T)
gseaplot2(gse_ALL, "GO:0043542",pvalue_table = T)
write.csv(gse_ALL@result,file = "gse_all_GO.csv")
#by KEGG
gse_KEGG <- gseKEGG(genelist_sort ,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
gseaplot2(gse_KEGG, "hsa04380",pvalue_table = T)
write.csv(gse_KEGG@result,file = "gse_all_KEGG.csv")

#GO enrich
gene <- DEG$symbol
gene.df<-bitr(gene, fromType = "SYMBOL",
              toType = c("ENTREZID"),
              OrgDb = org.Hs.eg.db)
ego_cc<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_bp<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
dotplot(ego_cc,
        showCategory = 30,
        title="The GO_CC enrichment analysis of all DEGs ",
        color = 'pvalue')
ego_cc@result$qvalue = 0
ego_cc@result$p.adjust = 0

dotplot(ego_bp,
        showCategory = 30,
        title="The GO_BP enrichment analysis of all DEGs ",
        color = 'pvalue')
ego_bp@result$qvalue = 0
ego_bp@result$p.adjust = 0

dotplot(ego_mf,
        showCategory = 30,
        title="The GO_MF enrichment analysis of all DEGs ",
        color = 'pvalue')
ego_mf@result$qvalue = 0
ego_mf@result$p.adjust = 0

#KEGG enrich
kk<-enrichKEGG(gene      =gene.df$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.9)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk@result$qvalue = 0
kk@result$p.adjust = 0
barplot(kk,showCategory = 30, title="The KEGG enrichment analysis of all DEGs",color = 'pvalue')
cnetplot(kk,foldChange=gene, circular = TRUE, colorEdge = TRUE)
emapplot(kk,foldChange=gene, circular = TRUE, colorEdge = TRUE,color = "pvalue")

FC <- DEG$logFC
names(FC) <- DEG$symbol
class(FC)
FC
select_FC <- FC[c('METTL3','METTL14','FOXO1','FOXO3','ALDOA','PGK1','PYGL','ENO1','PGAM1','HK1','PFKFB3','PGM2','PCNA','RPE','KRAS','ETF1','LDHA','CDK1','PPIA','IL13RA2','CXCR4','PKM','AURKA','IFRD1','HNRNPC','YWHAQ','ODC1','VBP1','KPNA2','SRM','CCT5','XRCC6','RAN','GLO1','TALDO1','GNPDA1','PSMC4','CCL2')]
select <- t(scale(t(data1[c('METTL3','METTL14','FOXO1','FOXO3','ALDOA','PGK1','PYGL','ENO1','PGAM1','HK1','PFKFB3','PGM2','PCNA','RPE','KRAS','ETF1','LDHA','CDK1','PPIA','IL13RA2','CXCR4','PKM','AURKA','IFRD1','HNRNPC','YWHAQ','ODC1','VBP1','KPNA2','SRM','CCT5','XRCC6','RAN','GLO1','TALDO1','GNPDA1','PSMC4','CCL2'),])))
select[1:4,1:4]
group <- data.frame(group=group_list)
rownames(group)=colnames(select)
pheatmap(select,show_colnames =F,show_rownames = T, 
         cluster_cols = F,cluster_rows = F,
         annotation_col=group,
         border_color = NA,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

samples <- rep(c('cirrhosis', 'normal'), c(4, 4))
heat <- Heatmap(dat, 
                col = colorRampPalette(c('navy', 'white', 'firebrick3'))(100), #定义热图由低值到高值的渐变颜色
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
                show_row_names = T,  #not to show gene names
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('control' = '#00DAE0', 'treat' = '#FF9289')),  # #定义样本分组的颜色
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 6),
                cluster_rows = F,
                cluster_columns = F)
heat

