#read data
data <- read.csv("gene_fpkm(1).csv",header = T)
colnames(data) <- c("gene_id",	"hu_con_4",	"hu_con_5",	"lung_1",	"lung_2",	"gene_name","gene_chr",	"gene_start",	"gene_end",	"gene_strand",	"gene_length",	"gene_biotype",	"gene_description",	"tf_family")
data <- filter(data,!duplicated(data$gene_name))

#mean data according to duplicated rows
data1 <- data[,2:10]
data1 <- data1 %>%
  group_by(gene_name) %>%
  summarise_all(mean)

data1 <- as.data.frame(data1)
rownames(data)=data[,6]
data <- data[,-1]
datafilt <- data[,1:4]

#fpkm to tpm
expMatrix <- datafilt
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

#differential analysis
#强制比较顺序
group_list=c(rep('Normal',2),rep('Tumor',2))
group_list <- factor(group_list,levels = c("Normal","Tumor"),ordered = F)
#data correction
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#to test if data needed transfer
exprSet <- log2(exprSet+1)
dat <- exprSet
####output######
exprSet <- as.data.frame(exprSet)
exprSet$id_ref <- rownames(exprSet)
exprSet <- exprSet[,c(5,1:4)]
write.table(exprSet,file = "exprset.txt",row.names = F)
write.csv(exprSet,file = "exprset.csv",row.names = F)
###############
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
#differential expression genes
if(T){
  logFC_t=1.5
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
#GSEA
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
#GO
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
#KEGG
gse_KEGG <- gseKEGG(genelist_sort ,
                organism     = 'hsa',
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.9,
                verbose      = FALSE)
gseaplot2(gse_KEGG, "hsa04380",pvalue_table = T)
write.csv(gse_KEGG@result,file = "gse_all_KEGG.csv")
