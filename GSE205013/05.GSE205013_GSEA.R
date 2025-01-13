library(data.table)  # 用于高效处理大数据集的库
library(dplyr)       # 数据操作和转换的库
#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)

expr<-read.csv("X12.LYPLA1_35_all_genes_with_logfc.csv")
head(expr)
geneList<-as.numeric(expr$logfoldchanges)
names(geneList)<-expr$names
head(geneList)
geneList<-sort(geneList,decreasing = T)#Named num
#准备gmt文件
geneset <- read.gmt("moffitt.Hs.symbols.gmt")  
#geneset$term = str_remove(geneset$term,"HALLMARK_")
table(geneset$term)
# gsea
egmt <- GSEA(geneList, TERM2GENE=geneset,pvalueCutoff = 1)

ridgeplot(egmt)+
  scale_y_discrete(labels=function(x) str_wrap(str_replace_all(x,"_"," "), width=30))
library(enrichplot)
n=2
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])
n=1
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])




expr<-read.csv("X12.ZDHHC7_35_all_genes_with_logfc.csv")
head(expr)
geneList<-as.numeric(expr$logfoldchanges)
names(geneList)<-expr$names
head(geneList)
geneList<-sort(geneList,decreasing = T)#Named num
#准备gmt文件
geneset <- read.gmt("moffitt.Hs.symbols.gmt")  
#geneset$term = str_remove(geneset$term,"HALLMARK_")
table(geneset$term)
# gsea
egmt <- GSEA(geneList, TERM2GENE=geneset,pvalueCutoff = 1)
ridgeplot(egmt)+
  scale_y_discrete(labels=function(x) str_wrap(str_replace_all(x,"_"," "), width=30))
library(enrichplot)
n=2
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])
n=1
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])



expr<-read.csv("X12.PTPN11_35_all_genes_with_logfc.csv")
head(expr)
geneList<-as.numeric(expr$logfoldchanges)
names(geneList)<-expr$names
head(geneList)
geneList<-sort(geneList,decreasing = T)#Named num
#准备gmt文件
geneset <- read.gmt("moffitt.Hs.symbols.gmt")  
#geneset$term = str_remove(geneset$term,"HALLMARK_")
table(geneset$term)
# gsea
egmt <- GSEA(geneList, TERM2GENE=geneset,pvalueCutoff = 1)
ridgeplot(egmt)+
  scale_y_discrete(labels=function(x) str_wrap(str_replace_all(x,"_"," "), width=30))
library(enrichplot)
n=2
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])
n=1
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])



expr<-read.csv("X12.ZDHHC20_35_all_genes_with_logfc.csv")
head(expr)
geneList<-as.numeric(expr$logfoldchanges)
names(geneList)<-expr$names
head(geneList)
geneList<-sort(geneList,decreasing = T)#Named num
#准备gmt文件
geneset <- read.gmt("moffitt.Hs.symbols.gmt")  
#geneset$term = str_remove(geneset$term,"HALLMARK_")
table(geneset$term)
# gsea
egmt <- GSEA(geneList, TERM2GENE=geneset,pvalueCutoff = 1)
ridgeplot(egmt)+
  scale_y_discrete(labels=function(x) str_wrap(str_replace_all(x,"_"," "), width=30))
library(enrichplot)
n=2
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])
n=1
gseaplot2(egmt, geneSetID = n, title = egmt$Description[n])

