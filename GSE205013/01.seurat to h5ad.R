#预备PYTHON分析
library(stringr)
setwd("F:/scRNAdata/GSE205013/GSE205013_RAW")
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste(str_split(y[1],'_',simplify = T)[,2],collapse = '')
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"features.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})
folders=list.files('./')
folders

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
setwd('F:/scRNAdata/GSE205013/')
save(sceList,file = '0.total_sceList.Rdata')

load('0.total_sceList.Rdata')
pbmc <- merge(sceList[[1]], 
              y=c(sceList[[2]], sceList[[3]],sceList[[4]], sceList[[5]],sceList[[6]], sceList[[7]],
              sceList[[8]],sceList[[9]],sceList[[10]], sceList[[11]],sceList[[12]], sceList[[13]],
              sceList[[14]],sceList[[15]],sceList[[16]], sceList[[17]],sceList[[18]], sceList[[19]],
              sceList[[20]]),
              add.cell.ids = folders,project = "untreated")
pbmc

#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)

# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
sceList=0
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)

# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total.v3.h5seurat", overwrite = TRUE)

# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total.v3.h5seurat", dest = "h5ad", overwrite = TRUE)


library(SeuratDisk)
library(Seurat)
#拆分
load('0.total_sceList.Rdata')
setwd("F:/scRNAdata/GSE205013/GSE205013_RAW")
folders=list.files('./')
setwd('F:/scRNAdata/GSE205013/')
####merge1######
folder01<-folders[1:5]
pbmc <- merge(sceList[[1]], 
              y=c(sceList[[2]], sceList[[3]],sceList[[4]], sceList[[5]]),
              add.cell.ids = folder01,project = "untreated")
pbmc
#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)
# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
sceList=0
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)
# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total1-5.v3.h5seurat", overwrite = TRUE)
# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total1-5.v3.h5seurat", dest = "h5ad", overwrite = TRUE)


####merge2######
folder02<-folders[6:10]
pbmc <- merge(sceList[[6]], 
              y=c(sceList[[7]], sceList[[8]],sceList[[9]], sceList[[10]]),
              add.cell.ids = folder02,project = "untreated")
pbmc
#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)
# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)
# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total6-10.v3.h5seurat", overwrite = TRUE)
# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total6-10.v3.h5seurat", dest = "h5ad", overwrite = TRUE)

####merge3######
folder03<-folders[11:12]
pbmc <- merge(sceList[[11]], 
              y=c(sceList[[12]]),
              add.cell.ids = folder03,project = "untreated")
pbmc
#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)
# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)
# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total11-12.v3.h5seurat", overwrite = TRUE)
# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total11-12.v3.h5seurat", dest = "h5ad", overwrite = TRUE)


####merge4######
folder04<-folders[13:17]
pbmc <- merge(sceList[[16]], 
              y=c(sceList[[17]], sceList[[14]],sceList[[15]],sceList[[13]]),
              add.cell.ids = folder04,project = "untreated")
pbmc
#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)
# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)
# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total13-17.v3.h5seurat", overwrite = TRUE)
# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total13-17.v3.h5seurat", dest = "h5ad", overwrite = TRUE)



####merge5######
folder05<-folders[18:20]
pbmc <- merge(sceList[[18]], 
              y=c(sceList[[19]], sceList[[20]]),
              add.cell.ids = folder05,project = "untreated")
pbmc
#add source
pbmc$source<-'pancreas'
pbmc$source[pbmc$orig.ident%in%c('P01','P02','P11','P16','P18','P21','P24','P25','P27')]<-'liver'
table(pbmc$orig.ident,pbmc$source)

#add stage
pbmc$stage<-'resectable'
pbmc$stage[pbmc$orig.ident%in%c('P13','P22')]<-'borderline'
pbmc$stage[pbmc$orig.ident%in%c('P16','P20','P09','P01','P25','P24','P27','P21','P26','P02','P17','P11','P18')]<-'metastatic'
table(pbmc$orig.ident,pbmc$stage)
# 设置Seurat对象的assay版本为v3
options(Seurat.object.assay.version = "v3")
#help: https://github.com/satijalab/seurat/issues/8304
DefaultAssay(pbmc) <- "RNA"
##Join layers
New<- JoinLayers(pbmc, overwrite = TRUE)
# 将Seurat对象转换为Seurat v3的assay结构
object_v3 <- CreateSeuratObject(
  counts = GetAssayData(New, assay="RNA", slot = "counts"),
  meta.data = New@meta.data
)
# 保存转换后的Seurat对象为.h5seurat格式
SaveH5Seurat(object_v3, filename="pbmc_total18-20.v3.h5seurat", overwrite = TRUE)
# 将.h5seurat文件转换为h5ad格式
Convert("pbmc_total18-20.v3.h5seurat", dest = "h5ad", overwrite = TRUE)
