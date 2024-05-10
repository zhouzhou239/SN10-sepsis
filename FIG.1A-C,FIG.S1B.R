#载入程序包
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(DOSE)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(harmony)
library(pheatmap)
library(hdf5r)
library(scater)
library(ggsci)
library(dplyr)
library(Matrix)
library(cowplot)
if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
setwd("D:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP")
##1 单细胞数据读入
#方法一
# 读入数据,使用目录向量合并
obj_dir <- c("A1","A2","A3","A4","A5","A6","R1","R2","R3","R4","R5","R6")#读入样本文件的相对路径
names(obj_dir) = c("A1","A2","A3","A4","A5","A6","R1","R2","R3","R4","R5","R6")#设置样本名，用于后续在seurat对象中区分样本,filtered h5 文件用  Read10X_h5() 函数读取， 3个文件 用  Read10X() 函数读取,分别创建 seurat 对象后，merge seurat 对象就可以了.

counts <- Read10X(data.dir = obj_dir)#读入表达矩阵
scRNA0 = CreateSeuratObject(counts, min.cells=1)#依据表达矩阵创建seurat对象

counts1 <-Read10X_h5("HC1.h5")
counts2 <-Read10X_h5("HC2.h5")
counts3 <-Read10X_h5("HC3.h5")
scRNA2 = CreateSeuratObject(counts1, project = "HC1", min.cells=1)#依据表达矩阵创建seurat对象
scRNA3 = CreateSeuratObject(counts2, project = "HC2", min.cells=1)#依据表达矩阵创建seurat对象
scRNA4 = CreateSeuratObject(counts3, project = "HC3", min.cells=1)#依据表达矩阵创建seurat对象
#合并不同来源的数据
scRNA1 <-merge(scRNA0,c(scRNA2,scRNA3,scRNA4))
dim(scRNA1)
table(scRNA1@meta.data$orig.ident)
#保存读入的rds
saveRDS(scRNA1, "Step1_merged_sample.rds")
#生成seurat对象list
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")

##2 单细胞数据质控
#本项目中的质控标准为：保留基因数大于 200 、UMI数大于 1000、log10GenesPerUMI大于 0.7、
#线粒体UMI占比低于 20%、血红细胞基因占比低于 5%的细胞作为高质量细胞
#然后使用DoubletFinder软件进行双细胞去除，进行下游分析。
 ###2.1 标准质控
#质控
for (i in 1:15) 
{
  j = i+1 
  if (i  == 1 ) {
    tem <- merge(scRNAlist[[i]], y=scRNAlist[[j]])
  }
  if (i == 15) {
    print("MIssion finished")
  }
  if (i != 1  && i !=15 ) {
    tem <- merge(tem,scRNAlist[[j]])
    
  }

}#15个样本依次整合
scRNA <- tem
Idents(scRNA) <- 'orig.ident'
table(Idents(scRNA))#查看不同orig.ident的细胞数

# 将每个细胞每个UMI的基因数目添加到元数据中
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)
# 计算线粒体比率
scRNA$mitoRatio <- PercentageFeatureSet(object =scRNA, pattern = "^MT-")#因为人的线粒体基因是MT-开头的，所以可以通过这种方式匹配，如果是小鼠的话，是以“mt-”开头的。其他特殊物种可以直接提供线粒体基因的列表
scRNA$mitoRatio <- scRNA@meta.data$mitoRatio / 100
#计算红细胞基因比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA$percentHB<-PercentageFeatureSet(scRNA, features=HB.genes) #此处就是根据提供的红细胞基因列表计算红细胞基因表达的比例

# 创建元数据数据框
metadata <- scRNA@meta.data
# 为元数据添加细胞ID
metadata$cells <- rownames(metadata)
metadata$sample<-metadata$orig.ident
# 重命名列
metadata <- metadata %>%
        dplyr::rename(nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

# 将元数据添加回Seurat对象中
scRNA@meta.data <- metadata
                           
# 任何时候都要创建.RData对象保存进度
dir.create("data")
save(scRNA, file="data/merged_filtered_seurat_scRNA.RData")

# 使用选择的阈值筛掉低质量读写 - 
filtered_scRNA <- subset(x = scRNA, 
                         subset= (nUMI >= 1000) & 
                           (nGene >= 200) & 
                           (log10GenesPerUMI >= 0.70) & 
                           (mitoRatio <= 0.20)& (percentHB <=0.05)) 
table(Idents(filtered_scRNA))#查看不同orig.ident的细胞数

# 提取计数
counts <- GetAssayData(object = filtered_scRNA, slot = "counts")

# 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
nonzero <- counts > 0

# 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
keep_genes <- Matrix::rowSums(nonzero) >= 10

# 仅保留那些在10个以上细胞中表达的基因
filtered_counts <- counts[keep_genes, ]

# 重新赋值给经过过滤的Seurat对象
filtered_scRNA <- CreateSeuratObject(filtered_counts, meta.data = filtered_scRNA@meta.data)

save(filtered_scRNA, file="data/filtered_scRNA.RData")

###2.2 去除双细胞
# remove doublets
load("data/filtered_scRNA.RData")
scRNA1<-filtered_scRNA
obj = SplitObject(scRNA1,split.by = "orig.ident")
obj_rm=list()#创建空list用于存放去除双细胞以后的seurat对象
doublets_plot = list()#创建空list用于存放双细胞的分布情况
pc.num = 1:30#设置的经验PC维数
#dir.create("SingleCell_QC")
RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
  ){
    ## 寻找最优pK值
    sweep.res.list <- paramSweep_v3(object, PCs = pc.num, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
    bcmvn <- find.pK(sweep.stats)#求出最大bcmvn值
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()#最大bcmvn值对应的pk值为最优
    ## 排除不能检出的同源doublets，优化期望的doublets数量
    homotypic.prop <- modelHomotypic(object$seurat_clusters)   # 最好提供celltype，因为cluster和celltype之间存在多对一的情况
    nExp_poi <- round(doublet.rate*ncol(object)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    seu.scored <- doubletFinder_v3(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    # 选出储存双细胞预测结果的列
    cname <-colnames(seu.scored[[]])
    DF<-cname[grep('^DF',cname)]
    seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
    
    # 去除双细胞
    seu.removed <- subset(seu.scored, subset = doublet != 1)
    p1 <- DimPlot(seu.scored, group.by = DF)
    res.list <- list("plot"=p1, "obj"=seu.removed)
    return(res.list)
  }
#对每一个样本进行标准聚类的操作
for( i in names(obj)){
    obj[[i]] <- NormalizeData(obj[[i]])
    obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
    obj[[i]] <- ScaleData(obj[[i]])
    obj[[i]] <- RunPCA(obj[[i]])
    obj[[i]] <- RunUMAP(obj[[i]], dims = 1:20)
    obj[[i]] <- FindNeighbors(obj[[i]], dims = pc.num) %>% FindClusters(resolution = 0.5)#此处选择了一个较小的resolution，避免over-clustering的情况
    tmp <- RemoveDoublets(obj[[i]], doublet.rate=ncol(obj[[i]])*8*1e-6,pc.num=pc.num)# 1000细胞对应的doublets rate是0.8%， 如果细胞数大于1000，可以根据每1000个细胞增加0.8%的方式进行估算和替换
    obj_rm[[i]] <- tmp$obj
    doublets_plot[[i]] <- tmp$plot
  }

for (i in 1:15) 
{
  j = i+1 
  if (i  == 1 ) {
    tem <- merge(obj_rm[[i]], y=obj_rm[[j]])
  }
  if (i == 15) {
    print("MIssion finished")
  }
  if (i != 1  && i !=15 ) {
    tem <- merge(tem,obj_rm[[j]])
    
  }

}
scRNA <- tem
Idents(scRNA) <- 'orig.ident'
saveRDS(scRNA,"data/Step2_After_QC.rds")

## 3 数据归一化与标准化
#数据归一化与标准化
scRNA <- NormalizeData(scRNA)#归一化
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst")#选择高变基因
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))#标准化

#再次检查seurat对象的结构
str(scRNA)

##4 数据批次矫正与聚类
#查看pca结果
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))#pca降维
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")#每个样本在pca降维结果中的分布情况
plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") #选取前30维绘制elbowplot
plot3 <- plot1+plot2
dir.create("Clustering3")
ggsave("Clustering/pca2.png", plot = plot3, width = 8, height = 4)
plot3
saveRDS(scRNA,"data/Step2_After_SCALE.rds")
##细胞聚类
pc.num=1:5
scRNA <- FindNeighbors(scRNA, dims = pc.num) #KNN + SNN
scRNA <- FindClusters(scRNA, resolution = 0.3)#Louvain
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)#将聚类结果另存为data.frame
write.csv(cell_cluster,'Clustering3/cell_cluster12.csv',row.names = F, quote = F)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)#UMAP二次降维
#group_by_cluster
plot3 = DimPlot(scRNA, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
ggsave("Clustering3/UMAP_PC5.png",plot3)
plot3

#查看每个样本在UMAP降维图中的分布
plot4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')

#查看每个样本中的cluster组成情况
plot5 = DimPlot(scRNA, reduction = "umap", split.by='orig.ident')
ggsave("Clustering3/UMAP_cluster12_sample1.png", plot = plot4, width = 6, height = 5)
ggsave("Clustering3/UMAP_cluster12_sample.png", plot = plot5, width = 25, height = 5)
saveRDS(scRNA,"Step3_Clustering12.rds")

###4.1 批次效应矫正（MNN）
#去批次后进行聚类
scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")#将合并的seurat对象重新拆分
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))#对每一个样本进行归一化
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))#对每一个样本寻找高变基因
scRNA_mnn <- RunFastMNN(object.list = scRNAlist)#fastMNN去批次
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)#对去完批次以及合并后的对象重新进行高变基因的选取
scRNA_mnn1 <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = 1:12)#UMAP降维
scRNA_mnn1 <- FindNeighbors(scRNA_mnn1, reduction = "mnn", dims = 1:12)#SNN + KNN
scRNA_mnn1 <- FindClusters(scRNA_mnn1, resolution = 0.3)#Louvain
scRNA_mnn1@meta.data$group ="pat,HC"
scRNA_mnn1@meta.data[which(scRNA_mnn1@meta.data$orig.ident %in% c("HC1","HC2","HC3") ), "group"]="HC"
scRNA_mnn1@meta.data[which(scRNA_mnn1@meta.data$orig.ident %in% c("A1","A3","A2","A4","A5","A6","R1","R3","R2","R4","R5","R6") ), "group"]="pat"
p1 <- DimPlot(scRNA_mnn, group.by = "orig.ident", pt.size=0.1) + 
  ggtitle("Integrated by fastMNN")#去批次后的可视化
p2 <- DimPlot(scRNA, group.by="orig.ident", pt.size=0.1) + 
  ggtitle("No integrated")#去批次前的可视化
p = p1 + p2 + plot_layout(guides='collect')
#比较去批次前后聚类结果的变化
p
saveRDS(scRNA_mnn1,"Step4_MNN_cluster12.rds")

p3 <- DimPlot(scRNA_mnn1,  pt.size=0.1) #去批次后重降维聚类结果的可视化
p4 <- DimPlot(scRNA_mnn1, split.by="group", pt.size=0.1) #去批次后重降维聚类结果的分样本可视化展示
#p5 = p3 + p4 + plot_layout(guides='collect') #合并两者结果，将legend统一放到右边
#比较各个样本内cluster的组成情况
#p5
ggsave("Clustering3/fastMNN1.pdf", plot = p3, width = 6, height =5)
ggsave("Clustering3/fastMNN2.pdf", plot = p4, width = 11, height =5)

##5 marker基因鉴定
scRNA3<-scRNA_mnn1
#鉴定 marker 基因
dir.create("Marker")
# 鉴定各 cluster 的 Marker 基因，选定的cluster为实验组，剩余的所有cluster为对照组
all.markers = FindAllMarkers(scRNA3, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25, 
                                     only.pos = TRUE)
head(all.markers)
#簇中基因平均表达相对于所有其他簇中平均表达的最小log2倍变化。默认值为 0.25。
# 根据 FoldChange 进行排序选取每群细胞的 top10 Marker 基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#此处的n可以替换为任意数字，如果只想看top1，就把10替换成1

# 保存 Marker 基因
write.table(all.markers, 
            "Marker/all_Markers_of_each_clusters_pc8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Marker/top10_Markers_of_each_clusters_pc8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
# 不同组之间进行比较
#先对样本进行分组
scRNA3@meta.data$group ="AA1,AA2,RR1,RR2,HC"## 增加一列标签 
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A1","A3") ), "group"]="AA1" ## A1 ,A3代表2个样本名。
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A2","A4","A5","A6") ), "group"]="AA2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R1","R3") ), "group"]="RR1"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R2","R4","R5","R6") ), "group"]="RR2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("HC1","HC2","HC3") ), "group"]="HC"

scRNA3@meta.data$group3 ="AA1,AA2,RR1,RR2,HC"## 增加一列标签 
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A1","A2","A3") ), "group3"]="AA1" ## A1 ,A3代表2个样本名。
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A4","A5","A6") ), "group3"]="AA2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R1","R2","R3") ), "group3"]="RR1"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R4","R5","R6") ), "group3"]="RR2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("HC1","HC2","HC3") ), "group3"]="HC"
## 默认识别的信息改成 组信息
p5<-DimPlot(scRNA3,group.by = 'seurat_clusters', split.by = "group3")
ggsave('Integrate2_umap3.png', p5, width=25, height=4)
Idents(scRNA3) <- 'group3'
AA.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = c("AA1", "AA2"), #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(AA.HC.markers, 
            "Marker/AA.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = c("RR1", "RR2"), #这里可以替换成选择的实验组
                                ident.2 ="HC" , #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR.HC.markers, 
            "Marker/RR.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
AA1.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = "AA1", #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(AA1.HC.markers, 
            "Marker/AA1.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
AA2.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = "AA2", #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(AA2.HC.markers, 
            "Marker/AA2.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR1.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = "RR1", #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR1.HC.markers, 
            "Marker/RR1.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR2.HC.markers <- FindMarkers(scRNA3, 
                                ident.1 = "RR2", #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR2.HC.markers, 
            "Marker/RR2.HC.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR1.AA1.markers <- FindMarkers(scRNA3, 
                                ident.1 = "RR1", #这里可以替换成选择的实验组
                                ident.2 = "AA1", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR1.AA1.markers, 
            "Marker/RR1.AA1.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR2.AA2.markers <- FindMarkers(scRNA3, 
                                ident.1 = "RR2", #这里可以替换成选择的实验组
                                ident.2 = "AA2", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR2.AA2.markers, 
            "Marker/RR2.AA2.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
AA2.AA1.markers <- FindMarkers(scRNA3, 
                                ident.1 = "AA2", #这里可以替换成选择的实验组
                                ident.2 = "AA1", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(AA2.AA1.markers, 
            "Marker/AA2.AA1.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
RR2.RR1.markers <- FindMarkers(scRNA3, 
                                ident.1 = "RR2", #这里可以替换成选择的实验组
                                ident.2 = "RR1", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(RR2.RR1.markers, 
            "Marker/RR2.RR1.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")

# 绘制热图
scRNA3 <- ScaleData(scRNA3, features = row.names(scRNA3))#对所有的基因进行标准化，因为差异基因可能出现在高变基因以外
#saveRDS(scRNA3,"Step5_marker11_MNN.rds")
heatmap_plot = DoHeatmap(object = scRNA3, 
                         features = as.vector(top10$gene), #对每个cluster的top10 marker进行可视化
                         group.by = "celltype", #按组进行分列
                         group.bar = T, #绘制colorbar提示cluster的位置
                         size = 2) +
  theme(axis.text.y = element_text(size = 4))

  # 保存绘图结果
ggsave("Marker/top10_marker_of_each_celltype_heatmap_PC8.pdf", width = 12, height = 12,
       plot = heatmap_plot)
ggsave("Marker/top10_marker_of_each_celltype_heatmap_PC8.png", width = 12, height = 12,
       plot = heatmap_plot)
write.table(all.markers, 
            "Marker/all_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Marker/top10_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
##6 差异基因鉴定
#将数据集的idents从clusters切换回orig.ident(各样本)
#Idents(scRNA3) <- 'orig.ident'（只有两个样本时，无需分组）
#先对样本进行分组
scRNA3@meta.data$group ="AA1,AA2,RR1,RR2,HC"## 增加一列标签 
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A1","A3") ), "group"]="AA1" ## A1 ,A3代表2个样本名。
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A2","A4","A5","A6") ), "group"]="AA2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R1","R3") ), "group"]="RR1"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R2","R4","R5","R6") ), "group"]="RR2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("HC1","HC2","HC3") ), "group"]="HC"
## 默认识别的信息改成 组信息
Idents(scRNA3) <- 'group'
# 鉴定各样本的差异基因
Diff_exp = FindAllMarkers(scRNA3, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
head(Diff_exp)
top10 = Diff_exp %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#选择每个样本前10个差异基因进行可视化
dir.create("Diffexp3")
write.table(Diff_exp, 
            "Diffexp3/group_diff_group_PC8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Diffexp3/group_top10_diff_group_PC8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

heatmap_plot = DoHeatmap(object = scRNA3, 
                         features = as.vector(top10$gene), 
                         group.by = "group3", #按样本进行分列
                         group.bar = T, 
                         size = 2) +
  theme(axis.text.y = element_text(size = 4))
  # 保存绘图结果
ggsave("Diffexp3/top10_marker_of_each_group_heatmap_PC8.pdf",
       plot = heatmap_plot)
ggsave("Diffexp3/top10_marker_of_each_group_heatmap_PC8.png", 
       plot = heatmap_plot)
heatmap_plot
##7 富集分析
#富集分析
library(ReactomePA)
genes_symbol <- as.character(Diff_exp$gene)
#转换ID
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
dir.create("enrichment")
#GO
ego <- enrichGO(gene = id,#使用ENTREZID进行富集分析
                OrgDb = org.Hs.eg.db,#选择人的数据库进行富集分析
                ont = "BP",#选择其中实现的分子功能条目进行富集分析，#subontology这里选"MF“,也可以分别选"BP","CC"或"ALL"；
                pAdjustMethod = "BH",#选择BH作为p值矫正的方法
                pvalueCutoff = 0.05,#p值的阈值
                qvalueCutoff = 0.05,#q值的阈值
                readable = TRUE)
test<-ego[ego@result[["ID"]]]
write.table(test, 
            "enrichment/group3_BP.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
#选择指定条目的进行后续分析
nfkb_pathway <-  ego[ego@result[["ID"]]  %in% c("GO:0051092","GO:0043122", "GO:0038061", "GO:0007249", 
                                                 "GO:0043123","GO:1901224","GO:1901222") ]

Top_ten_nfkb_pathway <-  ego[ego@result[["ID"]]  %in% c("GO:0002181",
"GO:0006119",
"GO:0009060",
"GO:0045333",
"GO:0046034",
"GO:0016032",
"GO:0019058",
"GO:0019646",
"GO:0015980",
"GO:0001819",
"GO:0051092","GO:0043122", "GO:0038061", "GO:0007249", "GO:0043123","GO:1901224","GO:1901222") ]
#或者nfkb_pathway <- subset(ego,ego@result[["ID"]]  %in% c("GO:0043122", "GO:0051092", "GO:0007249", "GO:0043123"))
ego@result <- Top_ten_nfkb_pathway

GO_dot <- dotplot(ego,showCategory = 17)
GO_bar <- barplot(ego,showCategory = 17)
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("enrichment/GO_results_pc8.pdf", plot=res_plot, width = 14,height = 12)
ggsave("enrichment/GO_results_pc8.png", plot=res_plot, width = 14,height = 12)

GO_dot

GO_bar

ego_all <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",#对所有三个分类进行富集分析
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
#替换某特定的值
ego1<-ego_all
test<-ego1[ego1@result[["ID"]]]
test[175,test[175,]=="BP"]<-"NF kappaB"
test[228,test[228,]=="BP"]<-"NF kappaB"
test[241,test[241,]=="BP"]<-"NF kappaB"
test[279,test[279,]=="BP"]<-"NF kappaB"
test[317,test[317,]=="BP"]<-"NF kappaB"
test[476,test[476,]=="BP"]<-"NF kappaB"
test[504,test[504,]=="BP"]<-"NF kappaB"
ego_all@result <- test
#画图
GO_dot <-dotplot(ego_all,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") #对三个分类进行分割展示

GO_bar <-barplot(ego_all,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
ggsave("enrichment3/GO_results_dot.pdf", plot=GO_dot, width = 10,height = 17)
ggsave("enrichment3/GO_results_dot.png", plot=GO_dot, width = 10,height = 17)
ggsave("enrichment3/GO_results_BAR.pdf", plot=GO_bar, width = 10,height =17)
ggsave("enrichment3/GO_results_BAR.png", plot=GO_bar, width = 10,height = 17)

##8 细胞类型鉴定
dir.create("SingleR3")
Idents(scRNA3) <- 'seurat_clusters'
# 读入参考数据集
load("HumanPrimaryCellAtlas_hpca.se_human.RData")
library(SingleR)
library(celldex)

norm_count = GetAssayData(scRNA3, 
                          layer="data") #归一化后的矩阵slot="data"
pred<- SingleR(test = norm_count, 
               ref = hpca.se, 
               labels = hpca.se$label.main)
#将注释结果添加到seurat对象里
scRNA3$singleR_celltype = pred$labels

# 对细胞类型注释结果进行可视化
celltype_plot = DimPlot(scRNA3, 
                        reduction = "umap", 
                        group.by = "singleR_celltype")
#ggsave("SingleR/celltype_plot.pdf", celltype_plot)
#ggsave("SingleR/celltype_plot.png", celltype_plot)
celltype_plot

# 通过热图了解每个细胞对应细胞类型的打分情况
score_heatmap <-  plotScoreHeatmap(pred, clusters = scRNA3@meta.data$seurat_clusters, order.by = "clusters")
ggsave("SingleR3/score_heatmap.pdf", score_heatmap, width = 8, height = 10)
ggsave("SingleR3/score_heatmap.png", score_heatmap, width = 8, height = 10)
score_heatmap

# 通过热图了解每个cluster对应细胞类型的打分情况
tab <- table(cluster=scRNA3@meta.data$seurat_clusters, label=pred$labels) 
cluster_type <- pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 
ggsave("SingleR3/cluster_type.pdf", cluster_type)
ggsave("SingleR3/cluster_type.png", cluster_type)
cluster_type

#手动注释
#select_genes_Fibroblast<- c('ACTA2')
#select_genes_Endothelium<- c('PTPRB','PECAM1')
select_genes_Epithelium<- c('TPPP3','KRT18','KRT5','KRT4','KRT13','CAPS','SCGB1A1','SCGB3A1','CRCT1','SPRR3','ECM1','ASS1')
#select_genes_Mast<- c('KIT','TPSAB1')
select_genes_Neutrophil<- c('HCAR3','PROK2', 'PI3')
select_genes_Dendritic<- c('CD1C')
select_genes_pDendritic<- c('IL3RA','GZMB','SERPINF1','ITM2C')
select_genes_Monocyte<- c('FCN1','CD300E','TNIP3')
select_genes_Macrophage<- c('CD68','FABP4','CFD','RBP4','CST3','C1QB','C1QA')
#select_genes_Bcell<- c('MS4A1','CD79A','CD79B','FCRL5')
#select_genes_Plasma<- c('MZB1','JCHAIN')
#select_genes_Proliferative_signal<- c('MKI67','TOP2A','STMN1')
#select_genes_NK_NKT<- c('FGFBP2','GNLY','KLRD1')
select_genes_Tcell<- c('CD3D','CD3E','KLRD1')
#select_genes_NK_Tcell<- c('CD3D','CD3E','FGFBP2','GNLY','KLRD1')

#p1 <-DotPlot(scRNA3, features = select_genes_Fibroblast ) + coord_flip()
#p2 <-DotPlot(scRNA3, features = select_genes_Endothelium ) + coord_flip()
p3 <-DotPlot(scRNA3, features = select_genes_Epithelium ) + coord_flip()
p4 <-DotPlot(scRNA3, features = select_genes_Macrophage ) + coord_flip()
#p5 <-DotPlot(scRNA3, features = select_genes_Mast) + coord_flip()
p6 <-DotPlot(scRNA3, features = select_genes_Neutrophil ) + coord_flip()
p7 <-DotPlot(scRNA3, features = select_genes_Dendritic ) + coord_flip()
p8 <-DotPlot(scRNA3, features = select_genes_pDendritic ) + coord_flip()
p9 <-DotPlot(scRNA3, features = select_genes_Monocyte) + coord_flip()
#p10 <-DotPlot(scRNA3, features = select_genes_Bcell ) + coord_flip()
#p11 <-DotPlot(scRNA3, features = select_genes_Plasma ) + coord_flip()
#p12 <-DotPlot(scRNA3, features = select_genes_Proliferative_signal ) + coord_flip()
#p13 <-DotPlot(scRNA3, features = select_genes_NK_Tcell ) + coord_flip()
p14 <-DotPlot(scRNA3, features = select_genes_Tcell ) + coord_flip()
#p15 <-DotPlot(scRNA3, features = select_genes_NK_NKT ) + coord_flip()

#ggsave("SingleR2/select_genes_Fibroblast.png", p1, width = 8, height = 10)
#ggsave("SingleR2/select_genes_Endothelium.png", p2, width = 8, height = 10)
ggsave("SingleR3/select_genes_Epithelium .png", p3, width = 8, height = 10)
ggsave("SingleR3/select_genes_Macrophage.png", p4, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Mast.png", p5, width = 8, height = 10)
ggsave("SingleR3/select_genes_Neutrophil.png", p6, width = 8, height = 10)
#ggsave("SingleR3/select_genes_cDendritic.png", p7, width = 8, height = 10)
#ggsave("SingleR3/select_genes_pDendritic.png", p8, width = 8, height = 10)
ggsave("SingleR3/select_genes_Monocyte.png", p9, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Bcell.png", p10, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Plasma.png", p11, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Proliferative_signal.png", p12, width = 8, height = 10)
#ggsave("SingleR3/select_genes_NK_Tell.png", p13, width = 8, height = 10)
ggsave("SingleR3/select_genes_Tcell.png", p14, width = 8, height = 10)
#ggsave("SingleR3/select_genes_NK_NKT.png", p14, width = 8, height = 10)


# 将singleR得到的细胞类型加入到对象的meta.data中，后续用于高级分析
scRNA3<- RenameIdents(scRNA3, '0'='Macrophage','1'='Neutrophils','2'='Neutrophils','3'='Epithelial_cells','4'='Monocyte','5'='Macrophage','6'='Macrophage','7'='Epithelial_cells','8'='T_cells','9'='Macrophage','10'='Epithelial_cells')
scRNA3[["celltype"]] <- Idents(scRNA3)
head(scRNA3@meta.data)
saveRDS(scRNA3, "Step5_CellTyping.rds")
#scRNA3<- RenameIdents(scRNA3, 'Macrophage'='myeloid_cells','Neutrophils'='myeloid_cells','Monocyte'='myeloid_cells','Epithelial_cells'='Epithelial_cells','T_cells'='T_cells')


# 对细胞类型注释结果进行可视化
celltype_plot = DimPlot(scRNA3, 
                        reduction = "umap", 
                        group.by = "celltype")
ggsave("SingleR3/celltype_plot_z.pdf", celltype_plot)
ggsave("SingleR3/celltype_plot_z.png", celltype_plot)
celltype_plot

p1 <- DimPlot(scRNA3, group.by = "seurat_clusters")
p2 <- DimPlot(scRNA3, group.by = "celltype2")
p3 <- DimPlot(scRNA3, group.by = "celltype")
#查看cluster与细胞类型的对应关系
p4 = p1+p3+p2
ggsave("SingleR3/celltype_contrast_z.pdf",p4,width = 18, height = 6)
ggsave("SingleR3/celltype_contrast_z.png",p4, width = 18, height = 6)
p4

# 通过热图了解每个cluster对应细胞类型的打分情况
tab <- table(cluster=scRNA3@meta.data$seurat_clusters, label=scRNA3@meta.data$celltype) 
cluster_type <- pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 
ggsave("SingleR3/cluster_type_z.pdf", cluster_type)
ggsave("SingleR3/cluster_type_z.png", cluster_type)
cluster_type