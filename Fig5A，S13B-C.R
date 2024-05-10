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
library(org.Mm.eg.db)
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

setwd("E:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP/data/mouse/mouse")
##1 单细胞数据读入
scRNA<-readRDS("Step6_Ctrl_M_IFX_SN10.rds")
scRNA= UpdateSeuratObject(object = scRNA)
## 3 数据归一化与标准化
#数据归一化与标准化
scRNA <- NormalizeData(scRNA)#归一化
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst")#选择高变基因
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))#标准化

#再次检查seurat对象的结构
str(scRNA)
dir.create("Clustering")
#查看pca结果
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))#pca降维
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")#每个样本在pca降维结果中的分布情况
plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") #选取前30维绘制elbowplot
plot3 <- plot1+plot2
ggsave("Clustering/pca2.png", plot = plot3, width = 8, height = 4)
plot3
Idents(scRNA) <- 'seurat_clusters'
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)#将聚类结果另存为data.frame
write.csv(cell_cluster,'Clustering/cell_cluster1.csv',row.names = F, quote = F)
#group_by_cluster
plot3 = DimPlot(scRNA, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
ggsave("Clustering/UMAP_PC19.pdf",plot3)
plot3

#查看每个样本在UMAP降维图中的分布
plot4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')

#查看每个样本中的cluster组成情况
plot5 = DimPlot(scRNA, reduction = "umap", split.by='orig.ident')
ggsave("Clustering/UMAP_cluster19_sample1.pdf", plot = plot4, width = 6, height = 5)
ggsave("Clustering/UMAP_cluster19_sample.pdf", plot = plot5, width = 25, height = 5)

##5 marker基因鉴定
scRNA3<-scRNA
#先对样本进行分组
scRNA3@meta.data$group ="Ctrl,MOD,IFX,SN10,SN61"## 增加一列标签 
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("Ctrl1","Ctrl2") ), "group"]="Ctrl" ## A1 ,A3代表2个样本名。
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("MOD1","MOD2") ), "group"]="M"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("IFX1","IFX2") ), "group"]="IFX"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("SN101","SN102") ), "group"]="SN10"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("SN611","SN612") ), "group"]="SN61"

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
## 默认识别的信息改成 组信息
p5<-DimPlot(scRNA3,group.by = 'seurat_clusters', split.by = "group")
ggsave('Clustering/Integrate3_UMAP.pdf', p5, width=25, height=4)
Idents(scRNA3) <- 'group'
M.Ctrl.markers <- FindMarkers(scRNA3, 
                                ident.1 = "M", #这里可以替换成选择的实验组
                                ident.2 = "Ctrl", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(M.Ctrl.markers, 
            "Marker/M.Ctrl.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
IFX.Ctrl.markers <- FindMarkers(scRNA3, 
                                ident.1 = "IFX", #这里可以替换成选择的实验组
                                ident.2 =  "Ctrl", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(IFX.Ctrl.markers, 
            "Marker/IFX.Ctrl.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
SN10.Ctrl.markers <- FindMarkers(scRNA3, 
                                ident.1 = "SN10", #这里可以替换成选择的实验组
                                ident.2 = "Ctrl", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(SN10.Ctrl.markers, 
            "Marker/SN10.Ctrl.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
SN61.Ctrl.markers <- FindMarkers(scRNA3, 
                                ident.1 = "SN61", #这里可以替换成选择的实验组
                                ident.2 = "Ctrl", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(SN61.Ctrl.markers, 
            "Marker/SN61.Ctrl.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
M.SN61.markers <- FindMarkers(scRNA3, 
                                ident.1 = "M", #这里可以替换成选择的实验组
                                ident.2 = "SN61", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(M.SN61.markers, 
            "Marker/M.SN61.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
M.SN10.markers <- FindMarkers(scRNA3, 
                                ident.1 = "M", #这里可以替换成选择的实验组
                                ident.2 = "SN10", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(M.SN10.markers, 
            "Marker/M.SN10.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
M.IFX.markers <- FindMarkers(scRNA3, 
                                ident.1 = "M", #这里可以替换成选择的实验组
                                ident.2 = "IFX", #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(M.IFX.markers, 
            "Marker/M.IFX.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
SN10.IFX.markers <- FindMarkers(scRNA3, 
                                ident.1 = "SN10", #这里可以替换成选择的实验组
                                ident.2 ="IFX" , #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(SN10.IFX.markers, 
            "Marker/SN10.IFX.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
SN61.IFX.markers <- FindMarkers(scRNA3, 
                                ident.1 = "SN61", #这里可以替换成选择的实验组
                                ident.2 ="IFX" , #这里可以替换成选择的对照组
                                min.pct = 0.25)
write.table(SN61.IFX.markers, 
            "Marker/SN61.IFX.markers.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")

# 绘制热图
scRNA3 <- ScaleData(scRNA3,features = row.names(scRNA3))#对所有的基因进行标准化，因为差异基因可能出现在高变基因以外
heatmap_plot = DoHeatmap(object = scRNA3, 
                         features = as.vector(top10$gene), #对每个cluster的top10 marker进行可视化
                         group.by = "seurat_clusters", #按组进行分列
                         group.bar = T, #绘制colorbar提示cluster的位置
                         size = 2) +
  theme(axis.text.y = element_text(size = 4))

  # 保存绘图结果
ggsave("Marker/top10_marker_of_each_cluster20_heatmap.pdf", width = 12, height = 12,
       plot = heatmap_plot)
ggsave("Marker/top10_marker_of_each_cluster20_heatmap.png", width = 12, height = 12,
       plot = heatmap_plot)

##6 差异基因鉴定
#将数据集的idents从clusters切换回orig.ident(各样本)
#Idents(scRNA3) <- 'orig.ident'（只有两个样本时，无需分组）
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
            "Diffexp3/group11_diff_PC8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Diffexp3/group11_top10-diff_PC8.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

heatmap_plot = DoHeatmap(object = scRNA3, 
                         features = as.vector(top10$gene), 
                         group.by = "group", #按样本进行分列
                         group.bar = T, 
                         size = 2) +
  theme(axis.text.y = element_text(size = 4))
  # 保存绘图结果
ggsave("Diffexp3/top10_marker_of_each_cluster20_heatmap_PC8.pdf",
       plot = heatmap_plot)
ggsave("Diffexp3/top10_marker_of_each_cluster20_heatmap_PC8.png", 
       plot = heatmap_plot)
heatmap_plot
##7 富集分析
#富集分析
library(ReactomePA)
genes_symbol <- as.character(Diff_exp$gene)
#转换ID
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
id = as.character(eg[,2])
dir.create("enrichment3")
#GO
ego <- enrichGO(gene = id,#使用ENTREZID进行富集分析
                OrgDb = org.Mm.eg.db,#选择人的数据库进行富集分析
                ont = "BP",#选择其中实现的分子功能条目进行富集分析，#subontology这里选"MF“,也可以分别选"BP","CC"或"ALL"；
                pAdjustMethod = "BH",#选择BH作为p值矫正的方法
                pvalueCutoff = 0.05,#p值的阈值
                qvalueCutoff = 0.05,#q值的阈值
                readable = TRUE)
test<-ego[ego@result[["ID"]]]
write.table(test, 
            "enrichment3/group15_BP.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

GO_dot <- dotplot(ego,showCategory = 10)
GO_bar <- barplot(ego,showCategory = 10)
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("enrichment3/GO_results_BP_TOP10.pdf", plot=res_plot, width = 14,height = 12)
ggsave("enrichment3/GO_results_BP_TOP10.png", plot=res_plot, width = 14,height = 12)

GO_dot

GO_bar

ego_all <- enrichGO(gene = id,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",#对所有三个分类进行富集分析
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
test<-ego_all[ego_all@result[["ID"]]]
write.table(test, 
            "enrichment3/group15_all.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
#画图
GO_dot <-dotplot(ego_all,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") #对三个分类进行分割展示
GO_bar <-barplot(ego_all,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
ggsave("enrichment3/GO_results_dot.pdf", plot=GO_dot, width = 10,height = 20)
ggsave("enrichment3/GO_results_dot.png", plot=GO_dot, width = 10,height = 20)
ggsave("enrichment3/GO_results_BAR.pdf", plot=GO_bar, width = 10,height =20)
ggsave("enrichment3/GO_results_BAR.png", plot=GO_bar, width = 10,height = 20)


#选择指定条目的进行后续分析
ego<-ego_all
nfkb_pathway <-  ego[ego@result[["ID"]]  %in% c("GO:0007249",
"GO:0043122",
"GO:0051092",
"GO:0043123",
"GO:0038061",
"GO:1901222","GO:1901223",
"GO:1901224",
"GO:0043124",
"GO:0032088",
"GO:0007252",
"GO:0051059")]
TNF_pathway<-  ego[ego@result[["ID"]]  %in% c(
"GO:1903555","GO:1903556",
"GO:0032680","GO:0032640",
"GO:1903557",
"GO:0032760",
"GO:0034612",
"GO:0071356","GO:0071706",
"GO:0032720","GO:0033209",
"GO:1903556",
"GO:0010803",
"GO:1903265")]

#或者nfkb_pathway <- subset(ego,ego@result[["ID"]]  %in% c("GO:0043122", "GO:0051092", "GO:0007249", "GO:0043123"))
ego1<-ego
ego2<-ego
ego1@result <- nfkb_pathway
ego2@result <- TNF_pathway
GO_dot1 <- dotplot(ego1,showCategory = 17)
GO_bar1 <- barplot(ego1,showCategory = 17)
res_plot1 <- CombinePlots(list(GO_dot1,GO_bar1), nrow=1)
ggsave("enrichment3/GO_results_nfkb_pathway.pdf", plot=res_plot1, width = 14,height = 12)
ggsave("enrichment3/GO_results_nfkb_pathway.png", plot=res_plot1, width = 14,height = 12)
GO_dot2 <- dotplot(ego2,showCategory = 17)
GO_bar2 <- barplot(ego2,showCategory = 17)
res_plot2 <- CombinePlots(list(GO_dot2,GO_bar2), nrow=1)
ggsave("enrichment3/GO_results_TNF_pathway.pdf", plot=res_plot2, width = 14,height = 12)
ggsave("enrichment3/GO_results_TNF_pathway.png", plot=res_plot2, width = 14,height = 12)

##8 细胞类型鉴定
dir.create("SingleR3")
Idents(scRNA3) <- 'seurat_clusters'
# 读入参考数据集
Mouse.se=MouseRNAseqData()
library(SingleR)
library(celldex)

norm_count = GetAssayData(scRNA3, 
                          slot="data") #归一化后的矩阵
pred<- SingleR(test = norm_count, 
               ref = Mouse.se, 
               labels = Mouse.se$label.main)
#将注释结果添加到seurat对象里
scRNA3$singleR_celltype = pred$labels

# 对细胞类型注释结果进行可视化
celltype_plot = DimPlot(scRNA3, 
                        reduction = "umap", 
                        group.by = "singleR_celltype")
ggsave("SingleR3/celltype_plot1.pdf", celltype_plot)
ggsave("SingleR3/celltype_plot1.png", celltype_plot)
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
#select_genes_Fibroblast<- c('Acta2')
#select_genes_Endothelium<- c('Ptprb','Pecam1')
select_genes_Epithelium<- c('Tppp3','Krt18','Krt13','Caps','Scgb1a1','Scgb3a1','Ecm1')
#select_genes_Mast<- c('Kit','Tpsab1')
select_genes_Neutrophil<- c('Prok2','S100a9', 'S100a8', 'Retnlg')
select_genes_Dendritic<- c('Cd1c','Clec9a','Xcr1','Wdfy4','Cpvl','Clec10a')
select_genes_pDendritic<- c('Il3ra','Gzmb','Serpinf1','Itm2c', 'Tcf4', 'Irf7', 'Clic3', 'Derl3')
#select_genes_Monocyte<- c('Tnip3')
select_genes_Macrophage<- c('Spp1','Cd68','Fabp4','C1qb','C1qa')
#select_genes_Bcell<- c('Ms4a1','Cd79a','Cd79b','Fcrl5')
#select_genes_Plasma<- c('Mzb1','Jchain')
#select_genes_Proliferative_signal<- c('Mki67','Top2a','Stmn1')
#select_genes_NK<- c('Fgfbp2','Gnly','Klrd1')
select_genes_Tcell<- c('Cd3d','Cd3e','Cd3g')
#select_genes_NK_Tcell<- c('Cd3d','Cd3e','Cd3g','Fgfbp2','Gnly','Klrd1')
select_genes<-c('Spp1','Cd68','Fabp4','Lgmn','Prok2','S100a9', 'S100a8', 'Retnlg','Cd3d','Cd3e','Cd3g',
                'Ccr7','Fscn1','Krt18','Scgb1a1','Scgb3a1')


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
#p15 <-DotPlot(scRNA3, features = select_genes_NK ) + coord_flip()

#ggsave("SingleR3/select_genes_Fibroblast.png", p1, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Endothelium.png", p2, width = 8, height = 10)
ggsave("SingleR3/select_genes_Epithelium2.png", p3, width = 8, height = 10)
ggsave("SingleR3/select_genes_Macrophage2.png", p4, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Mast.png", p5, width = 8, height = 10)
ggsave("SingleR3/select_genes_Neutrophil2.png", p6, width = 8, height = 10)
ggsave("SingleR3/select_genes_cDendritic.png", p7, width = 8, height = 10)
ggsave("SingleR3/select_genes_pDendritic2.png", p8, width = 8, height = 10)
ggsave("SingleR3/select_genes_Monocyte2.png", p9, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Bcell.png", p10, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Plasma.png", p11, width = 8, height = 10)
#ggsave("SingleR3/select_genes_Proliferative_signal.png", p12, width = 8, height = 10)
#ggsave("SingleR3/select_genes_NK_Tell.png", p13, width = 8, height = 10)
ggsave("SingleR3/select_genes_Tcell2.png", p14, width = 8, height = 10)
#ggsave("SingleR3/select_genes_NK.png", p15, width = 8, height = 10)


# 将singleR得到的细胞类型加入到对象的meta.data中，后续用于高级分析
scRNA3<- RenameIdents(scRNA3, '0'='Macrophage','1'='Macrophage','18'='Neutrophils','2'='Macrophage','3'='Macrophage','4'='Epithelial_cells',
                              '13'='DC','5'='Macrophage','6'='Epithelial_cells','7'='Macrophage','8'='Macrophage','9'='Macrophage',
                               '10'='T_cells','11'='Macrophage','12'='Macrophage','14'='Macrophage','17'='Macrophage','15'='Epithelial_cells',
                               '16'='Epithelial_cells','19'='Epithelial_cells')
scRNA3[["celltype"]] <- Idents(scRNA3)
head(scRNA3@meta.data)
saveRDS(scRNA3, "Step5_CellTyping_new.rds")

# 对细胞类型注释结果进行可视化
celltype_plot = DimPlot(scRNA3, 
                        reduction = "umap", 
                        group.by = "celltype1")
ggsave("SingleR3/celltype_plot_z1.pdf", celltype_plot)
ggsave("SingleR3/celltype_plot_z1.png", celltype_plot)
celltype_plot

p2 <- DimPlot(scRNA3, group.by = "celltype",split.by='group')
p3 <- DimPlot(scRNA3, group.by = 'group')
p1 <- DimPlot(scRNA3, group.by = "celltype")
p4 <- DimPlot(scRNA3, group.by = "seurat_clusters")
#查看cluster与细胞类型的对应关系
ggsave("SingleR3/celltype_contrast_z3.pdf",p3,width = 5, height = 4)
ggsave("SingleR3/celltype_contrast_z1.pdf",p1, width = 5, height = 4)
ggsave("SingleR3/celltype_contrast_z2.pdf",p2, width = 16, height =4)
ggsave("SingleR3/celltype_contrast_z1.pdf",p4, width = 5, height = 4)

# 通过热图了解每个cluster对应细胞类型的打分情况
tab <- table(cluster=scRNA3@meta.data$seurat_clusters, label=scRNA3@meta.data$celltype) 
cluster_type <- pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 
ggsave("SingleR3/cluster_type_z0.pdf", cluster_type)
ggsave("SingleR3/cluster_type_z0.png", cluster_type)
cluster_type