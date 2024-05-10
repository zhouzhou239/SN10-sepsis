library(ClusterGVis)
library(org.Hs.eg.db)
library(ggplot2)
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
setwd("E:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP/data/human")
#导入单细胞上游分析结果文件
scRNA3<-readRDS("Step5_CellTyping.rds")
Idents(scRNA3) <- 'celltype'
#寻找marker基因
scRNA3.markers.all <- FindAllMarkers(scRNA3,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

#挑选表达量前十的marker基因

scRNA3.markers <- scRNA3.markers.all %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 20, wt = avg_log2FC)

#利用prepareDataFromscRNA函数准备数据,showAverage 参数设为 TRUE 则表示对 基因细胞亚群一样的细胞取均值进行绘图,否则就是所有细胞进行绘图,默认使用 seurat 对象的 RNA assay 的 data 数据。

st.data <- prepareDataFromscRNA(object = scRNA3,
                                diffData = scRNA3.markers.all,
                                showAverage = TRUE)

#对每个cluster进行富集分析，这里采用GO富集分析，有需要的可以选择kegg富集分析

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314)

#挑选需要展示的marker基因
markGenes = unique(scRNA3.markers$gene)[sample(1:length(unique(scRNA3.markers$gene)),50,
                                             replace = F)]

#绘制cluster基因表达折线图

P1=visCluster(object = st.data,
           plot.type = "line")
ggsave("QUSHI/select_genes_BP.pdf", P1, width = 9, height =5)

#绘制热图
pdf('QUSHI/sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
plot.type = "heatmap",
column_names_rot = 45,
markGenes = markGenes,
cluster.order = c(1:9))
dev.off()

#绘制热图并添加富集注释和分组折线图
pdf('QUSHI/sc2_GO_BP.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
plot.type = "both",
column_names_rot = 45,
show_row_dend = F,
markGenes = markGenes,
markGenes.side = "left",
annoTerm.data = enrich,
line.side = "left",
cluster.order = c(1:9),
go.col = rep(jjAnno::useMyCol("stallion",n = 10),each = 5),
add.bar = T)
dev.off()
#KEGG
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "KEGG",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314)

#绘制热图并添加富集注释和分组折线图
pdf('QUSHI/FIG.S1C.sc2_KEGG3.pdf',height = 14,width = 14,onefile = F)
visCluster(object = st.data,
plot.type = "both",
column_names_rot = 45,
show_row_dend = F,
markGenes = markGenes,
markGenes.side = "left",
annoTerm.data = enrich,
line.side = "left",
cluster.order = c(1:9),
go.col = rep(jjAnno::useMyCol("stallion",n = 10),each = 5),
add.bar = T)
dev.off()
