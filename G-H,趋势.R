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
scRNA1<-readRDS("Step5_HC_AA1_RR1.rds")
scRNA2<-readRDS("Step5_HC_AA2_RR2.rds")
Idents(scRNA1) <- 'group3'
Idents(scRNA2) <- 'group3'
#G寻找marker基因
scRNA1.markers.all <- FindAllMarkers(scRNA1,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

#挑选表达量前十的marker基因

scRNA1.markers <- scRNA1.markers.all %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 20, wt = avg_log2FC)

#利用prepareDataFromscRNA函数准备数据,showAverage 参数设为 TRUE 则表示对 基因细胞亚群一样的细胞取均值进行绘图,否则就是所有细胞进行绘图,默认使用 seurat 对象的 RNA assay 的 data 数据。

st.data <- prepareDataFromscRNA(object = scRNA1,
                                diffData = scRNA1.markers.all,
                                showAverage = TRUE)

#对每个cluster进行富集分析，这里采用GO富集分析，有需要的可以选择kegg富集分析

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,topn = 50,
                        seed = 5201314)
write.table(enrich, 
            "QUSHI/GO_BP_terms_groupAA1.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
#挑选需要展示的marker基因
markGenes = unique(scRNA1.markers$gene)[sample(1:length(unique(scRNA1.markers$gene)),40,
                                             replace = F)]

#绘制cluster基因表达折线图

P1=visCluster(object = st.data,ncol = 1,
           plot.type = "line",sample.order=c("HC","AA1","RR1"))
ggsave("QUSHI/select_genes_BP_groupAA1.RR1.pdf", P1, width = 5, height =9)


#H寻找marker基因
scRNA2.markers.all <- FindAllMarkers(scRNA2,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

#挑选表达量前十的marker基因

scRNA2.markers <- scRNA2.markers.all %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 20, wt = avg_log2FC)

#利用prepareDataFromscRNA函数准备数据,showAverage 参数设为 TRUE 则表示对 基因细胞亚群一样的细胞取均值进行绘图,否则就是所有细胞进行绘图,默认使用 seurat 对象的 RNA assay 的 data 数据。

st.data <- prepareDataFromscRNA(object = scRNA2,
                                diffData = scRNA2.markers.all,
                                showAverage = TRUE)

#对每个cluster进行富集分析，这里采用GO富集分析，有需要的可以选择kegg富集分析

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,topn = 50,
                        seed = 5201314)
write.table(enrich, 
            "QUSHI/GO_BP_terms_groupAA2.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

#挑选需要展示的marker基因
markGenes = unique(scRNA2.markers$gene)[sample(1:length(unique(scRNA2.markers$gene)),40,
                                             replace = F)]

#绘制cluster基因表达折线图

P1=visCluster(object = st.data,ncol = 1,
           plot.type = "line",sample.order=c("HC","AA2","RR2"))
ggsave("QUSHI/select_genes_BP_groupAA2.pdf", P1, width = 5, height =9)

