setwd("E:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP/data/human")
scRNA3<-readRDS("Step6_Epithelial_cells.rds")
scRNA3@meta.data$group3 ="AA1,AA2,RR1,RR2,HC"## 增加一列标签 
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A1","A2","A3") ), "group3"]="AA1" ## A1 ,A3代表2个样本名。
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("A4","A5","A6") ), "group3"]="AA2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R1","R2","R3") ), "group3"]="RR1"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("R4","R5","R6") ), "group3"]="RR2"
scRNA3@meta.data[which(scRNA3@meta.data$orig.ident %in% c("HC1","HC2","HC3") ), "group3"]="HC"
## 默认识别的信息改成 组信息
Idents(scRNA3) <- 'group3'
Patients_Healthy<- FindMarkers(scRNA3, 
                                ident.1 = c("AA1", "AA2","RR1", "RR2"), #这里可以替换成选择的实验组
                                ident.2 = "HC", #这里可以替换成选择的对照组
                                min.pct = 0.25)

#GESA新可视化#
library(ggsci)
library(gggsea)
library(GSVA)
#读取基因集（本示例使用的 MSigDB 数据库中的 hallmark）
hallmark <- clusterProfiler::read.gmt("h.all.v2023.2.Hs.symbols.gmt")  
#读取基因列表（包含基因名称以及 log2 转化后的 Fold Change 信息）
Patients_Healthy <- Patients_Healthy %>%
  rownames_to_column("gene")
geneList = Patients_Healthy$avg_log2FC 
names(geneList) = Patients_Healthy$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
gsea.re1<- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = hallmark,  #基因集
    pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH')  #指定 p 值校正方法
     
#GSEA(geneList, exponent = 1, nPerm = 1000, minGSSize = 10,
#       maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE,
#       TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
##exponent: weight of each step
##nPerm置换检验的次数，默认为1000
##minGSSize富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤掉，默认为10
##maxGSSize富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤掉，默认为500
##verbose是否输出提示信息
##seed是否使结果具有可重复性
##by选择使用的统计学方法，默认为fgsea
#输出
write.table(gsea.re1, 'gsea_Patients_Healthy_hallmark_Epithelial_cells.xls', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<=0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
library(fgsea)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
  nproc = 0,#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_Patients_Healthy_Epithelial_cells.RData')
#画图
p1=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p1, filename = 'GSEA_Patients_Healthy_Epithelial_cells.pdf', width =6, height =3.2)