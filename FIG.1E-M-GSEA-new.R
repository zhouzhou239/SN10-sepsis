library(forcats)
library(ggstance)
##arrange去给NES的绝对值从大到小排序
Patients_Healthy <- Patients_Healthy %>%
  rownames_to_column("gene")

geneList = Patients_Healthy$avg_log2FC 
names(geneList) = Patients_Healthy$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_Patients_Healthy_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'Patients_Healthy_all.pdf', width =6, height =4)


AA1.HC.markers <- AA1.HC.markers %>%
  rownames_to_column("gene")

geneList = AA1.HC.markers$avg_log2FC 
names(geneList) = AA1.HC.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_AA1.HC.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'AA1.HC.markers_all.pdf', width =6, height =4)


AA2.HC.markers <- AA2.HC.markers %>%
  rownames_to_column("gene")

geneList = AA2.HC.markers$avg_log2FC 
names(geneList) = AA2.HC.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_AA2.HC.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'AA2.HC.markers_all.pdf', width =6, height =4)

RR1.HC.markers <- RR1.HC.markers %>%
  rownames_to_column("gene")

geneList = RR1.HC.markers$avg_log2FC 
names(geneList) = RR1.HC.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_RR1.HC.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'RR1.HC.markers_all.pdf', width =6, height =4)

RR2.HC.markers <- RR2.HC.markers %>%
  rownames_to_column("gene")

geneList = RR2.HC.markers$avg_log2FC 
names(geneList) = RR2.HC.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_RR2.HC.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'RR2.HC.markers_all.pdf', width =6, height =4)

RR1.AA1.markers <- RR1.AA1.markers %>%
  rownames_to_column("gene")

geneList = RR1.AA1.markers$avg_log2FC 
names(geneList) = RR1.AA1.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_RR1.AA1.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'RR1.AA1.markers_all.pdf', width =6, height =4)


RR2.AA2.markers <- RR2.AA2.markers %>%
  rownames_to_column("gene")

geneList = RR2.AA2.markers$avg_log2FC 
names(geneList) = RR2.AA2.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_RR2.AA2.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'RR2.AA2.markers_all.pdf', width =6, height =4)

RR2.RR1.markers <- RR2.RR1.markers %>%
  rownames_to_column("gene")

geneList = RR2.RR1.markers$avg_log2FC 
names(geneList) = RR2.RR1.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_RR2.RR1.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'RR2.RR1.markers_all.pdf', width =6, height =4)

AA2.AA1.markers <- AA2.AA1.markers %>%
  rownames_to_column("gene")

geneList = AA2.AA1.markers$avg_log2FC 
names(geneList) = AA2.AA1.markers$gene
geneList = sort(geneList,decreasing = T)
geneList[1:10]
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 0.5)
egmt1 <- egmt[egmt@result$p.adjust<=0.05,]
egmt@result <- egmt1

#结果转化
y=data.frame(egmt)
head(y)
write.table(y, 
            "GSEA_AA2.AA1.markers_all.xls", 
            col.names = T, 
            row.names = T, 
            sep = "\t")
#GSEA结果可视化
#气泡图，展示geneset被激活还是抑制
p1=dotplot(egmt,split=".sign")+facet_grid(~.sign)+ theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8))
ggsave(p1, filename = 'AA2.AA1.markers_all.pdf', width =6, height =4)



#GESA新可视化#
library(ggsci)
library(gggsea)
library(GSVA)
#读取基因集（本示例使用的 MSigDB 数据库中的 hallmark）
hallmark <- clusterProfiler::read.gmt("h.all.v2023.2.Hs.symbols.gmt")  
#读取基因列表（包含基因名称以及 log2 转化后的 Fold Change 信息）
AA1.HC.markers <- AA1.HC.markers %>%
  rownames_to_column("gene")
geneList = AA1.HC.markers$avg_log2FC 
names(geneList) = AA1.HC.markers$gene
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
write.table(gsea.re1, 'gsea_AA1.HC.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_AA1.HC.markers.RData')
#画图
p1=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p1, filename = 'GSEA_AA1.HC.markers.pdf', width =6, height =1.2)

AA2.HC.markers <- AA2.HC.markers %>%
  rownames_to_column("gene")
geneList = AA2.HC.markers$avg_log2FC 
names(geneList) = AA2.HC.markers$gene
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
write.table(gsea.re1, 'gsea_AA2.HC.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_AA2.HC.markers.RData')
#画图
p2=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p2, filename = 'GSEA_AA2.HC.markers.pdf', width =6, height =1.6)

RR1.HC.markers <- RR1.HC.markers %>%
  rownames_to_column("gene")
geneList = RR1.HC.markers$avg_log2FC 
names(geneList) = RR1.HC.markers$gene
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
write.table(gsea.re1, 'gsea_RR1.HC.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_RR1.HC.markers.RData')
#画图
p3=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p3, filename = 'GSEA_RR1.HC.markers.pdf', width =6, height =1.8)

RR2.HC.markers <- RR2.HC.markers %>%
  rownames_to_column("gene")
geneList = RR2.HC.markers$avg_log2FC 
names(geneList) = RR2.HC.markers$gene
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
write.table(gsea.re1, 'gsea_RR2.HC.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_RR2.HC.markers.RData')
#画图
p4=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p4, filename = 'GSEA_RR2.HC.markers.pdf', width =6, height =2.6)

RR2.RR1.markers <- RR2.RR1.markers %>%
  rownames_to_column("gene")
geneList = RR2.RR1.markers$avg_log2FC 
names(geneList) = RR2.RR1.markers$gene
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
write.table(gsea.re1, 'gsea_RR2.RR1.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_RR2.RR1.markers.RData')
#画图
p5=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p5, filename = 'GSEA_RR2.RR1.markers.pdf', width =6, height =1.4)

AA2.AA1.markers <- AA2.AA1.markers %>%
  rownames_to_column("gene")
geneList = AA2.AA1.markers$avg_log2FC 
names(geneList) = AA2.AA1.markers$gene
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
write.table(gsea.re1, 'gsea_AA2.AA1.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_AA2.AA1.markers.RData')
#画图
p6=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p6, filename = 'GSEA_AA2.AA1.markers.pdf', width =6, height =1.4)


RR1.AA1.markers <- RR1.AA1.markers %>%
  rownames_to_column("gene")
geneList = RR1.AA1.markers$avg_log2FC 
names(geneList) = RR1.AA1.markers$gene
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
write.table(gsea.re1, 'gsea_RR1.AA1.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_RR1.AA1.markers.RData')
#画图
p7=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p7, filename = 'GSEA_RR1.AA1.markers.pdf', width =6, height =1.6)

RR2.AA2.markers <- RR2.AA2.markers %>%
  rownames_to_column("gene")
geneList = RR2.AA2.markers$avg_log2FC 
names(geneList) = RR2.AA2.markers$gene
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
write.table(gsea.re1, 'gsea_RR2.AA2.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_RR2.AA2.markers.RData')
#画图
p8=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p8, filename = 'GSEA_RR2.AA2.markers.pdf', width =6, height =1.8)


pat.HC.markers <- pat.HC.markers %>%
  rownames_to_column("gene")
geneList = pat.HC.markers$avg_log2FC 
names(geneList) = pat.HC.markers$gene
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
write.table(gsea.re1, 'gsea_pat.HC.markers_hallmark.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
#fgsea分析
## 这里去掉了基因集前缀
##hallmark$ont <- str_remove(hallmark$term,"HALLMARK_")
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
##hallmark.list <- hallmark %>% split(.$ont) %>% lapply( "[[", 2)
gsea.re2 <-fgseaMultilevel(
  pathways = hallmark.list,
 stats = geneList,
  sampleSize = 101,
  minSize = 1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
  maxSize = 10000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
  eps = 1e-50,
  scoreType = c("std", "pos", "neg"),
# nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
  gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
  BPPARAM = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
  nPermSimple = 1000,#置换检验的次数
  absEps = NULL)

#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05&NES>0,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea_pat.HC.markers.RData')
#画图
p9=plotGseaTable(hallmark.list[g2$pathway],
              geneList, pathwayLabelStyle =list(size=8), headerLabelStyle =list(size=8),
              valueStyle =list(size=8),
              axisLabelStyle = list(size=8),
              gsea.re2,gseaParam = 1,
              colwidths = c(0.5,0.2,0.1,0.1,0.1))

ggsave(p9, filename = 'GSEA_pat.HC.markers.pdf', width =6, height =2.4)

