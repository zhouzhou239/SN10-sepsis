deg<- read_excel("E:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP/data/mouse/mouse/Marker/M.Ctrl.markers.xls", 
    col_types = c("text", "numeric", "numeric", 
        "numeric", "numeric", "numeric"))
deg$symbol<-deg$Names
deg$log2FoldChange<-deg$avg_log2FC
#一、ENTREZID转换
library(clusterProfiler)
library(org.Hs.eg.db)      #Hs代表人类
s2e <- bitr(unique(deg$symbol), 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
dim(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
#输入数据准备,拿出up基因和down基因的ENTREZID，组成diff基因
gene_up = deg[deg$log2FoldChange >0,'ENTREZID'] 
gene_down = deg[deg$log2FoldChange <0,'ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']

#####超几何分布来判断是否富集
#上调基因的通路富集
gene_select1 <- sample(gene_up$ENTREZID)
 kk.up <- enricher(gene = gene_select1,
                      TERM2GENE = com_ano[c('pathway_id', 'KEGG.ID')], 
                      TERM2NAME = com_ano[c('pathway_id', 'pathway_name')],
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH', 
                      qvalueCutoff = 1, 
                      maxGSSize = 500)
head(kk.up)[,1:6]
#下调基因的通路富集
gene_select2 <- sample(gene_down$ENTREZID)
kk.down <- enricher(gene = gene_select2,
                      TERM2GENE = com_ano[c('pathway_id', 'KEGG.ID')], 
                      TERM2NAME = com_ano[c('pathway_id', 'pathway_name')],
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH', 
                      qvalueCutoff = 1, 
                      maxGSSize = 500)
head(kk.down )[,1:6]
#接着将得到的结果加上一列分组信息，以备后续作图使用
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)

down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)
#####绘图
#自定义绘图函数keggplot
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10Pvalue") +
    coord_flip() + 
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
  return(g_kegg)
}
##每次使用直接将代码复制放在project里面，用之前
#source("functions.R")
##之后直接使用函数即可,绘制图像仅需要一句代码即可
g_kegg = kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave("M_vs._Ctrl/KEGG_results_up_down.png", plot=g_kegg, width = 12,height = 30)
ggsave("M_vs._Ctrl/KEGG_results__up_down.pdf", plot=g_kegg, width = 12,height = 30)
#绘图
kk.up <- DOSE::setReadable(kk.up, 
                                         OrgDb="org.Hs.eg.db", #人类的注释包是org.Hs.eg.db，小鼠是org.Mm.eg.db
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
kk.down <- DOSE::setReadable(kk.down, 
                                         OrgDb="org.Hs.eg.db", #人类的注释包是org.Hs.eg.db，小鼠是org.Mm.eg.db
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
#Error in UseMethod("rescale") : "rescale"没有适用于"AsIs"目标对象的方法,解决方案
dropAsis <- function(x){
    cls <- class(x)
    structure(x, class = setdiff(cls, "AsIs"))
  }
rescale.AsIs <- function(x, ...){
  # 自定义dropAsis方法
  dropAsis <- function(x){
    cls <- class(x)
    structure(x, class = setdiff(cls, "AsIs"))
  }
  # 调用本来的rescale方法
  scales:::rescale(dropAsis(x), ...)
}
#在Rstudio界面中，点击stop退出debug状态。
#然后将ggplot_build.ggplot函数退出debug模式（记得再运行一下刚才的rescale.AsIs函数的定义）。
undebug(ggplot2:::ggplot_build.ggplot)
stopifnot(exists("rescale.AsIs")) 
#选特定信号通路
kk.up1_pathway_M_vs._Ctrl <- kk.up[kk.up@result[["ID"]]  %in% c("hsa04064",
"hsa04668",
"hsa04062",
"hsa05163",
"hsa05167",
"hsa04061",
"hsa04657",
"hsa04621")]
kk.up1_pathway_M_vs._Ctrl$group<-"M_vs._Ctrl"


kk.up1_pathway<-rbind(kk.up1_pathway_M_vs._Ctrl,kk.up1_pathway_AA1_vs._HC,kk.up1_pathway_AA2_vs._HC,kk.up1_pathway_RR1_vs._HC,kk.up1_pathway_RR2_vs._HC)
kk.up1_pathway<-na.omit(kk.up1_pathway)

A<-kk.up1_pathway
write.table(A, 
            "Figure2/diff_group_2.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
A$Description <- as.factor(A$Description,levels=c())
A$Description <- fct_inorder(A$Description)
p1=ggplot(A, aes(group, Description)) +
    geom_point(aes(color=p.adjust, size=GeneRatio))+theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=30,hjust = 1,vjust=1))+
    scale_color_gradient2(low='red',mid ='grey',high='blue',midpoint = 0.012)+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
ggsave("Fig.M_vs._Ctrl.KEGG_results_dotplot_groups (Up).pdf", plot=p1, width = 7,height = 5)

#画图参考网址https://zhuanlan.zhihu.com/p/445134223