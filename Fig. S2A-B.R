#KEGG上调和下调
#.首先加载所需要的包
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(massdatabase)
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method",'auto')
library(readxl)
library(tidyr)
#加载自己的数据
#Patients_Healthy差异基因
 deg<- read_excel("D:/单细胞测序原始数据和报告-欧易生物-6个人+动物/singleP/data/Diffexp3/Patients_Healthy.xlsx", 
    col_types = c("text", "numeric", "numeric", 
        "numeric", "numeric", "numeric"))
 deg$symbol<-deg$Names
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
#输入数据整理
gene_up = deg[deg$log2FoldChange >0,'ENTREZID'] 
gene_down = deg[deg$log2FoldChange <0,'ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#KEGG富集分析
## 下载tgo的通路代谢物关系（别人下载）
download_kegg_pathway(path = "kegg_tgo_pathway",
                    sleep = 1,
                    organism = "hsa")#mmu是小鼠，hsa是人

## 然后读取所下载的数据，如果文件在data文件夹的mmu里面有kegg_tgo_pathway文件，目录在data按照下面读取文件
hsa_data <- 
  read_kegg_pathway(path = "kegg_tgo_pathway")
#，目录在mmu按照下面读取文件
#hsa_data <- read_kegg_pathway()
class(hsa_data)

## 转变数据库形式
kegg_pathway_database <-
  convert_kegg2metpath(data = hsa_data, path = ".")

##查看有多少通路
length(unlist(kegg_pathway_database@pathway_name))

path_num <- c(1:8,353:356)  ## 排除没有对应基因的几个通路,hsa为c(1:8),mmu为c(1:9)
## 提取1到356的通路(hsa是356，mmu是104)
result1 <- NULL
for (i in 1:356) {
  a <-  hsa_data[[i]]$pathway_id  # 提取id
  b <- hsa_data[[i]]$pathway_name # 提取名称
  e <- hsa_data[[i]]$gene_list # 提取gene对应关系
  
  c <- hsa_data[[i]]$describtion
  if ( i %in% path_num){
    d <-  'None;None'
  }else {
    d <- hsa_data[[i]]$pathway_class  # 提取前两层级注释
    Pathway_re <- e%>% mutate(pathway_id = a,
                              pathway_name = b,
                              pathway_class= d) %>%
      separate(pathway_class,into =c("Pathway1","Pathway2"),sep = ';',remove = FALSE)  ##前两层分开
    result1 <- rbind(result1,Pathway_re)}
}
write.table(result1, 'hsa_KEGG_path.csv', sep = '\t',row.names = FALSE,quote = FALSE)


result2 <- NULL
com_num <- c(1:8,353:356)    ## 排除没有对应化合物的几个通路
for (i in 1:356) {
  a <-  hsa_data[[i]]$pathway_id  # 提取id
  b <- hsa_data[[i]]$pathway_name # 提取名称
  f <- hsa_data[[i]]$compound_list # 提取化合物对应关系
  # c <- hsa_data[[i]]$describtion
  if ( i %in% com_num){
    d <-  'None;None'
  }else {
    d <- hsa_data[[i]]$pathway_class
    Compound_re <- f %>% mutate(pathway_id = a,
                                pathway_name = b,
                                pathway_class= d)%>%
      separate(pathway_class,into =c("Pathway1","Pathway2"),sep = ';',remove = FALSE)
    result2 <- rbind(result2,Compound_re)
  }
}
write.table(result2, 'tgo_KEGG_com.csv', sep = '\t',row.names = FALSE,quote = FALSE)



library(clusterProfiler)
##制作背景文件
com_ano <- result1 %>% dplyr::select(KEGG.ID,pathway_id,pathway_name)

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
ggsave("PATIENTS_VS_HC/KEGG_results_up_down.png", plot=g_kegg, width = 12,height = 30)
ggsave("PATIENTS_VS_HC/KEGG_results__up_down.pdf", plot=g_kegg, width = 12,height = 30)
#绘图
kk.up <- DOSE::setReadable(kk.up, 
                                         OrgDb="org.Hs.eg.db", #人类的注释包是org.Hs.eg.db，小鼠是org.Mm.eg.db
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
kk.down <- DOSE::setReadable(kk.down, 
                                         OrgDb="org.Hs.eg.db", #人类的注释包是org.Hs.eg.db，小鼠是org.Mm.eg.db
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
write.csv(as.data.frame(kk.up@result),file="PATIENTS_VS_HC/KEGG_enrichment_kk.up.csv")
write.csv(as.data.frame(kk.down@result),file="PATIENTS_VS_HC/KEGG_enrichment_kk.down.csv")
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
library(enrichplot)
library(ggnewscale)
x2 <- pairwise_termsim(kk.up)
x3 <- pairwise_termsim(kk.down)
p1=dotplot(kk.up,title="Group Patients vs.Healthy (Up)")+dotplot(kk.down,title="Group Patients vs.Healthy (Down)")+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))  #富集气泡图 
p2_1=cnetplot(kk.up,foldchange =log2(foldchange), 
                 node_label="all",# categ##circular = TRUE,#
                   colorEdge = TRUE)#网络图展示富集功能和基因的包含关系
p2_2=cnetplot(kk.down,foldchange =log2(foldchange), 
                 node_label="all",# categ##circular = TRUE,#
                   colorEdge = TRUE)
p3_1=emapplot(x2,foldchange =log2(foldchange), showCategory = 10,
                 node_label="all",# categ##circular = TRUE,#
                   colorEdge = TRUE) #网络图展示各富集功能之间共有基因关系
p3_2=emapplot(x3,foldchange =log2(foldchange), showCategory = 10,
                 node_label="all",# categ##circular = TRUE,#
                   colorEdge = TRUE) 
p4_1=heatplot(kk.up, showCategory = 10)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))  #富集气泡图 
p4_2=heatplot(kk.down, showCategory = 10)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))  #富集气泡图 
p5=barplot(kk.up,showCategory=10,title="Group Patients vs.Healthy (Up)")+barplot(kk.down,showCategory=10,title="Group Patients vs.Healthy (Down)")+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))  #富集气泡图 
ggsave("PATIENTS_VS_HC/KEGG_results_dotplot.png", plot=p1, width = 20,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_dotplot.pdf", plot=p1, width = 20,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up.png", plot=p2_1, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up.pdf", plot=p2_1, width = 24,height =14)
ggsave("PATIENTS_VS_HC/KEGG_results_emapplot_up.png", plot=p3_1, width = 15,height = 15)
ggsave("PATIENTS_VS_HC/KEGG_results_emapplot_up.pdf", plot=p3_1, width = 15,height = 15)
ggsave("PATIENTS_VS_HC/KEGG_results_heatplot_up.png", plot=p4_1, width = 20,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_heatplot_up.pdf", plot=p4_1, width = 20,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down.png", plot=p2_2, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down.pdf", plot=p2_2, width = 24,height =14)
ggsave("PATIENTS_VS_HC/KEGG_results_emapplot_down.png", plot=p3_2, width = 15,height = 15)
ggsave("PATIENTS_VS_HC/KEGG_results_emapplot_down.pdf", plot=p3_2, width = 15,height = 15)
ggsave("PATIENTS_VS_HC/KEGG_results_heatplot_down.png", plot=p4_2, width = 30,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_heatplot_down.pdf", plot=p4_2, width = 30,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_barplot.png", plot=p5, width = 20,height = 7)
ggsave("PATIENTS_VS_HC/KEGG_results_barplot.pdf", plot=p5, width = 20,height = 7)
#hilight.params选择展示的通路，如category = "Cell cycle"，其他通路变灰。
#layout选择布局模式，可选：'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'。
y<-c("NF-kappa B signaling pathway","TNF signaling pathway")
p6 <- cnetplot(kk.up,layout = "star", showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
p7 <- cnetplot(kk.up, showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up1.png", plot=p6, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up1.pdf", plot=p6, width = 24,height =14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up2.png", plot=p7, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up2.pdf", plot=p7, width = 24,height =14)

y<-c("NF-kappa B signaling pathway","TNF signaling pathway","Chemokine signaling pathway","Viral protein interaction with cytokine and cytokine receptor",
"Human cytomegalovirus infection","Cytokine-cytokine receptor interaction")
p8 <- cnetplot(kk.up, showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up3.png", plot=p8, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_up3.pdf", plot=p8, width = 24,height =14)
#down
y<-c("NF-kappa B signaling pathway","TNF signaling pathway")
p6 <- cnetplot(kk.down,layout = "star", showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
p7 <- cnetplot(kk.down, showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down1.png", plot=p6, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down1.pdf", plot=p6, width = 24,height =14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down2.png", plot=p7, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down2.pdf", plot=p7, width = 24,height =14)

y<-c("NF-kappa B signaling pathway","TNF signaling pathway","Chemokine signaling pathway","Viral protein interaction with cytokine and cytokine receptor",
"Human cytomegalovirus infection","Cytokine-cytokine receptor interaction")
p8 <- cnetplot(kk.down, showCategory = y,
              color.params = list(foldChange = FC, edge = F),
              cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 1),
              hilight.params = list(category = y, alpha_hilight = 1, alpha_no_hilight = 0.3),
              node_label ="all",colorEdge = TRUE)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down3.png", plot=p8, width = 24,height = 14)
ggsave("PATIENTS_VS_HC/KEGG_results_cnetplot_down3.pdf", plot=p8, width = 24,height =14)
#选特定信号通路
kk.up1_pathway <- kk.up[kk.up@result[["ID"]]  %in% c("hsa04064",
"hsa04668",
"hsa04062",
"hsa05163",
"hsa05167",
"hsa04061",
"hsa04657",
"hsa04621")]
kk.up1<-kk.up
kk.up1@result <- kk.up1_pathway
p1=dotplot(kk.up1,title="Group PATIENTS vs.HC (Up)")+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))  #富集气泡图 
p5=barplot(kk.up1,showCategory=10,title="Group PATIENTS vs.HC (Up)")+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
ggsave("KEGG_results_dotplot_PATIENTS_VS_HC (Up).pdf", plot=p1, width = 10,height = 7)
ggsave("KEGG_results_barplot_PATIENTS_VS_HC (Up).pdf", plot=p5, width = 10,height = 7)


#其他类似或者#

kk.up1_pathway<-rbind(kk.up1_pathway_Pat_vs._HC,kk.up1_pathway_AA1_vs._HC,kk.up1_pathway_AA2_vs._HC,kk.up1_pathway_RR1_vs._HC,kk.up1_pathway_RR2_vs._HC)
kk.up1_pathway<-na.omit(kk.up1_pathway)

A<-kk.up1_pathway
write.table(A, 
            "Fig. S2A-diff_group_1.xlsx", 
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
ggsave("Fig.S2A.KEGG_results_dotplot_groups (Up).pdf", plot=p1, width = 7,height = 5)

kk.up2_pathway<-rbind(kk.up1_pathway_AA2_vs._AA1,kk.up1_pathway_RR1_vs._AA1,kk.up1_pathway_RR2_vs._AA2,kk.up1_pathway_RR2_vs._RR1)
kk.up2_pathway<-na.omit(kk.up1_pathway)

A<-kk.up2_pathway
write.table(A, 
            "Fig. S2B-diff_group_2.1.xlsx", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
A$Description <- as.factor(A$Description,levels=c())
A$Description <- fct_inorder(A$Description)
p2=ggplot(A, aes(group, Description)) +
    geom_point(aes(color=p.adjust, size=GeneRatio))+theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=30,hjust = 1,vjust=1))+
    scale_color_gradient2(low='red',mid ='grey',high='blue',midpoint = 0.012)+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
ggsave("Fig.S2B.KEGG_results_dotplot_groups (Up).pdf", plot=p2, width = 7,height = 5)