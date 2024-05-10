
#1、堆叠柱状图.这是比较普通也是最常用的细胞比例可视化方法。
#这种图在做微生物菌群的研究中非常常见。具体的思路是计算各个样本中细胞群的比例，形成数据框之后转化为长数据，用ggplot绘制即可。
dir.create("SingleR3/FOUR_GROUP")
table(Neut3$group)#查看各组细胞数
prop.table(table(Idents(Neut3)))
Idents(Neut3) <- 'celltype'
table(Idents(Neut3), Neut3$group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Neut3), Neut3$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,linewidth = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample group',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))
ggsave("SingleR3/FOUR_GROUP/celltype_prop1.pdf",  width = 12,height = 5)
ggsave("SingleR3/FOUR_GROUP/celltype_prop1.png", width = 12,height = 5)
#2、批量统计图
table(Neut3$orig.ident)#查看各组细胞数
prop.table(table(Idents(Neut3)))
table(Idents(Neut3), Neut3$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Neut3), Neut3$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

###添加分组信息
sample <- c("Ctrl_1","Ctrl_2","M_1","M_2","IFX_1","IFX_2","SN10_1","SN10_2")
group3 <- c("Ctrl","Ctrl","M","M","IFX","IFX","SN10","SN10")
samples <- data.frame(sample, group3)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group3']#R添加列

###作图展示
pplist = list()
sce_groups = c('Macrophage','Neutrophils','Monocyte','Epithelial_cells','T_cells')
library(ggplot2)
library(dplyr)
library(ggpubr)

#dplyr包经常和其他包的函数有冲突，需要选择一下优先级
library(conflicted)
conflict_prefer("filter","dplyr")
conflict_prefer("select","dplyr")
conflict_scout()

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =1)
  
  ###组间t检验分析，两者比较
  #labely = max(cellper_$percent)
 # compare_means(percent ~ group,  data = cellper_)
 # my_comparisons <- list( c("HC","patA","patR1","patR2"))
  #pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  #pplist[[group_]] = pp1

###组间多组比较
  labely = max(cellper_$percent)
  #Global_test
compare_means(percent ~ group,data = cellper_, method ="anova")
# Default method = "kruskal.test” for multiple groups
ggplot(cellper_, x = "group",y = "percent",color = "group")+stat_compare_means()
# Change method to anova
#ggplot(cellper_, x = "group",y = "percent",color = "group")+stat_compare_means(method = "anova")
# Perorm pairwise comparisons
compare_means(percent ~ group,  data = cellper_)
# Visualize: Specify the comparisons you want
my_comparisons<- list(c("M","Ctrl"),c("IFX","Ctrl"),c("SN10","Ctrl"),c("SN61","Ctrl"),c("IFX","M"),c("SN10","M"))
ggplot(cellper_, x = "group",y = "percent",color = "group")
pp1 =pp1 +stat_compare_means(comparisons =my_comparisons) # Add pairwise comparisons p-value
#如果t检验：pp1 =pp1 +stat_compare_means(comparisons =my_comparisons,method = "t.test",lable="p.signif") 
#如果需要指定标签显示的Y轴位置，可使用参数label.y
#ggplot(cellper_, x = "group",y = "percent",color = "group")+
#stat_compare_means(comparisons = my_comparisons, label.y = c(0.95,1,1.05,1.1,1.15,1.2))+
#stat_compare_means(label.y = 1.25)
pplist[[group_]] = pp1
}

library(cowplot)
p1=plot_grid(pplist[['Macrophage']]+stat_compare_means(label.y = 1.05),# Add global p-value
          pplist[['Neutrophils']]+stat_compare_means(label.y = 0.013),
          pplist[['Monocyte']]+stat_compare_means(label.y = 0.033),
          pplist[['Epithelial_cells']]+stat_compare_means(label.y =0.28),
          pplist[['T_cells']]+stat_compare_means(label.y = 0.056))
#如果方法分析，plot_grid(pplist[['Macrophage']]+stat_compare_means(method = "anova",label.y = 1.6),# Add global p-value
          #pplist[['Neutrophils']]+stat_compare_means(method = "anova",label.y = 1.6),
          #pplist[['Monocyte']]+stat_compare_means(method = "anova",label.y = 1.2),
          #pplist[['Epithelial_cells']]+stat_compare_means(method = "anova",label.y = 1.1),
          #pplist[['T_cells']]+stat_compare_means(method = "anova",label.y = 0.25))
ggsave("SingleR3/FOUR_GROUP/celltype_CPMPARE1.pdf",p1,width = 12, height = 8)
ggsave("SingleR3/FOUR_GROUP/celltype_CPMPARE1.png",p1, width = 12, height = 8)
#条图
#采用ggplot2绘制误差线默认是上下两个方向均绘出，但有时对于柱状图只显示一个方向的误差线效果更好。想要实现这一目的，可以修改geom_errorbar的ymax/ymin的参数（但会多显示一条直线）或者geom_errorbar在geom_bar之前（前提要求误差线比柱子要短）。下面提供一个彻底解决该问题的方法---新增geom_uperrorbar函数。
library(ggplot2)
#' @export
#' @rdname geom_linerange
geom_uperrorbar <- function(mapping = NULL, data = NULL,
                          stat = "identity", position = "identity",
                          ...,
                          na.rm = FALSE,
                          orientation = NA,
                          show.legend = NA,
                          inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
  default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
    alpha = NA),

  draw_key = draw_key_path,

  required_aes = c("x|y", "ymin|xmin", "ymax|xmax"),

  setup_params = function(data, params) {
    GeomLinerange$setup_params(data, params)
  },

  extra_params = c("na.rm", "orientation"),

  setup_data = function(data, params) {
    data$flipped_aes <- params$flipped_aes
    data <- flip_data(data, params$flipped_aes)
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)
    data <- transform(data,
      xmin = x - width / 2, xmax = x + width / 2, width = NULL
    )
    flip_data(data, params$flipped_aes)
  },

  draw_panel = function(data, panel_params, coord, width = NULL, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    #x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x,    NA, data$xmin, data$xmax))
    #y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$ymin, NA, data$ymin, data$ymin))
    sel <- data$y < 0 
    data$ymax[sel] <- data$ymin[sel]
    x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x))
    y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y))
    data <- new_data_frame(list(
      x = x,
      y = y,
      colour = rep(data$colour, each = 5),
      alpha = rep(data$alpha, each = 5),
      size = rep(data$size, each = 5),
      linetype = rep(data$linetype, each = 5),
      group = rep(1:(nrow(data)), each = 5),
      row.names = 1:(nrow(data) * 5)
    ))
    data <- flip_data(data, flipped_aes)
    GeomPath$draw_panel(data, panel_params, coord)
  }
)

new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }
  
  class(x) <- "data.frame"
  
  attr(x, "row.names") <- .set_row_names(n)
  x
}

#具体作图1
table(Neut3$group)#查看各组细胞数
prop.table(table(Idents(Neut3)))
Idents(Neut3) <- 'celltype'
table(Idents(Neut3), Neut3$group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Neut3), Neut3$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
colnames(Cellratio)<-c('celltype','group','Percentage')
Cellratio$Percentage<-Cellratio$Percentage*100#为了做图好看

table(Neut3$orig.ident)#查看各组细胞数
prop.table(table(Idents(Neut3)))
table(Idents(Neut3), Neut3$orig.ident)#各组不同细胞群细胞数
Cellratio1 <- prop.table(table(Idents(Neut3), Neut3$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio1 <- data.frame(Cellratio1)
library(reshape2)
cellper <- dcast(Cellratio1,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
###添加分组信息
sample <- c("Ctrl_1","Ctrl_2","M_1","M_2","IFX_1","IFX_2","SN10_1","SN10_2")
group3 <- c("Ctrl","Ctrl","M","M","IFX","IFX","SN10","SN10")
samples <- data.frame(sample, group3)#创建数据框
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group3']
ID<-c(1,2,3,4,5,6,7,8)
samples <- data.frame(sample, ID)
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$ID <- samples[rownames(cellper),'ID']#R添加列
#横向变纵向
ntt <- c('Macrophage','Neutrophils','Monocyte','Epithelial_cells','T_cells')##
cellper %>%
  mutate(
    ID = as.numeric(ID)) %>%
  filter(!is.na(ID)) %>%
  select(ID,all_of(ntt),group,sample) -> tmp_cv
melt(tmp_cv, id.vars = c("ID","group","sample"),
     measure.vars = ntt,
     variable.name = "celltype", 
     value.name = "Percentage") -> cellper2
cellper2$Percentage<-cellper2$Percentage*100#为了做图好看
# 获取各组各个celltype的平均值
df_mean <- aggregate(cellper2$Percentage, list(group=cellper2$group,
                                               celltype=cellper2$celltype), mean, na.rm=T)
# 获取各组各个celltype的标准差
df_sd <- aggregate(cellper2$Percentage, by=list(group=cellper2$group,
                        celltype=cellper2$celltype), FUN=sd, na.rm=T)
# 合并mean和sd
colnames(df_mean)[3] <- "mean"
colnames(df_sd)[3] <- "sd"
df_stat <- merge(df_mean, df_sd, by=c("group", "celltype"))
str(df_stat)
# 合并数据
Cellratio2<-merge(Cellratio,df_stat,by=c("group","celltype"))

library("ggpubr")
dodge <- position_dodge(width=0.8)
pp1=ggplot(Cellratio2,aes(x=celltype,y=Percentage,fill=group))+
    geom_bar(stat = 'identity', 
             #柱状图位置并排:
             position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
             width = 0.8,      #设置柱子宽度,使变量之间分开
             color='black')+        
    labs(x='Celltype',y='Percentage')+ 
    geom_errorbar(aes(x=celltype,ymin =Percentage, ymax =Percentage + sd),stat = 'identity',
                position =dodge,width=0.4)+
    theme_bw(base_size = 18)+  
    theme(axis.text = element_text(colour = 'black'))
###组间多组比较
d1<-compare_means(Percentage~group,group.by ='celltype', data = cellper2,method = "t.test")

# Visualize: Specify the comparisons you want
#my_comparisons<- list(c("HC","AA1"),c("HC","AA2"),c("HC","RR1"),c("HC","RR2"),c("AA1","AA2"),c("AA1","RR1"),c("AA2","RR2"),c("RR1","RR2"))
##具体作图2
table(Neut3$group)#查看各组细胞数
prop.table(table(Idents(Neut3)))
Idents(Neut3) <- 'celltype'
table(Idents(Neut3), Neut3$group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Neut3), Neut3$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
colnames(Cellratio)<-c('celltype','group','Percentage')
Cellratio$Percentage<-Cellratio$Percentage*100#为了做图好看

table(Neut3$orig.ident)#查看各组细胞数
prop.table(table(Idents(Neut3)))
table(Idents(Neut3), Neut3$orig.ident)#各组不同细胞群细胞数
Cellratio1 <- prop.table(table(Idents(Neut3), Neut3$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio1 <- data.frame(Cellratio1)
library(reshape2)
cellper <- dcast(Cellratio1,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
###添加分组信息
sample <- c("Ctrl_1","Ctrl_2","M_1","M_2","IFX_1","IFX_2","SN10_1","SN10_2")
group3 <- c("Ctrl","Ctrl","M","M","IFX","IFX","SN10","SN10")
samples <- data.frame(sample, group3)#创建数据框
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group3']
ID<-c(1,2,3,4,5,6,7,8)
samples <- data.frame(sample, ID)
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$ID <- samples[rownames(cellper),'ID']#R添加列
#横向变纵向
ntt <- c('Macrophage','Neutrophils','Monocyte','Epithelial_cells','T_cells')##
cellper %>%
  mutate(
    ID = as.numeric(ID)) %>%
  filter(!is.na(ID)) %>%
  select(ID,all_of(ntt),group,sample) -> tmp_cv
melt(tmp_cv, id.vars = c("ID","group","sample"),
     measure.vars = ntt,
     variable.name = "celltype", 
     value.name = "Percentage") -> cellper2
cellper2$Percentage<-cellper2$Percentage*100#为了做图好看
tg <-cellper2
tg$celltype = as.factor (tg$celltype)
p1<- ggplot(tg, aes(x = celltype, y = Percentage, fill = group)) +
    geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge()) +
    stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black",
                 width = 0.25,position = position_dodge( .9))
p1=p1 +  ylim(0, 105)#设置y轴范围
ggsave("SingleR3/FOUR_GROUP/celltype_CPMPARE3.pdf",p1,width = 12, height = 8)
ggsave("SingleR3/FOUR_GROUP/celltype_CPMPARE3.png",p1, width = 12, height = 8)
#方差齐性检验，利用car下的Levene检验
library(car); library(purrr)
tg %>% split(.$celltype) %>% map(~leveneTest(Percentage ~ group, data = .x, center = mean))
#不同dose下，Levene检验p值，这李表明各组总体方差不同(var.equal = FALSE)。
d1<-compare_means(Percentage~group,group.by ='celltype', data = cellper2,method = "t.test")
d2<-compare_means(Percentage~group,group.by ='celltype', data = cellper2,method = "kruskal.test")
anno1=0.047
anno2=0.002
anno3=0.034
anno4=0.037
p2=p1+geom_signif(annotation=formatC(anno1, digits=2),textsize=4,
               y_position=100, xmin=0.66, xmax=1.12, 
               tip_length = c(0.01, 0.01))
p3=p2+geom_signif(annotation=formatC(anno2, digits=2),textsize=4,
                y_position=94, xmin=0.89, xmax=1.35, 
               tip_length = c(0.01, 0.01))
p4=p3+geom_signif(annotation=formatC(anno3, digits=2),textsize=4,
               y_position=25, xmin=2.7, xmax=3.1, 
               tip_length = c(0.01, 0.01))
p5=p4+geom_signif(annotation=formatC(anno4, digits=2),textsize=4,
                                  y_position=19, xmin=2.9, xmax=3.3, 
                                 tip_length = c(0.01, 0.01))
p7=p5+ylim(0,105)
p7=p7+ theme(axis.text.x= element_text(  face="bold",   colour="black", size=10))+
 theme(axis.text.y= element_text(colour="black", size=10))+ 
 theme( axis.line.y = element_line(  colour = "black",   
                                        linewidth =0.5,   linetype = "solid",  lineend="round"  ))+
theme(panel.background = element_rect(fill = "white"))+
theme(axis.title.x=element_text(vjust=1,face="bold",  
                                      size=12),  # X axis title
            axis.title.y=element_text(size=12,face="bold",
                                      color = "black"))
ggsave("SingleR3/FOUR_GROUP/celltype_t_test2.png",width = 14, height = 8, dpi = 300)
ggsave("SingleR3/FOUR_GROUP/celltype_t_test2.pdf", width =14, height = 8, dpi = 300)
anno8="*"
anno9="**"
p2=p1+geom_signif(annotation=formatC(anno8, digits=2),textsize=4,
               y_position=100, xmin=0.66, xmax=1.12, 
               tip_length = c(0.01, 0.01))
p3=p2+geom_signif(annotation=formatC(anno9, digits=2),textsize=4,
                y_position=94, xmin=0.89, xmax=1.35, 
               tip_length = c(0.01, 0.01))
p4=p3+geom_signif(annotation=formatC(anno8, digits=2),textsize=4,
               y_position=25, xmin=2.7, xmax=3.1, 
               tip_length = c(0.01, 0.01))
p5=p4+geom_signif(annotation=formatC(anno8, digits=2),textsize=4,
                                  y_position=19, xmin=2.9, xmax=3.3, 
                                 tip_length = c(0.01, 0.01))
p7=p5+ylim(0,105)
p7=p7+ theme(axis.text.x= element_text(  face="bold",   colour="black", size=10))+
 theme(axis.text.y= element_text(colour="black", size=10))+ 
 theme( axis.line.y = element_line(  colour = "black",   
                                        linewidth =0.5,   linetype = "solid",  lineend="round"  ))+
theme(panel.background = element_rect(fill = "white"))+
theme(axis.title.x=element_text(vjust=1,face="bold",  
                                      size=12),  # X axis title
            axis.title.y=element_text(size=12,face="bold",
                                      color = "black"))
ggsave("SingleR3/FOUR_GROUP/celltype_t_test2.png",width = 14, height = 8, dpi = 300)
ggsave("SingleR3/FOUR_GROUP/celltype_t_test2.pdf", width =14, height = 8, dpi = 300)