# devtools::install_github("Hy4m/linkET", force = TRUE)
#https://mp.weixin.qq.com/s/q6nkujgTYlbOQpkENjyyxA
rm(list = ls())
library(org.Hs.eg.db)
library(reactome.db)
library(xlsx)
library(ReactomePA)
library(stringr)
library(BiocGenerics)
library(clusterProfiler)
library(enrichplot)
library(future)
library(future.apply)
library(linkET)
library(easyTCGA)
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(ggridges)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
setwd("E://YUE实验（新）/10.单基因富集分析/")
load(file = "TCGA-PAAD_mrna_expr_tpm.rdata")
expr <- log2(mrna_expr_tpm+1)
dim(expr)
#设定物种
GO_database <- 'org.Hs.eg.db ' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#hsa 人；mmu 小鼠
KEGG_database <- 'hsa' #
#提取一下单基因表达量
RPN1_expr <- expr["PML",]
RPN1_expr[,1:4]

##################第一种相关性分析
#自定义一个函数
cormatrixes <- function(x,y){
  tem <- linkET::correlate(t(x),t(y),engine = "WGCNA")
  tem1 <- as.data.frame(tem[[1]])
  tem1 <- cbind(rownames(tem1),tem1)
  tem1_long <- reshape2::melt(tem1,value.name = "correlation")
  tem2 <- as.data.frame(tem[[2]])
  tem2 <- cbind(rownames(tem2),tem2)
  tem2_long <- reshape2::melt(tem2,value.name = "pvalue")
  result <- cbind(tem1_long,tem2_long$pvalue)
  names(result) <- c("v1","v2","correlation","pvalue")
  return(result)
}
# 确保两个两个矩阵样本顺序一样
identical(colnames(expr),colnames(NLN_expr))
## [1] TRUE
cor_res <- cormatrixes(NLN_expr,expr)
head(cor_res)
cor_res <- na.omit(cor_res)
suppressMessages(library(dplyr))

NLN_related_mrna <- cor_res %>% 
  filter(correlation > 0.6, pvalue < 0.05) %>% 
  distinct(v2) %>% 
  pull(v2)

length(NLN_related_mrna)
head(NLN_related_mrna)
#如果用相关性分析
#把所有基因跟单个基因的相关性系数当做LogFC
#接https://www.jingege.wang/2020/11/08/%e5%9f%ba%e4%ba%8e%e5%8d%95%e5%9f%ba%e5%9b%a0%e6%89%b9%e9%87%8f%e7%9b%b8%e5%85%b3%e6%80%a7%e5%88%86%e6%9e%90%e7%9a%84gsea/

############第二种根据表达量分组
#根据表达量分组
set.seed(810)
sample_group <- ifelse(expr["PML",] > median(t(expr["PML",])), "high","low")
sample_group <- factor(sample_group, levels = c("low","high"))
library(easyTCGA)
deg_res <- diff_analysis(exprset = expr, 
                         group = sample_group,
                         is_count = F,
                         logFC_cut = 0.585, #绘制火山图后改回0.585
                         pvalue_cut = 0.05,
                         adjpvalue_cut = 0.05,
                         save = F
)
#做火山图先把log改0
# 提取limma的结果
deg_limma <- deg_res$deg_limma
head(deg_limma)
####先绘制一个火山图吧~
#https://mp.weixin.qq.com/s/jGDax2AD1iNdOpo2eQrvEQ
library(tidyverse)
res=deg_limma
names(res)
res <- res |> rownames_to_column("SYMBOL") |> arrange(desc(logFC))
##给差异基因打标签，logFC > 1且 padj < 0.05认为是上调基因，logFC < -1且 padj < 0.05认为是下调基因
df <- res |>  
  mutate(significant = case_when(logFC > 0.585 & adj.P.Val < 0.05 ~ "Up",
                                 abs(logFC) < 0.585 | adj.P.Val > 0.05 ~ "None",
                                 logFC < -0.585 & adj.P.Val < 0.05 ~ "Down"))
df$significant |> as.factor()

head(df)
#loading packages---------------------------
library(ggrepel)
library(ggfun)
library(grid)
names(df)
df <- df %>%
  dplyr::rename(padj = adj.P.Val)
#指定显示上调前5名的基因，和下调前5名的基因
p <- ggplot(data = df) + 
  geom_point(aes(x = logFC, y = -log10(padj), 
                 color = logFC,
                 size = -log10(padj))) + 
  geom_text_repel(data =  df %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::arrange(desc(-log10(padj))) %>%
                    dplyr::slice(1:5) %>%
                    dplyr::filter(significant == "Up"),
                  aes(x = logFC, y = -log10(padj), label = SYMBOL),
                  nudge_x = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  direction = "y",
                  hjust = "left",
                  max.overlaps = 200
  )+
  geom_text_repel(data =  df %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::filter(significant != "Up") %>%
                    dplyr::arrange(desc(-log10(padj))) %>%
                    dplyr::slice(1:5) %>%
                    dplyr::filter(significant == "Down"),
                  aes(x = logFC, y = -log10(padj), label = SYMBOL),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "left",
                  max.overlaps = 200
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 2) +  # Updated lines
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  #ggtitle(label = "Volcano Plot") + 
  xlim(c(-2, 2)) +  # Change range to -2 and 2
  ylim(c(-1, 20)) + # Changed range to -1 and 30
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)
  ) + 
  annotate(geom = "text", x = 1.5, y = 0, label = "p = 0.05", size = 5) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = -1, 
    xmax = -0.585,
    ymin = 18,
    ymax = 18
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = -1, 
    xmax = -0.585,
    ymin = 18,
    ymax = 18
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = 1, 
    xmax = 0.585,
    ymin = 18,
    ymax = 18
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = 1, 
    xmax = 0.585,
    ymin = 18,
    ymax = 18
  )
print(p)

#############################
########绘制完火山图重新改limma值
#设置种子
set.seed(810)
suppressMessages(library(clusterProfiler))
#gene ID转换
deg_entrezid <- bitr(deg_limma$genesymbol,fromType = "SYMBOL"
                     ,toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

####GSEA分析####
dd=deg_limma
dd$mRNAs <- rownames(dd)
#制作genelist
gene <- dd$mRNAs## 转换
library(clusterProfiler)
gene = bitr(gene, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=dd$logFC,
                      SYMBOL = dd$mRNAs)
gene_df <- merge(gene_df,gene,by='SYMBOL')
genelist <- gene_df
names(genelist)
colnames(genelist)[colnames(genelist) == "ENTREZID"] <- "gene"
genelist <- genelist %>% 
  dplyr::mutate(gene = as.character(gene))
## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt('h.all.v2023.2.Hs.entrez.gmt')
# 需要网络
GSEA <- GSEA(geneList,TERM2GENE =hallmarks)

###GSEA可视化
#GSEA峰峦图:
ridgeplot(GSEA,core_enrichment = TRUE,
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE)+
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ggtitle("NLN")  +  # Add title
  scale_fill_gradientn(colors = c("rosybrown1", "white", "#CDFF3E")) 
####美化峰峦图####
gsearesult <- GSEA@result
names(gsearesult)
gsearesult <- gsearesult %>% 
  dplyr::mutate(log10P = -log10(p.adjust)) %>%
  separate_rows(core_enrichment, sep = "/")
gsearesult <- gsearesult %>%
  left_join(genelist, by = c("core_enrichment" = "gene")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
custom_colors <- colorRampPalette(c("#8075ad", "#f5edf0","#f5da73", "#ffb800"))(100)
# 使用 iris 数据集，根据不同物种绘制峰峦图
p <- ggplot(gsearesult, aes(x = logFC, y = Description, fill = log10P)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  labs(x = "Log2FoldChange", y = "") +
  scale_fill_gradientn(colors = custom_colors) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p
ggsave("PML-GSEA-ridge.pdf", plot = p, height = 5.5, width = 7.5) #height = 5.5, width = 7.5
custom_colors1 <- colorRampPalette(c("#008bd0", "#eeeeee", "#fc4e00"))(100)
p1 <- ggplot(gsearesult, aes(x = 1, y = Description)) +
  geom_point(aes(size = abs(NES), color = NES)) + # 气泡大小和颜色映射 NES 数据
  scale_color_gradientn(colors = custom_colors1, name = "NES", limits = c(-4, 2)) + # 颜色渐变
  scale_size_continuous(range = c(2, 6)) + # 气泡大小范围
  labs(x = NULL, y = NULL) + # 去除 x 和 y 轴标签
  theme_minimal() + # 简洁主题
  theme(
    panel.background = element_blank(), # 无背景
    panel.grid = element_blank(),
    panel.border = element_blank(), # 无边框
    axis.ticks = element_blank(), # 去除刻度
    axis.text.x = element_blank() # 去除 x 轴标签
  )
p1
ggsave("PML-GSEA-rigde_point.pdf", p1, width = 5,height = 6)
# KEGG
#分析
KEGG<-enrichKEGG(deg_entrezid$ENTREZID,#KEGG富集分析
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
eKEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")##转回gene symbol
dt <- eKEGG@result
dt <- dt %>%
  filter(p.adjust < 0.05)
dotplot(eKEGG, showCategory = 30,orderBy = "GeneRatio")
write.csv(dt,file="PML-KEGG.csv",quote=F,row.names=F) #保存富集结果
ego=read.csv("PML-KEGG.csv", header = T)
#读取kegg富集结果文件
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
#创建一个包含所需列的数据框go，其中包括分类、ID、术语、相关基因、和调整后的p值等信息。
#读取基因的logFC文件
id.fc <- gene_df
genelist <- data.frame(ID = id.fc$SYMBOL, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
#read.table函数从文件中读取包含基因和logFC值的信息，并将其存储在名为genelist的数据框中。然后将genelist的行名设置为基因ID。
circ <- circle_dat(go, genelist)
#用circle_dat函数创建一个用于圆圈图的数据对象。
termNum = 10 #限定term数目
geneNum = nrow(genelist) #限定蛋白数目
chord <- chord_dat(circ, genelist[1:geneNum,],
                   go$Term[1:termNum]
                   )


# 获取 Set3 调色板，并根据需要重复颜色以适应数据行数
termCol <- brewer.pal(12, "Set3")
termCol <- rep(termCol, length.out = 14)  # 扩展调色板以适应 14 行数据
pdf(file="PML-KEGG.pdf",width = 11.5,height = 9)
#使用GOCluster函数创建一个用于树状图的数据对象，不同的参数设置调整图的外观和布局。
# 确保 termCol 被正确传递到 GOCluster 函数中
names(circ)
# 对logFC列进行排序
circ <- circ[order(circ$logFC), ]
GOCluster(circ,
          go$Term[1:termNum],
          lfc.space = 0.2,  # 倍数跟树间的空隙大小
          lfc.width = 2,    # 变化倍数的圆圈宽度
          term.col = termCol[1:termNum],  # 自定义 term 的颜色
          term.space = 0.2,  # 倍数跟 term 间的空隙大小
          term.width = 4)    # 富集 term 的圆圈宽度
dev.off()

####GO####
GO <- enrichGO(deg_entrezid$ENTREZID,
                   OrgDb = "org.Hs.eg.db",
                   ont = "ALL",
                   readable = T
)
##GO-条图##
dt <- GO@result[,c("ONTOLOGY","Description","Count","p.adjust")]
# 添加新的列 -Log10(p.adjust)
dt$`-Log10(p.adjust)` <- -log10(dt$p.adjust)
# # 提取 ONTOLOGY 为 BP 的前 20 行
# dt_BP <- dt[dt$ONTOLOGY == "BP", ][1:20, ]
# #初步绘制柱状图
# names(dt_BP)
# p1 <- ggplot(dt_BP)+
#   geom_col(aes(x = reorder(Description,-`-Log10(p.adjust)`),
#                y = Count),
#            color='black',   #柱子边框颜色
#            width=0.6,       #柱子宽度
#            fill='#A59ACA')+ #柱子填充色 文献的颜色：#d5a478
#   labs(x=NULL,y='Count',
#        title = 'GO-BP Terms TOP20')+
#   theme_test(base_size = 15)+ #主题基本大小
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.text = element_text(color = 'black',face = 'bold',size=8),
#         plot.margin = margin(1,0.5,0.5,2.5,'cm'),
#         panel.border = element_rect(size = 1),
#         axis.title = element_text(face = 'bold',size = 12),
#         plot.title = element_text(face = 'bold',
#                                   size=15,hjust = 0.5))
# p1
# #双y轴绘制
# p2 <- p1+
#   scale_y_continuous(expand = c(0,0),limits = c(0,80), #count最大值+20
#                      sec.axis = sec_axis(~./5,
#                                          name = '-Log10(p.adjust)',
#                                          breaks = seq(0,20,2)))+ #-log10的数值
#   geom_line(aes(x= reorder(Description,-`-Log10(p.adjust)`),
#                 y=`-Log10(p.adjust)`*5, #中位数
#                 group=1),
#             linetype=3,cex=0.6)+
#   geom_point(aes(x= reorder(Description,-`-Log10(p.adjust)`),
#                  y=`-Log10(p.adjust)`*5),
#              color = "black", fill = '#589c47', shape = 21, size=3.5)
# p2
# p3 <- p2+
#   geom_text(aes(x= reorder(Description,-`-Log10(p.adjust)`),
#                 y=Count,
#                 label=Count),
#             vjust=-0.5,size=3.5,fontface='bold')
# p3
# #添加图例
# #这个图例没有一处不需要改！！！！
# #都要改！！！！！！！！！！！！
# #先定义一个用来画框的数据框：
# df <- data.frame(a=c(2.5,2.5,18.5,18.5),
#                  b=c(80,71,71,80))
# p4 <- p3+
#   geom_line(data = df,aes(a,b),cex=0.5)+
#   geom_rect(aes(xmin=4,xmax=5.4,ymin=74,ymax=78),
#             fill='#A59ACA',color='black')+ #文献中颜色 #d5a478
#   annotate(geom='text',x=7,y=76,label='Count',
#            fontface='bold',size=4)+
#   annotate('segment',x=10.9,xend = 12.6,y=76,yend = 76,
#            linetype=3,cex=0.5)+
#   annotate(geom='point', x=11.7,y=76,
#            color = "black", fill = '#589c47', shape = 21, size = 3)+
#   annotate('text',x=15.4,y=76,label='-Log10(p.adjust)',
#            fontface='bold',size=4)
# p4


# 对数据框按-Log10(p.adjust)列进行降序排序
dt_sorted <- dt[order(dt$`-Log10(p.adjust)`, decreasing = TRUE), ]
# 提取前20行的数据
top <- head(dt_sorted, 20)
#初步绘制柱状图
names(top)
# 为每种ONTOLOGY分配颜色
ontology_colors <- c("BP" = '#A59ACA', "CC" = 'rosybrown3', "MF" = '#d5a478')

p5 <- ggplot(top) +
  geom_col(aes(x = reorder(Description, -`-Log10(p.adjust)`),
               y = Count,
               fill = ONTOLOGY),   # 根据Group填充颜色
           color = 'black',        # 柱子边框颜色
           width = 0.6) +          # 柱子宽度
  scale_fill_manual(values = ontology_colors) +  # 设置自定义颜色
  labs(x = NULL, y = 'Count') +
  theme_test(base_size = 15) +     # 主题基本大小
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black', face = 'bold', size = 8),
        plot.margin = margin(1, 0.5, 0.5, 2.5, 'cm'),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = 'bold', size = 12),
        plot.title = element_text(face = 'bold', size = 15, hjust = 0.5),
        legend.position = "left",    # 图例位置
        legend.key.size = unit(0.5, "cm"),  # 图例大小
        legend.text = element_text(size = 8)) +  # 图例文本大小
  guides(fill = guide_legend(title = NULL))  # 移除图例标题

p5
#双y轴绘制
p6 <- p5+
  scale_y_continuous(expand = c(0,0),limits = c(0,80), #count最大值+20 #原80
                     sec.axis = sec_axis(~./5,
                                         name = '-Log10(p.adjust)',
                                         breaks = seq(0,80,2)))+ #-log10的数值
  geom_line(aes(x= reorder(Description,-`-Log10(p.adjust)`),
                y=`-Log10(p.adjust)`*5, #中位数
                group=1),
            linetype=3,cex=0.6)+
  geom_point(aes(x= reorder(Description,-`-Log10(p.adjust)`),
                 y=`-Log10(p.adjust)`*5),
             color = "black", fill = '#589c47', shape = 21, size=3.5)
p6
p7 <- p6+
  geom_text(aes(x= reorder(Description,-`-Log10(p.adjust)`),
                y=Count,
                label=Count),
            vjust=-0.5,size=3.5,fontface='bold')
p7
#添加图例
#这个图例没有一处不需要改！！！！
#都要改！！！！！！！！！！！！
#先定义一个用来画框的数据框：
df <- data.frame(a=c(2.5,2.5,18.5,18.5), 
                 b=c(80,71,71,80))
p8 <- p7+
  geom_line(data = df,aes(a,b),cex=0.5)+
  geom_rect(aes(xmin=4,xmax=4.4,ymin=74,ymax=78),
            fill='#A59ACA',color='black')+ #文献中颜色 #d5a478
  geom_rect(aes(xmin=4.5,xmax=4.9,ymin=74,ymax=78),
            fill='rosybrown3',color='black')+ 
  geom_rect(aes(xmin=5,xmax=5.4,ymin=74,ymax=78),
            fill='#d5a478',color='black')+
  annotate(geom='text',x=7,y=76,label='Count',
           fontface='bold',size=4)+
  annotate('segment',x=10.9,xend = 12.6,y=76,yend = 76,
           linetype=3,cex=0.5)+
  annotate(geom='point', x=11.7,y=76,
           color = "black", fill = '#589c47', shape = 21, size = 3)+
  annotate('text',x=14.5,y=76,label='-Log10(p.adjust)',
           fontface='bold',size=4)
p8
########RPN1
#添加图例
#这个图例没有一处不需要改！！！！
#都要改！！！！！！！！！！！！
#先定义一个用来画框的数据框：
df <- data.frame(a=c(2.5,2.5,18.5,18.5), 
                 b=c(125,110,110,125))
p8 <- p7+
  geom_line(data = df,aes(a,b),cex=0.5)+
  geom_rect(aes(xmin=4,xmax=4.4,ymin=115,ymax=120),
            fill='#A59ACA',color='black')+ #文献中颜色 #d5a478
  geom_rect(aes(xmin=4.5,xmax=4.9,ymin=115,ymax=120),
            fill='rosybrown3',color='black')+
  #geom_rect(aes(xmin=5,xmax=5.4,ymin=230,ymax=248),
            #fill='#d5a478',color='black')+
  annotate(geom='text',x=6.5,y=118,label='Count',
           fontface='bold',size=4)+
  annotate('segment',x=10.9,xend = 12.6,y=118,yend = 118,
           linetype=3,cex=0.5)+
  annotate(geom='point', x=11.7,y=118,
           color = "black", fill = '#589c47', shape = 21, size = 3)+
  annotate('text',x=14.5,y=118,label='-Log10(p.adjust)',
           fontface='bold',size=4)
p8

#GO条目的气泡图
#https://blog.csdn.net/sweet_yemen/article/details/125496149
#选择所有GOterm
dotplot(GO, showCategory = 30,orderBy = "GeneRatio")
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names=F) #保存富集结果
ego=read.table("GO.txt", header = T,sep="\t",check.names=F)
#读取kegg富集结果文件
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
#创建一个包含所需列的数据框go，其中包括分类、ID、术语、相关基因、和调整后的p值等信息。
#读取基因的logFC文件
id.fc <- gene_df
genelist <- data.frame(ID = id.fc$SYMBOL, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
#read.table函数从文件中读取包含基因和logFC值的信息，并将其存储在名为genelist的数据框中。然后将genelist的行名设置为基因ID。
circ <- circle_dat(go, genelist)
#使用泡泡图展现更多维度的数据
GOBubble(circ, labels = 5)#label的意思是在绘图时仅对-log（p.adj）>3的term进行标注

#按照富集类型将泡泡图分开绘制、添加标题、更改颜色并且指定阈值
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)
#title指定标题，colour指定颜色，display指定按不同类型分开绘制。

####GO-新弦图####
GO_df <- GO@result
names(GO_df)
# [1] "ONTOLOGY"       "ID"             "Description"   
# [4] "GeneRatio"      "BgRatio"        "RichFactor"    
# [7] "FoldEnrichment" "zScore"         "pvalue"        
# [10] "p.adjust"       "qvalue"         "geneID"        
# [13] "Count"   
table(GO_df$ONTOLOGY)
# BP  CC  MF 
# 372  84  64 
# 提取 BP、CC、MF 的前 7
GO_df <- GO_df %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 7) %>%
  dplyr::ungroup()
write.csv(GO_df,"PML-go分类前7result.csv",row.names = F)
GO_ID <- GO_df %>% 
  dplyr::select(ID,geneID,p.adjust) %>%
  tidyr::separate_rows(geneID, sep = "/") 
DEG <- deg_limma
DEG <- DEG %>%
  dplyr::rename(geneID = genesymbol)
plot_df <- GO_ID %>%
  dplyr::left_join(DEG, by = c("geneID" = "geneID")) %>%
  dplyr::mutate(`-log10(p.adjust)` = -log10(p.adjust))
normalize <- function(x, new_min = 1, new_max = 2) {
  (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
}
GO_ID <- plot_df %>% 
  dplyr::select(ID) %>%
  dplyr::distinct(ID) %>%
  dplyr::rename(Chr = ID) %>%
  dplyr::mutate(Start = 1,
                End = 3)
plot_df <- plot_df %>%
  dplyr::mutate(Change = case_when(
    logFC > 0.585 ~ "Upregulated",      
    logFC < -0.585 ~ "Downregulated",   
    TRUE ~ "Unchanged"             # 其他情况表示无变化
  ))
names(plot_df)
# [1] "ID"               "geneID"           "p.adjust"        
# [4] "logFC"            "AveExpr"          "t"               
# [7] "P.Value"          "adj.P.Val"        "B"               
# [10] "-log10(p.adjust)""Change" 
Gene_df <- plot_df %>%
  dplyr::select(ID, geneID, logFC, Change) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(location = 1:n()) %>%
  dplyr::mutate(location = normalize(location, 1.1, 2.9)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(Chr = ID,
                Start = location) %>%
  dplyr::mutate(End = Start) %>%
  dplyr::select(1, 5, 6, 3)
Gene_df2 <- plot_df %>%
  dplyr::select(ID, `-log10(p.adjust)`) %>%
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::mutate(Start = 1.5,
                End = 2.5) %>%
  dplyr::select(1,3,4,2)
#绘图
pdf(file = "PML-GO新圈图.pdf",
    height = 8,
    width = 9)
circos.par("start.degree" = 180)

circos.genomicInitialize(GO_ID, 
                         plotType = c("labels"), 
                         axis.labels.cex = 0.8 * par("cex"),
                         labels.cex = 1 * par("cex"),
                         track.height = 0.05)

circos.genomicTrackPlotRegion(
  GO_ID, 
  track.height = 0.1, 
  stack = TRUE, 
  bg.border = NA,
  track.margin = c(0, 0), 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "#a8ddb5", border = "black", ...)
  })

circos.genomicTrack(
  Gene_df, 
  track.height = 0.25, 
  bg.col = "#f0f0f0",
  bg.border = NA,
  panel.fun = function(region, value, ...){
    for (i in c(-4, -1.8, 0, 2.8, 4)) {
      circos.lines(seq(1, 3, 0.05),
                   rep(i, 41),
                   col = "#000000",
                   lwd = 0.15,
                   lty = 2)
    }
    circos.yaxis(labels.cex = 0.5, 
                 lwd = 0.1,
                 tick.length = convert_x(0.2, "mm"))
    circos.genomicPoints(
      region, 
      value, 
      col = "#000000",
      pch = 21, 
      cex = 1,
      bg = ifelse(value > 0, "#fa9fb5", "#2ca25f"))
  }
)

col_fun2 = colorRamp2(breaks = c(0, 7.5, 15), colors = c("#fde0dd", "#fa9fb5", "#c51b8a"))

circos.genomicTrack(
  Gene_df2, 
  track.height = 0.35, 
  bg.col = NA,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data("sector.index")  
    circos.genomicRect(region, value, 
                       col = col_fun2(value[[1]]), 
                       border = NA, 
                       ytop.column = 1, 
                       ybottom = 0,
                       ...) 
  }
)

legend1 <- Legend(
  at = c(1, 2), 
  labels = c("Up regulated", "Down regulated"), 
  title = "log10(p.adjust)", 
  type = "points", pch = NA, 
  background = c("#fa9fb5", "#2ca25f"))

legend2 <- Legend(
  col_fun = col_fun2,
  title = "-log10(p.adjust)",
  direction = "horizontal")

pushViewport(viewport(x = 0.1, y = 0.2))
grid.draw(legend1)
y_coord <- 0.2
upViewport()

pushViewport(viewport(x = 0.9, y = 0.2))
grid.draw(legend2)
y_coord <- 0.2
upViewport()

circos.clear()

dev.off()
