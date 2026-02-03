#https://mp.weixin.qq.com/s/RRUKsHEbg9YiMmZoqTfoIQ
####处理数据####
rm(list = ls())
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggpubr)
library(tibble)
#remotes::install_github('jorvlan/raincloudplots')
library(raincloudplots)
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(tidyverse)
setwd("E://YUE实验（新）/11.免疫浸润1/")
load("TCGA-PAAD_mrna_expr_tpm.rdata")
tpm <- mrna_expr_tpm
k=as.numeric(stringr::str_sub(colnames(tpm),14,15)) < 10
tpm=tpm[,k] # 只保留肿瘤样本colnames(tpm)=stringr::str_sub(colnames(tpm),1,16)tpm=tpm[,!duplicated(colnames(tpm))]
colnames(tpm)=stringr::str_sub(colnames(tpm),1,12)
tpm=tpm[,!duplicated(colnames(tpm))]
tpm=tpm[!rowSums(tpm==0)==ncol(tpm),] # 删除全是0的基因
write.table(tpm, file = "tpm.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
save(tpm,file="PAAD_tpm.Rdata")

####immunedeconv####
rm(list=ls())
load("PAAD_tpm.Rdata")
library(immunedeconv)
####epic####
epic=deconvolute(tpm,'epic')

####estimate####
estimate=deconvolute_estimate(tpm)
write.csv(estimate,"estimate.csv",row.names = T)

es=as.data.frame(t(estimate))
names(es)
# [1] "StromalScore"  "ImmuneScore"   "ESTIMATEScore"
# [4] "TumorPurity"  
group <- read.csv("group.csv",row.names = 1)
# 检查两个数据框的行名是否一致
common_rownames <- intersect(rownames(es), rownames(group))
# 保留两个数据框中具有相同行名的部分
es <- es[common_rownames, ]
group <- group[common_rownames, ]
# 按行名合并
data <- cbind(es, group)
#estimate-StromalScore
# 提取 "high" 和 "low" 组的分数，假设使用 StromalScore 列
data_high <- data[data$group == "high", "StromalScore"]  # "high" 组
data_low <- data[data$group == "low", "StromalScore"]    # "low" 组
# 使用 data_1x1 函数进行转换
dt1 <- data_1x1(
  array_1 = data_high,  # 分组1 ("high")
  array_2 = data_low,   # 分组2 ("low")
  jit_distance = 0.09,  # 抖动点间的距离
  jit_seed = 123        # 设置随机种子
)
# 查看转换后的数据
head(dt1)
#1.3 配对云雨图绘制：
head(dt1) #对两组中的相同id进行配对；
colnames(dt1)
wilcox_test_result <- wilcox.test(data_high, data_low)
p_value <- wilcox_test_result$p.value
# 格式化 p 值
p_value_formatted <- paste("p =", format(p_value, digits = 3))
p1 <- raincloud_1x1_repmes(
  data = dt1,
  colors = (c('#6F99AD','#BF5960')),
  fills = (c("#6F99AD","#BF5960")),
  size = 1.6,
  alpha = 0.5,
  line_color = 'grey',
  line_alpha = 0.8,
  align_clouds = FALSE
) +
  scale_x_continuous(breaks = c(0.6, 2.4),
                     labels = c("low", "high"),
                     limits = c(0, 3)) +
  labs(x = "", y = "StromalScore") +
  theme_classic() +
  # 调整x轴标签字体大小
  theme(axis.text.x = element_text(size = 14)) +  # 调大x轴标签的字体大小
  # 添加 p 值注释
  annotate("text", x = 1.5, y = max(dt1$y, na.rm = TRUE) * 1.05, 
           label = p_value_formatted, vjust = -1) +
  # 添加横线
  annotate("segment", x = 0.6, xend = 2.4, 
           y = max(dt1$y, na.rm = TRUE) * 1.03, 
           yend = max(dt1$y, na.rm = TRUE) * 1.03)

p1
ggsave("estimate-StromalScore.pdf", plot = p1, width = 6, height = 8)
#estimate-ImmuneScore
data_high <- data[data$group == "high", "ImmuneScore"]  # "high" 组
data_low <- data[data$group == "low", "ImmuneScore"]    # "low" 组
# 使用 data_1x1 函数进行转换
dt1 <- data_1x1(
  array_1 = data_high,  # 分组1 ("high")
  array_2 = data_low,   # 分组2 ("low")
  jit_distance = 0.09,  # 抖动点间的距离
  jit_seed = 123        # 设置随机种子
)
# 查看转换后的数据
head(dt1)
#1.3 配对云雨图绘制：
head(dt1) #对两组中的相同id进行配对；
colnames(dt1)
wilcox_test_result <- wilcox.test(data_high, data_low)
p_value <- wilcox_test_result$p.value
# 格式化 p 值
p_value_formatted <- paste("p =", format(p_value, digits = 3))
p2 <- raincloud_1x1_repmes(
  data = dt1,
  colors = (c('#6F99AD','#BF5960')),
  fills = (c("#6F99AD","#BF5960")),
  size = 1.6,
  alpha = 0.5,
  line_color = 'grey',
  line_alpha = 0.8,
  align_clouds = FALSE
) +
  scale_x_continuous(breaks = c(0.6, 2.4),
                     labels = c("low", "high"),
                     limits = c(0, 3)) +
  labs(x = "", y = "ImmuneScore") +
  theme_classic() +
  # 调整x轴标签字体大小
  theme(axis.text.x = element_text(size = 14)) +  # 调大x轴标签的字体大小
  # 添加 p 值注释
  annotate("text", x = 1.5, y = max(dt1$y, na.rm = TRUE) * 1.05, 
           label = p_value_formatted, vjust = -1) +
  # 添加横线
  annotate("segment", x = 0.6, xend = 2.4, 
           y = max(dt1$y, na.rm = TRUE) * 1.03, 
           yend = max(dt1$y, na.rm = TRUE) * 1.03)

p2
ggsave("estimate-ImmuneScore.pdf", plot = p2, width = 6, height = 8)
#estimate-ESTIMATEScore
data_high <- data[data$group == "high", "ESTIMATEScore"]  # "high" 组
data_low <- data[data$group == "low", "ESTIMATEScore"]    # "low" 组
# 使用 data_1x1 函数进行转换
dt1 <- data_1x1(
  array_1 = data_high,  # 分组1 ("high")
  array_2 = data_low,   # 分组2 ("low")
  jit_distance = 0.09,  # 抖动点间的距离
  jit_seed = 123        # 设置随机种子
)
# 查看转换后的数据
head(dt1)
#1.3 配对云雨图绘制：
head(dt1) #对两组中的相同id进行配对；
colnames(dt1)
wilcox_test_result <- wilcox.test(data_high, data_low)
p_value <- wilcox_test_result$p.value
# 格式化 p 值
p_value_formatted <- paste("p =", format(p_value, digits = 3))
p3 <- raincloud_1x1_repmes(
  data = dt1,
  colors = (c('#6F99AD','#BF5960')),
  fills = (c("#6F99AD","#BF5960")),
  size = 1.6,
  alpha = 0.5,
  line_color = 'grey',
  line_alpha = 0.8,
  align_clouds = FALSE
) +
  scale_x_continuous(breaks = c(0.6, 2.4),
                     labels = c("low", "high"),
                     limits = c(0, 3)) +
  labs(x = "", y = "ESTIMATEScore") +
  theme_classic() +
  # 调整x轴标签字体大小
  theme(axis.text.x = element_text(size = 14)) +  # 调大x轴标签的字体大小
  # 添加 p 值注释
  annotate("text", x = 1.5, y = max(dt1$y, na.rm = TRUE) * 1.05, 
           label = p_value_formatted, vjust = -1) +
  # 添加横线
  annotate("segment", x = 0.6, xend = 2.4, 
           y = max(dt1$y, na.rm = TRUE) * 1.03, 
           yend = max(dt1$y, na.rm = TRUE) * 1.03)

p3
ggsave("estimate-ESTIMATEScore.pdf", plot = p3, width = 6, height = 8)
#estimate-TumorPurity
data_high <- data[data$group == "high", "TumorPurity"]  # "high" 组
data_low <- data[data$group == "low", "TumorPurity"]    # "low" 组
# 使用 data_1x1 函数进行转换
dt1 <- data_1x1(
  array_1 = data_high,  # 分组1 ("high")
  array_2 = data_low,   # 分组2 ("low")
  jit_distance = 0.09,  # 抖动点间的距离
  jit_seed = 123        # 设置随机种子
)
# 查看转换后的数据
head(dt1)
#1.3 配对云雨图绘制：
head(dt1) #对两组中的相同id进行配对；
colnames(dt1)
wilcox_test_result <- wilcox.test(data_high, data_low)
p_value <- wilcox_test_result$p.value
# 格式化 p 值
p_value_formatted <- paste("p =", format(p_value, digits = 3))
p4 <- raincloud_1x1_repmes(
  data = dt1,
  colors = (c('#6F99AD','#BF5960')),
  fills = (c("#6F99AD","#BF5960")),
  size = 1.6,
  alpha = 0.5,
  line_color = 'grey',
  line_alpha = 0.8,
  align_clouds = FALSE
) +
  scale_x_continuous(breaks = c(0.6, 2.4),
                     labels = c("low", "high"),
                     limits = c(0, 3)) +
  labs(x = "", y = "TumorPurity") +
  theme_classic() +
  # 调整x轴标签字体大小
  theme(axis.text.x = element_text(size = 14)) +  # 调大x轴标签的字体大小
  # 添加 p 值注释
  annotate("text", x = 1.5, y = max(dt1$y, na.rm = TRUE) * 1.05, 
           label = p_value_formatted, vjust = -1) +
  # 添加横线
  annotate("segment", x = 0.6, xend = 2.4, 
           y = max(dt1$y, na.rm = TRUE) * 1.03, 
           yend = max(dt1$y, na.rm = TRUE) * 1.03)

p4
ggsave("estimate-TumorPurity.pdf", plot = p4, width = 6, height = 8)


####quantiseq####
quantiseq=deconvolute(tpm,'quantiseq')
write.csv(quantiseq,"quantiseq.csv",row.names = T)

#可视化
p5 <- quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(quantiseq))) +
  theme(axis.text.y = element_blank(),  # 移除 Y 轴的文字
        axis.ticks.y = element_blank())  # 移除 Y 轴的刻度
p5
# 用 ggsave 保存成 PNG 文件
ggsave("quantiseq_plot.pdf", plot = p5, width = 8, height = 12, units = "in", dpi = 300)


####xcell####
xcell=deconvolute(tpm,'xcell')
write.csv(xcell,"xcell.csv",row.names = T)
xc = as.data.frame(xcell)
rownames(xc) <- xc[,1]
xc <- xc[,-1]
#标准化
for(i in colnames(xc)) { 
  xc[,i] <- (xc[,i] - min(xc[,i])) / (max(xc[,i]) - min(xc[,i]))
}
####分组可视化
library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)
xc=t(xc)
xc=as.data.frame(xc)
xc$X <- rownames(xc)
# Replace dots with hyphens in the specific column
xc$X <- gsub("\\.", "-", xc$X)
fen<-read.csv("group.csv")
xc<-merge(fen,xc,"X")
mydatassgsea<-melt(
  xc,
  id.vars=c("X","group"),
  variable.name="immunecell",
  value.name="tpm"
)

ylabname <- paste("immunecell", "expression")
mydatassgsea=as.data.frame(mydatassgsea)
colnames(mydatassgsea) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value(秩和检验)
pvalues <- sapply(mydatassgsea$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydatassgsea, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydatassgsea$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', ''))
mydatassgsea<-mydatassgsea[,-1]
# 画box plot
xcell <- ggplot(mydatassgsea, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5, width = 0.6, position = position_dodge(width = 0.75)) + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") + 
  theme(axis.text.x = element_text(colour="black", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),  # 移除x轴标题
        axis.title.y = element_blank()) +  # 移除y轴标题
  geom_text(aes(x=gene, y=max(mydatassgsea$tpm) * 1.1, label = sigcode), data=pv, inherit.aes=F)

xcell
# 保存为PDF并指定宽度和高度
ggsave("xcell.pdf", plot = xcell, width = 15, height = 7)


####mcp_counter####
mcp_counter=deconvolute(tpm,'mcp_counter')

#timer
indications=c(rep("PAAD",ncol(tpm)))
timer=deconvolute(tpm,'timer',indications = indications)

####cibersort####
library(devtools)
#获得token
#https://blog.csdn.net/weixin_57360787/article/details/139073650
# usethis::create_github_token()
# Sys.setenv(GITHUB_PAT = "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")
# Sys.getenv("GITHUB_PAT")
# options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#perm置换次数=1000，QN分位数归一化=TRUE
TRUEresults <- cibersort(sig_matrix, "tpm.txt",perm = 1000,QN = T)
TRUEresults
results=TRUEresults[,1:22] #最后3列不需要
results=t(results)
results <- as.data.frame(results)  # 先将矩阵转换为数据框
results$cell_type <- rownames(results)  # 将 rownames 添加为新列
rownames(results) <- NULL  # 移除行名
# 使用 dplyr::select 将 cell_type 列设为第一列
library(dplyr)
results <- results %>% dplyr::select(cell_type, everything())
write.csv(results,'cibersort.csv',row.names = F)
cibersort <- results


####cibersort_abs####
library(IOBR)
ciber_res<-deconvo_tme(eset = tpm, method = "cibersort_abs", 
                       arrays = TRUE, perm = 1000)
ciber_res
# 删除 ciber_res 的后四列
cibersort_abs <- ciber_res %>%
  select(-((ncol(ciber_res)-3):ncol(ciber_res)))
cibersort_abs=t(cibersort_abs)
colnames(cibersort_abs) <- cibersort_abs[1,]
cibersort_abs=cibersort_abs[-1,]
cibersort_abs=as.data.frame(cibersort_abs)
# 将行名作为新列 "cell_type" 添加到数据框中
cibersort_abs$cell_type <- rownames(cibersort_abs)
# 重置行名为默认数值索引
rownames(cibersort_abs) <- NULL
# 将 cell_type 列移动到第一列
cibersort_abs <- cibersort_abs[, c("cell_type", setdiff(names(cibersort_abs), "cell_type"))]

save(tpm, quantiseq, xcell, epic, mcp_counter, timer, cibersort, cibersort_abs, file = "浸润结果.Rdata")

#清一下数据
rm(list = ls())
gc()
load("浸润结果.Rdata")
####合并结果-处理####
quantiseq$cell_type=paste(quantiseq$cell_type,'QUANTISEQ',sep = '_')
xcell$cell_type=paste(xcell$cell_type,'XCELL',sep = '_')
##这里注意xcell后三行是分数，看保存的时候想保存什么
#xcell1 <- xcell[-((nrow(xcell)-2):nrow(xcell)), ]
epic$cell_type=paste(epic$cell_type,'EPIC',sep = '_')
mcp_counter$cell_type=paste(mcp_counter$cell_type,'MCPCOUNTER',sep = '_')
timer$cell_type=paste(timer$cell_type,'TIMER',sep = '_')
cibersort$cell_type=paste(cibersort$cell_type,'CIBERSORT',sep = '_')
# # 替换所有 _CIBERSORT 为 _CIBERSORT-ABS
# cibersort_abs$cell_type <- gsub("_CIBERSORT", "_CIBERSORT-ABS", cibersort_abs$cell_type)
# cibersort_abs[] <- lapply(cibersort_abs, as.character)
# str(cibersort)
# # 保留第一列为字符类型，其他列转换为数值类型
# cibersort_abs <- cibersort_abs %>%
#   mutate(across(-cell_type, as.numeric))
# write.csv(cibersort_abs,"cibersort_abs.csv",row.names = F)
# cibersort <- read.csv("cibersort.csv")
# # 修改列名，从第二列开始将 "." 替换为 "-"
# colnames(cibersort)[2:ncol(cibersort)] <- gsub("\\.", "-", colnames(cibersort)[2:ncol(cibersort)])
####合并结果####
cell_ratio=dplyr::bind_rows(quantiseq,xcell,epic,mcp_counter,timer,cibersort,cibersort_abs)%>%  
  tibble::column_to_rownames("cell_type")
#write.csv(cell_ratio, "cell_ratio未标准化.csv", row.names = TRUE)

#####cell_ratio与riskscore的相关性####
riskscore <- read.csv("riskscore.csv",header = T,row.names = 1)
immucell=t(cell_ratio)%>%  
  as.data.frame() %>%  
  tibble::add_column(ID=stringr::str_sub(rownames(.),1,12))%>%  
  dplyr::filter(!duplicated(ID)) %>% # remove duplicated samples randomly  
  tibble::remove_rownames(.) %>% 
  tibble::column_to_rownames("ID")%>%  
  dplyr::filter(rownames(.) %in% rownames(riskscore))
identical(sort(rownames(immucell)), sort(rownames(riskscore)))
#与riskscore相关性分析
library(psych)
cor <-psych::corr.test(riskscore$riskScore, immucell, method = 'spearman',adjust="none")
cmt <-t(as.data.frame(cor$r))%>%  
  as.data.frame()%>%  
  tibble::add_column(rownames(.))%>%  
  tibble::remove_rownames(.)
colnames(cmt)=c("Correlation Coefficient","Immune cell")
#为了画图好看，需要固定Immune cell顺序
df=cmt
# 定义 Software 对应的颜色映射
software_colors <- c('QUANTISEQ' = "#f27767", 'XCELL' = "#bc9826", 
                     'EPIC' = "#53b449", 'MCPCOUNTER' = "#23b892", 
                     'TIMER' = "#1cb6ec", 'CIBERSORT' = "#d269a4",
                     'CIBERSORT-ABS' = "lightpink2")


df$Software <- sub(".*_", "", df$`Immune cell`)
df$Immune_cell_label <- factor(df$`Immune cell`, levels = unique(df$`Immune cell`))

# 创建 y_cols 映射，确保基于 Software 后缀进行颜色匹配
y_cols <- software_colors[df$Software]

p7 <- ggplot(df, aes(x = `Correlation Coefficient`, y = Immune_cell_label)) +  
  geom_point(aes(color = Software), shape = 16, size = 5) +  
  scale_color_manual(values = software_colors) +  # 手动设置 Software 对应的颜色
  theme_bw() +  # 使用白色背景，并保留网格
  theme(axis.text.y = element_text(size = 10, color = y_cols),  # 根据 Software 分配的 y 轴颜色
        axis.text.x = element_text(size = 10, color = 'black'),        
        axis.title.x = element_text(size = 12),        
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA),  # 设置黑色边框
        panel.grid.major = element_line(color = "gray80"),  # 设置主网格线颜色
        panel.grid.minor = element_blank()) +  # 移除次要网格线
  labs(y = "Immune Cell", x = "Correlation Coefficient")

p7
# 保存图片
ggsave("免疫浸润与riskscore相关性.pdf", plot = p7, width = 10, height = 20, dpi = 300)

#####分组cibersort####
group <- read.csv("group.csv")
a=cibersort
head(a)
# 转置数据框
a <- t(a)
# 将第一行设置为列名
colnames(a) <- a[1, ]
# 删除第一行（因为已经用作列名）
a <- a[-1, ]
# 将行名转换为一列并命名为 'ID'
a <- data.frame(X = rownames(a), a, row.names = NULL)
# 将结果转化为数据框，以确保所有列的数据类型一致
a <- as.data.frame(a)
head(a)
# 假设 a 的 X 列包含要处理的字符串
a$X <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", a$X)
fen<-read.csv("group.csv")
a<-merge(fen,a,"X")
library(reshape2)

mydata1<-melt(
  a,
  id.vars=c("X","group"),
  variable.name="immunecell",
  value.name="tpm"
)
# 去掉后缀_CIBERSORT
mydata1$immunecell <- gsub("_CIBERSORT", "", mydata1$immunecell)
library(ggpubr)
library(ggplot2)

ylabname <- paste("immunecell", "expression")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value(秩和检验)
pvalues <- sapply(mydata1$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydata1, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydata1$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', ''))
mydata1<-mydata1[,-1]
# 画box plot
mydata1$tpm <- as.numeric(as.character(mydata1$tpm))
p.box <- ggplot(mydata1, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5) + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  
  # 移除轴标题
  theme(axis.text.x = element_text(colour="black", size = 9,
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),  # 移除x轴标题
        axis.title.y = element_blank()) +  # 移除y轴标题
  
  geom_text(aes(x=gene, y=max(mydata1$tpm) * 1.1,
                label = pv$sigcode),
            data=pv, 
            inherit.aes=F)

p.box

####cibersort套图####
cib_exist=t(cibersort)
# 假设 cibersort 是你的数据框
colnames(cib_exist) <- cib_exist[1, ]  # 将第一行设为列名
cib_exist <- cib_exist[-1, ]  # 删除第一行
cib_exist = cib_exist[,which(apply(cib_exist, 2, var)!=0)]#删除0方差列的细胞
cib_exist <- as.data.frame(cib_exist)
cib_exist$X <- rownames(cib_exist)
# 4. 如果需要，可以重新设置行名为 NULL
rownames(cib_exist) <- NULL
group <- read.csv("group.csv")
meltdf_cib_exist <- merge(group,cib_exist,by = "X")

library(reshape2)
# 将数据转换为长格式，X 列改为 sample，保留 group 列
meltdf_cib_exist <- melt(meltdf_cib_exist, id.vars = c("X", "group"), 
                variable.name = "variable", value.name = "value")
# 将 X 列重命名为 sample
names(meltdf_cib_exist)[names(meltdf_cib_exist) == "X"] <- "sample"
# 使用 colnames 修改列名
#colnames(meltdf_cib_exist)[colnames(meltdf_cib_exist) == "group"] <- "status1"
# 使用 sub 函数删除 variable 列中每个值的最后的 _CIBERSORT
meltdf_cib_exist$variable <- sub("_CIBERSORT$", "", meltdf_cib_exist$variable)
plot_tme = function(meltdf) {
  
  p_stack_df = ggplot(meltdf, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity") + 
    theme_test() + 
    scale_fill_manual(values = cols) +
    theme(legend.position = 'bottom') +
    theme(plot.title = element_text(size = rel(2), hjust = 0.5), 
          axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +  # Set y limits to 1
    facet_nested(. ~ group, drop = TRUE, scale = "free", space = "free", switch = "x",
                 strip = strip_nested(background_x = elem_list_rect(fill = c('firebrick','steelblue')), 
                                      by_layer_x = FALSE))
  
  p_box_df = ggplot(meltdf, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(width = 0.4, lwd = 0.2, color = 'black', outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    theme(legend.position = 'none') +
    scale_fill_manual(values = cols)
  
  p_box_comp = ggplot(meltdf, aes(x = variable, y = value, fill = group, color = group)) +
    geom_boxplot(width = 0.5, lwd = 0.2, color = 'black', outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    theme_bw() +
    theme(legend.position = 'top') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    stat_compare_means(aes(group = group, label = after_stat(p.signif)), hide.ns = FALSE,
                       method = "t.test") +
    scale_fill_manual(values = c('firebrick','steelblue'))
  
  p_heat = ggplot(meltdf, aes(x = sample, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = c('white', 'steelblue', 'firebrick')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill = NA, color = 'black')) +
    labs(fill = "Z score", x = NULL, y = NULL)
  
  plot_comb = (p_stack_df) / (p_box_df | p_box_comp) / p_heat
  return(plot_comb)
}


library(IOBR)
library(ggh4x)
cols <- IOBR::palettes(category = "random", palette = 4,
                       show_col = F,
                       show_message = F)
plot_exist_cib = plot_tme(meltdf_cib_exist)

print(plot_exist_cib)

ggsave("cibersort套图.pdf", plot = plot_exist_cib, width = 14, height = 14, units = "in")

#热图
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(gplots)
cell_ratio <- read.csv("所有免疫浸润结果.csv",row.names = 1)
data_total <- t(cell_ratio)
data_total <- as.data.frame(data_total)
# 替换行名中的.为-，并保留到-符号的前三位
new_rownames <- sapply(strsplit(gsub("\\.", "-", rownames(data_total)), "-"), function(x) {
  if (length(x) > 3) {
    return(paste(x[1:3], collapse = "-"))
  } else {
    return(paste(x, collapse = "-"))
  }
})
# 确保行名唯一
new_rownames <- make.unique(new_rownames)
# 设置新的行名
rownames(data_total) <- new_rownames
head(data_total)[1:5,1:5]
data_total=t(data_total)
#scale默认对列进行均一化，这里我们是要对不同算法结果，也就是行，进行均一化，因此需要转置
tme_combine2 <- scale(data_total,
                      center = TRUE, scale = TRUE)

pheatmap(tme_combine2,
         border_color = NA,
         color = colorRampPalette(c("darkblue", "white","darkred"))(100),
         cluster_rows = FALSE, # 行不聚类 
         cluster_cols = FALSE, # 列不聚类
         show_rownames = T, # 不显示行名
         show_colnames = FALSE)



####ssgsea-29####
library("GSVA")
library("limma")
library("GSEABase")
library("xlsx")
library("pheatmap")
library("ggpubr")
library("edgeR")
library("limma")
library("data.table")
library("TCGAbiolinks")
library("clusterProfiler")
library("org.Hs.eg.db")
library("tidyr")
library(survival)
library(survminer)
library(cowplot)
exp_data1 = log2(tpm+1)
# 将 tpm 保存为 txt 格式
write.table(exp_data1, file = "tpm.txt", sep = "\t", col.names = TRUE, quote = FALSE)

geneSet=getGmt("ssgsea/immune_geneset.txt", geneIdType=SymbolIdentifier())#读取我们所提供的免疫细胞基因文件
expr_paad.new <- exp_data1

####以上注意是改了患者格式
# 计算富集得分
class(expr_paad.new)
#[1] "matrix" "array" 
expr_paad.new <- as.matrix(expr_paad.new)
#gsvapar <- gsvaParam(expr_paad.new, geneSet) 
gsvapar <- ssgseaParam(
  exprData = expr_paad.new,
  geneSets = geneSet,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE)
#开始计算
gsva_data <- gsva(gsvapar)
#标准化
for(i in colnames(gsva_data)) { # 修正了 'conames' 为 'colnames'
  gsva_data[,i] <- (gsva_data[,i] - min(gsva_data[,i])) / (max(gsva_data[,i]) - min(gsva_data[,i]))
}

gsva_data=rbind(id=colnames(gsva_data),gsva_data)
gsva_data[1:4,1:4]  # 结果展示
write.table(gsva_data,file="ssGSEA_Score.txt",sep="\t",quote=F,col.names=F)#此处的输出文件即为ssGSEA富集得分文件

#热图绘图展示更直观
#添加样本注释(tumor和normal分组）：
#group_list <- ifelse(str_sub(colnames(res),14,15) < 10,"tumor","normal")
group_list<-read.csv("group.csv")
annotation <- data.frame(group_list)
annotation<-annotation[order(annotation$group),]
#annotation<-annotation[order(annotation$X),]
rownames(annotation) <- annotation$X
annotation$X<-NULL
gsva_data <- gsva_data[, rownames(annotation)]
head(annotation)
#热图
# 定义 group 的颜色
group_colors <- list(group = c("low" = "blue", "high" = "red"))
# 绘制热图
pheatmap(gsva_data,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation,
         annotation_colors = group_colors,
         fontsize = 10)

#箱线图
#宽数据转换为长数据(ggplot2绘图格式)：
library(reshape2)
dt <- as.data.frame(gsva_data)
#宽数据转换为长数据(ggplot2绘图格式)：
dt <- dt %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)
#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)
#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)
dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))
length(unique(dtt$cell_type))

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")
#配色挑选
# 加载paletteer包
library(paletteer)
d_palettes <- palettes_d_names #查看离散型色板(paletteer包)
col <- paletteer_d("khroma::discreterainbow",n=29)
#重新绘图(代码相同)：
p8 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p8

####分组可视化
library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)
ssgsea <- gsva_data
ssgsea=t(ssgsea)
ssgsea=as.data.frame(ssgsea)
ssgsea$X <- rownames(ssgsea)
# Replace dots with hyphens in the specific column
ssgsea$X <- gsub("\\.", "-", ssgsea$X)
fen<-read.csv("group.csv")
ssgsea<-merge(fen,ssgsea,"X")
mydatassgsea<-melt(
  ssgsea,
  id.vars=c("X","group"),
  variable.name="immunecell",
  value.name="tpm"
)
library(ggpubr)
library(ggplot2)

ylabname <- paste("immunecell", "expression")
mydatassgsea=as.data.frame(mydatassgsea)
colnames(mydatassgsea) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value(秩和检验)
pvalues <- sapply(mydatassgsea$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydatassgsea, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydatassgsea$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', ''))
mydatassgsea<-mydatassgsea[,-1]
# 画box plot
p.box <- ggplot(mydatassgsea, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5, width = 0.6, position = position_dodge(width = 0.75)) + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") + 
  theme(axis.text.x = element_text(colour="black", size = 11, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),  # 移除x轴标题
        axis.title.y = element_blank()) +  # 移除y轴标题
  geom_text(aes(x=gene, y=max(mydatassgsea$tpm) * 1.1, label = sigcode), data=pv, inherit.aes=F)

p.box

####ssgsea-28####
library(GSVA)
exp_data1 = log2(tpm+1)
load("ssgsea/TISIDB肿瘤浸润淋巴细胞基因集.rdata")
geneSet=tisidb_cell
expr_paad.new <- exp_data1

####以上注意是改了患者格式
# 计算富集得分
class(expr_paad.new)
#[1] "matrix" "array" 
expr_paad.new <- as.matrix(expr_paad.new)
#gsvapar <- gsvaParam(expr_paad.new, geneSet) 
gsvapar <- ssgseaParam(
  exprData = expr_paad.new,
  geneSets = geneSet,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE)
#开始计算
gsva_data <- gsva(gsvapar)
#标准化
for(i in colnames(gsva_data)) { # 修正了 'conames' 为 'colnames'
  gsva_data[,i] <- (gsva_data[,i] - min(gsva_data[,i])) / (max(gsva_data[,i]) - min(gsva_data[,i]))
}

gsva_data=rbind(id=colnames(gsva_data),gsva_data)
gsva_data[1:4,1:4]  # 结果展示
write.table(gsva_data,file="ssGSEA28-Score.txt",sep="\t",quote=F,col.names=F)#此处的输出文件即为ssGSEA富集得分文件

#热图绘图展示更直观
#添加样本注释(tumor和normal分组）：
#group_list <- ifelse(str_sub(colnames(res),14,15) < 10,"tumor","normal")
group_list<-read.csv("group.csv")
annotation <- data.frame(group_list)
annotation<-annotation[order(annotation$group),]
#annotation<-annotation[order(annotation$X),]
rownames(annotation) <- annotation$X
annotation$X<-NULL
gsva_data <- gsva_data[, rownames(annotation)]
head(annotation)
#热图
# 定义 group 的颜色
group_colors <- list(group = c("low" = "blue", "high" = "red"))
# 绘制热图
pheatmap(gsva_data,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation,
         annotation_colors = group_colors,
         fontsize = 10)

#箱线图
#宽数据转换为长数据(ggplot2绘图格式)：
library(reshape2)
dt <- as.data.frame(gsva_data)
#宽数据转换为长数据(ggplot2绘图格式)：
dt <- dt %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)
#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)
#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)
dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))
length(unique(dtt$cell_type))

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")
#配色挑选
# 使用新的色板 khroma::discreterainbow，n=28
col <- paletteer_d("ggsci::default_igv", n = 28)
# 绘制带有新色板的ggplot
ggplot(dtt, aes(cell_type, proportion)) +
  geom_violin(width = 2.0, aes(color = cell_type)) +
  geom_boxplot(width = 0.2, fill = "black") + 
  theme_bw() + 
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_color_manual(values = col, name = NULL)
####分组可视化
ssgsea <- gsva_data
ssgsea=t(ssgsea)
ssgsea=as.data.frame(ssgsea)
ssgsea$X <- rownames(ssgsea)
# Replace dots with hyphens in the specific column
ssgsea$X <- gsub("\\.", "-", ssgsea$X)
fen<-read.csv("group.csv")
ssgsea<-merge(fen,ssgsea,"X")
mydatassgsea<-melt(
  ssgsea,
  id.vars=c("X","group"),
  variable.name="immunecell",
  value.name="tpm"
)
library(ggpubr)
library(ggplot2)

ylabname <- paste("immunecell", "expression")
mydatassgsea=as.data.frame(mydatassgsea)
colnames(mydatassgsea) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value(秩和检验)
pvalues <- sapply(mydatassgsea$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydatassgsea, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydatassgsea$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', ''))
mydatassgsea<-mydatassgsea[,-1]
# 画box plot
p.box <- ggplot(mydatassgsea, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5, width = 0.6, position = position_dodge(width = 0.75)) + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set2") + 
  scale_color_brewer(palette = "Set2") + 
  theme(axis.text.x = element_text(colour="black", size = 9, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),  # 移除x轴标题
        axis.title.y = element_blank()) +  # 移除y轴标题
  geom_text(aes(x=gene, y=max(mydatassgsea$tpm) * 1.1, label = sigcode), data=pv, inherit.aes=F)

p.box

####单基因mantel热图####
library(linkET)
library(ggplot2)
library(dplyr)
# 假设 tpm 是一个以基因为行名的数据框
gene <- subset(tpm, rownames(tpm) %in% c("PML", "RPN1"))
colnames(gene) <- sub("^(([^-]*-){2}[^-]*).*", "\\1", colnames(gene))
gene=t(gene)
gene=as.data.frame(gene)
riskscore <- read.csv("riskscore.csv")
rownames(riskscore) <- riskscore$X
riskscore$X <- NULL
# 将 gene 和 riskscore 的行名转换为一列
gene$SampleID <- rownames(gene)
riskscore$SampleID <- rownames(riskscore)

# 按照 SampleID 合并两个数据框
cli <- merge(gene, riskscore, by = "SampleID")
# 如果需要，合并后可以将 SampleID 列重新设置为行名
rownames(cli) <- cli$SampleID
cli$SampleID <- NULL  # 删除 SampleID 列
cibersort <- read.csv("cibersort.csv")
rownames(cibersort) <- cibersort$cell_type
cibersort$cell_type <- NULL
cibersort=t(cibersort)
cibersort=as.data.frame(cibersort)
# 第一步：将行名中的 '.' 替换为 '-'
rownames(cibersort) <- gsub("\\.", "-", rownames(cibersort))
# 找到两个数据框中相同的样本（行名）
common_samples <- intersect(rownames(cibersort), rownames(cli))

# 只保留相同样本的行
cibersort <- cibersort[common_samples, , drop = FALSE]
cli <- cli[common_samples, , drop = FALSE]
mantel <- mantel_test(cli, cibersort,
                      spec_select = list(RPN1 = 1,
                                         PML = 1,                                         
                                         riskScore = 3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))) %>%
  mutate(
    size = cut(r, breaks = c(-Inf, 0.5, 0.6, 0.7, Inf),
               labels = c("<0.5", "0.5~0.6", "0.6~0.7", ">=0.7"),
               right = FALSE),
    color = ifelse(p < 0.05, r, NA)
  )
#> `mantel_test()` using 'bray' dist method for 'spec'.
#> #> `mantel_test()` using 'euclidean' dist method for 'env'.
library(ggnewscale)
library(RColorBrewer)
#这个是小方形带星号的
# p10 <- qcorrplot(correlate(cibersort),           
#                  grid_col = "grey50",          
#                  grid_size = 0.5,          
#                  type = "upper",           
#                  diag = FALSE) +  
#   geom_square() +  
#   geom_mark(size = 4,            
#             only_mark = T,            
#             sig_level = c(0.05, 0.01, 0.001),            
#             sig_thres = 0.05,            
#             colour = 'black') +  
#   geom_couple(data = mantel,              
#               aes(color = pd, size = rd),                
#               label.size = 3.88,              
#               label.family = "",              
#               label.fontface = 1,              
#               nudge_x = 0.2,              
#               curvature = nice_curvature(by = "from")) +      
#   scale_fill_gradient2(low = "#6A5ACD", 
#                        mid = "#FFDDC1", 
#                        high = "violetred3", 
#                        midpoint = 0, 
#                        limits = c(-1, 1)) +  
#   scale_size_manual(values = c(0.5, 1.5, 3)) +  
#   scale_color_manual(values = color_pal(4, alpha = 0.6)) +    
#   guides(size = guide_legend(title = "Mantel's r",                                                            
#                              order = 2,                                                            
#                              keyheight = unit(0.5, "cm")),                    
#          colour = guide_legend(title = "Mantel's p",                                                                 order = 1,                                                                
#                                keyheight = unit(0.5, "cm")),                   
#          fill = guide_colorbar(title = "Pearson's r",                                
#                                keyheight = unit(2.2, "cm"),                               
#                                keywidth = unit(0.5, "cm"),                               
#                                order = 3)) +   
#   theme(legend.box.spacing = unit(0, "pt"))
# 
# p10

p10 <- qcorrplot(correlate(cibersort), method = "spearman", type = "upper", diag = FALSE, grid_col = NA) +
  geom_tile(aes(fill = r), color = "snow2", size = 0.2) +  # 使用 geom_tile 生成方形块，去掉圆点，保留细黑色边框
  #geom_point(aes(size = abs(r), fill = r), shape = 21) + 
  scale_fill_gradient2(low = "#6A5ACD", mid = "#FFDDC1", high = "violetred3", midpoint = 0, limits = c(-1, 1)) +  # 颜色渐变
  new_scale("size") +
  geom_couple(data = mantel, aes(color = pd, size = rd,
                                 linetype = ifelse(r > 0, "positive", "negative")), label.size = 3.88, label.family = "",
              label.fontface = 1, nudge_x = 0.2, curvature = nice_curvature(by = "from")) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#f46d43","#762a83","#CCCCCC99")) +
  scale_linetype_manual(values = c("positive" = "solid", "negative" = "dashed"), 
                        labels = c("positive" = "Positive", "negative" = "Negative")) +  
  guides(size = guide_legend(title = "Mantel's r", order = 2, keyheight = unit(0.5, "cm")),
         colour = guide_legend(title = "Mantel's p", order = 1, keyheight = unit(0.5, "cm")),
         fill = guide_colorbar(title = "Spearman's r", keyheight = unit(2.2, "cm"), keywidth = unit(0.5, "cm"), order = 3)) +
  theme(legend.box.spacing = unit(0, "pt"))
p10
ggsave(p10, file = "cibersort.mantel.pdf", width = 10, height = 10)
#单基因与ssgseamantel
ssgsea <- read.table("ssGSEA_Score.txt",header = T,row.names = 1)
ssgsea = t(ssgsea)
rownames(ssgsea) <- gsub("\\.", "-", rownames(ssgsea))
ssgsea=as.data.frame(ssgsea)
# 找到两个数据框中相同的样本（行名）
common_samples <- intersect(rownames(ssgsea), rownames(cli))

# 只保留相同样本的行
ssgsea <- ssgsea[common_samples, , drop = FALSE]
cli <- cli[common_samples, , drop = FALSE]
mantel <- mantel_test(cli, ssgsea, 
                      spec_select = list(RPN1 = 1,
                                         PML = 1,
                                         riskScore = 3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))) %>%
  mutate(
    size = cut(r, breaks = c(-Inf, 0.5, 0.6, 0.7, Inf),
               labels = c("<0.5", "0.5~0.6", "0.6~0.7", ">=0.7"),
               right = FALSE),
    color = ifelse(p < 0.05, r, NA)
  )
#一个复现，我不喜欢
# p11<- qcorrplot(correlate(ssgsea), type = "lower", diag=FALSE, grid_col = NA) +      
#   geom_point(shape=21, size=6, fill = NA, stroke = 0.35, color = "black") +      
#   geom_point(aes(size=abs(r), fill=r), shape=21, stroke = 0.35, color = "black") +      
#   scale_size(range = c(1, 8), guide = "none") + new_scale("size") +      
#   geom_couple(data = mantel, aes(color = pd, size = rd), label.size = 3.88, label.family = "", 
#               label.fontface = 1, nudge_x = 0.2, curvature = nice_curvature(by = "from")) +      
#   scale_fill_gradientn(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4), 
#                        colors = rev(brewer.pal(11, "Spectral"))) +  
#   scale_size_manual(values = c(0.5, 1.5, 3)) +  
#   scale_color_manual(values = color_pal(4, alpha = 0.6)) +    
#   guides(size = guide_legend(title = "Mantel's r", order = 2, keyheight = unit(0.5, "cm")), 
#          colour = guide_legend(title = "Mantel's p", order = 1, keyheight = unit(0.5, "cm")), 
#          fill = guide_colorbar(title = "Pearson's r", keyheight = unit(2.2, "cm"), keywidth = unit(0.5, "cm"), order = 3)) +   
#   theme(legend.box.spacing = unit(0, "pt"))
# p11
p11 <- qcorrplot(correlate(ssgsea), method = "spearman", type = "lower", diag = FALSE, grid_col = NA) +
  geom_tile(aes(fill = r), color = "snow2", size = 0.2) +  # 使用 geom_tile 生成方形块，去掉圆点，保留细黑色边框
  #geom_point(aes(size = abs(r), fill = r), shape = 21) + 
  scale_fill_gradient2(low = "#6A5ACD", mid = "#FFDDC1", high = "violetred3", midpoint = 0, limits = c(-1, 1)) +  # 颜色渐变
  new_scale("size") +
  geom_couple(data = mantel, aes(color = pd, size = rd,
                                 linetype = ifelse(r > 0, "positive", "negative")), label.size = 3.88, label.family = "",
              label.fontface = 1, nudge_x = 0.2, curvature = nice_curvature(by = "from")) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#f46d43","#762a83","#CCCCCC99")) +
  scale_linetype_manual(values = c("positive" = "solid", "negative" = "dashed"), 
                        labels = c("positive" = "Positive", "negative" = "Negative")) +  
  guides(size = guide_legend(title = "Mantel's r", order = 2, keyheight = unit(0.5, "cm")),
         colour = guide_legend(title = "Mantel's p", order = 1, keyheight = unit(0.5, "cm")),
         fill = guide_colorbar(title = "Spearman's r", keyheight = unit(2.2, "cm"), keywidth = unit(0.5, "cm"), order = 3)) +
  theme(legend.box.spacing = unit(0, "pt"))
p11
ggsave(p11, file = "ssgsea.mantel.pdf", width = 10, height = 10)
#将11和12合并
# library(cowplot)
# ggdraw()+
#   draw_plot(p10,scale = 0.6,x=0.2,y=0.11)+
#   draw_plot(p11,scale = 0.6,x=-0.21,y=-0.14)

####panel热图####
#对cell_ratio进行标准化
range(cell_ratio)
cell_ratio2 <- scale(cell_ratio,
                     center = TRUE, scale = TRUE)
group <- read.csv("group.csv",row.names = 1)
cell_ratio1 <- as.data.frame(t(cell_ratio2))
# Step 2: 将行名转换为一列以便进行合并
cell_ratio1$ID <- rownames(cell_ratio1)
group$ID <- rownames(group)
# Step 3: 按照 'ID' 列进行合并
tme_mid <- merge(cell_ratio1, group, by = "ID", all = TRUE)
# Step 4: 如果需要可以将 'ID' 列重新设置为行名
rownames(tme_mid) <- tme_mid$ID
tme_mid$ID <- NULL
tme_mid <- unique(tme_mid)
result_tme_comb = as.data.frame(t(tme_mid))
# 将 'group' 行设为第一行
result_tme_comb <- result_tme_comb[c("group", setdiff(rownames(result_tme_comb), "group")), ]

heat_data=result_tme_comb[-c(1),]

library(dplyr)

gene_group <- data.frame(cell = rownames(heat_data)) %>%
  mutate(cell_group = tolower(sub(".*_", "", cell)))

for(i in names(heat_data)){heat_data[,i]<-as.numeric(heat_data[,i])}

heat_data = heat_data%>%
  mutate(cell_group=gene_group$cell_group)%>%
  dplyr::select(cell_group,everything())

str(heat_data)
score <- read.csv("riskscore.csv")
cli <- read.csv("数字版clinical.csv")
cli <- merge(cli,score,by = "X")
#cli <- read.csv("预后数据.csv")
colnames(cli)[colnames(cli) == "X"] <- "sample"
table(cli$N)
cli[] <- lapply(cli, function(x) {
  x <- gsub("X", "NA", x)        # 替换含有X的值为"NA"
  x[is.na(x)] <- "NA"            # 将缺失值NA替换为"NA"
  return(x)
})
cli <- cli %>%
  # 修改 status 列为 Alive 和 Dead
  mutate(status = factor(status, levels = c('0', '1'), labels = c('Alive', 'Dead'))) %>%
  
  # 修改 Age 列为 <=65 和 >65
  mutate(Age = ifelse(Age <= 65, "<=65", ">65")) %>%
  
  # 修改 Gender 列为 Female 和 Male (假设 0 为 Female, 1 为 Male，若不确定可以调整)
  mutate(Gender = ifelse(Gender == 0, "Female", "Male")) 
  
  # 修改 Grade 列，加上 G
  #mutate(Grade = paste0("G", Grade)) %>%
  
  # 修改 Stage 列，加上 StageI 等
  #mutate(Stage = paste0("Stage", Stage)) %>%
  
  # 修改 T 列，将 2, 3, 4 加上 T
  #mutate(T = ifelse(T %in% c(2, 3, 4), paste0("T", T), T)) %>%
  
  # 修改 M 列，将 0 和 1 加上 M
  #mutate(M = paste0("M", M)) %>%
  
  # 修改 N 列，将 0 和 1 加上 N
  #mutate(N = paste0("N", N))
str(cli)
cli$time <- as.numeric(cli$time)
result_tidy<-
  heat_data %>%
  #这里将行名命名为新的列
  as_tibble(rownames = 'cell') %>%
  #这里选择去掉哪些列来进行归一化
  mutate_at(vars(-colnames(gene_group)), scale) %>%
  #宽数据转长数据
  pivot_longer(cols = -colnames(gene_group), 
               names_to = "sample", 
               values_to = "Value")%>%
  left_join(cli,by='sample')
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
pal_anno1<-pal_nejm()(7)

pal_plot1<-
  #调色板指定范围
  #seq范围越小颜色越深
  circlize::colorRamp2(seq(-3, 3, length.out = 3), 
                       # rev(RColorBrewer::brewer.pal(9, "RdBu"))
                       c('steelblue','white','firebrick')
  )
result_tidy <- na.omit(result_tidy)
result_tidy1=result_tidy
table(result_tidy$cell_group)
result_tidy1 = result_tidy %>%
  mutate(cell = gsub('_QUANTISEQ', '', cell)) %>%
  mutate(cell = gsub('_XCELL', ' ', cell)) %>%
  mutate(cell = gsub('_EPIC', '  ', cell)) %>%
  mutate(cell = gsub('_MCPCOUNTER', '   ', cell)) %>%
  mutate(cell = gsub('_TIMER', '    ', cell)) %>%
  mutate(cell = gsub('_CIBERSORT-ABS', '     ', cell)) %>%  # 移到这里
  mutate(cell = gsub('_CIBERSORT', '      ', cell))
g<-read.csv("group.csv")
library(dplyr)
# 将文件 g 和 result_tidy1 按照 g 文件的 X 列和 result_tidy1 的 sample 列合并
result_tidy1 <- result_tidy1 %>%
  left_join(g %>% select(X, group), by = c("sample" = "X"))
#获取 result_tidy1 中的样本名称
samples_in_result <- unique(result_tidy1$sample)
#过滤 tme_mid 中不存在于 result_tidy1 中的样本
tme_mid <- tme_mid[rownames(tme_mid) %in% samples_in_result, ]

# 确保按照 result_tidy1 数据中的顺序显示 cell
result_tidy1$cell <- factor(result_tidy1$cell, levels = unique(result_tidy1$cell))

plot_heatmap1 <- result_tidy1 %>%
  group_by(`cell_group`, `group`) %>%
  heatmap(
    .row = cell, 
    .column = sample, 
    .value = Value,
    .scale = "row",
    palette_grouping = list(pal_anno1, c("#66C2A5", "steelblue")),
    palette_value = pal_plot1,
    cluster_rows = F,
    border = T,
    row_names_gp = gpar(fontsize = 8),
    show_column_names = F,   # 隐藏样本名称
    show_row_names = T,      # 显示 cell 名称
    column_order = rownames(tme_mid),
    column_title = NULL,     # 移除样本标签
    row_title = NULL         # 移除 cell 标签
  )

plot_heatmap1

col1<-c('#fee5d9','#fcae91','#fb6a4a','#cb181d')
col2<-c('#f2f0f7','#cbc9e2','#9e9ac8','#6a51a3')
# 确保 riskScore 列是数值型
result_tidy1$riskScore <- as.numeric(result_tidy1$riskScore)
plot_heatmap = plot_heatmap1%>%
  add_tile(`Age`,palette = col1)%>%
  add_tile(`Gender`,palette = c('pink','skyblue'))%>%
  add_tile(`Grade`,palette = col1)%>%
  add_tile(`Stage`,palette = col1)%>%
  add_tile(`T`,palette = col2)%>%
  add_tile(`M`,palette = col2)%>%
  add_tile(`N`,palette = col2)%>%
  add_tile(`status`,palette=c('slategray2','pink3'))%>%
  add_line(`riskScore`)%>%
  add_point('time') #gpar（）设置颜色

plot_heatmap
pdf('免疫热图含NA.pdf',width = 13,height = 15)
print(plot_heatmap)
dev.off()

####单基因分组免疫浸润####
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)
data = log2(tpm+1)
data=t(data)
data=as.data.frame(data)
names(data)
# 保留基因 列
data <- data[, "RPN1", drop = FALSE]
data=avereps(data)
#记得基因改
gene=colnames(data)[1]
data=as.data.frame(data)
data$gene=ifelse(data[,gene]>median(data[,gene]), "high", "low")
gsva_data <- read.table("ssGSEA_Score.txt", header = TRUE, sep = "\t", row.names = 1)
immune=t(gsva_data)
rownames(immune) <- gsub("\\.", "-", rownames(immune))
#cibersort
gsva_data <- cibersort
row.names(gsva_data) <- gsva_data$cell_type
gsva_data$cell_type <- NULL
# 删除最后一个 "_" 符号及其后的内容
rownames(gsva_data) <- sub("_[^_]*$", "", rownames(gsva_data))
immune=t(gsva_data)
#cibersort结

# 确保数据已被正确读取并转换为数据框
data <- as.data.frame(data)
immune <- as.data.frame(immune)
# 添加行名作为新列
data$rownames <- rownames(data)
immune$rownames <- rownames(immune)
# 使用 merge 函数进行合并
rt <- merge(immune, data, by = "rownames")

# 将行名列重新设置为行名，并移除行名列
rownames(rt) <- rt$rownames
rt$rownames <- NULL
data=rt[,-(ncol(rt)-1)] #将最后一列high low删除
data=melt(data,id.vars=c("gene")) #整理成了三列 第一列highlow,细胞，value
colnames(data)=c("gene", "Immune", "Expression")
group=levels(factor(data$gene))
data$gene=factor(data$gene, levels=c("low","high"))
#想绘制相关性散点图，接1365行
#bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=c("#0066FF","#D20A13")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(
    aes(group = gene),
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c( "***", "**", "*", "")
    ),
    label = "p.signif"
  )
#保存哦
pdf(file="NLNssgsea差异.pdf", width=15, height=7)
print(boxplot)
dev.off()
ylabname <- paste("Immune", "Expression")
data=as.data.frame(data)
# 计算p value(秩和检验)
pvalues <- sapply(data$Immune, function(x) {
  res <- wilcox.test(as.numeric(Expression) ~ gene, data = subset(data, Immune == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = data$Immune, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', ''))
# 自定义颜色映射
custom_colors <- c('high' = '#B43970', 'low' = '#66BC98')
# 绘制 boxplot
p.box <- ggplot(data, aes(x = Immune, y = Expression, color = gene, fill = gene)) +
  geom_boxplot(alpha = 0.6, width = 0.7, position = position_dodge(width = 0.75)) + 
  theme_classic() + 
  scale_fill_manual(values = custom_colors) +  # 修改填充颜色
  scale_color_manual(values = custom_colors) +  # 修改边框颜色
  theme(
    axis.text.x = element_text(colour = "black", size = 11, angle = 90, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),  # 移除 x 轴标题
    axis.title.y = element_blank()   # 移除 y 轴标题
  ) +
  geom_text(aes(x = gene, y = max(data$Expression) * 1.1, label = sigcode), 
            data = pv, inherit.aes = F)
# 显示图形
pdf(file="RPN1-cibersort差异.pdf", width=15, height=7)
p.box
dev.off()


####单基因与cell_ratio相关性####
####cell_ratio与riskscore的相关性####
data = log2(tpm+1)
data=t(data)
data=as.data.frame(data)
names(data)
# 保留基因 列
data <- data[, "RPN1", drop = FALSE]
immucell=t(cell_ratio)%>%  
  as.data.frame() %>%  
  tibble::add_column(ID=stringr::str_sub(rownames(.),1,12))%>%  
  dplyr::filter(!duplicated(ID)) %>% # remove duplicated samples randomly  
  tibble::remove_rownames(.) %>% 
  tibble::column_to_rownames("ID")%>%  
  dplyr::filter(rownames(.) %in% rownames(data))
identical(sort(rownames(immucell)), sort(rownames(data)))
#与基因相关性分析
library(psych)
cor <-psych::corr.test(data$RPN1, immucell, method = 'spearman',adjust="none")
cmt <-t(as.data.frame(cor$r))%>%  
  as.data.frame()%>%  
  tibble::add_column(rownames(.))%>%  
  tibble::remove_rownames(.)
#这是长图
colnames(cmt)=c("Correlation Coefficient","Immune cell")
#不要星星的话到这里就行了
#下面是加星星
# 添加 p.adj 列
cmt$`p.adj` <- t(as.data.frame(cor$p.adj))

# 创建显著性标记列
cmt$Significance <- cut(cmt$p.adj, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                        labels = c("***", "**", "*", ""))

# 固定 Immune cell 顺序
df <- cmt

# 定义 Software 对应的颜色映射
software_colors <- c('QUANTISEQ' = "#f27767", 'XCELL' = "#bc9826", 
                     'EPIC' = "#53b449", 'MCPCOUNTER' = "#23b892", 
                     'TIMER' = "#1cb6ec", 'CIBERSORT' = "#d269a4",
                     'CIBERSORT-ABS' = "lightpink2")

df$Software <- sub(".*_", "", df$`Immune cell`)
df$Immune_cell_label <- factor(df$`Immune cell`, levels = unique(df$`Immune cell`))

# 创建 y_cols 映射，确保基于 Software 后缀进行颜色匹配
y_cols <- software_colors[df$Software]
# 绘图
p12 <- ggplot(df, aes(x = `Correlation Coefficient`, y = Immune_cell_label)) +  
  geom_point(aes(color = Software), shape = 16, size = 5) +  
  scale_color_manual(values = software_colors) +  # 手动设置 Software 对应的颜色
  theme_bw() +  # 使用白色背景，并保留网格
  theme(axis.text.y = element_text(size = 10, color = y_cols),  # 根据 Software 分配的 y 轴颜色
        axis.text.x = element_text(size = 10, color = 'black'),        
        axis.title.x = element_text(size = 12),        
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA),  # 设置黑色边框
        panel.grid.major = element_line(color = "gray80"),  # 设置主网格线颜色
        panel.grid.minor = element_blank()) +  # 移除次要网格线
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1),  # 设置从-0.5到0.5，每隔0.1显示一个刻度
                     labels = scales::number_format(accuracy = 0.1)) +  # 保留一位小数
  coord_cartesian(xlim = c(-0.5, 0.5)) +  # 强制设置 x 轴的显示范围
  labs(y = "Immune Cell", x = "Correlation Coefficient")

# 添加显著性标记
p12 <- p12 + geom_text(aes(label = Significance, x = 0.495),  # 将标记放置在右侧
                       hjust = 0, size = 5)

p12
ggsave("免疫浸润与RPN1相关性.pdf", plot = p12, width = 10, height = 20, dpi = 300)
####有相关性的细胞与基因的雷达图####
#install.packages("fmsb")
library(fmsb)
dd=df
dd <- dd[grepl("\\*", dd$Significance), ]
ddcor=dd
dd=t(dd)#将数据进行行列转置
rownames(dd)
# 将第一行（行名为 "Immune cell"）设置为列名
colnames(dd) <- as.character(unlist(dd["Immune cell", ]))
dd=dd["Correlation Coefficient",,drop=F]#保留相关性系数这一行，其他行全部删除
# 保存原来的列名
col_names <- colnames(dd)
# 将数据转换为数值类型
dd <- as.data.frame(lapply(dd, as.numeric))
# 恢复原来的列名
colnames(dd) <- col_names
maxValue=ceiling(max(abs(dd))*10)/10
dd=rbind(rep(maxValue,ncol(dd)),rep(-maxValue,ncol(dd)),dd)
colors="red"
corStat=ddcor
colnames(dd)=paste0(colnames(dd),corStat$Significance)
radarchart( dd, axistype=1 ,
            pcol = "#88419d", pfcol = scales::alpha("#fb9a99", 0.5),
            plwd=2 ,
            plty=1,
            cglcol="grey", cglty=1, 
            caxislabels=seq(-maxValue,maxValue,maxValue/2),
            centerzero = TRUE, # 强制中心点对齐
            cglwd=1.2, axislabcol="blue",vlcex=0.5 )
####基因与免疫浸润相关性蝴蝶图####
colnames(cmt)=c("correlation","label")
cmt1 <- cmt %>% 
  filter(correlation < -0.2 | correlation > 0.2)
df=cmt1
df$category <- sub(".*_", "", df$label)
# 删除最后一个 `_` 后面的内容
df$label <- sub("_(?!.*_).*", "", df$label, perl = TRUE)
library(ggnewscale)
library(tidyverse)
library(ggforce)
# 中心点位置
center <- data.frame(x = 0, y = 0, label = "RPN1")
colorMap <- structure(
  c("#c5523f", '#bc9822', '#7DA88E', '#23b892',"#1cb6ec","#82649c",'lightpink2'),
  names = c('QUANTISEQ', 'XCELL', 'EPIC', 'MCPCOUNTER',
            'TIMER','CIBERSORT','CIBERSORT-ABS')
)
theme_no_axis <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank()
)
names(df)
added_data <- df %>%
  mutate(x = ifelse(correlation < 0, -1, 1)) %>%
  group_by(category) %>%
  arrange(label, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(y = n():1, y = y - mean(y))
adjust <- 0.5
bezier.data <- dplyr::select(added_data, x, y) %>%
  tibble::rowid_to_column('group') %>%
  mutate(type = ifelse(x < 0, 'neg', 'pos')) %>%
  group_by(group) %>%
  group_modify(function(data, group) {
    if (data$x[1] < 0) {
      data.frame(
        x = c(0, -adjust, data$x + adjust, data$x),
        y = c(0, 0, data$y, data$y),
        type = data$type
      )
    } else {
      data.frame(
        x = c(0, adjust, data$x - adjust, data$x),
        y = c(0, 0, data$y, data$y),
        type = data$type
      )
    }
  })
# 绘图
p <- ggplot(added_data, aes(x = x, y = y)) +
  # 曲线连接
  geom_bezier(aes(group = group, 
                  colour = type),
              data = bezier.data, 
              linewidth = 0.4,
              show.legend = FALSE) +
  scale_colour_manual(values = c('#2B393B', '#7F2749')) +
  new_scale_colour() +
  # 中心点
  geom_point(data = center, size = 8, color = "#2B393B") +
  # 外部点
  geom_point(aes(size = abs(correlation), fill = category),
             shape = 21) +
  # 添加点标签
  geom_text(aes(label = label, colour = category), 
            hjust = ifelse(added_data$x > 0, -0, 1), vjust = 0.5,
            nudge_x = ifelse(added_data$x > 0, 0.2, -0.2),
            size = 4, alpha = 0.8,
            show.legend = FALSE) +
  # 添加中心点标签
  geom_text(aes(label = label), 
            data = center, 
            fontface = "italic", size = 6, hjust = 0.5, vjust = 0.5,
            nudge_x = 1.5, colour = "#2B393B") +
  # 调整点的大小和颜色
  scale_colour_manual(values = colorMap) +
  scale_fill_manual(values = colorMap) +
  lims(x = c(-3, 3)) +
  # 美化图形
  theme_no_axis

p
# | label: plot-complex
low_pos <- min(added_data$y)
p1 <- p + annotate(geom = 'segment', 
                   x = -.5, y = low_pos - 1, 
                   xend = -2.5, yend = low_pos - 1,
                   arrow = arrow(angle = 20, length = unit(0.03, 'npc'), 
                                 type = 'closed')) +
  annotate(geom = 'text', x = -2.3, y = low_pos - 2, 
           label = 'Negative correlation',
           hjust = 0) +
  annotate(geom = 'segment', 
           x = .5, y = low_pos - 1, 
           xend = 2.5, yend = low_pos - 1,
           arrow = arrow(angle = 20, length = unit(0.03, 'npc'), 
                         type = 'closed')) +
  annotate(geom = 'text', x = 2.3, y = low_pos - 2, 
           label = 'Positive correlation',
           hjust = 1)
p1
pdf("RPN1-correlation.pdf", width = 8, height = 8)  # 设置 PDF 文件宽度和高度
print(p1)
dev.off()

####基因与estimate相关性####
library("corrplot")
library(linkET)
estimate <- read.csv("estimate.csv")
rownames(estimate)<-estimate$X
estimate$X<-NULL
estimate=t(estimate)
riskscore<-read.csv("riskscore.csv")
rownames(riskscore)<-riskscore$X
riskscore$X<-NULL
gene <- c("RPN1","PML")
# 提取基因表达量
tp <- tpm[gene, ]
tp<-t(as.data.frame(tp))
gg<-merge(tp,riskscore,by = "row.names")
rownames(gg)<-gg$Row.names
gg$Row.names<-NULL
# 假设你的评分数据为 estimate
rownames(estimate) <- gsub("\\.", "-", rownames(estimate))
# 取共同患者
common_patients <- intersect(rownames(gg), rownames(estimate))
# 对 gg 和 estimate 进行子集化
gg <- gg[common_patients, ]
estimate <- estimate[common_patients, ]
identical(rownames(gg),rownames(estimate))
cor_res <- correlate(gg, estimate,method = "spearman")
qcorrplot(cor_res) +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
# 先整理下数据
df_r <- cor_res$r %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1, names_to = "cell_type", values_to = "correlation")
df_p <- cor_res$p %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1, names_to = "cell_type", values_to = "pvalue")
df_cor <- df_r %>%
  left_join(df_p) %>%
  mutate(stars = cut(pvalue, breaks = c(-Inf, 0.05, 0.01, 0.001, Inf), right = F, labels = c("***", "**", "*", " ")))

head(df_cor)
library(ggplot2)
# 设置 gene 列的顺序为 NLN、RPN1、riskScore
df_cor$gene <- factor(df_cor$gene, levels = c("RPN1", "PML", "riskScore"))
# 创建图形对象
p9 <- ggplot(df_cor, aes(gene, cell_type)) +  
  geom_tile(aes(fill = correlation)) +
  geom_text(aes(label = stars), color = "black", size = 4) +
  scale_fill_gradient2(high = '#6A5ACD', low = '#32CD32', mid = 'white',  
                       limit = c(-1, 1), name = paste0("*    p < 0.05", "\n\n", 
                                                       "**  p < 0.01", "\n\n", 
                                                       "*** p < 0.001", "\n\n", 
                                                       "Correlation")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
p9
p9 <- ggplot(df_cor, aes(gene, cell_type)) +  
  geom_point(aes(size = abs(correlation), fill = correlation), shape = 21, color = "black") + # 气泡带边框
  geom_text(aes(label = stars), color = "black", size = 4) +  # 添加显著性标记
  scale_size(range = c(3, 10)) +  # 控制气泡大小范围
  scale_fill_gradient2(high = '#6A5ACD', low = '#32CD32', mid = 'white',  
                       limit = c(-1, 1), name = paste0("*    p < 0.05", "\n\n", 
                                                       "**  p < 0.01", "\n\n", 
                                                       "*** p < 0.001", "\n\n", 
                                                       "Correlation")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # 给整个图添加黑色边框
p9
# 保存为PDF
ggsave("基因与ssgsea相关性气泡图.pdf", plot = p9, device = "pdf", width = 6, height = 10)
