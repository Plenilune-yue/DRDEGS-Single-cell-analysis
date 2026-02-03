rm(list = ls())
setwd("E://YUE实验（新）/12.基因组/")
library(TCGAbiolinks)
library(dplyr)
#SNP数据
# query <- GDCquery(
#   project = "TCGA-PAAD", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# GDCdownload(query)
# GDCprepare(query, save = T,save.filename = "TCGA-PAAD_SNP.Rdata")
library(maftools)
load(file = "TCGA-PAAD_SNP.Rdata")
maf.paad <- data
write.csv(maf.paad,"paad-maf-data.csv")
class(maf.paad)
dim(maf.paad)
## [1] 24849   140
maf.paad[1:10,1:10]
#临床数据
library(magrittr)
clin <- read.csv("clinical.csv")
names(clin)
clin <- clin[,c("Id","time","status","T","N","M","Stage","Gender","risk","riskScore")]
names(clin)<-c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender","risk","riskScore")
laml <- read.maf(
  maf.paad,clinicalData = clin,isTCGA = TRUE)
# 展示样本层面的信息
getSampleSummary(laml)
# 展示基因层面的信息
getGeneSummary(laml)
# 展示临床相关信息
getClinicalData(laml)
# 获取 MAF 文件的所有列字段
getFields(laml)
# 保存为文件
write.mafSummary(maf = laml, basename = 'paad.maf')
laml
#7.1 Plotting MAF summary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#成对的Fisher精确检验分析突变的互斥性和共现性
par(oma = c(3, 4, 5, 1))
output <- somaticInteractions(maf=laml, top=25, pvalue=c(0.05, 0.01)) # 使用变量捕获函数输出
output # 第一次
output # 第二次
write.table(output$pairs, file="somaticInteractions.pairwise.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(output$gene_sets, file="somaticInteractions.comet.tsv", quote=FALSE, row.names=FALSE, sep="\t")
laml@data[["TumorVAF"]] <- laml@data$t_alt_count / laml@data$t_depth
plotVaf(maf = laml, vafCol = 'TumorVAF')
#7.2.1 Drawing oncoplots
oncoplot(maf=laml, top=30, borderCol=NULL,draw_titv = TRUE,
         fontSize = 0.5,
         gene_mar = 8)
#我要进行高低组分析
# 从临床数据中提取性别对应的"Tumor_Sample_Barcode"
clin.high <- subset(clin, risk=="high")$Tumor_Sample_Barcode
clin.high <- clin.high[clin.high != "TCGA-IB-7651"]
clin.low <- subset(clin, risk=="low")$Tumor_Sample_Barcode
# 使用subsetMaf构建高低组的MAF对象
paad.high <- subsetMaf(maf=laml, tsb=clin.high, isTCGA=TRUE)
paad.low <- subsetMaf(maf=laml, tsb=clin.low, isTCGA=TRUE)
paad.high@clinical.data$group <- "high"
paad.low@clinical.data$group <- "low"
#新版的瀑布图
#高组
pdf("high突变瀑布图.pdf",width = 8,height = 5)
oncoplot(maf = paad.high,
         top = 20,
         fontSize = 0.6,
         showTumorSampleBarcodes = F,
         clinicalFeatures = 'group',draw_titv = T,borderCol=NULL)
dev.off()
#low组
pdf("low突变瀑布图.pdf",width = 8,height = 5)
oncoplot(maf = paad.low,
         top = 20,
         fontSize = 0.6,
         showTumorSampleBarcodes = F,
         clinicalFeatures = 'group',
         draw_titv = T,
         borderCol = NULL,
         annotationColor = list(group = c("low" = "blue")))
dev.off()
#绘制一下瀑布图
library(maftools)
library(ComplexHeatmap)
library(circlize)
oncoplot(maf=laml, clinicalFeatures = 'risk',
         annotationDat = clin,borderCol=NULL)
# 设置高低风险组的颜色
annotationColors <- list(
  risk = c(low = "blue", high = "red")
)
# 使用 oncoplot 绘制图并添加分组信息
oncoplot(
  maf = paad.low,
  clinicalFeatures = 'risk',
  annotationDat = clin,
  annotationColor = annotationColors,
  borderCol = NULL,
  top = 20,
  fontSize = 0.5,
  gene_mar = 6
)
#转换（transition）和颠换（transversion）的统计和可视化
paad.titv <- titv(maf=paad.high, plot=FALSE, useSyn=TRUE)
plotTiTv(res=paad.titv)
# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=paad.high, m2=paad.low, m1Name="high", m2Name="low", minMut=5)
#fvsm <- mafCompare(m1=paad.high, m2=paad.low, m1Name="high", m2Name="low", minMut = 2, useCNV =FALSE)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="high_vs_low.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#绘制一个森林图
#挺不好看的，我是为了更好地分组
forestPlot(mafCompareRes = fvsm, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
#2. 比较两个队列的oncoplot
#使用coOncoplot并排绘制两个队列的oncoplot：
# 获取突变最频繁的前 20 个基因
high_gene_summary <- getGeneSummary(paad.high)
low_gene_summary <- getGeneSummary(paad.low)
# 获取前 20 个突变最频繁的基因
top_genes_high <- head(high_gene_summary$Hugo_Symbol, 20)
top_genes_low <- head(low_gene_summary$Hugo_Symbol, 20)

# 合并两个列表，并去重
top_genes <- unique(c(top_genes_high, top_genes_low))

# 绘制共突变图
coOncoplot(
  m1 = paad.high, 
  m2 = paad.low, 
  m1Name = "high", 
  m2Name = "low", 
  borderCol = NULL, 
  genes = top_genes,
  geneNamefont = 0.5,
  removeNonMutated = TRUE
)
#上面是我改的，下面这个是原代码
# 绘制共突变图
#143\144是10.29新增
coOncoplot(m1=paad.high, m2=paad.low, m1Name = 'high', m2Name = 'low', genes = g, removeNonMutated = TRUE)
coBarplot(m1=paad.high, m2=paad.low, m1Name = 'high', m2Name = 'low')

coOncoplot(m1=paad.high, m2=paad.low, m1Name="high", m2Name="low",borderCol=NULL,removeNonMutated = TRUE)
#9.5.3 Co-bar plots
coBarplot(m1 = paad.high, m2 = paad.low, m1Name = "high", m2Name = "low",genes = top_genes)
# 绘制共突变条形图
coBarplot(
  m1 = paad.high, 
  m2 = paad.low, 
  m1Name = "high", 
  m2Name = "low", 
  geneSize = 0.4,
  genes = top_genes,
  normalize = TRUE,
  titleSize = 1,
  showPct = TRUE,
  pctSize = 0.8,
  axisSize = 0.8,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  geneMar = 3
)
##计算每个样本TMB
maf = tmb(maf = laml,
          captureSize = 50,
          logScale = TRUE)   #绘制
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
rownames(maf) <- maf$sample
write.csv(maf,"tmb_res.csv")
#与肿瘤比较
par(oma = c(3, 4, 5, 1))
maf.TMB <- tcgaCompare(
  maf = laml, cohortName = 'PAAD-maf', 
  logscale = TRUE, capture_size = 50
)
####突变互斥####
#exclusive/co-occurance event analysis on top 10 mutated genes.
Interact <- somaticInteractions(maf = paad.low, top = 20, pvalue = c(0.05, 0.1))
#提取P值结果
Interact$gene_sets

#下面，我们来计算每个样本的 TMB 值
get_TMB <- function(file) {
  # 需要用到的列
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # 删除这些突变类型
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # 读取文件
  df <- read_csv(file, col_select = use_cols)
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # 计算 VAF
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}
#读取刚才下载 MAF 文件时保存的整理后的 CSV 文件
library(readr)
paad.csv <- get_TMB('paad-maf-data.csv')
max(paad.csv$TMB)
#[1] 425.5333
#现在，先删除该样本，然后根据 TMB 中位值将样本分为高低两组
paad.tmb <- paad.csv %>% 
  filter(TMB != max(TMB)) %>%
  mutate(
    label = if_else(TMB >= median(TMB), "high", "low"), 
    .before = "MaxVAF"
  ) %>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12))

#然后对两类样本进行生存分析
# 生存分析
library(survival)
library(survminer)
data.surv <- clin %>% filter(!is.na(time)) %>%
  mutate(
    time = time / 365, 
    status = if_else(status == "Alive", 0, 1)
  ) %>% 
  inner_join(paad.tmb, by = "Tumor_Sample_Barcode") %>%
  dplyr::select(c("Tumor_Sample_Barcode", "time", "status", "label"))

km_fit <- survfit(Surv(time, status) ~ label, data = data.surv)

# 绘制生存曲线并添加阴影
ggsurvplot(
  km_fit, data = data.surv,
  pval = TRUE, surv.median.line = "hv",
  legend.labs = c("TMB-High", "TMB-Low"),
  legend.title = "Group",
  title = "Overall survival",
  xlab = "Time(years)",
  risk.table = TRUE,
  palette = c("#8D5A85", "#FFCC4F"),
  legend = "right",  # 将图例放置在右侧
  conf.int = TRUE  # 添加置信区间阴影
)

#结合模型的high和low
# 创建组合变量
# 合并数据
clin1 <- clin
clin1$time <- clin1$time /365
tmb <- data.surv[,c("Tumor_Sample_Barcode","label")]
# 创建组合变量
# 合并数据
merged_data <- merge(clin1, tmb, by = "Tumor_Sample_Barcode")
merged_data$group <- paste(merged_data$label, merged_data$risk, sep = "-")

# 创建生存对象
km_fit <- survfit(Surv(time, status) ~ group, data = merged_data)

# 绘制生存曲线并添加置信区间
surv_plot <- ggsurvplot(
  km_fit, data = merged_data,
  pval = TRUE, surv.median.line = "hv",
  legend.labs = c("High-TMB-High-risk", "High-TMB-Low-risk", "Low-TMB-High-risk", "Low-TMB-Low-risk"),
  legend.title = "Group",
  title = "Overall survival",
  xlab = "Time(years)",
  risk.table = TRUE,
  palette = c("darkorange", "dodgerblue", "forestgreen", "darkorchid"),  # 修改调色板
  conf.int = TRUE,  # 添加置信区间阴影
  ggtheme = theme_minimal(),  # 使用简约主题
  legend = "right"  # 将图例放置在右侧
)
print(surv_plot)
# #处理一下临床信息，去掉 NA 或 不明确的样本
# data.cate <- inner_join(paad.tmb, clin, by = "Tumor_Sample_Barcode") %>%
#   filter_at(vars(T, N, M), all_vars(!endsWith(., "X") & !is.na(.))) %>%
#   mutate(
#     N = gsub("(N\\d).*", "\\1", N, perl = TRUE),
#     stage = if_else(startsWith(stage, c("Stage III", "Stage IV")), "Stage III-IV", "Stage I-II")
#   )
# ggboxplot(data.cate, x = "risk", y = "TMB",
#           fill = "risk") +
#   stat_compare_means(comparisons = list(
#     c("high", "low")
#   )) +
#   stat_compare_means(label.y = 4.5)
# #####换一种
# png(paste0(gene,"_risk.png"),width = 1200,height = 1200,res = 300)
# # 计算p值
# p_value <- wilcox.test(TMB ~ risk, data = data.cate)$p.value
# p_value_label <- sprintf("p = %.3f", p_value) # 格式化p值到小数点后三位
# 
# # 绘图
# # 修改数据框中risk列的因子水平
# data.cate$risk <- factor(data.cate$risk, levels = c("high", "low"), labels = c("High", "Low"))
# 
# # 重新绘制图形
# ggplot(data=data.cate, aes(x=risk, y=TMB, colour=risk)) +
#   geom_violin(alpha=0.8, scale='width', trim=TRUE) +
#   geom_boxplot(aes(fill=risk), alpha=0.5, size=1.5, width=0.3) +
#   geom_jitter(alpha=0.3, size=3) +
#   scale_fill_manual(limits=c("High", "Low"), values=c("#FB8C82", "#B0C8CA")) +
#   scale_color_manual(limits=c("High", "Low"), values=c("#FB8C82", "#B0C8CA")) +
#   geom_signif(comparisons=list(c("High", "Low")),
#               map_signif_level=FALSE,
#               annotations=p_value_label,
#               tip_length=c(0,0),
#               y_position=c(4),
#               size=1,
#               textsize=4,
#               test="wilcox.test",
#               color="black") +
#   theme_bw() +
#   guides(fill=guide_legend(title="Risk"),
#          color=guide_legend(title="Risk")) +
#   labs(x="", y="TMB")
# 
# dev.off()

output <- somaticInteractions(maf=paad.low, top=20, pvalue=c(0.05, 0.01)) # 使用变量捕获函数输出




####TMB与riskScore相关散点图####
names(paad.csv)
tmb=paad.tmb
tmb=data.frame(tmb)
# 提取两列
tmb_selected <- tmb[, c("Tumor_Sample_Barcode", "TMB")]
rs <- clin[,c("Tumor_Sample_Barcode","riskScore","risk")]
cor.data <- merge(tmb_selected,rs,by="Tumor_Sample_Barcode")
table(cor.data$risk)
library(corrplot)
library(ggplot2)
library(ggpubr)
#两两相关性散点图 ggplot
#4.1 计算两两相关性
#取子集：
names(cor.data)
riskScore <- cor.data[,"riskScore"]
TMB <- cor.data[,"TMB"]
#计算相关系数，method可选"pearson", "kendall", "spearman":
cor.test (riskScore, TMB, method="pearson")
#ggplot绘图
p1 <- ggplot(cor.data, aes(x = riskScore, y = TMB))
p2 <- p1 + geom_point() 
p3 <- p2 + geom_smooth(method="lm")
p3
#换一种
#两两相关性散点图 ggscatter
ggscatter(cor.data, x = "riskScore", y = "TMB", 
          add = "reg.line",conf.int = TRUE, 
          fill = "lightgray")
#更改ggscatter绘图参数，添加相关系数及p值等
ggscatter(cor.data, x = "riskScore", y = "TMB",  
          color = "red3",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          cor.coef.coord = c(25, 3.5),
          cor.coef.size = 8)
####tmb与riskscore相关性美丽版####
library(ggplot2)
library(ggside)
library(cols4all)
df = cor.data
#散点图绘制：
p <- ggplot(data = df, 
            mapping = aes(x = riskScore, y = TMB, 
                          color = risk, fill = risk)) + 
  geom_point(size = 2) + #添加散点图
  geom_smooth(method = 'lm', formula = y ~ x, se = T) #添加拟合曲线
p
#自定义主题和配色：
c4a_gui()
#mycol <- c4a('light',2)
mycol <- c("#EE8866", "#77AADD") # 将颜色对调
mytheme <- theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

p1 <- p +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  mytheme
p1
#在散点图上方添加密度分布曲线：
p2 <- p1 + 
  geom_xsidedensity(aes(y = after_stat(density)), 
                    alpha = 0.3,
                    position = "stack")
p2 #子图和主图共享x轴，独立y轴
#密度曲线+箱线图：
p7 <- p1 + 
  geom_xsidedensity(aes(y = after_stat(density)), 
                    alpha = 0.3,
                    position = "stack") +
  geom_ysideboxplot(aes(x = riskScore), 
                    orientation = "x", alpha = 0.3) + 
  scale_ysidex_discrete() +
  theme_ggside_classic() +
  theme(ggside.panel.scale = 0.3,
        ggside.axis.text.x = element_text(angle = 45,hjust = 1))
p7
# 计算整体的相关系数和p值
cor_test <- cor.test(df$riskScore, df$TMB, method = "spearman")

# 提取相关系数和p值
cor_value <- round(cor_test$estimate, 3)
p_value <- ifelse(cor_test$p.value < 0.001, "<0.001", signif(cor_test$p.value, 3))

# 绘制图形并手动添加相关系数和p值的文本标签
p7 <- p1 + 
  geom_xsidedensity(aes(y = after_stat(density)), 
                    alpha = 0.3,
                    position = "stack") +
  geom_ysideboxplot(aes(x = riskScore), 
                    orientation = "x", alpha = 0.3) + 
  scale_ysidex_discrete() +
  theme_ggside_classic() +
  theme(ggside.panel.scale = 0.3,
        ggside.axis.text.x = element_text(angle = 45, hjust = 1)) +
  # 手动添加相关系数和p值的文本标签
  annotate("text", x = 3, y = 3, 
           label = paste("R =", cor_value, ", p =", p_value), 
           hjust = 0, size = 5, color = "black")

p7
# 计算high和low组的TMB的p值
wilcox_test <- wilcox.test(TMB ~ risk, data = df)
p_value_risk <- ifelse(wilcox_test$p.value < 0.001, "<0.001", signif(wilcox_test$p.value, 3))
# 在p7的基础上，添加高低组TMB差异的p值并标注
p8 <- p7 + 
  ggside::geom_ysidetext(
    aes(x = max(riskScore) * 0.5),  # Adjust horizontal position
    y = 3.2,  # Adjust vertical position
    label = paste("p =", p_value_risk), 
    size = 4,
    color = "black"
  )
p8


##########################################################################
#TMB先到这
#临床富集分析
clin_enrich <- clinicalEnrichment(maf=laml, clinicalFeature="risk")
write.table(clin_enrich$pairwise_comparision, file="clin_enrich_pair.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(clin_enrich$groupwise_comparision, file="clin_enrich_group.tsv", quote=FALSE, row.names=FALSE, sep="\t")
plotEnrichmentResults(enrich_res=clin_enrich, pVal=0.001)

# 肿瘤异质性和MATH score
# 1. 肿瘤样品的异质性
# 指定Tumor_Sample_Barcodes为"TCGA-55-7283"
vafclust <- inferHeterogeneity(maf=laml, tsb="TCGA-55-7283")
# 保存结果到文件
write.table(vafclust$clusterData, file="vaf_clustdata.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(vafclust$clusterMeans, file="vaf_clustmean.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#输出结果中各突变基因VAF的数据储存在$clusterData中，聚类数量和平均VAF储存在$clusterMeans中。用plotClusters函数展示聚类结果：
plotClusters(clusters=vafclust)

#循环出每一个患者的math
# 设定结果保存目录
result_dir <- "result"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}
# 获取所有患者的样本条码
patients <- unique(laml@data$Tumor_Sample_Barcode)
# 创建空的列表用于存储所有患者的结果
all_cluster_data <- list()
all_cluster_means <- list()
# 循环处理每个患者
for (patient in patients) {
  # 执行异质性推断
  vafclust <- inferHeterogeneity(maf=laml, tsb=patient)
  
  # 添加患者信息到结果数据中
  cluster_data <- vafclust$clusterData
  cluster_data$Patient <- patient
  all_cluster_data[[patient]] <- cluster_data
  
  cluster_means <- vafclust$clusterMeans
  cluster_means$Patient <- patient
  all_cluster_means[[patient]] <- cluster_means
  
  # 生成并保存聚类图
  plot_file <- file.path(result_dir, paste0(patient, "_cluster_plot.png"))
  png(filename = plot_file)
  tryCatch({
    plotClusters(clusters=vafclust)
  }, error = function(e) {
    cat("Error in plotClusters for patient:", patient, "\n")
  })
  dev.off()
}
library(data.table) # 加载data.table包来使用rbindlist
# 合并所有患者的结果
all_cluster_data_combined <- rbindlist(all_cluster_data, fill=TRUE)
all_cluster_means_combined <- rbindlist(all_cluster_means, fill=TRUE)
# 保存所有患者的结果到文件
write.csv(all_cluster_data_combined,file=file.path(result_dir,"all_vaf_clustdata.csv"))
write.csv(all_cluster_means_combined, file=file.path(result_dir, "all_vaf_clustmean.csv"))
cat("All patients have been processed and results saved in the 'result' directory.\n")
#TMB与高低患者math的组间差异
#提取tmb
math<-all_cluster_data_combined[,c("Tumor_Sample_Barcode","MATH")]
math <- unique(math)
mathl<-clin[,c("Tumor_Sample_Barcode","risk")]
math<-merge(math,mathl,by="Tumor_Sample_Barcode")
ggboxplot(math, x = "risk", y = "MATH",
          fill = "risk") +
  stat_compare_means(comparisons = list(
    c("high", "low")
  )) +
  stat_compare_means(label.y = 4.5)

#高低组的GISTIC_2.0分析
rm(list = ls())
setwd("E://YUE实验（新）/12.基因组/")
# 读取数据
segment <- read.table("segment_file.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marker <- read.table("marker_file.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clin <- read.csv("clinical.csv")
names(clin)
clin <- clin[,c("Id","risk")]
names(clin)<-c("Sample","risk")
# 提取 Sample 列，保留 '-' 分隔符之前的前三部分
segment$Sample <- sapply(segment$Sample, function(x) {
  parts <- unlist(strsplit(x, "-"))  # 将字符串按照 '-' 分隔
  paste(parts[1:3], collapse = "-")  # 保留前三段并重新组合
})
data <- merge(segment,clin,by="Sample")
high <- data[data$risk == "high", ]
low <- data[data$risk == "low", ]
# 删除最后一列
high <- high[, -ncol(high)]
low <- low[, -ncol(low)]
write.table(high, file = "high.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(low, file = "low.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#然后上传GISTIC平台