#install_github("ZJU-UoE-CCW-LAB/scCDC",force = TRUE)
suppressMessages({
  library(splatter)
  library(Seurat)
  library(msigdbr)
  library(harmony)
  library(patchwork)
  library(singleseqgset)
  library(heatmap3)
  library(GEOquery)
  library(Seurat) 
  library(scater)
  library(SingleCellExperiment)
  library(DoubletFinder)
  library(viridisLite)
  library(viridis)
  library(clustree)
  library(MAST)
  library(scCDC)
  library(Matrix)
  library(destiny)
  library(ggbeeswarm)
  library(ggthemes)
  library(dplyr)
  library(tidyr)
  library(Seurat)
  library(stringr)
  library(colorRamps)
  library(colorRamps)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(singleseqgset)
  library(devtools)
})
setwd("~/data/GSE212966_RAW")
####清洗数据####
#p1
PDAC1 <- "~/data/GSE212966_RAW/PDAC1"  # 你的数据存放路径
expression_matrix_p1 <- ReadMtx(
  mtx = file.path(PDAC1,"GSM6567157_PDAC1_matrix.mtx.gz"), features = file.path(PDAC1,"GSM6567157_PDAC1_genes.tsv.gz"),
  cells = file.path(PDAC1,"GSM6567157_PDAC1_barcodes.tsv.gz")
)
seurat_object.p1 <- CreateSeuratObject(counts = expression_matrix_p1)
seurat_object.p1[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p1, pattern = "^MT-")
seurat_object.p1 <- NormalizeData(seurat_object.p1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p1 <- FindVariableFeatures(seurat_object.p1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p1)
seurat_object.p1 <- ScaleData(seurat_object.p1, features = all.genes)
seurat_object.p1 <- RunPCA(seurat_object.p1, features = VariableFeatures(object = seurat_object.p1))
seurat_object.p1 <- RunUMAP(seurat_object.p1, dims = 1:20)
seurat_object.p1 <- FindNeighbors(seurat_object.p1, dims = 1:20)
seurat_object.p1 <- FindClusters(seurat_object.p1, resolution = 0.5)
seurat_object.p1@meta.data$orig.ident <- "PDAC1"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p1, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
#双细胞比例计算
annotations <- seurat_object.p1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p1)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p1$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p1 <- doubletFinder(seurat_object.p1, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
seurat_object.p1@meta.data$DF.classifications_0.25_0.22_441#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p1, reduction = "umap", group.by = "DF.classifications_0.25_0.22_441")#看一下双细胞和单细胞
#去除双细胞
pbmc_singlet <- subset(seurat_object.p1, cells= rownames(seurat_object.p1@meta.data[seurat_object.p1@meta.data$DF.classifications_0.25_0.22_441=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
PDAC1 <- pbmc_singlet
saveRDS(PDAC1, file = "result/PDAC1.rds")

#p2
PDAC2 <- "~/data/GSE212966_RAW/PDAC2"  # 你的数据存放路径
expression_matrix_p2 <- ReadMtx(
  mtx = file.path(PDAC2,"GSM6567159_PDAC2_matrix.mtx.gz"), features = file.path(PDAC2,"GSM6567159_PDAC2_genes.tsv.gz"),
  cells = file.path(PDAC2,"GSM6567159_PDAC2_barcodes.tsv.gz")
)
seurat_object.p2 <- CreateSeuratObject(counts = expression_matrix_p2)
seurat_object.p2[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p2, pattern = "^MT-")
seurat_object.p2 <- NormalizeData(seurat_object.p2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p2 <- FindVariableFeatures(seurat_object.p2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p2)
seurat_object.p2 <- ScaleData(seurat_object.p2, features = all.genes)
seurat_object.p2 <- RunPCA(seurat_object.p2, features = VariableFeatures(object = seurat_object.p2))
seurat_object.p2 <- RunUMAP(seurat_object.p2, dims = 1:20)
seurat_object.p2 <- FindNeighbors(seurat_object.p2, dims = 1:20)
seurat_object.p2 <- FindClusters(seurat_object.p2, resolution = 0.5)
seurat_object.p2@meta.data$orig.ident <- "PDAC2"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p2, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
#双细胞比例计算
annotations <- seurat_object.p2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p2)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p2$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p2 <- doubletFinder(seurat_object.p2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.p2@meta.data$
seurat_object.p2@meta.data$DF.classifications_0.25_0.26_298#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p2, reduction = "umap", group.by = "DF.classifications_0.25_0.26_298")#看一下双细胞和单细胞
#去除双细胞
pbmc_singlet <- subset(seurat_object.p2, cells= rownames(seurat_object.p2@meta.data[seurat_object.p2@meta.data$DF.classifications_0.25_0.26_298=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
PDAC2 <- pbmc_singlet
saveRDS(PDAC2, file = "result/PDAC2.rds")

#p3
PDAC3 <- "~/data/GSE212966_RAW/PDAC3"  # 你的数据存放路径
expression_matrix_p3 <- ReadMtx(
  mtx = file.path(PDAC3,"GSM6567160_PDAC3_matrix.mtx.gz"), features = file.path(PDAC3,"GSM6567160_PDAC3_genes.tsv.gz"),
  cells = file.path(PDAC3,"GSM6567160_PDAC3_barcodes.tsv.gz")
)
seurat_object.p3 <- CreateSeuratObject(counts = expression_matrix_p3)
seurat_object.p3[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p3, pattern = "^MT-")
seurat_object.p3 <- NormalizeData(seurat_object.p3, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p3 <- FindVariableFeatures(seurat_object.p3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p3)
seurat_object.p3 <- ScaleData(seurat_object.p3, features = all.genes)
seurat_object.p3 <- RunPCA(seurat_object.p3, features = VariableFeatures(object = seurat_object.p3))
seurat_object.p3 <- RunUMAP(seurat_object.p3, dims = 1:20)
seurat_object.p3 <- FindNeighbors(seurat_object.p3, dims = 1:20)
seurat_object.p3 <- FindClusters(seurat_object.p3, resolution = 0.5)
seurat_object.p3@meta.data$orig.ident <- "PDAC3"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p3, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
#双细胞比例计算
annotations <- seurat_object.p3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p3)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p3$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p3 <- doubletFinder(seurat_object.p3, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.p3@meta.data$   #看DF值
seurat_object.p3@meta.data$DF.classifications_0.25_0.08_311#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p3, reduction = "umap", group.by = "DF.classifications_0.25_0.08_311")#看一下双细胞和单细胞,DF这里每次都要改
#去除双细胞
pbmc_singlet <- subset(seurat_object.p3, cells= rownames(seurat_object.p3@meta.data[seurat_object.p3@meta.data$DF.classifications_0.25_0.08_311=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
PDAC3 <- pbmc_singlet
saveRDS(PDAC3, file = "result/PDAC3.rds")

#p4
PDAC4 <- "~/data/GSE212966_RAW/PDAC4"  # 你的数据存放路径
expression_matrix_p4 <- ReadMtx(
  mtx = file.path(PDAC4,"GSM6567161_PDAC4_matrix.mtx.gz"), features = file.path(PDAC4,"GSM6567161_PDAC4_genes.tsv.gz"),
  cells = file.path(PDAC4,"GSM6567161_PDAC4_barcodes.tsv.gz")
)
seurat_object.p4 <- CreateSeuratObject(counts = expression_matrix_p4)
seurat_object.p4[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p4, pattern = "^MT-")
seurat_object.p4 <- NormalizeData(seurat_object.p4, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p4 <- FindVariableFeatures(seurat_object.p4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p4)
seurat_object.p4 <- ScaleData(seurat_object.p4, features = all.genes)
seurat_object.p4 <- RunPCA(seurat_object.p4, features = VariableFeatures(object = seurat_object.p4))
seurat_object.p4 <- RunUMAP(seurat_object.p4, dims = 1:20)
seurat_object.p4 <- FindNeighbors(seurat_object.p4, dims = 1:20)
seurat_object.p4 <- FindClusters(seurat_object.p4, resolution = 0.5)
seurat_object.p4@meta.data$orig.ident <- "PDAC4"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p4, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
#双细胞比例计算
annotations <- seurat_object.p4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p4)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p4$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p4 <- doubletFinder(seurat_object.p4, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.p4@meta.data$   #看DF值
seurat_object.p4@meta.data$DF.classifications_0.25_0.09_267#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p4, reduction = "umap", group.by = "DF.classifications_0.25_0.09_267")#看一下双细胞和单细胞,DF这里每次都要改
#去除双细胞
pbmc_singlet <- subset(seurat_object.p4, cells= rownames(seurat_object.p4@meta.data[seurat_object.p4@meta.data$DF.classifications_0.25_0.09_267=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
PDAC4 <- pbmc_singlet
saveRDS(PDAC4, file = "result/PDAC4.rds")

#p5
PDAC5 <- "~/data/GSE212966_RAW/PDAC5"  # 你的数据存放路径
expression_matrix_p5 <- ReadMtx(
  mtx = file.path(PDAC5,"GSM6567163_PDAC5_matrix.mtx.gz"), features = file.path(PDAC5,"GSM6567163_PDAC5_genes.tsv.gz"),
  cells = file.path(PDAC5,"GSM6567163_PDAC5_barcodes.tsv.gz")
)
seurat_object.p5 <- CreateSeuratObject(counts = expression_matrix_p5)
seurat_object.p5[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p5, pattern = "^MT-")
seurat_object.p5 <- NormalizeData(seurat_object.p5, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p5 <- FindVariableFeatures(seurat_object.p5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p5)
seurat_object.p5 <- ScaleData(seurat_object.p5, features = all.genes)
seurat_object.p5 <- RunPCA(seurat_object.p5, features = VariableFeatures(object = seurat_object.p5))
seurat_object.p5 <- RunUMAP(seurat_object.p5, dims = 1:20)
seurat_object.p5 <- FindNeighbors(seurat_object.p5, dims = 1:20)
seurat_object.p5 <- FindClusters(seurat_object.p5, resolution = 0.5)
seurat_object.p5@meta.data$orig.ident <- "PDAC5"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p5, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
#双细胞比例计算
annotations <- seurat_object.p5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p5)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p5$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p5 <- doubletFinder(seurat_object.p5, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.p5@meta.data$   #看DF值
seurat_object.p5@meta.data$DF.classifications_0.25_0.09_167#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p5, reduction = "umap", group.by = "DF.classifications_0.25_0.09_167")#看一下双细胞和单细胞,DF这里每次都要改
#去除双细胞
pbmc_singlet <- subset(seurat_object.p5, cells= rownames(seurat_object.p5@meta.data[seurat_object.p5@meta.data$DF.classifications_0.25_0.09_167=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
  PDAC5 <- pbmc_singlet
saveRDS(PDAC5, file = "result/PDAC5.rds")

#p6
PDAC6 <- "~/data/GSE212966_RAW/PDAC6"  # 你的数据存放路径
expression_matrix_p6 <- ReadMtx(
  mtx = file.path(PDAC6,"GSM6567164_PDAC6_matrix.mtx.gz"), features = file.path(PDAC6,"GSM6567164_PDAC6_genes.tsv.gz"),
  cells = file.path(PDAC6,"GSM6567164_PDAC6_barcodes.tsv.gz")
)
seurat_object.p6 <- CreateSeuratObject(counts = expression_matrix_p6)
seurat_object.p6[["percent.mt"]] <- PercentageFeatureSet(seurat_object.p6, pattern = "^MT-")
seurat_object.p6 <- NormalizeData(seurat_object.p6, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.p6 <- FindVariableFeatures(seurat_object.p6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.p6)
seurat_object.p6 <- ScaleData(seurat_object.p6, features = all.genes)
seurat_object.p6 <- RunPCA(seurat_object.p6, features = VariableFeatures(object = seurat_object.p6))
seurat_object.p6 <- RunUMAP(seurat_object.p6, dims = 1:20)
seurat_object.p6 <- FindNeighbors(seurat_object.p6, dims = 1:20)
seurat_object.p6 <- FindClusters(seurat_object.p6, resolution = 0.5)
seurat_object.p6@meta.data$orig.ident <- "PDAC6"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.p6, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
##最佳参数：
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
#双细胞比例计算
annotations <- seurat_object.p6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.p6)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.p6$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.p6 <- doubletFinder(seurat_object.p6, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.p6@meta.data$   #看DF值
seurat_object.p6@meta.data$DF.classifications_0.25_0.3_209#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.p6, reduction = "umap", group.by = "DF.classifications_0.25_0.3_209")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.p6@meta.data$DF.classifications_0.25_0.3_209)
#去除双细胞
pbmc_singlet <- subset(seurat_object.p6, cells= rownames(seurat_object.p6@meta.data[seurat_object.p6@meta.data$DF.classifications_0.25_0.3_209=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
PDAC6 <- pbmc_singlet
saveRDS(PDAC6, file = "result/PDAC6.rds")

#a1
ADJ1 <- "~/data/GSE212966_RAW/ADJ1"  # 你的数据存放路径
expression_matrix_a1 <- ReadMtx(
  mtx = file.path(ADJ1,"GSM6567165_ADJ1_matrix.mtx.gz"), features = file.path(ADJ1,"GSM6567165_ADJ1_genes.tsv.gz"),
  cells = file.path(ADJ1,"GSM6567165_ADJ1_barcodes.tsv.gz")
)
seurat_object.a1 <- CreateSeuratObject(counts = expression_matrix_a1)
seurat_object.a1[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a1, pattern = "^MT-")
seurat_object.a1 <- NormalizeData(seurat_object.a1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a1 <- FindVariableFeatures(seurat_object.a1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a1)
seurat_object.a1 <- ScaleData(seurat_object.a1, features = all.genes)
seurat_object.a1 <- RunPCA(seurat_object.a1, features = VariableFeatures(object = seurat_object.a1))
seurat_object.a1 <- RunUMAP(seurat_object.a1, dims = 1:20)
seurat_object.a1 <- FindNeighbors(seurat_object.a1, dims = 1:20)
seurat_object.a1 <- FindClusters(seurat_object.a1, resolution = 0.5)
seurat_object.a1@meta.data$orig.ident <- "ADJ1"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a1, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a1)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a1$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a1 <- doubletFinder(seurat_object.a1, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a1@meta.data$   #看DF值
seurat_object.a1@meta.data$DF.classifications_0.25_0.09_278#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a1, reduction = "umap", group.by = "DF.classifications_0.25_0.09_278")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a1@meta.data$DF.classifications_0.25_0.09_278)
pbmc_singlet <- subset(seurat_object.a1, cells= rownames(seurat_object.a1@meta.data[seurat_object.a1@meta.data$DF.classifications_0.25_0.09_278=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
pbmc_singlet@meta.data$
ADJ1 <- pbmc_singlet
saveRDS(ADJ1, file = "result/ADJ1.rds")

#a2
ADJ2 <- "~/data/GSE212966_RAW/ADJ2"  # 你的数据存放路径
expression_matrix_a2 <- ReadMtx(
  mtx = file.path(ADJ2,"GSM6567166_ADJ2_matrix.mtx.gz"), features = file.path(ADJ2,"GSM6567166_ADJ2_genes.tsv.gz"),
  cells = file.path(ADJ2,"GSM6567166_ADJ2_barcodes.tsv.gz")
)
seurat_object.a2 <- CreateSeuratObject(counts = expression_matrix_a2)
seurat_object.a2[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a2, pattern = "^MT-")
seurat_object.a2 <- NormalizeData(seurat_object.a2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a2 <- FindVariableFeatures(seurat_object.a2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a2)
seurat_object.a2 <- ScaleData(seurat_object.a2, features = all.genes)
seurat_object.a2 <- RunPCA(seurat_object.a2, features = VariableFeatures(object = seurat_object.a2))
seurat_object.a2 <- RunUMAP(seurat_object.a2, dims = 1:20)
seurat_object.a2 <- FindNeighbors(seurat_object.a2, dims = 1:20)
seurat_object.a2 <- FindClusters(seurat_object.a2, resolution = 0.5)
seurat_object.a2@meta.data$orig.ident <- "ADJ2"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a2, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a2)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a2$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a2 <- doubletFinder(seurat_object.a2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a2@meta.data$   #看DF值
seurat_object.a2@meta.data$DF.classifications_0.25_0.17_265#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a2, reduction = "umap", group.by = "DF.classifications_0.25_0.17_265")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a2@meta.data$DF.classifications_0.25_0.17_265)
pbmc_singlet <- subset(seurat_object.a2, cells= rownames(seurat_object.a2@meta.data[seurat_object.a2@meta.data$DF.classifications_0.25_0.17_265=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
ADJ2 <- pbmc_singlet
saveRDS(ADJ2, file = "result/ADJ2.rds")

#a3
ADJ3 <- "~/data/GSE212966_RAW/ADJ3"  # 你的数据存放路径
expression_matrix_a3 <- ReadMtx(
  mtx = file.path(ADJ3,"GSM6567167_ADJ3_matrix.mtx.gz"), features = file.path(ADJ3,"GSM6567167_ADJ3_genes.tsv.gz"),
  cells = file.path(ADJ3,"GSM6567167_ADJ3_barcodes.tsv.gz")
)
seurat_object.a3 <- CreateSeuratObject(counts = expression_matrix_a3)
seurat_object.a3[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a3, pattern = "^MT-")
seurat_object.a3 <- NormalizeData(seurat_object.a3, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a3 <- FindVariableFeatures(seurat_object.a3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a3)
seurat_object.a3 <- ScaleData(seurat_object.a3, features = all.genes)
seurat_object.a3 <- RunPCA(seurat_object.a3, features = VariableFeatures(object = seurat_object.a3))
seurat_object.a3 <- RunUMAP(seurat_object.a3, dims = 1:20)
seurat_object.a3 <- FindNeighbors(seurat_object.a3, dims = 1:20)
seurat_object.a3 <- FindClusters(seurat_object.a3, resolution = 0.5)
seurat_object.a3@meta.data$orig.ident <- "ADJ3"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a3, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a3)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a3$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a3 <- doubletFinder(seurat_object.a3, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a3@meta.data$   #看DF值
seurat_object.a3@meta.data$DF.classifications_0.25_0.005_74#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a3, reduction = "umap", group.by = "DF.classifications_0.25_0.005_74")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a3@meta.data$DF.classifications_0.25_0.005_74)
pbmc_singlet <- subset(seurat_object.a3, cells= rownames(seurat_object.a3@meta.data[seurat_object.a3@meta.data$DF.classifications_0.25_0.005_74=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
ADJ3 <- pbmc_singlet
saveRDS(ADJ3, file = "result/ADJ3.rds")

#a4
ADJ4 <- "~/data/GSE212966_RAW/ADJ4"  # 你的数据存放路径
expression_matrix_a4 <- ReadMtx(
  mtx = file.path(ADJ4,"GSM6567169_ADJ4_matrix.mtx.gz"), features = file.path(ADJ4,"GSM6567169_ADJ4_genes.tsv.gz"),
  cells = file.path(ADJ4,"GSM6567169_ADJ4_barcodes.tsv.gz")
)
seurat_object.a4 <- CreateSeuratObject(counts = expression_matrix_a4)
seurat_object.a4[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a4, pattern = "^MT-")
seurat_object.a4 <- NormalizeData(seurat_object.a4, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a4 <- FindVariableFeatures(seurat_object.a4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a4)
seurat_object.a4 <- ScaleData(seurat_object.a4, features = all.genes)
seurat_object.a4 <- RunPCA(seurat_object.a4, features = VariableFeatures(object = seurat_object.a4))
seurat_object.a4 <- RunUMAP(seurat_object.a4, dims = 1:20)
seurat_object.a4 <- FindNeighbors(seurat_object.a4, dims = 1:20)
seurat_object.a4 <- FindClusters(seurat_object.a4, resolution = 0.5)
seurat_object.a4@meta.data$orig.ident <- "ADJ4"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a4, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a4)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a4$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a4 <- doubletFinder(seurat_object.a4, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a4@meta.data$   #看DF值
seurat_object.a4@meta.data$DF.classifications_0.25_0.16_240#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a4, reduction = "umap", group.by = "DF.classifications_0.25_0.16_240")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a4@meta.data$DF.classifications_0.25_0.16_240)
pbmc_singlet <- subset(seurat_object.a4, cells= rownames(seurat_object.a4@meta.data[seurat_object.a4@meta.data$DF.classifications_0.25_0.16_240=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
ADJ4 <- pbmc_singlet
saveRDS(ADJ4, file = "result/ADJ4.rds")

#a5
ADJ5 <- "~/data/GSE212966_RAW/ADJ5"  # 你的数据存放路径
expression_matrix_a5 <- ReadMtx(
  mtx = file.path(ADJ5,"GSM6567170_ADJ5_matrix.mtx.gz"), features = file.path(ADJ5,"GSM6567170_ADJ5_genes.tsv.gz"),
  cells = file.path(ADJ5,"GSM6567170_ADJ5_barcodes.tsv.gz")
)
seurat_object.a5 <- CreateSeuratObject(counts = expression_matrix_a5)
seurat_object.a5[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a5, pattern = "^MT-")
seurat_object.a5 <- NormalizeData(seurat_object.a5, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a5 <- FindVariableFeatures(seurat_object.a5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a5)
seurat_object.a5 <- ScaleData(seurat_object.a5, features = all.genes)
seurat_object.a5 <- RunPCA(seurat_object.a5, features = VariableFeatures(object = seurat_object.a5))
seurat_object.a5 <- RunUMAP(seurat_object.a5, dims = 1:20)
seurat_object.a5 <- FindNeighbors(seurat_object.a5, dims = 1:20)
seurat_object.a5 <- FindClusters(seurat_object.a5, resolution = 0.5)
seurat_object.a5@meta.data$orig.ident <- "ADJ5"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a5, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a5)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a5$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a5 <- doubletFinder(seurat_object.a5, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a5@meta.data$   #看DF值
seurat_object.a5@meta.data$DF.classifications_0.25_0.005_163#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a5, reduction = "umap", group.by = "DF.classifications_0.25_0.005_163")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a5@meta.data$DF.classifications_0.25_0.005_163)
pbmc_singlet <- subset(seurat_object.a5, cells= rownames(seurat_object.a5@meta.data[seurat_object.a5@meta.data$DF.classifications_0.25_0.005_163=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
ADJ5 <- pbmc_singlet
saveRDS(ADJ5, file = "result/ADJ5.rds")

#a6
ADJ6 <- "~/data/GSE212966_RAW/ADJ6"  # 你的数据存放路径
expression_matrix_a6 <- ReadMtx(
  mtx = file.path(ADJ6,"GSM6567171_ADJ6_matrix.mtx.gz"), features = file.path(ADJ6,"GSM6567171_ADJ6_genes.tsv.gz"),
  cells = file.path(ADJ6,"GSM6567171_ADJ6_barcodes.tsv.gz")
)
seurat_object.a6 <- CreateSeuratObject(counts = expression_matrix_a6)
seurat_object.a6[["percent.mt"]] <- PercentageFeatureSet(seurat_object.a6, pattern = "^MT-")
seurat_object.a6 <- NormalizeData(seurat_object.a6, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object.a6 <- FindVariableFeatures(seurat_object.a6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object.a6)
seurat_object.a6 <- ScaleData(seurat_object.a6, features = all.genes)
seurat_object.a6 <- RunPCA(seurat_object.a6, features = VariableFeatures(object = seurat_object.a6))
seurat_object.a6 <- RunUMAP(seurat_object.a6, dims = 1:20)
seurat_object.a6 <- FindNeighbors(seurat_object.a6, dims = 1:20)
seurat_object.a6 <- FindClusters(seurat_object.a6, resolution = 0.5)
seurat_object.a6@meta.data$orig.ident <- "ADJ6"
#确定参数，寻找最优pk值
sweep.res.list_pbmc <- paramSweep(seurat_object.a6, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc) #可以看到最佳参数的点
mpK<-as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
mpK
annotations <- seurat_object.a6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(seurat_object.a6)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
nExp_poi <- round(DoubletRate*length(seurat_object.a6$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object.a6 <- doubletFinder(seurat_object.a6, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)#502个
#seurat_object.a6@meta.data$   #看DF值
seurat_object.a6@meta.data$DF.classifications_0.25_0.27_274#看一下上一步跑出来的DF.数值
DimPlot(seurat_object.a6, reduction = "umap", group.by = "DF.classifications_0.25_0.27_274")#看一下双细胞和单细胞,DF这里每次都要改
table(seurat_object.a6@meta.data$DF.classifications_0.25_0.27_274)
pbmc_singlet <- subset(seurat_object.a6, cells= rownames(seurat_object.a6@meta.data[seurat_object.a6@meta.data$DF.classifications_0.25_0.27_274=="Singlet",]))
DimPlot(pbmc_singlet, reduction = "umap",label = TRUE,pt.size = 1.5)
ADJ6 <- pbmc_singlet
saveRDS(ADJ6, file = "result/ADJ6.rds")

pdac.combined <- merge(PDAC1, y = c(PDAC2, PDAC3,PDAC4,PDAC5,PDAC6), add.cell.ids = c("PDAC1", "PDAC2", "PDAC3","PDAC4","PDAC5","PDAC6"), project = "PDAC")
pdac.combined@meta.data$tissue <- "PDAC"
adj.combined <- merge(ADJ1, y = c(ADJ2, ADJ3,ADJ4,ADJ5,ADJ6), add.cell.ids = c("ADJ1", "ADJ2", "ADJ3","ADJ4","ADJ5","ADJ6"), project = "ADJ")
adj.combined@meta.data$tissue <- "ADJ"
combined <- merge(pdac.combined,adj.combined)
saveRDS(combined, file="result/GSE212966_after_doublefinder.rds" )

####去RNA污染####
library(decontX)#devtools::install_github("campbio/decontX")
library(tidydr)
sce <- combined
LayerData(sce, assay = "RNA", layer = "counts")
sce <- JoinLayers(sce)
dim(sce[["RNA"]]$counts )
#[1] 36601 68087
colnames(sce@meta.data)
sce1 <- sce #留一个sce
saveRDS(sce1, file="result/sce1.rds" )
sce <- readRDS("result/sce1.rds")
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
count.feature.ls <- sce@meta.data[, c("nCount_RNA", "nFeature_RNA")]
count.feature.ls %<>% map(log10) %>% map(~c(10^(mean(.x) + 3*sd(.x)), 10^(mean(.x) - 3*sd(.x))))
sce <- subset(sce, subset = nFeature_RNA > 200 &
                     nFeature_RNA < count.feature.ls[[2]][1] &
                     nCount_RNA < count.feature.ls[[1]][1] &
                     percent.mt < 20)
####整合去批次-查看批次####
sce <- NormalizeData(sce,normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
#降维时的一些参数教程：https://mp.weixin.qq.com/s/2DCGx7GwDKnzub313ay2tQ
sce <- FindVariableFeatures(sce, 
                                 selection.method = "vst", 
                                 nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce),npcs = 110, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:40, umap.method = "uwot", 
                       metric = "cosine", reduction = "harmony", 
                       n.neighbors=30, min.dist=0.01,spread=1)
p1 <- DimPlot(sce, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sce, label = TRUE, group.by = "tissue")
####整合数据####
pdac.combined.raw <- sce
pdac.combined.raw[["percent.mt"]] <- PercentageFeatureSet(pdac.combined.raw, pattern = "^MT-")
pdac.combined.raw <- PercentageFeatureSet(pdac.combined.raw, pattern = "^KRT", col.name = "percent.krt")
pdac.combined.raw <- PercentageFeatureSet(pdac.combined.raw, pattern = "MGP", col.name = "percent.MGP")
pdac.combined.raw <- SCTransform(pdac.combined.raw , vars.to.regress = c("percent.mt", "percent.krt", "percent.MGP"), verbose = TRUE)
pbmc <- pdac.combined.raw
pbmc <- RunPCA(pbmc, eatures = VariableFeatures(object = sce),npcs = 110, verbose = FALSE)
pdf("elbowplot.pdf",width = 8,height = 8)
ElbowPlot(pbmc, ndims = 50)
dev.off()
pbmc <- RunUMAP(pbmc, dims = 1:13, verbose = FALSE,
                n.neighbors=50, min.dist=0.01,spread=1)
pbmc <- FindNeighbors(pbmc, dims = 1:13, verbose = FALSE)
pbmc <- FindClusters(pbmc,verbose = FALSE)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pbmc, label = TRUE)
combined_plot <-  p1 |p2
pdf(file = "result/1.降维/combined_umap.pdf", width = 14, height = 6) # 宽度增加以适应两个图
print(combined_plot)
dev.off()
####自动注释(⚠️不使用)####
#SingleR https://zhuanlan.zhihu.com/p/448269413
#remotes::install_github("LTLA/celldex")
library(SingleR)
library(Seurat)
library(pheatmap)
library(celldex)
library(ggplot2)
library(plyr)
library(cowplot)
library(ggh4x)
library(viridis)
# ref <- BlueprintEncodeData() ##下载注释数据库-网址https://cloud.tencent.com/developer/article/1820060
# save(ref,file = 'BlueprintEncodeData.Rdata')
# load("HumanPrimaryCellAtlas_hpca.se_human.RData")
load("BlueprintEncodeData.Rdata")
#进行SingleR注释
hms_cluster <- pbmc
hms_cluster_for_SingleR <- GetAssayData(hms_cluster, slot="data") ##获取标准化矩阵
hms_cluster.hesc <- SingleR(test = hms_cluster_for_SingleR, ref = ref, labels = ref$label.main) #
hms_cluster.hesc
hms_cluster.hesc$labels
#seurat 和 SingleR的table表
a <- table(hms_cluster.hesc$labels, hms_cluster@meta.data$seurat_clusters)# pbmc@meta.data  pbmc的meta文件，包含了seurat的聚类结果
write.csv(a, "result/1.降维/SingleR注释.csv")
#保存亚群名到hms_cluster
hms_cluster@meta.data$labels <-hms_cluster.hesc$labels
pdf("result/1.降维/hms_cluster_bySingleR.pdf",width = 21,height = 8)
DimPlot(hms_cluster, group.by = c("seurat_clusters", "labels"),reduction = "umap",label = T)
dev.off()
saveRDS(hms_cluster, file = "result/1.降维/hms_cluster_bySingleR.rds")
####手动注释####
#hms_cluster <- pbmc  
genes <- list(
  "Fibroblast" = c("COL1A1", "COL1A2", "COL3A1", "DCN"),
  "B cell" = c("MS4A1", "CD79A","LINC01781"),
  "T cell" = c("IL7R", "GZMK","TIGIT"),
  "Pancreatic ductal cell" = c("CLDN4", "KRT19", "KRT7"),
  "Macrophage" = c("CD68", "CD14", "CD163"),
  "Endothelial cell" = c("VWF",  "PECAM1", "CLDN5"),
  "NK cells" = c("GNLY", "GZMB", "KLRD1"),
  "Acinar cell" = c("PRSS1","CTRB1","PRSS2"),
  "Neitrophil" = c("G0S2", "S100A9", "S100A8"),
  "Pancreatic stellate cell" = c("ADIRF", "RGS5", "PDGFRB"),
  "Mast cell" = c("CPA3", "TPSAB1","TPSB2"),
  "Endocrine cell" = c("CHGA", "INS", "CHGB"),
  "Plasma cell" = c("IGHA1", "JCHAIN", "MZB1"),
  "Schwann cell" = c("CDH19", "LGI4","S100B")
)
p <- DotPlot(hms_cluster,
        features = genes,
        cols = c("#ffffff", "#448444")
) +
  RotatedAxis() + # 来自Seurat
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
  )
p$data$feature.groups2 <- factor(p$data$feature.groups, 
                                 levels = c("Fibroblast","B cell","T cell",
                                            "Pancreatic ductal cell","Macrophage","Endothelial cell","NK cells","Acinar cell",
                                            "Neitrophil","Pancreatic stellate cell","Mast cell","Endocrine cell",
                                            "Plasma cell","Schwann cell"))
color_map <- c("Fibroblast" = "#FFFFB3",
               "B cell" = "#B3DE69",
               "T cell" = "#BEBADA",
               "Pancreatic ductal cell" = "#a44e89",
               "Macrophage" = "#8DD3C7",
               "Endothelial cell" = "#FB8072",
               "NK cells" = "#CCEBC5",
               "Acinar cell" = "#FCCDE5",
               "Neitrophil" = "#80B1D3",
               "Pancreatic stellate cell" = "#FDB462",
               "Mast cell" = "#33A02C",
               "Endocrine cell" = '#e4b8a5',  # Not specified, will default or be omitted
               "Plasma cell" = "#1F78B4",
               "Schwann cell" = "#D9D9D9")
strip <- strip_themed(
  background_x = elem_list_rect(fill = color_map[levels(p$data$feature.groups2)])
)
p$data %>% 
  ggplot(aes(x = features.plot,
             y = id)) + 
  geom_point(aes(size = pct.exp, 
                 color = avg.exp.scaled)) + 
  facet_wrap2(~feature.groups2, 
              scales = "free_x", 
              strip = strip, 
              nrow = 1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5, 
                                   vjust = 0.3, 
                                   color = "black"),
        axis.title = element_blank(),
        strip.background = element_rect(color = "white"),
        axis.text.y = element_blank()) + 
  scale_color_gradient(low = "#ffffff",
                       high = "#448444", 
                       name = "avg.exp") -> p

p
df <- data.frame(x = 0, y = levels(hms_cluster), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y)
pl <- ggplot(df, aes(x, y, color = factor(y))) +
  geom_point(size = 6, shape = 15, show.legend = F) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +
  theme(
    plot.margin = margin(r=0),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
pl
vv <- plot_grid(pl, p, align = "h", axis="bt", rel_widths = c(1.5, 35))
vv
ggsave("result/1.降维/注释DotPlot.pdf", plot = vv, width = 15, height = 10)
ann.ids <- c("T cell",                  #0    
             "T cell",                  #1    
             "B cell",                  #2    
             "Fibroblast",             #3
             "T cell",              #4
             "T cell",        #5
             "Macrophage",                  #6
             "Fibroblast",            #7
             "Endothelial cell",            #8
             "Macrophage",      #9
             "Acinar cell",           #10
             "Pancreatic ductal cell" ,                 #11
             "T cell",     #12
             "Pancreatic ductal cell",       #13
             "Neitrophil",       #14
             "Pancreatic stellate cell",    #15
             "Pancreatic ductal cell",    #16
             "T cell",               #17
             "Acinar cell",                 #18
             "Fibroblast",                #19
             "Acinar cell",                #20
             "Mast cell",        #21
             "Neitrophil",      #22
             "Pancreatic ductal cell",  #23
             "Fibroblast",       #24
             "Plasma cell",     #25
             "NK cell",      #26
             "Fibroblast",#27
             "Schwann cell",#28
             "Macrophage"#29
)
seuratidens=mapvalues(Idents(hms_cluster), from = levels(Idents(hms_cluster)), to = ann.ids)
Idents(hms_cluster)=seuratidens
hms_cluster$cellType=Idents(hms_cluster)
table(hms_cluster$cellType)
# T cell                   B cell               Fibroblast
# 20400                     4420                     8699
# Macrophage         Endothelial cell              Acinar cell
# 5055                     2359                     4454
# Pancreatic ductal cell               Neitrophil Pancreatic stellate cell
# 5697                     2335                     1498
# Mast cell              Plasma cell                  NK cell
# 965                      495                      412
# Schwann cell
# 224
p1 <- DimPlot(hms_cluster, reduction = "umap", group.by = "seurat_clusters",label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
combined_plot <-  p1 |p2
pdf(file = "result/1.降维/combined_umap.pdf", width = 14, height = 6) # 宽度增加以适应两个图
print(combined_plot)
dev.off()
saveRDS(hms_cluster, file = "result/1.降维/注释后seurat.rds")
####细胞代谢分析####
# install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
# BiocManager::install("AUCell")
# BiocManager::install("GSEABase",force = T)
# BiocManager::install("GSVA",force = T)
# devtools::install_github("YosefLab/VISION@v2.1.0") #有一个curl修复还有好多依赖包
# devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA <- hms_cluster
sc.metabolism.SeuratV5 <- function (obj, method = "VISION", imputation = F, ncores = 2, 
                                    metabolism.type = "KEGG") 
{
  countexp <- GetAssayData(obj, layer='counts')
  countexp <- data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  if (metabolism.type == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (metabolism.type == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  }
  cat("Start quantify the metabolism activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  }
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  obj@assays$METABOLISM$score <- signature_exp
  obj
}
Idents(scRNA) <- "cellType"
res <-sc.metabolism.SeuratV5(obj = scRNA,
                             method = "AUCell", # VISION、AUCell、ssgsea和gsva
                             imputation =F, ncores = 2, 
                             metabolism.type = "KEGG") # KEGG和REACTOME
# 需要修改DimPlot.metabolism中的UMAP大小写
DimPlot.metabolismV5 <- function (obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, 
                                  size = 1) 
{
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  if (dimention.reduction.type == "umap") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, 
    ]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = umap_1, 
                                                y = umap_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
    ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
                                                y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  plot
}
# check一下有哪些pathway
pathways <- res@assays$METABOLISM$score
head(rownames(pathways))
# [1] "Glycolysis / Gluconeogenesis"            
# [2] "Citrate cycle (TCA cycle)"               
# [3] "Pentose phosphate pathway"               
# [4] "Pentose and glucuronate interconversions"
# [5] "Fructose and mannose metabolism"         
# [6] "Galactose metabolism"
# Dimplot
DimPlot.metabolismV5(obj = res, 
                     pathway = "Glycolysis / Gluconeogenesis", 
                     dimention.reduction.type = "umap",  
                     dimention.reduction.run = F, size = 1)
# DotPlot
input.pathway <- rownames(res@assays[["METABOLISM"]][["score"]])[1:30]
pdf("result/1.降维/metabolism-dotplot.pdf", width = 8, height = 8)
DotPlot.metabolism(obj = res,
                   pathway = input.pathway,
                   phenotype = "cellType", # 这个参数需按需修改
                   norm = "y")
dev.off()
# BoxPlot
# 由于开发者默认吧obj定义为countexp.Seurat，所以还需要重新命名一下
countexp.Seurat <- res
pdf("result/1.降维/metabolism-boxplot.pdf", width = 20, height = 30)
BoxPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway, 
                   phenotype = "cellType", #这个参数需按需修改
                   ncol = 4)
dev.off()
#单细胞代谢水平计算
sce_Metal_exp = countexp.Seurat
sce_Metal_exp$celltype = sce_Metal_exp$cellType
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$celltype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.celltype),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG
c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:5]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]
#绘制单细胞代谢热图
rownames(c_k_d) = c_k_d$KEGG
pdf("result/1.降维/metabolism-pheatmap.pdf", width = 6, height = 8)
pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
dev.off()
####细胞比例####
#https://mp.weixin.qq.com/s/K-yvh9sSZcy3yQPHHaZxSA
phe <- hms_cluster
phe$cellType <- as.character(Idents(phe))
phe <- phe@meta.data
colnames(phe)
table(phe$orig.ident)
table(phe$cellType)
#取出需要的数据
group_dotplot_df <- phe %>%
  group_by(tissue,cellType) %>%
  tally() %>%                                 
  group_by(tissue) %>%                         
  mutate(percentage = n / sum(n)) 
stage_dotplot_df <- phe %>%
  group_by(orig.ident, cellType) %>%
  tally() %>%                                 
  group_by(orig.ident) %>%                         
  mutate(percentage = n / sum(n)) 
colnames(stage_dotplot_df)[1] <- "tissue"
group_dotplot_df_all <- rbind(group_dotplot_df, stage_dotplot_df)#合并全部数据
head(group_dotplot_df_all)
#将group列转换为因子，并设置级别
group_dotplot_df_all$tissue <- factor(group_dotplot_df_all$tissue, 
                                     levels = c("ADJ", "PDAC","ADJ1", "ADJ2", "ADJ3","ADJ4","ADJ5","ADJ6",
                                                "PDAC1", "PDAC2", "PDAC3","PDAC4","PDAC5","PDAC6"))
col <- c("Fibroblast" = "#FFFFB3",
         "B cell" = "#B3DE69",
         "T cell" = "#BEBADA",
         "Pancreatic ductal cell" = "#a44e89",
         "Macrophage" = "#8DD3C7",
         "Endothelial cell" = "#FB8072",
         "NK cells" = "#CCEBC5",
         "Acinar cell" = "#FCCDE5",
         "Neitrophil" = "#80B1D3",
         "Pancreatic stellate cell" = "#FDB462",
         "Mast cell" = "#33A02C",
         "Endocrine cell" = '#e4b8a5',  # Not specified, will default or be omitted
         "Plasma cell" = "#1F78B4",
         "Schwann cell" = "#D9D9D9")
#气泡图
dotplot <- ggplot(group_dotplot_df_all, aes(x = tissue, y = cellType, size = percentage, color = cellType)) + 
  geom_point(alpha = 0.7) +
  scale_color_manual(values = col) +  # 自定义颜色
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),  # 调整图例标题字体大小
    legend.text = element_text(size = 8)    # 调整图例标签字体大小
  ) +
  labs(y = "Cell proportion", x = NULL, size = "Cell proportion", color = "Clusters")
dotplot
ggsave('result/1.降维/细胞比例dotplot.pdf',width = 10,height = 7)
#X轴堆叠柱状图
barplot_x <- ggplot(group_dotplot_df_all, aes(x = tissue, y = n, fill = cellType)) +
  geom_bar(stat = "identity") +  # 使用原始数值
  scale_fill_manual(values = col) +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )+
  labs(y = "Number of cells", x = NULL)
barplot_x
ggsave('result/1.降维/细胞比例barplot_x.pdf',width=7, height=6)
#y轴条形图
celltype_sum <- as.data.frame(table(phe$cellType))
colnames(celltype_sum) <- c("cellType", "total_n")
head(celltype_sum)
barplot_y <- ggplot(celltype_sum, aes(x = cellType, y = total_n, fill = cellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col) +
  labs(title = "Cell type composition", x = "Cell Type", y = "Number of cells") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
barplot_y
ggsave('result/1.降维/细胞比例barplot_y.pdf',width=7, height=5)
##细胞比例条图
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(scales)
library(readr)
library(gplots)
library(tibble)
library(grid)
library(rlang)
library(plotrix)
library(ggsci)
table(hms_cluster$cellType)
sample_table <- as.data.frame(table(hms_cluster@meta.data$orig.ident, 
                                    hms_cluster@meta.data$cellType))
names(sample_table) <- c("Samples", "celltype", "cellNumber")
sample_table$Samples <- factor(sample_table$Samples, 
                               levels = c("ADJ", "PDAC","ADJ1", "ADJ2", "ADJ3","ADJ4","ADJ5","ADJ6",
                                          "PDAC1", "PDAC2", "PDAC3","PDAC4","PDAC5","PDAC6"))
col <- c("Fibroblast" = "#FFFFB3",
               "B cell" = "#B3DE69",
               "T cell" = "#BEBADA",
               "Pancreatic ductal cell" = "#a44e89",
               "Macrophage" = "#8DD3C7",
               "Endothelial cell" = "#FB8072",
               "NK cells" = "#CCEBC5",
               "Acinar cell" = "#FCCDE5",
               "Neitrophil" = "#80B1D3",
               "Pancreatic stellate cell" = "#FDB462",
               "Mast cell" = "#33A02C",
               "Endocrine cell" = '#e4b8a5',  # Not specified, will default or be omitted
               "Plasma cell" = "#1F78B4",
               "Schwann cell" = "#D9D9D9")
# 横轴为样本, 纵轴为细胞比例
p <- ggplot(sample_table, aes(x = Samples, weight = cellNumber, fill = celltype)) +  
  geom_bar(position = "fill", width = 0.7, size = 0.5, colour = '#222222') +  
  scale_fill_manual(values = col) +   
  theme(panel.grid = element_blank(),        
        panel.background = element_rect(fill = "transparent", colour = NA),        
        axis.line.x = element_line(colour = "black"),        
        axis.line.y = element_line(colour = "black"),        
        plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5, size = 16)) + 
  labs(y = "Percentage") + 
  RotatedAxis()
pdf("result/1.降维/细胞比例条图.pdf", width = 8, height = 6)
p
dev.off()
####差异基因-富集分析####
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(scRNAtoolVis)
library(ClusterGVis)
library(org.Hs.eg.db)
seurat_merge_v5 <- hms_cluster
table(seurat_merge_v5$cellType)
Idents(seurat_merge_v5) <- seurat_merge_v5$cellType
seurat_merge_v5 <- PrepSCTFindMarkers(seurat_merge_v5)# 调整 SCT 数据以便于 marker 检测
library(future)
plan("multicore", workers = 8)  # 设置 8 核心并行，按自己电脑调整核心数
options(future.globals.maxSize = 20 * 1024^3)  # 允许最大内存，比如 20GB
pbmc.markers.all <- Seurat::FindAllMarkers(
  seurat_merge_v5,
  only.pos = FALSE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)
plan("sequential")  # 恢复到单线程模式
write.csv(pbmc.markers.all,file = "result/2.富集分析/FindAllMarkers.csv")
col <- c(
  "Macrophage" = "#8DD3C7", 
  "Fibroblast" = "#FFFFB3",
  "T cell" = "#BEBADA",
  "Pancreatic ductal cell" = "#a44e89",
  "Endothelial cell" = "#FB8072",
  "Neitrophil" = "#80B1D3",
  "Pancreatic stellate cell" = "#FDB462",
  "B cell" = "#B3DE69",
  "Acinar cell" = "#FCCDE5",
  "NK cells" = "#CCEBC5",
  "Plasma cell" = "#1F78B4",
  "Mast cell" = "#33A02C",
  "Schwann cell" ="#D9D9D9"
)
mygene <- c('RPN1','PML') #用myMarkers参数展示指定基因
p <- jjVolcano(diffData = pbmc.markers.all, 
          tile.col = col, aesCol = c('purple','orange'),
          myMarkers = mygene, polar = T)
p
ggsave(filename = "result/2.富集分析/cluster火山图.pdf",   # 文件名
       plot = p,            # 图形对象
       width = 10,                      # 宽度（单位：英寸）
       height = 10,                     # 高度（单位：英寸）
       dpi = 300)                      # 分辨率（适用于位图格式）
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 4, wt = avg_log2FC)
source("result/2.富集分析/prepareDataFromscRNA.R")
st.data <- prepareDataFromscRNA(object = seurat_merge_v5,
                                diffData = pbmc.markers,
                                showAverage = TRUE)
#loop plot
lapply(unique(pbmc.markers$cluster), function(x){
  tmp <- pbmc.markers |> dplyr::filter(cluster == x)
  
  # plot
  p <- Seurat::FeaturePlot(object = pbmc,
                           features = tmp$gene,
                           ncol = 4)
  
  return(p)
}) -> gglist
# assign names
names(gglist) <- paste("C",1:13,sep = "")
source("result/2.富集分析/visCluster.R")
source("result/2.富集分析/enrichCluster.R")
# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
head(enrich,3)
#            group                                         Description       pvalue    ratio
# GO:0002573    C1                   myeloid leukocyte differentiation 0.0003941646 66.66667
# GO:0050870    C1            positive regulation of T cell activation 0.0005222343 66.66667
# GO:1903039    C1 positive regulation of leukocyte cell-cell adhesion 0.0006265618 66.66667
# barplot
palette = c(
  "Grays", "Light Grays",
  "Blues 2", "Blues 3",
  "Purples 2", "Purples 3",
  "Reds 2", "Reds 3",
  "Greens 2", "Greens 3",
  "Oranges", "YlOrRd", "PuBu"
)
# loop
lapply(seq_along(unique(enrich$group)), function(x){
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description,levels = tmp$Description)
  
  # plot
  p <-
    ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue),y = Description,fill = -log10(pvalue)),
             width = 0.75) +
    geom_line(aes(x = log10(ratio),y = as.numeric(Description)),color = "grey50") +
    geom_point(aes(x = log10(ratio),y = Description),size = 3,color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~.,name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  return(p)
}) -> gglist

# assign names
names(gglist) <- paste("C",1:13,sep = "")
# insert bar plot
pdf('result/2.富集分析/sc_ggplot_go.pdf',height = 25,width = 16,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           line.side = "left",
           column_names_rot = 45,
           markGenes = pbmc.markers$gene,
           cluster.order = c(1:13),
           ggplot.panel.arg = c(5,0.5,16,"grey90",NA),
           gglist = gglist)
dev.off()
####plot1cell-umap####
library(plot1cell)
human_data <- readRDS("result/1.降维/注释后seurat.rds")
col <- c(
  "Acinar cell" = "#FCCDE5",
  "B cell" = "#B3DE69",
  "Endothelial cell" = "#FB8072",
  "Fibroblast" = "#FFFFB3",
  "Macrophage" = "#8DD3C7", 
  "Mast cell" = "#33A02C",
  "Neitrophil" = "#80B1D3",
  "NK cells" = "#CCEBC5",
  "Pancreatic ductal cell" = "#a44e89",
  "Pancreatic stellate cell" = "#FDB462",
  "Plasma cell" = "#1F78B4",
  "Schwann cell" ="#D9D9D9",
  "T cell" = "#BEBADA"
)

circ_data <- prepare_circlize_data(human_data, scale = 0.8)
cluster_colors<-rand_color(length(names(table(human_data$cellType))))
group_colors<-rand_color(length(names(table(human_data$tissue))))
rep_colors<-rand_color(length(names(table(human_data$orig.ident))))
pdf(file = "result/1.降维/plot1cell.pdf", width = 8, height = 8)  
plot_circlize(circ_data, do.label = T, pt.size = 0.5, 
              col.use = col ,bg.color = 'white', 
              kde2d.n = 500, repel = T, label.cex = 0.8)
add_track(circ_data, group = "tissue", colors = group_colors, track_num = 2) 
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) 
#添加legend
legend("topright", 
       legend = unique(human_data$tissue),
       col = group_colors,
       pch = 15,
       cex=0.5,
       pt.cex=2,
       box.lwd = 0,
       bg=NULL,
       y.intersp = 1,
       box.lty=0)
legend("topleft", 
       legend = unique(human_data$orig.ident),
       col = rep_colors,
       pch = 15,
       cex=0.5,
       pt.cex=2,
       box.lwd = 0,
       bg=NULL,
       y.intersp = 1.5,
       box.lty=0)
dev.off()
####ggcircle####
#https://mp.weixin.qq.com/s/wAm1OQ3q1yWeXd_xLzOZ0g
#devtools::install_github("junjunlab/ggcirclize")
library(SeuratData)
library(Seurat)
library(ggnewscale)
library(ggcirclize)
pbmc <- human_data
pbmc <- UpdateSeuratObject(object = pbmc)
meta <- pbmc@meta.data
meta <- cbind(pbmc@reductions$umap@cell.embeddings,meta)
meta$cell_id <- rownames(meta)
meta$x <- 1:nrow(meta)
# colnames(meta)
# [1] "UMAP_1"             "UMAP_2"             "orig.ident"         "nCount_RNA"
# [5] "nFeature_RNA"       "seurat_annotations" "percent.mt"         "RNA_snn_res.0.5"
# [9] "seurat_clusters"    "cell_id"            "x"
col <- c(
  "Acinar cell" = "#FCCDE5",
  "B cell" = "#B3DE69",
  "Endothelial cell" = "#FB8072",
  "Fibroblast" = "#FFFFB3",
  "Macrophage" = "#8DD3C7", 
  "Mast cell" = "#33A02C",
  "Neitrophil" = "#80B1D3",
  "NK cell" = "#CCEBC5",
  "Pancreatic ductal cell" = "#a44e89",
  "Pancreatic stellate cell" = "#FDB462",
  "Plasma cell" = "#1F78B4",
  "Schwann cell" ="#D9D9D9",
  "T cell" = "#BEBADA"
)
my_30_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(30)
circle <- ggcirclize(meta, aes(x = x, y = 1, end = 360,
                               sector = cellType, sector.bg.col = NA)) +
  geom_tracksector(aes(r0 = 0.9, r1 = 1,  # ❗ 这里前移填充圈
                       fill = cellType),
                   strip.label.space = 0.05, show.legend = F,
                   add.yaxis = F, add.xaxis = F) +
  scale_fill_manual(values = col) +  # ✅ 设置 cellType 填充色
  
  new_scale_color() +
  geom_trackpoint(aes(r0 = 0.85, r1 = 0.9,  # ❗ 这里前移 seurat_clusters
                      color = seurat_clusters),
                  strip.label = F,
                  add.yaxis = F, add.xaxis = F) +
  scale_color_manual(values = my_30_colors) +  # ✅ seurat_clusters 颜色
  
  new_scale_color() +
  geom_trackpoint(aes(r0 = 0.8, r1 = 0.85,  # ❗ 这里前移 percent.mt
                      color = percent.mt),
                  strip.label = F,
                  add.yaxis = F, add.xaxis = F) +
  scale_color_gradient(low = "white", high = "purple") +
  
  new_scale_color() +
  geom_cpoint(aes(x = umap_1, y = umap_2,
                  color = cellType, cluster = cellType),
              psize = 0.65) +
  scale_color_manual(values = col)

circle
ggsave("result/1.降维/ggciclize.pdf", plot = circle, width = 15, height = 10)
####3D降维####
library(plotly)
cell <- readRDS("result/1.降维/注释后seurat.rds")
names(cell@reductions) #查看降维方式
cell <- RunUMAP(cell, 
                reduction = "pca",
                dims=1:30,
                n.components=3,
                n.neighbors=10,
                min.dist=0.5)
# cell <- RunUMAP(cell, dims = 1:40, umap.method = "uwot", 
#         metric = "cosine", reduction = "pca", 
#         n.neighbors=30, min.dist=0.01,n.components=3,spread=1)
####三维umap####
umap <- cbind(as.data.frame(cell@reductions$umap@cell.embeddings),
              celltype = Idents(cell))
#"umap_1"   "umap_2"   "umap_3"   "umap_4"   "umap_5"   "celltype"     
label <- umap %>%
  group_by(celltype)%>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
head(umap)
head(label)
umap_3d <- as.data.frame(cell@reductions$umap@cell.embeddings)
head(umap_3d)
celltype <- Idents(cell)
table(celltype)
umap <- cbind(umap_3d,celltype) #数据框合并
head(umap)
p <- plot_ly(
  data = umap_3d,
  x = ~umap_1,
  y = ~umap_2,
  z = ~umap_3,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 3),
  color = ~celltype,
  colors = col
)
htmlwidgets::saveWidget(as_widget(p), "result/1.降维/celltype-umap.html")
####单基因三维umap中表达####
library(scCustomize)
gene <- FetchData(object = cell, 
                  vars = c("umap_1", "umap_2", "umap_3", "Expression"= 'RPN1'), slot = 'data')
p1 <- plot_ly(
  data = gene,
  x = ~umap_1, y = ~umap_2, z = ~umap_3,
  color = ~RPN1,    # 基于表达值着色
  opacity = 0.5,   # 透明度
  colors = colorRampPalette(c("grey", "red"))(256),  # 设置颜色映射：灰色到红色
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)  # 点大小
) %>%
  layout(title = 'RPN1')  # 图标题
htmlwidgets::saveWidget(as_widget(p1), "result/1.降维/RPN1-umap表达.html")
##平面
pdf(file = "result/1.降维/PML在celltype中表达.pdf",width =25,height = 7)
FeaturePlot_scCustom(seurat_object = cell, 
                     features = "PML", 
                     split.by = "cellType", 
                     num_columns = 7,
                     reduction = "umap")
dev.off()
#平面小提琴
plots <- VlnPlot(cell, features = "RPN1", split.by = "tissue", group.by = "cellType",
                 pt.size = 0, combine = FALSE)
pdf("result/1.降维/RPN1表达提琴图.pdf",width = 16,height = 4)
wrap_plots(plots = plots, ncol = 1)
dev.off()
FeaturePlot(cell, features = c("RPN1"),#, "PML"
           
            max.cutoff = 3,
            cols = c("grey", "red"))
#平面气泡
pdf("result/1.降维/两基因表达气泡图.pdf",width =6,height = 6)
DotPlot(cell, features = c("RPN1","PML"), 
        cols = c("#7FC97F", "#BEAED4"),  
        dot.scale = 8, group.by = "cellType") +
  RotatedAxis()
dev.off()
#x,y轴互换
p_list <- VlnPlot(cell, features = c("RPN1","PML"), split.by = "tissue", group.by = "cellType",
                  pt.size = 0, combine = FALSE)
p_list <- lapply(p_list, function(p) p + coord_flip())
CombinePlots(p_list)
pdf("result/1.降维/RPN1表达提琴图（换）.pdf",width = 4,height = 6)
CombinePlots(p_list)
dev.off()
####galaxy####
#https://mp.weixin.qq.com/s/CkYKNzXfLES2n3thkbK_Sg
# options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
# options("repos"=c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
# remotes::install_github("alserglab/mascarade")
pbmc<-readRDS("result/1.降维/注释后seurat.rds")
library(ggplot2)
library(ggSCvis)
library(dplyr)
library(grid)
library(mascarade)
#单基因
genes <- c("RPN1", "PML")
p <- ggscplot(object = pbmc,features = genes) +
  geom_density2d_filled(bins = 10,show.legend = F) +
  geom_scPoint(aes(color = value),size = 0.2,alpha = 0.08) +
  scale_color_gradient(low = "white",high = "#CC0033",name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_manual(values = colorRampPalette(c("white","#336633"))(10))
p
ggsave("result/3.galaxy/RPN1+PML-密度.pdf", plot = p, width = 10, height = 4.6)
galaxy <- ggscplot(object = pbmc,features = genes) +
  stat_density2d(geom = "raster",aes(fill = ..density..),
                 contour = F) +
  geom_scPoint(aes(color = value),size = 0.1,# ✅ 点大小调小
               alpha = 0.25       # ✅ 增加透明度（0~1之间）
               ) +
  scale_color_gradient(low = "white",high = "black",name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_viridis_c(option = "rocket",direction = 1,name = "cell density") +
  coord_cartesian(expand = F)
galaxy
ggsave("result/3.galaxy/RPN1+PML-galaxy.pdf", plot = galaxy, width = 10, height = 4.6)
#tissue
p1 <-
  ggscplot(object = pbmc) +
  stat_density2d(geom = "raster",aes(fill = ..density..),
                 contour = F,show.legend = F) + #需要p图的注释时show.legend = T
  geom_scPoint(color = "white",size = 0.1,alpha = 0.1) +
  facet_wrap(~tissue,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank()) +
  scale_fill_viridis_c(option = "magma",direction = 1) + #inferno
  coord_cartesian(expand = F)
p2 <-
  ggscplot(object = pbmc) +
  geom_scPoint(aes(color = cellType,
                   cluster = cellType),
               show.legend = F,
               label.gp = gpar(fontsize = 8,fontface = "bold.italic")) +
  scale_color_manual(values = col) +  # ✅ 添加颜色映射
  facet_wrap(~tissue,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  xlab("")
tissue <- cowplot::plot_grid(plotlist = list(p2,p1),ncol = 1)
tissue
ggsave("result/3.galaxy/ADJ-PDAC-galaxy.pdf", plot = tissue, width = 15, height = 10)
#all
p3 <-
  ggscplot(object = pbmc) +
  stat_density2d(geom = "raster",aes(fill = ..density..),
                 contour = F,show.legend = F) +
  geom_scPoint(color = "white",size = 0.1,alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank()) +
  scale_fill_viridis_c(option = "magma",direction = 1) + #inferno
  coord_cartesian(expand = F)
p4 <-
  ggscplot(object = pbmc) +
  geom_scPoint(aes(color = cellType,
                   cluster = cellType),
               show.legend = F,
               label.gp = gpar(fontsize = 8,fontface = "bold.italic")) +
  scale_color_manual(values = col) +  # ✅ 添加颜色映射
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  xlab("")
all <- cowplot::plot_grid(plotlist = list(p4,p3),ncol = 1)
all
ggsave("result/3.galaxy/ADJ+PDAC-galaxy.pdf", plot = all, width = 5, height = 10)

####拟时序####
#BiocManager::install("slingshot")
#BiocManager::install("BUSpaRse")
#devtools::install_github("zhanghao-njmu/SCP")
#教程参考：https://mp.weixin.qq.com/s/9CQVyhjk7yvt9dGVbqMhbA
library(SCP)
library(slingshot)  #https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
library(SingleCellExperiment)
library(RColorBrewer)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
#library(Matrix)
library(tradeSeq)
soj <-  readRDS("result/1.降维/注释后seurat.rds")
table(soj$cellType)
macrophage <- subset(soj, cellType == "Macrophage")
dim(macrophage)
#[1] 27792  5055
sce_myeloid <- NormalizeData(macrophage,normalization.method = "LogNormalize",
                             scale.factor = 1e4)
GetAssay(sce_myeloid,assay = "RNA")
sce_myeloid <- FindVariableFeatures(sce_myeloid,
                                    selection.method = "vst",nfeatures = 3000)
sce_myeloid <- ScaleData(sce_myeloid)
sce_myeloid <- RunPCA(object = sce_myeloid,dims = 1:110,pc.genes = VariableFeatures(sce_myeloid))
sce_myeloid <- RunHarmony(sce_myeloid, "orig.ident")
names(sce_myeloid@reductions)
sce_myeloid <- RunUMAP(sce_myeloid,  dims = 1:20, 
                         reduction = "harmony")
sce_myeloid
sce_myeloid <- FindNeighbors(sce_myeloid,dims = 1:20,reduction = "harmony")
sce_myeloid <- FindClusters(sce_myeloid, resolution = 0.5)
DimPlot(sce_myeloid,label = T,pt.size =2)
genes <- list(
  MKI67=c("MKI67","TOP2A","PCLAF","UBE2C","TK1"),
  HSP=c("HSPA6","SERPINH1","BAG3","HSPB1","HSPD1"),
  SPP1=c("SPP1","MARCO","FBP1","APOC1","LIPA"),
  IL1B=c("IL1B","IL1A","NLRP3","PTGS2","CCL3"),
  FOLR2=c("FOLR2","LYVE1","SELENOP","SLC40A1","MRC1"),
  MT=c("MT1H","MT1G","MT1X","MT1E","MT2A")
)
# 生成 DotPlot，并设置颜色为红黄渐变
p <- DotPlot(object = sce_myeloid, features = genes, dot.scale = 12, cols = c("yellow", "red")) + 
  theme_classic() +  # 保持框框样式
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(3, "line"),  # 调整图例条宽度
    legend.key.height = unit(0.5, "line")  # 调整图例条高度
  ) +
  guides(color = guide_colorbar(
    title = "Avg Expression",  # 自定义图例条标题
    frame.colour = "black",    # 图例条框框颜色
    ticks.colour = "black"
  ))
ann.ids <- c("IL1B",      #cluster0    
             "MT",      #cluster1    
             "SPP1",     #2    
             "SPP1",       #3
             "IL1B",             #4
             "SPP1",               #5
             "IL1B",                 #6
             "IL1B",            #7
             "HSP",   #8
             "IL1B",              #9
             "MKI67",              #10
             "HSP"                  #11
             )
seuratidens=mapvalues(Idents(sce_myeloid), from = levels(Idents(sce_myeloid)), to = ann.ids)
Idents(sce_myeloid)=seuratidens
sce_myeloid$cellType=Idents(sce_myeloid)
p <- DotPlot(object = sce_myeloid, features = genes, dot.scale = 12, cols = c("yellow", "red")) + 
  theme_classic() +  # 保持框框样式
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(3, "line"),  # 调整图例条宽度
    legend.key.height = unit(0.5, "line")  # 调整图例条高度
  ) +
  guides(color = guide_colorbar(
    title = "Avg Expression",  # 自定义图例条标题
    frame.colour = "black",    # 图例条框框颜色
    ticks.colour = "black"
  ))
DimPlot(sce_myeloid,label = T,pt.size =2)
RPN1 <- GetAssayData(sce_myeloid,layer = "data")["RPN1",]
threshold <- median(RPN1)
sce_myeloid$group <- ifelse(RPN1 > threshold, "high", "low")
table(sce_myeloid$group)
# high  low
# 1671 3384
head(sce_myeloid@meta.data)
sce_myeloid$cell.type=Idents(sce_myeloid)
subset_rnaslot=GetAssay(sce_myeloid,assay = "RNA")
subset_rnaslot=AddMetaData(sce_myeloid,metadata = sce_myeloid@meta.data)
head(subset_rnaslot@meta.data)
pdf("result/4.拟时序/巨噬细胞亚型.pdf",width = 5,height = 5)
DimPlot(subset_rnaslot,group.by ="cell.type",label = TRUE )
dev.off()
library(SCP)
CellDimPlot(subset_rnaslot, group.by = "cell.type", reduction = "UMAP", theme_use = "theme_blank")
pdf("result/4.拟时序/巨噬细胞亚型.pdf",width = 5,height = 5)
CellDimPlot(subset_rnaslot, group.by = "cell.type", reduction = "UMAP", 
            theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16))
dev.off()
saveRDS(subset_rnaslot, file = "result/4.拟时序/巨噬细胞亚型.rds")
#subset_rnaslot <- readRDS("result/4.拟时序/巨噬细胞亚型.rds")
#slingshot
subset_rnaslot.sce=as.SingleCellExperiment(subset_rnaslot)
head(subset_rnaslot@meta.data)
sds <- slingshot(subset_rnaslot.sce, #这里是SingleCellExperiment单细胞对象
                            reducedDim = 'UMAP', #降维方式
                            clusterLabels = subset_rnaslot.sce$cell.type,#celltype
                 start.clus= "MKi67",  # 选择起点
                 end.clus = NULL
                 ) 
colnames(colData(sds))
seu=subset_rnaslot
sds
summary(sds$slingPseudotime_1)
#点的
pdf("result/4.拟时序/slingshot-1.pdf",width = 5,height = 5)
plot(reducedDims(sds)$UMAP,col = pal[sds$cell.type],cex = 0.5, pch=16, asp=1)
lines(SlingshotDataSet(sds), lwd=2,col = 'black',type = 'lineages')
dev.off()
#线的
pdf("result/4.拟时序/slingshot-2.pdf",width = 5,height = 5)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
plot(reducedDims(sds)$UMAP, col = plotcol, pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sds), lwd=2, col='black')
dev.off()
sim <- sds
library(tradeSeq)
counts <- sim@assays@data$counts
crv <- SlingshotDataSet(sim)
set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # no more than 12
                   nGenes = 200, # 每个细胞纳入分析的基因数量，默认是500，这里为了节省示例计算时间
                   verbose = T)
set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
system.time({ 
  sce <- fitGAM(counts = counts, 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = 6, 
                verbose = FALSE)
})
table(rowData(sce)$tradeSeq$converged)
assoRes <- associationTest(sce)
head(assoRes)
startRes <- startVsEndTest(sce)
head(startRes)
# 取celltpye和配色信息
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
plotcol[is.na(plotcol)] <- "lightgrey" 
plotcol
coldata <- data.frame(celltype = sim@colData$cell.type,
                      plotcol = plotcol)
rownames(coldata) = colnames(sim)
# 把sce中的3000个细胞对应信息取出
filter_coldata <- coldata[colnames(sce),]
# 添加拟时序信息
filter_coldata$Pseudotime = sce$crv$pseudotime.Lineage1
# 提取 NLN 的表达量
top6_exp = sce@assays@data$counts["RPN1",  , drop = FALSE]
top6_exp = log2(top6_exp + 1) %>% t()
plt_data = cbind(filter_coldata, top6_exp)
colnames(plt_data)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors = getPalette(length(unique(plt_data$celltype)))

gene = 'RPN1'
p_test = ggscatter(data = plt_data,
                   x = 'Pseudotime',
                   y = gene,
                   color = 'celltype',
                   size = 0.6)+
  geom_smooth(se = F, color = 'orange')+
  theme_bw()+
  scale_color_manual(values = mycolors)+
  theme(legend.position = 'right')+
  guides(color = guide_legend(ncol = 1,  # 将图例分为两列显示
                              override.aes = list(size = 3)) )  # 调整图例展示的点大小
p_test
ggsave("result/4.拟时序/RPN1拟时序.pdf", p_test, width = 4, height = 3)
####cellchat-前期准备####
soj <- readRDS("result/1.降维/注释后seurat.rds")
macrophage <- subset(soj, cellType == "Macrophage")
macrophage <- NormalizeData(macrophage, normalization.method = "LogNormalize", scale.factor = 1e4)
# 按照RPN1表达将巨噬细胞分为high和low两组
RPN1 <- GetAssayData(macrophage, layer = "data")["RPN1", ]
threshold <- median(RPN1)
macrophage$RPN1_group <- ifelse(RPN1 > threshold, "high", "low")
# 创建新的cellType列，将巨噬细胞标记为RPN1-high或RPN1-low
macrophage$cellType <- paste0("Macrophage RPN1-", macrophage$RPN1_group)
# 提取非巨噬细胞
non_macrophage <- subset(soj, cellType != "Macrophage")
# 确保两个数据集有相同的列
macrophage$cellType <- as.factor(macrophage$cellType)
# 合并数据集
merged_soj <- merge(non_macrophage, macrophage)
table(merged_soj$cellType)
####cellchat####
#教程
#https://zhuanlan.zhihu.com/p/652379409
#https://blog.csdn.net/zfyyzhys/article/details/142408082
#https://mp.weixin.qq.com/s?__biz=Mzk0NDY1MTc5OQ==&mid=2247484634&idx=1&sn=cf37978c01e38e15788d1b21bcfdd5ec&chksm=c3202fa2f457a6b451d72ae4108f43efb30801192f79904a9bcc0966b36e055f3ee7b505941f&cur_album_id=3370419981043924993&scene=189#wechat_redirect
#https://mp.weixin.qq.com/s?__biz=Mzk0NDY1MTc5OQ==&mid=2247484634&idx=1&sn=cf37978c01e38e15788d1b21bcfdd5ec&chksm=c3202fa2f457a6b451d72ae4108f43efb30801192f79904a9bcc0966b36e055f3ee7b505941f&cur_album_id=3370419981043924993&scene=189#wechat_redirect
#https://mp.weixin.qq.com/s/d2LXP3joJDgQwDvjSob5vQ
#devtools::install_github("sqjin/CellChat")
library(Seurat)
library(CellChat) 
library(cowplot)
sce <- merged_soj
#创建CellChat 对象
# 设置默认的分析对象为 SCT
DefaultAssay(sce) <- 'SCT'
sce_input <- GetAssayData(sce, layer = 'data')
sce_meta <- sce@meta.data[,c("orig.ident","cellType")]
colnames(sce_meta) <- c("group", "labels")
identical(colnames(sce_input), rownames(sce_meta))
cellchat <- createCellChat(object = sce_input, meta = sce_meta, group.by = "labels")
levels(cellchat@idents)# 显示细胞标签的因子水平
groupSize <- as.numeric(table(cellchat@idents)) # 计算每个细胞组的细胞数量
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  # 即使使用整个数据库，这一步也是必要的
future::plan("multisession", workers = 4)  # 并行计算
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
df.net <- subsetCommunication(cellchat)
saveRDS(cellchat, file = "result/6.cellchat/cellchat.rds")
cellchat <- readRDS("result/6.cellchat/cellchat.rds")
pdf("result/6.cellchat/1.细胞通讯网络-细胞交互次数.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
dev.off()
pdf("result/6.cellchat/2.细胞通讯网络-细胞交互权重.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
dev.off()
#比较不同细胞通讯网络之间的权重
mat <- cellchat@net$weight
pdf("result/6.cellchat/3.不同细胞通讯网络权重.pdf")
par(mfrow = c(2,3),xpd= TRUE)
for(i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,edge.weight.max = max(mat),title.name= rownames(mat)[i])
}
dev.off()
#使用热图更详细地显示相互作用的不同数量或相互作用强度
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf("result/6.cellchat/4.细胞通讯网络-细胞交互热图.pdf",width = 8,height = 6)
gg1 + gg2
dev.off()
#系统分析细胞间通讯网络计算并可视化网络中心性得分
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

pdf("result/6.cellchat/5.发送者和接收者散点图.pdf")
netAnalysis_signalingRole_scatter(cellchat)
dev.off()
#识别对某些细胞群体出向或入向信号传递贡献最大的信号
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming")
pdf("result/6.cellchat/6.信号通路在细胞通讯中的作用.pdf",width = 12,height = 6)
ht1 + ht2
dev.off()

levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(5,6), 
                 targets.use = c(1:4,7:14), remove.isolate = FALSE)
ggsave("result/6.cellchat/7.macrophage-source.pdf",width = 6,height = 8)
#反过来
source_indices <- c(1:4, 7:14) #合并两个向量
netVisual_bubble(cellchat, sources.use = source_indices, 
                 targets.use = c(5,6), remove.isolate = FALSE)
ggsave("result/6.cellchat/7.macrophage-target.pdf",width = 6,height = 8)

cellchat@netP$pathways
pathways.show <- c("TGFb") 
vertex.receiver = seq(1,5) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver)
netAnalysis_contribution(cellchat, signaling = pathways.show)

####RPN1-high-low富集分析####
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(scRNAtoolVis)
library(ClusterGVis)
library(org.Hs.eg.db)
library(future)
seurat_merge_v5 <- merged_soj
table(seurat_merge_v5$cellType)
seurat_merge_v5 <- subset(seurat_merge_v5, 
                            subset = cellType %in% c("Macrophage RPN1-high", "Macrophage RPN1-low"))
table(seurat_merge_v5$cellType)
# Macrophage RPN1-high  Macrophage RPN1-low
# 1671                 3384
seurat_merge_v5 <- NormalizeData(seurat_merge_v5,normalization.method = "LogNormalize",
                             scale.factor = 1e4)
GetAssay(seurat_merge_v5,assay = "RNA")
seurat_merge_v5 <- FindVariableFeatures(seurat_merge_v5,
                                    selection.method = "vst",nfeatures = 3000)
seurat_merge_v5 <- ScaleData(seurat_merge_v5)
seurat_merge_v5 <- RunPCA(object = seurat_merge_v5,dims = 1:110,pc.genes = VariableFeatures(sce_myeloid))
seurat_merge_v5 <- RunHarmony(seurat_merge_v5, "orig.ident")
names(seurat_merge_v5@reductions)
seurat_merge_v5 <- RunUMAP(seurat_merge_v5,  dims = 1:20, 
                       reduction = "harmony")
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5,dims = 1:20,reduction = "harmony")
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.5)
DimPlot(seurat_merge_v5,label = T,pt.size =2)
Idents(seurat_merge_v5) <- seurat_merge_v5$cellType
seurat_merge_v5 <- PrepSCTFindMarkers(seurat_merge_v5)
plan("multicore", workers = 8)  # 设置 8 核心并行，按自己电脑调整核心数
options(future.globals.maxSize = 20 * 1024^3)  # 允许最大内存，比如 20GB
pbmc.markers.all <- Seurat::FindAllMarkers(
  seurat_merge_v5,
  only.pos = FALSE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)
plan("sequential")  # 恢复到单线程模式
col <- c(
  "Macrophage RPN1-high" = "#d68c86", 
  "Macrophage RPN1-low" = "#8baad3"
)
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 4, wt = avg_log2FC)
st.data <- prepareDataFromscRNA(object = seurat_merge_v5,
                                diffData = pbmc.markers,
                                showAverage = TRUE)
lapply(unique(pbmc.markers$cluster), function(x){
  tmp <- pbmc.markers |> dplyr::filter(cluster == x)
  
  # plot
  p <- Seurat::FeaturePlot(object = seurat_merge_v5,
                           features = tmp$gene,
                           ncol = 4)
  
  return(p)
}) -> gglist
names(gglist) <- paste("C",1:2,sep = "")
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
head(enrich,3)
# group                             Description       pvalue ratio
# GO:1905039    C1 carboxylic acid transmembrane transport 0.0003981575    50
# GO:1903825    C1    organic acid transmembrane transport 0.0004135812    50
# GO:0098739    C1           import across plasma membrane 0.0006929247    50
enrich.KEGG <- enrichCluster(object = st.data,
                             OrgDb = org.Hs.eg.db,
                             type = "KEGG",
                             organism = "hsa",
                             pvalueCutoff = 0.05,
                             topn = 5,
                             seed = 5201314)
head(enrich.KEGG,3)
palette = c("Purples2","Purples3","Blues2","Blues3","Reds2","Reds3","Greens2")
# loop
lapply(seq_along(unique(enrich$group)), function(x){
  # go plot
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description,levels = tmp$Description)
  
  # plot
  p <-
    ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue),y = Description,fill = -log10(pvalue)),
             width = 0.75) +
    geom_line(aes(x = log10(ratio),y = as.numeric(Description)),color = "grey50") +
    geom_point(aes(x = log10(ratio),y = Description),size = 3,color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~.,name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  # plot kegg
  tmp.kg <- enrich.KEGG |> dplyr::filter(group == unique(enrich.KEGG$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp.kg$Description <- factor(tmp.kg$Description,levels = tmp.kg$Description)
  
  # plot
  pk <-
    ggplot(tmp.kg) +
    geom_segment(aes(x = 0,xend = -log10(pvalue),y = Description,yend = Description),
                 lty = "dashed",linewidth = 0.75) +
    geom_point(aes(x = -log10(pvalue),y = Description,color = -log10(pvalue)),size = 5) +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    colorspace::scale_color_binned_sequential(palette = palette[x]) +
    ylab("") + xlab("-log10(pvalue)")
  
  # combine
  cb <- cowplot::plot_grid(plotlist = list(p,pk))
  
  return(cb)
}) -> gglist
# assign names
names(gglist) <- paste("C",1:2,sep = "")
# insert bar plot
pdf('result/6.cellchat/gokegg.pdf',height = 6,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           line.side = "left",
           column_names_rot = 45,
           markGenes = pbmc.markers$gene,
           cluster.order = c(1:2),
           ggplot.panel.arg = c(5,0.5,32,"grey90",NA),
           gglist = gglist)
dev.off()

####单基因boxplot####
##tumor normal
library(RColorBrewer)
library(ggsignif)
combined <- readRDS("result/1.降维/注释后seurat.rds")
tucom@assays$RNA$data[1:4,1:4]
# 提取 PDAC 样本
tucom <- combined[, combined@meta.data$tissue == "PDAC"]
tuc <- data.frame(
  expression = as.numeric(tucom@assays$RNA$data["PML", ]),
  type = "PDAC"
)
# 提取 ADJ 样本
nocom <- combined[, combined@meta.data$tissue == "ADJ"]
noc <- data.frame(
  expression = as.numeric(nocom@assays$RNA$data["PML", ]),
  type = "ADJ"
)
# 合并数据
subset_data <- rbind(tuc, noc)
head(subset_data)
display.brewer.all()
mycol <- brewer.pal(2,"Accent")
subset_data0 <- subset_data[subset_data$expression > 0, ]
subset_data0$log_expression <- log1p(subset_data0$expression)
wilcox_res <- wilcox.test(log_expression ~ type, data = subset_data0)#计算 p 值（Wilcoxon 检验）
p_value <- wilcox_res$p.value
#自定义 p 值标签
p_label <- ifelse(p_value < 0.001, "p < 0.001", 
                  ifelse(p_value < 0.05, "p < 0.05", 
                         paste0("p = ", format(p_value, digits = 3))))
p <- ggplot(subset_data0, aes(x = type, y = log_expression, fill = type)) +
  geom_violin(trim = TRUE, alpha = 0.8) +  # 小提琴图
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "white") +  # 箱线图
  geom_signif(comparisons = list(c("ADJ", "PDAC")), 
              annotations = p_label, 
              textsize = 5, step_increase = 0.1) +  # 直接标注 p 值
  scale_fill_manual(values = c("ADJ" = "#7FC97F", "PDAC" = "#BEAED4")) +  # 自定义颜色
  theme_bw() +
  ggtitle("PML")
pdf("result/5.箱线/PDAC_ADJ_PML.pdf",width = 5,height = 5)
p
dev.off()
####基因箱线图####
com <- readRDS("result/1.降维/注释后seurat.rds")
com@assays$RNA$data[1:4,1:4]
table(com@meta.data$cellType)
# T cell                   B cell               Fibroblast
# Macrophage         Endothelial cell              Acinar cell
# Pancreatic ductal cell               Neitrophil Pancreatic stellate cell
# Mast cell              Plasma cell                  NK cell
# Schwann cell
tcom <- com[,com@meta.data$cellType == "T cell"]
tcell <- as.data.frame(t(as.matrix(tcom@assays$RNA$data["PML", , drop = FALSE])))
tcell$type <- "T cell"
tcell[1:4,]

bcom <- com[,com@meta.data$cellType == "B cell"]
bcell <- as.data.frame(t(as.matrix(bcom@assays$RNA$data["PML", , drop = FALSE])))
class(bcell)
bcell$type <- "B cell"
bcell[1:4,]

fcom <- com[,com@meta.data$cellType == "Fibroblast"]
fcell <- as.data.frame(t(as.matrix(fcom@assays$RNA$data["PML", , drop = FALSE])))
class(fcell)
fcell$type <- "Fibroblast"
fcell[1:4,]

Mcom <- com[,com@meta.data$cellType == "Macrophage"]
Mcell <- as.data.frame(t(as.matrix(Mcom@assays$RNA$data["PML", , drop = FALSE])))
class(Mcell)
Mcell$type <- "Macrophage"
Mcell[1:4,]

Ecom <- com[,com@meta.data$cellType == "Endothelial cell"]
Ecell <- as.data.frame(t(as.matrix(Ecom@assays$RNA$data["PML", , drop = FALSE])))
class(Ecell)
Ecell$type <- "Endothelial cell"
Ecell[1:4,]

Acom <- com[,com@meta.data$cellType == "Acinar cell"]
Acell <- as.data.frame(t(as.matrix(Acom@assays$RNA$data["PML", , drop = FALSE])))
class(Acell)
Acell$type <- "Acinar cell"
Acell[1:4,]

Pcom <- com[,com@meta.data$cellType == "Pancreatic ductal cell"]
Pcell <- as.data.frame(t(as.matrix(Pcom@assays$RNA$data["PML", , drop = FALSE])))
class(Pcell)
Pcell$type <- "Pancreatic ductal cell"
Pcell[1:4,]

Ncom <- com[,com@meta.data$cellType == "Neitrophil"]
Ncell <- as.data.frame(t(as.matrix(Ncom@assays$RNA$data["PML", , drop = FALSE])))
class(Ncell)
Ncell$type <- "Neitrophil"
Ncell[1:4,]

pcom <- com[,com@meta.data$cellType == "Pancreatic stellate cell"]
pcell <- as.data.frame(t(as.matrix(pcom@assays$RNA$data["PML", , drop = FALSE])))
class(pcell)
pcell$type <- "Pancreatic stellate cell"
pcell[1:4,]

mcom <- com[,com@meta.data$cellType == "Mast cell"]
mcell <- as.data.frame(t(as.matrix(mcom@assays$RNA$data["PML", , drop = FALSE])))
class(mcell)
mcell$type <- "Mast cell"
mcell[1:4,]

ppcom <- com[,com@meta.data$cellType == "Plasma cell"]
ppcell <- as.data.frame(t(as.matrix(ppcom@assays$RNA$data["PML", , drop = FALSE])))
class(ppcell)
ppcell$type <- "Plasma cell"
ppcell[1:4,]

ncom <- com[,com@meta.data$cellType == "NK cell"]
ncell <- as.data.frame(t(as.matrix(ncom@assays$RNA$data["PML", , drop = FALSE])))
class(ncell)
ncell$type <- "NK cell"
ncell[1:4,]

scom <- com[,com@meta.data$cellType == "Schwann cell"]
scell <- as.data.frame(t(as.matrix(scom@assays$RNA$data["PML", , drop = FALSE])))
class(scell)
scell$type <- "Schwann cell"
scell[1:4,]

data <- rbind(tcell,bcell,fcell,Mcell,Ecell,Acell,Pcell,Ncell,pcell,mcell,ppcell,ncell,scell)
data1 <- data
data <- data[data$PML > 0, ]
data$PML <- log1p(data$PML)
library(ggstatsplot) # install.packages("ggstatsplot")
library(ggpubr)
col <- c(
  "Acinar cell" = "#FCCDE5",
  "B cell" = "#B3DE69",
  "Endothelial cell" = "#FB8072",
  "Fibroblast" = "#FFFFB3",
  "Macrophage" = "#8DD3C7", 
  "Mast cell" = "#33A02C",
  "Neitrophil" = "#80B1D3",
  "NK cell" = "#CCEBC5",
  "Pancreatic ductal cell" = "#a44e89",
  "Pancreatic stellate cell" = "#FDB462",
  "Plasma cell" = "#1F78B4",
  "Schwann cell" ="#D9D9D9",
  "T cell" = "#BEBADA"
)
pdf("result/5.箱线/PML-statsplot.pdf",width = 15,height = 12)
p <- ggbetweenstats(
  data = data,
  x = type,
  y = PML,
  title = "PML expression",
  messages = FALSE,
  palette = "Set3"
) 
p + scale_color_manual(values = col) + scale_fill_manual(values = col)
dev.off()

####巨噬细胞极化####
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(scRNAtoolVis)
library(ClusterGVis)
library(org.Hs.eg.db)
library(future)
soj <- readRDS("result/1.降维/注释后seurat.rds")
macrophage <- subset(soj, cellType == "Macrophage")
macrophage <- NormalizeData(macrophage, normalization.method = "LogNormalize", scale.factor = 1e4)
RPN1 <- GetAssayData(macrophage, layer = "data")["RPN1", ]
threshold <- median(RPN1)
macrophage$RPN1_group <- ifelse(RPN1 > threshold, "high", "low")
macrophage$cellType <- paste0("Macrophage RPN1-", macrophage$RPN1_group)
non_macrophage <- subset(soj, cellType != "Macrophage")
macrophage$cellType <- as.factor(macrophage$cellType)
merged_soj <- merge(non_macrophage, macrophage)
table(merged_soj$cellType)
seurat_merge_v5 <- merged_soj
table(seurat_merge_v5$cellType)
seurat_merge_v5 <- subset(seurat_merge_v5, 
                          subset = cellType %in% c("Macrophage RPN1-high","Macrophage RPN1-low"))
table(seurat_merge_v5$cellType)
# Macrophage RPN1-high  Macrophage RPN1-low
# 1671                 3384
seurat_merge_v5 <- NormalizeData(seurat_merge_v5,normalization.method = "LogNormalize",
                                 scale.factor = 1e4)
GetAssay(seurat_merge_v5,assay = "RNA")
seurat_merge_v5 <- FindVariableFeatures(seurat_merge_v5,
                                        selection.method = "vst",nfeatures = 3000)
seurat_merge_v5 <- ScaleData(seurat_merge_v5)
seurat_merge_v5 <- RunPCA(object = seurat_merge_v5,dims = 1:110,pc.genes = VariableFeatures(sce_myeloid))
seurat_merge_v5 <- RunHarmony(seurat_merge_v5, "orig.ident")
names(seurat_merge_v5@reductions)
seurat_merge_v5 <- RunUMAP(seurat_merge_v5,  dims = 1:20, 
                           reduction = "harmony")
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5,dims = 1:20,reduction = "harmony")
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.5)
DimPlot(seurat_merge_v5,group.by = "cellType",label = T,pt.size =2)
high <- subset(seurat_merge_v5, 
                          subset = cellType %in% c("Macrophage RPN1-high"))
high.markers <- FindAllMarkers(
  high, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.585, 
  p_val_adj = 0.05  # 设置校正后的 p 值小于 0.05
)
head(high.markers)
write.csv(high.markers,file = "result/7.巨噬细胞极化/RPN1-high-marker基因.csv")

####敲低####
# library(remotes)
# install_github('cailab-tamu/scTenifoldKnk')
library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)  
setwd("~/data/GSE212966_RAW/")
af <- readRDS("result/1.降维/注释后seurat.rds")
macrophage <- subset(af, cellType == "Macrophage")
countMatrix <- GetAssayData(macrophage, slot = "counts")
dim(countMatrix)
keep_genes <- rowSums(countMatrix > 0) >= 1000
count_filt <- countMatrix[keep_genes, ]
dim(count_filt)
"RPN1" %in% rownames(count_filt)
result <- scTenifoldKnk(countMatrix = count_filt, 
                        gKO = 'RPN1', #需要敲除的基因
                        #qc = TRUE,#是否进行QC
                        qc_mtThreshold = 0.1,#mt阈值
                        qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
)
write.csv(result$diffRegulation,'result/7.敲低/diffRegulation.csv')
top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 20)
ggplot(top_genes, aes(x=reorder(gene, FC), y=log10(FC))) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes (log10 FC)",
       x="Gene", y="log10(FC)")
sig <- subset(result$diffRegulation, p.adj < 0.05)
top20 <- head(sig[order(-sig$FC), ], 20)
p <- ggplot(top20, aes(x=reorder(gene, FC), y=log10(FC))) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes (p.adj < 0.05)",
       x="Gene", y="log10(FC)") +
  theme_minimal()
ggsave("result/7.敲低/敲除基因排名带RPN1.pdf",plot=p,width=5,height=5)
df1 <- subset(top20, gene != "RPN1")
p <- ggplot(df1, aes(x=reorder(gene, FC), y=log10(FC))) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes (p.adj < 0.05)",
       x="Gene", y="log10(FC)") +
  theme_minimal()
ggsave("result/7.敲低/敲除基因排名不带RPN1.pdf",plot=p,width=5,height=5)
df <- result$diffRegulation
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.05)
p2 <- ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value)",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()
ggsave("result/7.敲低/Z分数带RPN1.pdf",plot=p2,width=5,height=5)
df <- subset(result$diffRegulation, gene != "RPN1")
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.05)
p2 <- ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value)",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()
ggsave("result/7.敲低/Z分数不带RPN1.pdf",plot=p2,width=5,height=5)
library(openxlsx)          # 读取Excel文件
library(dplyr)             # 数据操作
library(clusterProfiler)   # GO富集分析
library(org.Hs.eg.db)      # 人类基因注释数据库
library(ggplot2)           # 绘图
genelist <- read.xlsx("result/7.敲低/genelist.xlsx", colNames = F) %>% pull(1)
go_res <- enrichGO(
  genelist,                # 输入的基因符号列表
  OrgDb = org.Hs.eg.db,    # 人类基因数据库
  keyType = "SYMBOL",      # 输入基因类型为符号（Symbol）
  ont = "ALL",             # 同时分析BP/CC/MF三类
  pvalueCutoff = 0.05,      # p值阈值（宽松）
  pAdjustMethod = "BH",    # p值校正方法（Benjamini-Hochberg）
  qvalueCutoff = 0.05,      # q值阈值（宽松）
  minGSSize = 3,           # 最小基因集大小
  maxGSSize = 500          # 最大基因集大小
)
goresult <- data.frame(go_res@result)  # 转换为数据框
write.table(goresult, "result/7.敲低/goresult.txt", sep="\t", row.names=F, quote=F)  # 保存表格
colos <- c('#7f4ea8','#009eff','#ffa61d')
classes <- c('BP','CC','MF')
for (i in 1:3) {
  
  class <- classes[i]
  colo <- colos[i]
  
  go_select <- goresult %>%
    filter(ONTOLOGY == class) %>%
    arrange(pvalue) %>% 
    head(10) %>%
    dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  go_select$GeneRatio <- sapply(
    strsplit(go_select$GeneRatio, "/"),
    function(x) as.numeric(x[1]) / as.numeric(x[2])
  )
  
  negline <- max(-log10(go_select$pvalue)) / 4
  gr_range <- range(go_select$GeneRatio)
  
  go_select$scaled_gr <- scales::rescale(
    go_select$GeneRatio,
    to = c(-negline*3/4, -negline*1/4),
    from = gr_range
  )
}
library(dplyr)
library(stringr)
p <- ggplot(go_select) +
  geom_path(aes(x = scaled_gr, y = Description, group = 1), # group = 1
            color = "gray40", 
            linewidth = 0.8,
            lineend = "round") +  
  geom_point(aes(x = scaled_gr, y = Description, size = Count),
             color = colo) +
  
  geom_bar(aes(x = -log10(pvalue), y = Description),
           stat = "identity",
           fill = paste0(colo,'50'),
           width = 0.8, 
           position = position_nudge(y = 0.2)) + 
  geom_text(
    aes(x = 0.05, y = Description, label = str_wrap(Description, width = 40)),
    hjust = 0, vjust = 1, size = 4.5, lineheight = 0.8, 
    nudge_y = 0.35, color = "black"
  )+
  
  scale_x_continuous(
    name = 'GeneRatio(%) and -Log10(pvalue)',
    limits = c(-negline, max(-log10(go_select$pvalue)) * 1.01),
    expand = c(0, 0)
  ) +
  scale_size_continuous(
    name = "Gene Count",
    range = c(3, 6)) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12), 
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 13), 
    axis.text.x.bottom = element_text(color = "black")
  ) +
  labs(x = NULL, title = class)

print(p)

ggsave(paste0('result/7.敲低/',class,'.pdf'), p, width = 5, height = 5)
