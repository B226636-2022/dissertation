library(Seurat)
library(scMerge)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(corrplot)
#function
coefficient_of_variation <- function(x) {
  cv <- sd(x) / mean(x)
  return(cv)
}

david_enrichment_visualization <- function(enrichment_file){
  enrichment_data <- read.csv(enrichment_file)
  ggplot(enrichment_data, aes(x = Term, y = -log10(PValue), fill = Fold.Enrichment)) +
    geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 1)) +
    facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +
    labs(x = "Term", y = "-log10(P-value)", title = "DAVID Enrichment Analysis") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),strip.text.x = element_text(size = 14)) +
    scale_fill_viridis_c(option = "plasma")  # 使用viridis颜色映射
}

#start to load data, the express matrix and the SEG list
single_es_tximport <- readRDS("~/dissertation/count/single_es_tximport.rds")
data("segList_ensemblGeneID", package = "scMerge")
SEG_list_sc <- segList_ensemblGeneID$mouse$mouse_scSEG
SEG_list_micro <- segList_ensemblGeneID$mouse$bulkMicroarrayHK
SEG_list_bulk <- segList_ensemblGeneID$mouse$bulkRNAseqHK
SEG_common<-intersect(intersect(SEG_list_bulk,SEG_list_micro),SEG_list_sc)

#start to get good cells by spike-in method
sce <- SingleCellExperiment(assays=list(counts=single_es_tximport$counts,abundance=single_es_tximport$abundance))
is.spike <- grepl("^ERCC", rownames(sce))
#This list of identifiers was obtained from ENSEMBL
mitogenes="ENSMUSG00000064336|ENSMUSG00000064337|ENSMUSG00000064338|ENSMUSG00000064339|ENSMUSG00000064340|ENSMUSG00000064341|ENSMUSG00000064342|ENSMUSG00000064343|ENSMUSG00000064344|ENSMUSG00000064345|ENSMUSG00000064346|ENSMUSG00000064347|ENSMUSG00000064348|ENSMUSG00000064349|ENSMUSG00000064350|ENSMUSG00000064351|ENSMUSG00000064352|ENSMUSG00000064353|ENSMUSG00000064354|ENSMUSG00000064355|ENSMUSG00000064356|ENSMUSG00000064357|ENSMUSG00000064358|ENSMUSG00000064359|ENSMUSG00000064360|ENSMUSG00000064361|ENSMUSG00000065947|ENSMUSG00000064363|ENSMUSG00000064364|ENSMUSG00000064365|ENSMUSG00000064366|ENSMUSG00000064367|ENSMUSG00000064368|ENSMUSG00000064369|ENSMUSG00000064370|ENSMUSG00000064371|ENSMUSG00000064372"
is.mito <- grepl(mitogenes, rownames(sce))
sce$stats <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
par(mfrow=c(1,2))
hist(sce$stats$sum/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$stats$detected, xlab="Number of detected genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
par(mfrow=c(1,2))
hist(sce$stats$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$stats$subsets_ERCC_percent, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
libsize.drop <- isOutlier(sce$stats$sum, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$stats$detected, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce$stats$subsets_Mt_percent, nmads=3, type="higher")
spike.drop <- isOutlier(sce$stats$subsets_ERCC_percent, nmads=3, type="higher")
#copy for later
sceX <-sce
sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
stable_cell<-colnames(sce)

express_matrix <- as.data.frame(single_es_tximport$counts)
seurat_obj_raw <- CreateSeuratObject(counts = express_matrix)
annotation1<-read.csv("anno2.csv",header = TRUE,row.names = "ERR")
seurat_obj_raw <- AddMetaData(seurat_obj_raw, metadata = annotation1)
seurat_obj_raw_fil <- seurat_obj_raw[,stable_cell]
seurat_obj_raw_fil <- NormalizeData(seurat_obj_raw_fil)
seurat_obj_raw_fil <- ScaleData(seurat_obj_raw_fil) 
seurat_obj_raw_fil <- FindVariableFeatures(seurat_obj_raw_fil, selection.method = "vst")
seurat_obj_raw_fil<- RunPCA(object = seurat_obj_raw_fil, features = rownames(seurat_obj_raw_fil))
p1<-DimPlot(seurat_obj_raw_fil,reduction= "pca",group.by = "Type")+ ggtitle("all")
seurat_seg_all_fil <- subset(seurat_obj_raw_fil, features=SEG_list_sc)
seurat_seg_all_fil <- RunPCA(seurat_seg_all_fil)
p2<-DimPlot(seurat_seg_all_fil, reduction = "pca",group.by = "Type")+ ggtitle("SEG_all")
pca_compare <- cowplot::plot_grid(p1, p2, ncol = 2)


#normalize the seurat_obj_raw and subset it for later use 
seurat_obj_raw <- NormalizeData(seurat_obj_raw)
seurat_obj_raw <- ScaleData(seurat_obj_raw)
seurat_obj_raw <- FindVariableFeatures(seurat_obj_raw, selection.method = "vst")
variable_gene <- VariableFeatures(seurat_obj_raw)
seurat_seg_all <- subset(seurat_obj_raw, features=SEG_list_sc)
seurat_seg_serum <- subset(seurat_seg_all, Type == "serum")
seurat_seg_2i <- subset(seurat_seg_all, Type == "2i")
seurat_seg_a2i <- subset(seurat_seg_all, Type == "a2i")


#start to draw the cv compare method for SEG and non-SEG
seg_cv <- apply(as.matrix(seurat_seg_all@assays$RNA@data), 1, coefficient_of_variation)
remaining_genes <- setdiff(rownames(seurat_obj_raw), SEG_list_sc)
reference_genes <- sample(remaining_genes, size = 722)
reference_expr <- as.data.frame(seurat_obj_raw@assays$RNA@data)
reference_expr<-reference_expr[reference_genes, ]
nonseg_cv <- apply(as.matrix(reference_expr), 1, coefficient_of_variation)
seg_plot <- data.frame(Gene = rownames(seurat_seg_all@assays$RNA@data), CV = seg_cv)
nonseg_plot <- data.frame(Gene = rownames(reference_expr),CV = nonseg_cv)
#also can use plot to finish this step
p1 <- ggplot(seg_plot, aes(x = Gene, y = CV)) +
  geom_point() +
  labs(title = "CV of SEG genes") +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  ylim(0, 30)
p2 <- ggplot(nonseg_plot, aes(x = Gene, y = CV)) +
  geom_point() +
  labs(title = "CV of other genes") +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  ylim(0, 30)
# 组合两个子图
cv_compare <- cowplot::plot_grid(p1, p2, ncol = 2)

#use SEG to explore outlying cells
##begin cluster
seurat_seg_all <- RunPCA(seurat_seg_all)
pca_co<-as.data.frame(Embeddings(seurat_seg_all, reduction = "pca"))
DimPlot(seurat_seg_all, reduction = "pca",group.by = "Type")+ggtitle("PCA for all cells(feature = SEG)")
#seurat_seg_all <- FindNeighbors(seurat_seg_all)
#seurat_seg_all <- FindClusters(seurat_seg_all, resolution = 0.15)
seurat_seg_all$seurat_clusters <- ifelse(pca_co$PC_1 > -6 & pca_co$PC_1 <= 4, "cluster0",
                             ifelse(pca_co$PC_1 > 4,"cluster1", "cluster2"))
seurat_seg_all@active.ident<- factor(seurat_seg_all$seurat_clusters)

DimPlot(seurat_seg_all, reduction = "pca",group.by = "cluster")+ ggtitle("Cluster for all cells")+
  geom_vline(xintercept = 4, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = -6, linetype = "dashed", color = "blue")

#find markers for each cluster
#markers <- FindMarkers(seurat_seg_all, ident.1 = "cluster0", ident.2 = "cluster1", min.pct = 0.25, logfc.threshold = 0.25)
markers <-FindAllMarkers(seurat_seg_all)
#markers <- FindMarkers(seurat_seg_all, ident.1 = 0, ident.2 = 1, min.pct = 0.25, logfc.threshold = 0.25)
#significant_markers <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 1 & markers$pct.1>0.9 & markers$pct.2>0.9, ]
significant_markers1 <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 1 & markers$pct.1>0.9 & markers$pct.2>0.9, ]
top5marker<-significant_markers1[significant_markers1["cluster"] == "cluster0",] %>% top_n(n = 5, wt = abs(avg_log2FC))
DotPlot(seurat_seg_all, features = unique(top5marker$gene), cluster.idents = T) + RotatedAxis()

#对marker基因进行功能注释
##load annotation result from david
david_enrichment_visualization("gocluster.csv")

#heatmap
#choose genes from sumo2
sumo2genes<-trimws(unlist(strsplit(as.character(enrichment_data[1,]$Genes), ",")))

# 创建一个biomaRt连接，选择老鼠的Ensembl数据库
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 使用getBM函数进行转换
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = sumo2genes,
                   mart = ensembl)
seurat_sumo2=subset(seurat_seg_all, features = sumo2genes)
#sumo2_matrix <- as.matrix(seurat_seg_all@assays$RNA@data[sumo2genes, ])
gene_name<-gene_info[match(sumo2genes, gene_info$ensembl_gene_id),]$external_gene_name
#rownames(sumo2_matrix)<-gene_name
#seurat_sumo2<-CreateSeuratObject(sumo2_matrix)
#seurat_sumo2$seurat_clusters <- seurat_seg_all$seurat_clusters
#seurat_sumo2@active.ident<- factor(seurat_sumo2$seurat_clusters)
DoHeatmap(seurat_seg_all,features = sumo2genes)+annotate("text", x = Inf+1, y = c(25:1),check_overlap = TRUE, label = gene_name, hjust = 0.2, vjust = 0.5, size = 2,parse = TRUE, fontface = "bold")
DoHeatmap(seurat_seg_all,features = sumo2genes, group.by = "Type")+annotate("text", x = Inf+1, y = c(25:1),check_overlap = TRUE, label = gene_name, hjust = 0.2, vjust = 0.5, size = 2,parse = TRUE)

#start to calculate the nono expression

nono_expression<-data.frame(counts_raw=as.data.frame(seurat_obj_raw[["RNA"]]@counts["ENSMUSG00000031311", ]),total_count=seurat_obj_raw@meta.data$nCount_RNA)
colnames(nono_expression) <- c("counts_raw","total_count")
nono_expression$cpm<-nono_expression$counts_raw/nono_expression$total_count*100000
barplot(nono_expression$cpm, main = "Nono Expression Level (Counts Per Million)",xlab = "Sample", ylab = "CPM", col = "skyblue")

gene_expression_normal <- GetAssayData(seurat_obj_raw, slot = "data")
correlation_matrix <- cor(t(as.data.frame(gene_expression_normal)), method = "spearman")
correlation_matrix <- as.data.frame(correlation_matrix)
diag(correlation_matrix) <- 0
target_gene_correlation <- as.data.frame(t(correlation_matrix["ENSMUSG00000031311", ]))
high_correlate<-rownames(target_gene_correlation %>% top_n(n=20, wt=abs(ENSMUSG00000031311)))
high_correlate[21]<-"ENSMUSG00000031311"
target_cor_matrix<-correlation_matrix[high_correlate,high_correlate]
target_cor_matrix<-as.matrix(target_cor_matrix)
rownames(target_cor_matrix)
gene_info_correlate <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = rownames(target_cor_matrix),
                   mart = ensembl)
gene_name_cor<-gene_info_correlate[match(rownames(target_cor_matrix), gene_info_correlate$ensembl_gene_id),]$external_gene_name
rownames(target_cor_matrix)<-gene_name_cor
colnames(target_cor_matrix)<-gene_name_cor
corrplot.mixed(target_cor_matrix,order = 'AOE',tl.cex=0.7)

#visulize david result
david_enrichment_visualization("cor_all.csv")

#sometimes slow, so first filter some zero feature, actually, i run on all feature
#nonzero <- as.matrix(express_matrix > 0)
#keepgenes<-Matrix::rowSums(nonzero)>=10
#filter_express<-rownames(express_matrix[keepgenes,])
seurat_obj_2i<-subset(seurat_obj_raw,Type == "2i")
#seurat_obj_2i<-subset(seurat_obj_2i,feature = filter_express)
gene_expression_normal <- GetAssayData(seurat_obj_2i, slot = "data")
gene_expression_normal <- as.data.frame(gene_expression_normal)
correlation_matrix_2i <- as.data.frame(cor(t(gene_expression_normal), method = "spearman"))
cor_2i_target<-as.data.frame(t(correlation_matrix_2i["ENSMUSG00000031311", ]))
high_correlate_2i<-rownames(cor_2i_target %>% top_n(n=21, wt=abs(ENSMUSG00000031311)))
target_2i_matrix<-as.matrix(correlation_matrix_2i[high_correlate_2i,high_correlate_2i])
gene_info_2i <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                             filters = "ensembl_gene_id",
                             values = rownames(target_2i_matrix),
                             mart = ensembl)
gene_info_2i[21,]<-c("ENSMUSG00000112508","CL57BL6")
gene_name_2i<-gene_info_2i[match(rownames(target_2i_matrix), gene_info_2i$ensembl_gene_id),]$external_gene_name
rownames(target_2i_matrix)<-gene_name_2i
colnames(target_2i_matrix)<-gene_name_2i
corrplot.mixed(target_2i_matrix,order = 'AOE',tl.cex=0.7)
#DAVID the final result
david_enrichment_visualization("cor_2i.csv")


#divide the cells depend on nono expression
#first calculate the nono's cpm distribution
min_threshold <- min(nono_expression$cpm)
max_threshold <- max(nono_expression$cpm)
num_steps <- 200  # 可以根据需要调整步数
# 初始化存储百分比结果的向量
percentages <- numeric(num_steps)
# 循环遍历不同的阈值
output=TRUE
for (i in 1:num_steps) {
  threshold <- min_threshold + (max_threshold - min_threshold) * (i - 1) / (num_steps - 1)
  above_threshold <- sum(nono_expression$cpm >= threshold)
  total <- nrow(nono_expression)
  percentages[i] <- 100 * (total - above_threshold) / total
  if ((percentages[i] > 90) && (output ==TRUE)){
    print(threshold) 
    print(percentages[i])
    output<-FALSE
  }
}
# 绘制可视化曲线
plot(seq(min_threshold, max_threshold, length.out = num_steps), percentages, 
     type = "l", xlab = "Threshold", ylab = "Percentage below threshold",
     main = "Percentage of Data Below Threshold", col = "blue", axes = FALSE)

# 添加坐标轴刻度
axis(1, at = seq(min_threshold, max_threshold, length.out = 5), 
     labels = seq(min_threshold, max_threshold, length.out = 5))
axis(2, at = seq(0, 100, by = 20), labels = seq(0, 100, by = 20))
# 添加坐标轴标题
mtext("Threshold", side = 1, line = 3)
mtext("Percentage below threshold", side = 2, line = 3)


high_subgroup <- nono_expression$cpm >= 83
low_subgroup <- nono_expression$cpm <= 23
#subset the cells according to the nono's expression
# 根据cpm值划分细胞为高、中和低亚群
seurat_obj_raw_nono<-seurat_obj_raw
seurat_seg_all_nono<-seurat_seg_all
seurat_obj_raw_nono@meta.data$Subgroup <- factor(ifelse(high_subgroup, "High",
                                                        ifelse(low_subgroup, "Low", "Medium")))
seurat_seg_all_nono@meta.data$Subgroup <- factor(ifelse(high_subgroup, "High",
                                                        ifelse(low_subgroup, "Low", "Medium")))
DoHeatmap(seurat_seg_all_nono,features = sumo2genes, group.by = "Subgroup")+annotate("text", x = Inf+1, y = c(25:1),check_overlap = TRUE, label = gene_name, hjust = 0.2, vjust = 0.5, size = 2,parse = TRUE, fontface = "bold")
DimPlot(seurat_seg_all_nono, reduction = "pca",group.by = "Subgroup")
seurat_obj_raw_nono@active.ident<- factor(seurat_obj_raw_nono$Subgroup)
try_marker<-FindAllMarkers(seurat_obj_raw_nono)
try_marker[try_marker$p_val_adj < 0.05 & abs(try_marker$avg_log2FC) > 1 & try_marker$pct.1>0.8 & try_marker$pct.2>0.8, ]

##to explore the nono expression in filter dataset
nono_expression_fil<-data.frame(counts_raw=as.data.frame(seurat_obj_raw_fil[["RNA"]]@counts["ENSMUSG00000031311", ]),total_count=seurat_obj_raw_fil@meta.data$nCount_RNA)
colnames(nono_expression_fil) <- c("counts_raw","total_count")
nono_expression_fil$cpm<-nono_expression_fil$counts_raw/nono_expression_fil$total_count*100000
barplot(nono_expression_fil$cpm, main = "Nono Expression Level (Counts Per Million)",xlab = "Sample", ylab = "CPM", col = "skyblue")

#divide the cells depend on nono expression
#first calculate the nono's cpm distribution
min_threshold <- min(nono_expression_fil$cpm)
max_threshold <- max(nono_expression_fil$cpm)
num_steps <- 200  # 可以根据需要调整步数
# 初始化存储百分比结果的向量
percentages <- numeric(num_steps)
# 循环遍历不同的阈值
output=TRUE
for (i in 1:num_steps) {
  threshold <- min_threshold + (max_threshold - min_threshold) * (i - 1) / (num_steps - 1)
  above_threshold <- sum(nono_expression$cpm >= threshold)
  total <- nrow(nono_expression)
  percentages[i] <- 100 * (total - above_threshold) / total
  if ((percentages[i] > 90) && (output ==TRUE)){
    print(threshold) 
    print(percentages[i])
    output<-FALSE
  }
}
# 绘制可视化曲线
plot(seq(min_threshold, max_threshold, length.out = num_steps), percentages, 
     type = "l", xlab = "Threshold", ylab = "Percentage below threshold",
     main = "Percentage of Data Below Threshold", col = "blue", axes = FALSE)

# 添加坐标轴刻度
axis(1, at = seq(min_threshold, max_threshold, length.out = 5), 
     labels = seq(min_threshold, max_threshold, length.out = 5))
axis(2, at = seq(0, 100, by = 20), labels = seq(0, 100, by = 20))
# 添加坐标轴标题
mtext("Threshold", side = 1, line = 3)
mtext("Percentage below threshold", side = 2, line = 3)


high_subgroup <- nono_expression_fil$cpm >= 83
low_subgroup <- nono_expression_fil$cpm <= 23
#subset the cells according to the nono's expression
# 根据cpm值划分细胞为高、中和低亚群
seurat_obj_raw_fil_nono <- seurat_obj_raw_fil
seurat_obj_raw_fil_nono@meta.data$Subgroup <- factor(ifelse(high_subgroup, "High",
                                                            ifelse(low_subgroup, "Low", "Medium")))
DoHeatmap(seurat_obj_raw_fil_nono,features = sumo2genes, group.by = "Subgroup")+annotate("text", x = Inf+1, y = c(25:1),check_overlap = TRUE, label = gene_name, hjust = 0.2, vjust = 0.5, size = 2,parse = TRUE)
DimPlot(seurat_obj_raw_fil_nono, reduction = "pca",group.by = "Subgroup")
seurat_obj_raw_fil_nono@active.ident<- factor(seurat_obj_raw_fil_nono$Subgroup)
try_marker<-FindAllMarkers(seurat_obj_raw_fil_nono)
try_marker[try_marker$p_val_adj < 0.05 & abs(try_marker$avg_log2FC) > 1 & try_marker$pct.1>0.8 & try_marker$pct.2>0.8, ]

#cor
gene_expression_normal_fil <- as.data.frame(GetAssayData(seurat_obj_raw_fil_nono, slot = "data"))
correlation_matrix_fil_nono <- as.data.frame(cor(t(gene_expression_normal_fil), method = "spearman"))
cor_fil_target<-as.data.frame(t(correlation_matrix_fil_nono["ENSMUSG00000031311", ]))
high_correlate_fil<-rownames(cor_fil_target %>% top_n(n=21, wt=abs(ENSMUSG00000031311)))
target_fil_matrix<-as.matrix(correlation_matrix_fil_nono[high_correlate_fil,high_correlate_fil])

gene_info_fil <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = rownames(target_fil_matrix),
                       mart = ensembl)
gene_name_fil<-gene_info_fil[match(rownames(target_fil_matrix), gene_info_fil$ensembl_gene_id),]$external_gene_name
rownames(target_fil_matrix)<-gene_name_fil
colnames(target_fil_matrix)<-gene_name_fil
corrplot.mixed(target_fil_matrix,order = 'AOE',tl.cex=0.7)

david_enrichment_visualization("cor_fil.csv")
