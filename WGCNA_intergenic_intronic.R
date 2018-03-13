#==========================================================================================
#
#    Count table preparation with only transcriptsin intergenic and intronic regions 
#
#==========================================================================================
################### Set working directory ###################
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_data/WGCNA')
## dir()
## Load packages
library(reshape2)
library(ggplot2)

##################### Read in data (counts table and SampleGroup) ########################
## Load in sample data
sample <- read.csv('SampleGroup.csv', header = T, row.names = 1, colClasses = 'factor')

## Load in count table based on TPM normalization
transcripts <- read.table('merge_tpm.txt', header = T, row.names = 1)
intergenic_intronic <- read.table('transcript_intergenic_intronic_bedtools.txt')
intergenic_intronic <- as.vector(as.matrix(intergenic_intronic))
colnames(transcripts) <- paste("DCD002", c(492:511),"SQ", sep = "")
transcripts <- round(transcripts,0) # round reads
summary(transcripts);dim(transcripts)
### Have a look at the transcripts distribution after TPM normalization
group_List <- sample$celltype
exprSet_L <- row.names(transcripts)    
exprSet_L <- cbind(exprSet_L, as.data.frame(transcripts)) # Must convert normalized_count into data.frame before cbind
exprSet_L <- melt(data = exprSet_L, id.vars = 'exprSet_L')
exprSet_L$group <- rep(group_List, each = nrow(transcripts))
#### boxplot
p <- ggplot(data = exprSet_L, aes(x = variable, y = value, fill = group))+ geom_boxplot()
p
p <- p +stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p
p <- p + theme_set(theme_set(theme_bw(base_size=20)))
p
p <- p + theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### violinplot
p <- ggplot(data = exprSet_L,aes(x = variable, y = value, fill = group))+geom_violin()+
  stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")+
  theme_set(theme_set(theme_bw(base_size=20)))+
  theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### histogram (Here we can see the negative binomial distribution of read counts of all genes)
p <- ggplot(exprSet_L, aes(value, fill = group))+geom_histogram(bins = 200)+facet_wrap(~variable, nrow = 4)
p
#### density plot
p <- ggplot(exprSet_L, aes(value, fill = group, col = group))+geom_density()+facet_wrap(~variable, nrow = 4)
p
p <- ggplot(exprSet_L, aes(value, col = group))+geom_density()
p
### Make 'intergenic_intronic_transcript' object with transcripts only in 
### intergenic and intronic regions
intergenic_intronic_transcript <- transcripts[row.names(transcripts) %in% intergenic_intronic,]
dim(intergenic_intronic_transcript)


## Load in count table based on NumRead counts and DESeq2 normalization
### load packages
library(DESeq2)
reads <- read.table('merge_expr.txt', header = T, row.names = 1)
colnames(reads) <- paste("DCD002", c(492:511),"SQ", sep = "")
reads <- round(reads,0) # round reads
### DESeq2 package and perform normalization
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = reads, colData = sample, design  = ~ celltype)
#### Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
#### Log transformation
rld <- rlog(dds, blind = F) # This step takes some time
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]
### After Size factor normalization and log transformation, we can check the expression matrix again and compare the 
### sequencing depth among different samples
exprSet_L <- row.names(rlogMat)    
exprSet_L <- cbind(exprSet_L, as.data.frame(rlogMat)) # Must convert normalized_count into data.frame before cbind
exprSet_L <- melt(data = exprSet_L, id.vars = 'exprSet_L')
exprSet_L$group <- rep(group_List, each = nrow(rlogMat))
summary(rlogMat); dim(rlogMat)
#### boxplot
p <- ggplot(data = exprSet_L, aes(x = variable, y = value, fill = group))+ geom_boxplot()
print(p)
p <- p +stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p
p <- p + theme_set(theme_set(theme_bw(base_size=20)))
p
p <- p + theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### violinplot
p <- ggplot(data = exprSet_L,aes(x = variable, y = value, fill = group))+geom_violin()+
  stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")+
  theme_set(theme_set(theme_bw(base_size=20)))+
  theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### histogram (Here we can see the negative binomial distribution of read counts of all genes)
p <- ggplot(exprSet_L, aes(value, fill = group))+geom_histogram(bins = 200)+facet_wrap(~variable, nrow = 4)
p
#### density plot
p <- ggplot(exprSet_L, aes(value, fill = group, col = group))+geom_density()+facet_wrap(~variable, nrow = 4)
p
p <- ggplot(exprSet_L, aes(value, col = group))+geom_density()
p


##################### Using scater for object creation and normalization #################
# Load packages
library('SingleCellExperiment')
library('scater')
library('SC3')
options(stringsAsFactors = F)
data <- SingleCellExperiment(assays = list(counts = as.matrix(transcripts)), colData = sample)
data <- calculateQCMetrics(data)
plotQC(data, type = 'highest-expression', n = 20)
data <- normaliseExprs(data, method = "RLE", return_log = T, return_norm_as_exprs = T)


#################### Biological Analysis -- Feature selection ##########################
library('M3Drop')
library('matrixStats')
plot(rowMeans(counts(data)), rowVars(counts(data)), xlab = 'Mean Expression', ylab = 'Variance')
plot(rowMeans(counts(data)), rowVars(counts(data)), log = 'xy', xlab = 'Mean Expression', ylab = 'Variance')
Brennecke_HVG <- BrenneckeGetVariableGenes(logcounts(data), fdr = 0.01, minBiolDisp = 0.5)
length(Brennecke_HVG)



##################### Make DESeq2 object and perform normalization ##########################
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, colData = sample, design  = ~ celltype)

# Calculation and normalization, ordering according to mad
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)   # normalized_counts would be used in the downstream analysis, including 'WGCNA'
# head(normalized_counts)
# normalized_counts_mad <- apply(normalized_counts, 1, var)
# normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]   
# head(normalized_counts)

# Log transformation
rld <- rlog(dds, blind = F)
rlogMat <- assay(rld)
# rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]




#################### Hierarchincal Clustering and Principal Component Analysis (PCA) ######################
# Selecting colors
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# Calculating the pearson correlation matrix and perform the sample hierarchical clustering
pearson_cor <- as.matrix(cor(rlogMat, method = 'pearson'))
head(pearson_cor)
hc <- hcluster(t(rlogMat), method="pearson")
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')

# Principal Component Analysis
## How to examine the loading scores to determine what variables have the larest effect on the graph
## Prcomp will give us three things: x, sdev, rotation
## x contains the prcincipal component for drawing a graph
## The prcomp() function calls the loading scores rotation, there are loading scores for each PC
head(rlogMat)
pca <- prcomp(t(rlogMat), scale. = T)   # remember to transpose the matrix
plot(pca$x[,1],pca$x[,2])

## loadings for PC1
loading_score1 <- pca$rotation[,1]
head(loading_score1)
gene_score1 <- abs(loading_score1)
gene_score_ranked1 <- sort(gene_score1, decreasing = T)
top_25_genes1 <- names(gene_score_ranked1[1:25])
top_25_genes1
pca$rotation[top_25_genes1,1]
## loadings for PC2
loading_score2 <- pca$rotation[,2]
head(loading_score2)
gene_score2 <- abs(loading_score2)
gene_score_ranked2 <- sort(gene_score2, decreasing = T)
top_25_genes2 <- names(gene_score_ranked2[1:25])
top_25_genes2
pca$rotation[top_25_genes2,1]
## Using plotPCA and ggbiplot to make PCA plot 
pca_data <- plotPCA(rld, intgroup = c('celltype'), returnData = T, ntop = 5000)    # loading for each principal component (compared with novel analysis with the same important loadings with novel transcripts)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=celltype))
p + geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
plotPCA(rld, intgroup=c("celltype"))
ggbiplot(pca ,choices = 1:2, obs.scale = T,labels = NULL, var.scale = T,groups = sample$celltype, ellipse = T, circle = T, var.axes = F, alpha = 0.6)+theme(legend.direction = 'horizontal', legend.position = 'top')+theme_bw()
ggcorrplot(pearson_cor, hc.order = T,lab = T)    




###################### Differential Gene Expression (DGE) ########################
sampleA <- 'beta'
sampleB <- 'acinal'
sampleC <- 'alpha'
sampleD <- 'delta'

####### Differential gene expression between beta cell and acinal cell #####
contrastV <- c('celltype', sampleA, sampleB)
res <- results(dds, contrast = contrastV) # extract results (beta vs acinal) from dds object
head(res)
##extract gene expression matrix of each type of cell from counts object
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
head(baseA)
baseMeanA <- as.data.frame(rowMeans(baseA))
colnames(baseMeanA) <- sampleA
head(baseMeanA)

baseB <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleB]
baseMeanB <- as.data.frame(rowMeans(baseB))
colnames(baseMeanB) <- sampleB
head(baseMeanB)

## make up the new res object containing average expression level of each type of cell
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1  # padj的缺失值赋值为1
res <- res[order(res$padj),] # 根据padj从小到大进行排序，排在上方的统计学越显著
head(res)
res$significance <- (res$padj<0.05) # 增加统计学是否有显著性一列

######## MA plot ########
plotMA(res[,c(4,5,10)], main = "MAplot")
com <- paste(sampleA,'_vs_', sampleB)
file_base <- paste('DESeq2',com,sep = '.')
write.table(as.data.frame(res),file = file_base,quote = F, row.names = F)
## Extract differential expressed gene (whole, upregulated gene, downregulated gene)
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleB,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)

######## Volcano plot ########
logFC <- res$log2FoldChange
FDR <- res$padj
head(res)
plot(logFC,log10(FDR)*(-1),col = ifelse(FDR<=0.01, ifelse(abs(logFC)>=1, 'red', 'black'), 'black'), xlab = 'logFC', ylab = '-1*log10(FDR)', main = "Volcano plot", pch = '.', ylim=c(1,100))
ggplot(data = res, aes(x=log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point()+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = -log10(0.05),lty=4, lwd=0.6,alpha=0.6)+
  geom_vline(xintercept = c(-1,1),lty=4, lwd=0.6,alpha=0.6)+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="Volcanoplot: beta_vs_acinal", x= 'log2(fold change)', y = '-log10(padj)')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(data=subset(res, -log10(padj) > 100), aes(label=ID), col= 'black', alpha = 0.8)

######### Heatmap ########
# Heatmap of beta cell and acinal cell
res_de_up_top50_id <- as.vector(head(res_de_up$ID,50))
res_de_down_top50_id <- as.vector(head(res_de_down$ID,50))
res_de_top100 <- c(res_de_up_top50_id,res_de_down_top50_id)
res_de_top100_expr <- normalized_counts[rownames(normalized_counts) %in% res_de_top100,c(2:6,1,16:20)]
head(res_de_top100_expr)
pheatmap(res_de_top100_expr,cluster_rows = T, scale = 'row',annotation_col = sample)


######## biomaRt ID conversion ###########
# This step is important since the enrichment analysis needs/ recognizes the entrezid rather than ensembl id
# Use biomaRt to concert gene ID (for enrichment analysis, entrezid should be used)
head(listMarts())
mart <- useMart('ensembl')
View(listDatasets(mart))  # select dataset "drerio_gene_ensembl"
ensembl <- useDataset('drerio_gene_ensembl', mart)


######## Enrichment analysis (GO, KEGG) using beta_higherthan_acinal data #########
## beta higher than acinal
entrezid_beta_higherthan_acinal<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(res_de_up),mart = ensembl)
head(entrezid_beta_higherthan_acinal)
dim(entrezid_beta_higherthan_acinal)
entrezid_beta_higherthan_acinal <- entrezid_beta_higherthan_acinal[is.na(entrezid_beta_higherthan_acinal$entrezgene) != T,]
summary(entrezid_beta_higherthan_acinal)
dim(res_de)
head(res_de)
geneList_foldchange_beta_higherthan_acinal <- res_de[res_de$ID %in% entrezid_beta_higherthan_acinal$ensembl_gene_id,c('ID','log2FoldChange')]
head(geneList_foldchange_beta_higherthan_acinal)
entrezid_beta_higherthan_acinal <- unique(entrezid_beta_higherthan_acinal$entrezgene)
### Biological Process
BP <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
head(result_BP)
?dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,120)
barplot(result_BP,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
barplot(result_BP,x = 'count', showCategory = 10)
plotGOgraph(result_BP)
### Molecular Function
MF <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
head(result_MF)
plotGOgraph(result_MF)
### Cellular Component
CC <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.045)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
head(result_CC)
plotGOgraph(result_CC)
### KEGG pathway
kk <- enrichKEGG(entrezid_beta_higherthan_acinal,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
head(result_kk)
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)
?cnetplot
result_kk$Description
head(geneList_foldchange_beta_higherthan_acinal)
# To get the pathway image: pv.out <- pathview(gene.data = geneList_foldchange_beta_higherthan_acinal$log2FoldChange,pathway.id = result_kk$ID, species = 'dre', kegg.native = T)



############ Enrichment analysis (GSEA) using differential expressed gene between beta-cell and acinal cell (upregulated and downregulated) #################
## Preparation of geneList for GSEA analysis
entrezid_DGE_beta_with_acinal<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(res_de),mart = ensembl)
head(entrezid_DGE_beta_with_acinal)
dim(entrezid_DGE_beta_with_acinal)
dim(res_de)
entrezid_DGE_beta_with_acinal <- entrezid_DGE_beta_with_acinal[is.na(entrezid_DGE_beta_with_acinal$entrezgene) != T,]
summary(entrezid_DGE_beta_with_acinal)
head(entrezid_DGE_beta_with_acinal)
geneList_foldchange_beta_with_acinal<- res_de[res_de$ID %in% entrezid_DGE_beta_with_acinal$ensembl_gene_id,c('ID','log2FoldChange')]
head(geneList_foldchange_beta_with_acinal)
dim(geneList_foldchange_beta_with_acinal)
summary(entrezid_DGE_beta_with_acinal)
summary(geneList_foldchange_beta_with_acinal)
names(geneList_foldchange_beta_with_acinal)= c("ensembl_gene_id", "log2FoldChange")
geneList_foldchange_beta_with_acinal <- merge(entrezid_DGE_beta_with_acinal,geneList_foldchange_beta_with_acinal,by = "ensembl_gene_id")
summary(geneList_foldchange_beta_with_acinal)
head(geneList_foldchange_beta_with_acinal)
geneList_foldchange_beta_with_acinal <- geneList_foldchange_beta_with_acinal[,-1]
geneList_foldchange_beta_with_acinal_2 <- geneList_foldchange_beta_with_acinal[,2]
names(geneList_foldchange_beta_with_acinal_2)<- geneList_foldchange_beta_with_acinal[,1]
geneList_foldchange_beta_with_acinal_2 <- sort(geneList_foldchange_beta_with_acinal_2,decreasing = T)
head(geneList_foldchange_beta_with_acinal_2)
## gsea analysis of GO--biological process
gsea_analysis_GO_BP <- gseGO(geneList_foldchange_beta_with_acinal_2,ont = "BP", OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_BP)
dim(gsea_analysis_GO_BP)
gseaplot(gsea_analysis_GO_BP,geneSetID = "GO:0006508")
dotplot(gsea_analysis_GO_BP, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_BP, showCategory = 30,fill = 'pvalue')
## gsea analysis of GO--molecular function
gsea_analysis_GO_MF <- gseGO(geneList_foldchange_beta_with_acinal_2, ont = 'MF', OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_MF)
dim(gsea_analysis_GO_MF)
gseaplot(gsea_analysis_GO_MF,geneSetID = 'GO:0046914')
dotplot(gsea_analysis_GO_MF, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_MF, showCategory = 30, fill = 'pvalue')
## gsea analysis of GO--cellular component
gsea_analysis_GO_CC <- gseGO(geneList_foldchange_beta_with_acinal_2, ont = 'CC', OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_CC)
gseaplot(gsea_analysis_GO_CC, geneSetID = 'GO:0005576')
dotplot(gsea_analysis_GO_CC, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_CC, showCategory = 30, fill = 'pvalue')
## gsea analysis of KEGG
gsea_analysis_KEGG <- gseKEGG(geneList_foldchange_beta_with_acinal_2,organism = 'dre')
dotplot(gsea_analysis_KEGG)
ridgeplot(gsea_analysis_KEGG)
pv.out <- pathview(gene.data = geneList_foldchange_beta_with_acinal_2,pathway.id = result_kk$ID, species = 'dre', kegg.native = T)



#======================================================================
#
#                                  WGCNA 
#
#======================================================================

############# STEP 1: Expression data and Phenotype data Preparation ##############
## Load package and basic setting
library(WGCNA)
options(stringsAsFactors = F)



####### read in the expression matrix ############
# counts/FPKM/RPKM/normalized counts are all can be used in WGCNA, the batch information is better to be taken into account
# head(normalized_counts)
# dim(normalized_counts)
# rownames(normalized_counts)
# colnames(normalized_counts)
# datExpr0 <- as.data.frame(t(normalized_counts))
# Instead of using NumReads from Salmon, it is better to use TPM 
# All you need to change is the expression matrix
# data1  <- read.table('merge_tpm.txt', header = T, row.names = 1)
# data1 <- normalized_counts
data1 <- intergenic_intronic_transcript   # For TPM
datExpr0 <- as.data.frame(t(data1)); dim(datExpr0)
# After transpose, have a look at your data
# dim(datExpr0)    
# head(names(datExpr0))
# head(rownames(datExpr0))

####### check missing value and filter #########
## check missing value ##
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ',')))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Remove samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ',')))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## filter using threshold ##
### normalized counts such as RPKM and log transformed counts are the input
### You can also use 'apply' function to select most variable genes (5000 genes or so) for downstream analysis
### Pick a threshold, normally genes for WGCNA should not be too large (< 20000)
n <- nrow(datExpr0); n
datExpr0[n+1,] <- apply(datExpr0[c(1:n),],2, function(x){log2(mean(x)+1)}); dim(datExpr0)
datExpr0[n+2,] <- apply(datExpr0[c(1:n),], 2, function(x){log2(sd(x)/mean(x)+1)}) ; dim(datExpr0)# 使用变异系数（coefficient of variance, CV)较正表达量高低对变异度的影响
datExpr1 <- as.data.frame(t(datExpr0))
names(datExpr1)[21] <- 'log2_mean';names(datExpr1)[22] <- 'log2_CV'; names(datExpr1)
head(datExpr1)[21]; head(datExpr1[22]); colnames(datExpr1)

### Use loess model to fit the curve in order to select the highly variable genes
p <- ggplot(datExpr1, aes(x = log2_mean, y = log2_CV))+ geom_point() + 
  geom_smooth(span = 0.2, method = 'loess', na.rm = T) + 
  geom_smooth(method = lm, col = 'red', na.rm = T) + 
  ylim(c(0.4,2.6)) +
  geom_vline(xintercept = seq(0,3,0.2), col = 'darkgreen', lty = 2) +
  theme_classic(); p
model_xlog2mean_ylog2CV <- loess(datExpr1$log2_CV ~ datExpr1$log2_mean, span = 0.2, method = 'loess')
summary(model_xlog2mean_ylog2CV)
prediction <- predict(object = model_xlog2mean_ylog2CV, data.frame(datExpr1$log2_mean), se = T)
head(prediction$fit)
datExpr0 <- datExpr1[datExpr1$log2_CV > prediction$fit & datExpr1$log2_mean > 2.5,1:20]; dim(datExpr0)






summary(as.vector(as.matrix(datExpr0[21,])))
hist(as.vector(as.matrix(datExpr0[21,],breaks = 1000)))
summary(as.vector(as.matrix(datExpr0)[22,]))
hist(as.vector(as.matrix(datExpr0)[22,]), breaks = 1000, main = "Histogram of CV (coefficient of Varince)")
datExpr0 <- datExpr0[1:n, datExpr0[n+2,]>1.5 & datExpr0[n+1,] >= 1 ] # 增加对表达量的限定，以减少噪音
dim(datExpr0)
filtered_normalized_counts <- datExpr0
head(filtered_normalized_counts)
dim(filtered_normalized_counts)

###### Sample Cluster (样本聚类) ######
## For sample clustering, we should used transformed count matrix. 
## datExpr0 has been transposed in the steps before, therefore rows are samples, columns are genes
## We use this transposed count matrix as input to do sample clustering
sampleTree <- hclust(dist(t(datExpr0)), method = 'average')
par(mfrow = c(1,1))
plot(sampleTree, main = "Sample clustering to detect outlier")

###### load trait data ######
traitData <- read.csv(file = 'trait_D.csv', header = T, row.names = 1, check.names = F)
head(traitData)
## convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(traitData, signed = F); traitColors;names(traitData)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(traitData), main = 'Sample dendrogram and trait heatmap')
###### save important files ########
save(datExpr0, file = 'normalized_forAnalysis.RData')
save(traitData, file = 'trait_forAnalysis.RData')

################ STEP 2: Network Construction #############
######## Select the best soft-thresholding power #########
## Choose a set of soft-thresholding powers
powers <- c(1:30)
## Call the Network Topological Analysis function
sft <- pickSoftThreshold(t(datExpr0), powerVector = powers, verbose = 5)  # This step takes some time.
par(mfrow = c(1,2), mar = c(6.5,8,3,3))
cex1 = 0.9
str(sft)  
# The output of 'sft' is a list object. powerEstimate is the estimated best power
# The second output of 'sft' is fitIndices, which is a matrix. The fifth column, 'mean.k' denote average connectivity.

## Scale-free topological fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]*sft$fitIndices[,2]), xlab = 'Soft Threshold (power)', 
     ylab = 'Scale free Topological Model Fit, signed R^2', type = 'n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]*sft$fitIndices[,2]), labels = powers, cex = cex1, col = 'red')
abline(h = 0.9, col = 'blue') 
# The blue line corresponds to using a R^2 cut-off of h
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n', main = paste('Mean Connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels = powers, cex = cex1, col = 'red')
## Choose the softPower
softPower <- sft$powerEstimate
softPower
Adjacency <- adjacency(t(datExpr0), power = softPower)
## Turn adjacency matrix into Topological matrix, this step takes some time ##
TOM <- TOMsimilarity(Adjacency)  
dissTOM <- 1-TOM
## Call the hierarchincal clustering function
geneTree <- hclust(as.dist(dissTOM), method = 'average')   # This step takes some time, to calculate the distance in gene pairs.
par(mfrow = c(1,1))
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity', labels = F, hang = 0.04)
# We like large modules, so we set the minimum module size relatively high
minModuleSize <- 8
# set cutHeight, otherwise the function will set a cutHeight automatically
# The next step may take some time
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = F, minClusterSize = minModuleSize)
table(dynamicMods)
## Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors,'Dynamic Tree Cut', dendroLabels = F, 
                    hang = 0.03, addGuide = T, guideHang = 0.05, main = 'Gene dendrogram and module colors')
## Calculate eigengene
MEList <- moduleEigengenes(t(datExpr0), dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')
par(mfrow = c(1,1))
plot(METree, main = 'Clustering of module eigengene', xlab = '', sub = '')
MEDissThres = 0.25 # set the threshold to make some branches together
abline(h = MEDissThres, col = 'red')
Merge <- mergeCloseModules(t(datExpr0), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- Merge$colors
table(mergedColors)
mergedMEs <- Merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'), 
                    dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05) # It takes time!!!
## Rename to module Colors
moduleColors <- mergedColors
colorOrder <- c('grey', standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
moduleLabels
MEs <- mergedMEs
# The next save step can take some time!
save(MEs, file = 'MEs_networkConstruction-stepBystep.RData')
save(TOM, file = 'TOM_networkConstruction-stepBystep.RData')
save(dissTOM, file = 'dissTOM_networkConstruction-stepBystep.RData')
save(moduleLabels, file = 'moduleLabels_networkConstruction-stepBystep.RData')
save(geneTree, file = 'geneTree_networkConstruction-stepBystep.RData')
save(sft, file = 'sft_networkConstruction-stepBystep.RData')
save(MEs, TOM, dissTOM, moduleLabels, geneTree, sft, file = 'networkConstruction-stepBystep.RData')


######################### Relate modules to external clinical traits ########################
nGene <- ncol(t(datExpr0))
nSample <- nrow(t(datExpr0))
moduleTraitCor <- cor(MEs, traitData, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSample)
textMatrix = paste(signif(moduleTraitCor,2),"\n(", signif(moduleTraitPvalue,1),")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mfrow = c(1,1))
par(mar = c(6,8.5,3,3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim <- c(-1,1),
               main = paste('Module-trait relationships'))
### Here we found that MEcyan and MEbrown modules are highly associated with beta-cell (celltype/cell function...)


##################### Visualizing the gene network #######################
nSelect <- 153
set.seed(10)
select <- sample(nGene, size = nSelect)
selectTOM <- dissTOM[select, select]
selectTree <- hclust(as.dist(selectTOM), method = 'average')
selectColors <- moduleColors[select]
plotDiss <- selectTOM^8
diag(plotDiss) <- NA
TOMplot(plotDiss, selectTree, selectColors, main = 'Network heatmap plot, select genes')

###################### Visualizing the gene network of eigengene ###################
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8)

##################### Module membership (MM) and Gene significance ######################
# Gene Significance, GS: 基因显著性参数，为非负数字信息，比如基于样本的临床信息(clinical traits)和基于每个基因的-log(p-value)等
# names (colors) of each module
modNames <- substring(names(MEs),3)
modNames
# 首先计算模块与基因的相关性矩阵
# MEs表示每个模块在样本里的值
geneModuleMembership <- as.data.frame(cor(t(datExpr0), MEs, use = 'p'))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSample))
head(MMPvalue)
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
# names of those traits
traitNames <- names(traitData)
geneTraitSignificance <- as.data.frame(cor(t(datExpr0), traitData, use = 'p'))
head(geneTraitSignificance)
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSample))
head(GSPvalue)
names(geneTraitSignificance) <- paste('GS.',traitNames,sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")

## Plot MM vs GS for each trait vs each module
for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){####进行这部分计算必须每个模块内基因数量大于2，由于前面设置了最小数量是30，这里可以不做这个判断，但是grey有可能会出现1个gene,它会导致代码运行的时候中断，故设置这一步
      
      #sizeGrWindow(7, 7)
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}


########### 提取指定模块的基因名 ###############
# Select module
module = 'black'
# Select module probes
probes <- colnames(datExpr0)
inModule <- probes[moduleColors == module]
modProbes_black <- probes[inModule]; length(modProbes_black)
modProbes_pink <- probes[inModule]; length(modProbes_pink)
modProbes_black

# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = ‘black’ 的模块进行分析
# 'black' module gene
module <- 'black'
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
black_module_index <- which(moduleColors == 'black') # the index of genes which belong to 'brown' module
length(colnames(datExpr0)[black_module_index])
length(rownames(filtered_normalized_counts)[black_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
black_module_transcriptName <- rownames(filtered_normalized_counts)[black_module_index]
black_module_transcriptName
# 'brown' module 有1095个基因


# 'cyan' module gene
module <- 'cyan'
column <- match(module, modNames)
moduleGenes <- moduleColors == module
cyan_module_index <- which(moduleColors == 'cyan')
cyan_module_transcriptName <- rownames(filtered_normalized_counts)[cyan_module_index]
cyan_module_transcriptName




#####
names(datExpr0)
probes = names(datExpr0)

################# export GS and MM ############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = "10_GS_and_MM.xls",sep="\t",row.names=F)


########################### Exporting to Cytoscape all one by one ##########################
# Select each module

for (mod in 1:nrow(table(moduleColors)))
{
  
  modules = names(table(moduleColors))[mod]
  # Select module probes
  probes = names(datExpr0)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-TPM-edges-", modules , ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-TPM-nodes-", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}



################# Select MEcyan and MEbrown modules for network analysis #############
############# check brown module #################
## load edges file 
beta_brown_module <- beta_brown_module[order(beta_brown_module$weight,decreasing = T),]
head(beta_brown_module)
summary(beta_brown_module)
## sort edges and identify which genes/transcripts have the highest number of edges
sort(table(beta_brown_module$fromNode), decreasing = T)
## select edges based on weight
beta_brown_module <- beta_brown_module[beta_brown_module$weight>0.1,]
dim(beta_brown_module)
## write tables with edges weight higher than 0.1 and usd this file as input of cytoscape for cytoHubba analysis
write.table(x = beta_brown_module,file = 'beta_brown_module.txt',row.names = T, sep = '\t')
## load csv file after cytoHubba analysis
cytoHubba <- read.csv('CytoHubba-brown-module-30genes-summary.csv', header = T, row.names = NULL)
head(cytoHubba)
cytoHubba <- cytoHubba[,-1]
cytoHubba_vector <- c(cytoHubba$MCC,cytoHubba$DNMC, cytoHubba$MNC, cytoHubba$Degree, cytoHubba$EPC, cytoHubba$BottleNeck, cytoHubba$EcCentricity, cytoHubba$Closeness, cytoHubba$Radiality, cytoHubba$Betweenness, cytoHubba$Stress, cytoHubba$ClusteringCoefficient)
## check which transcripts have the higher rate of hub genes predicted by 12 algorithms from cytoHubba
sort(table(cytoHubba_vector), decreasing = T)
## exist in >= 8 algorithms and check their coding potential according to CPAT and PLEK
## CPAT: TU92920, TU40817, TU45134, TU54701, TU18235, TU3649, TU40817, TU4167, TU48082, TU52640, 
##       TU70082, TU80891, TU81363, TU86851, TU88207, TU90967
## PLEK: TU92920, TU79222, TU45134, TU54701, TU7033, TU76925, TU11175, TU17541, TU18235, 
##       TU23207, TU25639, TU32167, TU3649, TU40511, TU43386, TU43858, TU45801, TU47952, 
##       TU48082, TU48752, TU49427, TU52640, TU54874, TU54799, TU64562, TU65179, TU66710, 
##       TU67121, TU68936, TU7033, TU76156, TU76925, TU77612, TU79222, TU80891, TU82490,
##       TU86851, TU88207, TU90967, TU91260, TU91902, TU91903, TU92981, TU93649, TU93920,
##       TU96550, TU96700
## Overlap: TU92920(10,kalrn, RhoGEF Kinase, Human & Mouse; full length sequence in zebrafish, kalrn; cytoplasm, nucleus; high expression in beta-cell; Sandberg low expression of kalrn; high peaks in ATAC/DNase-seq), 
##          TU45134(8,uncharacterized LOC, SNTG1, UGT8,Taf1, Human 2| Mouse 1; Full length DNA in zebrafish; cytoplasm, nucleus; high expression in beta-cell), 
##          TU54701(8,ELOVL4, JUN, uncharacterized LOC, GLS, ZNF621, HTN3, Csrp3, Cntn4, Serp2, Msantd3, Human | Mouse; ppdpfb in zebrafish; cytoplasm, nucleus; high expression in beta-cell)
##          TU80891(4, FAM, uncharacterized LOC, ALS2CR11, Zcchc4, Zfp704, uncharacterized LOC, Human | Mouse; zebrafish full length DNA, ptprnb transcript; cytoplasm, nucleus; high expression in beta-cell),
##          TU52640(3,GNA11, uncharacterized LOC, ZMYND8, Maml3, Mylk, Human | Mouse; zebrafish full length, gna11a; cytoplasm, nucleus; high expression in beta-cell), 
##          TU88207(2,---; cytoplasm, nucleus),
##          TU3649(1,---), TU48082(1, GABRB3, LITAF, Zdbf2, GM31214, Human | Mouse; cytoplasm, nucleus), 
##          TU90967(1, C9orf3, LOC; cytoplasm, nucleus),TU86851(1,ltgb6); cytoplasm, nucleus,  TU18235(1,ATOX1, MSI2, ACBD5, PHF14, HOOK3, Human & Mouse ATOX1, Human | Mouse; cytoplasm, nucleus),
cytoHubba_vector_brown_unique <- names(table(cytoHubba_vector))
write.table(x = cytoHubba_vector_unique, file = '/Users/mijiarui/RNA-seq/lncRNA_prediction/PLEK_prediction/cytoHubba_cyan_module_geneName.txt', row.names = F)


############ check cyan module ##############
beta_cyan_module <- read.table('CytoscapeInput-edges-cyan.txt', header = T, stringsAsFactors = F)
head(beta_cyan_module)
beta_cyan_module_attribute <- read.table('CytoscapeInput-nodes-cyan.txt', header = T, stringsAsFactors = F)
head(beta_cyan_module_attribute)
beta_cyan_module <- beta_cyan_module[order(beta_cyan_module$weight, decreasing = T),]
summary(beta_cyan_module)
beta_cyan_module <- beta_cyan_module[beta_cyan_module$weight>0.1,]
write.table(x = beta_cyan_module, file = "beta_cyan_module.txt", row.names = T, sep = '\t')
dim(beta_cyan_module)
cytoHubba <- read.csv('CytoHubba-cyan-module-30genes-summary.csv', header = T, row.names = NULL)
head(cytoHubba)
cytoHubba <- cytoHubba[,-1]
cytoHubba_vector <- c(cytoHubba$MCC,cytoHubba$DNMC, cytoHubba$MNC, cytoHubba$Degree, cytoHubba$EPC, cytoHubba$BottleNeck, cytoHubba$EcCentricity, cytoHubba$Closeness, cytoHubba$Radiality, cytoHubba$Betweenness, cytoHubba$Stress, cytoHubba$ClusteringCoefficient)
## check which transcripts have the higher rate of hub genes predicted by 12 algorithms from cytoHubba
sort(table(cytoHubba_vector), decreasing = T)
cytoHubba_vector_cyan_unique <- names(table(cytoHubba_vector))
cytoHubba_vector_cyan_unique
length(cytoHubba_vector_cyan_unique)
write.table(x = cytoHubba_vector_cyan_unique,file = '/Users/mijiarui/RNA-seq/lncRNA_prediction/CPAT_prediction/cytoHubba_brown_module_geneName.txt', row.names = F)
## exist in >= 8 algorithms and check their coding potential according to CPAT and PLEK
## CPAT: TU38615, TU11166, TU13810, TU14524, TU17004, TU17140, TU17617, TU24625, TU29445, TU31086, TU31813,
##       TU38615,TU43443, TU45170, TU45271, TU4839, TU51616, TU5187, TU53164, TU547, TU54953, TU5560,
##       TU64832, TU67592, TU68230, TU77233, TU78908, TU80024, TU80615, TU89114, TU8929, TU8930, TU97447, TU98365
##       
## PLEK: TU38615, TU83796, TU14524, TU16838, TU17004, TU17140, TU24625, TU2684, TU31086, TU35880, 
##       TU54264, TU55376, TU57001, TU60147, TU76597, TU77233, TU83796, TU89114, TU8930, TU96138
##       TU96830, TU97447, TU98365

## Overlap: TU38615 (8,LOC, Human | Mouse; zebrafish: gpt2l; nucleus; high expression in beta-cell), 
##          TU31086(5,----, Human | Mouse; Zebrafish --; cytoplasm, nucleus), 
##          TU17004(4,RAP1GAP, Human & Mouse; Zebrafish full length; cytoplasm, nucleus; high expression in beta-cell), 
##          TU17140(1,IKZF3, TNFRSF19, ACYP2, Pdkcc, Gm38907, Pot1b, Human | Mouse; Zebrafish full length; nucleus), 
##          TU24625(1,RAP1GAP, Human & Mouse; zebrafish full length; cytoplasm, nucleus, high expression in beta-cell),  TU77233(1,---; cytoplasm, nucleus), TU14524(1,---; cytoplasm, nucleus),
##          TU8930(1, SLC24A2, KCMF1,ZAK, ZFp160, Plekha6, Zkscan7, Human | Mouse; cytoplasm, nucleus), TU97447(1,LOC, Mouse; cytoplasm, nucleus), 
##          TU98365(1,PSMD2, NBPF7, IFNE, FAM111B, COL3A1, spata1, supt6, Human & Mouse (psmd2) ; cytoplasm, nucleus)



#====================================================================
#
#                     kalrn flanking sequence 
#
#====================================================================
#####   TU38615
mart <- useMart('ensembl') 
ensembl <- useDataset('drerio_gene_ensembl', mart)
kalrn_zebrafish_flanking <- getBM(attributes = 'entrezgene', filters = c('chromosome_name', 'start', 'end'), values = list(2, 3602294,3838363), mart = ensembl)
kalrn_zebrafish_flanking
kalrn_zebrafish_flanking$gene <- c('hat1', 'dlx1a', 'dlx2a', 'itga6a', 'crbrd1', 'dcaf17', 'mettl8','tlk1a', 'gad1a', 'sp5a', 'myo3b', 'kalrna', 'kcnj3a', 'galnt13', 'rprma', 'prpf40a', 'fmnl2a', 'nr4a3a', 'cytip')

ensembl <- useDataset('mmusculus_gene_ensembl', mart)
kalrn_mouse_flanking <- getBM(attributes = 'entrezgene', filters = c('chromosome_name', 'start', 'end'), values = list(16, 32969073,35514924), mart = ensembl)
kalrn_mouse_flanking
kalrn_mouse_flanking$gene <- c('lrch3', 'Iqcg', 'rpl35a', 'lmln', 'osbpl11', 'snx4', 'zpf148', 'slc12a8', 'heg1', 'muc13', 'itgb5', 'umps', 'kalrn', 'ropn1',
                               'ccdc14', 'mylk', 'hacd2', 'adcy5', 'sec22a', 'pdia5')
ensembl <- useDataset('hsapiens_gene_ensembl', mart)
kalrn_human_flanking <- getBM(attributes = 'entrezgene', filters = c('chromosome_name', 'start', 'end'), values = list(3, 123388215  ,125395949), mart = ensembl)
kalrn_human_flanking$gene <- c('adcy5', 'hacd2', 'mylk-as1', 'mylk', 'mylk-as2', 'ccdc14', 'ropn1', 'kalrn', 'mir5002', 'mir6083', 'umps', 'mir544b',
                               'itgb5', 'muc13', 'heg1','slc12a8','mir5092', 'znf148')












cytohubba_ivory_30genes <- read.csv('cytohubba_ivory_30genes.csv', header = T)
sort(table(as.vector(as.matrix(cytohubba_ivory_30genes))), decreasing = T)
# CPAT: TU89152(chr8:8,038,215-8,074,578), 
#       TU2818(chr1:47,776,158-47,795,909), 
#       TU42457(chr2:59,127,073-59,167,027),
#       TU31064(chr17:52,686,801-52,690,132, high expression, meis2a near), 
#       TU68230(chr3:58,237,791-58,243,202), 
#       TU7260(chr11:1,289,009-1,291,053), 
#       TU14147(chr13:7,763,511-7,763,901), 
#       TU31086(chr17:52,690,277-52,692,381, high expression, meis2a near), 
#       TU58500(chr24:21,520,391-21,520,978, intronic,lnx2 important for exocrine cell differentiation, near pdx1), 
#       TU92926(chr9:5,327,117-5,327,687), 
#       TU31063(chr17:52,686,807-52,690,132)
#       TU57170(chr23:46,222,963-46,228,837), 
#       TU63869(chr3:10,050,970-10,065,216), 