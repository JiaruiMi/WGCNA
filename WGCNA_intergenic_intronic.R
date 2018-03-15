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
library(ggrepel)

##################### Read in data (SampleGroup) ########################
## Load in sample data
sample <- read.csv('SampleGroup.csv', header = T, row.names = 1, colClasses = 'factor')



##################### Read in data (TPM) and data visualization ########################
## Load in count table based on TPM normalization
transcripts <- read.table('merge_tpm.txt', header = T, row.names = 1)
intergenic_intronic <- read.table('transcript_intergenic_intronic_bedtools.txt')
intergenic_intronic <- as.vector(as.matrix(intergenic_intronic))
colnames(transcripts) <- paste("DCD002", c(492:511),"SQ", sep = "")
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

##################### Read in data (NumReads) and data visualization ########################
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



#==============================================================================================
#
#     Scater object creation, Gene filtering and Exploratory data analysis (SC3, M3Drop)
#                     Get the genelist of highly variable genes
#
#==============================================================================================

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

#### Using Scater Highly Variable Genes, we can do clustering and PCA plot
##### Clustering with normalized_counts derived from DESeq2
library('gplots')
library('pheatmap')
library('amap')
library('RColorBrewer')
normalized_counts_HVG <- normalized_counts[rownames(normalized_counts) %in% Brennecke_HVG,]
pearson_cor <- as.matrix(cor(normalized_counts_HVG, method = 'pearson'))
head(pearson_cor)
hc <- hcluster(t(normalized_counts_HVG), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)

##### Clusterig with Log-transformed normalized counts derived from DESeq2
rlogMat_HVG <- rlogMat[rownames(rlogMat) %in% Brennecke_HVG,]
pearson_cor <- as.matrix(cor(rlogMat_HVG, method = 'pearson'))
hc <- hcluster(t(rlogMat_HVG), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)



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



#======================================================================================
#
#                                          WGCNA 
#
#======================================================================================
## So, for the WGCNA, we need to first select proper transcripts for downstream analysis. For, we need to
## choose the transcripts that are mapping to the non-exonic regions (intergenic and intronic regions); 
## Second, we need to filter the transcripts with less than 2 exons and less than 200 nt; Actually, you have
## to prepare your data according to the bed file and select the gene list in adcance. Then you prepare the
## count table with the qualified transcripts listed and do the down stream analysis.


############# STEP 0: Prepare your count table (TPM) for downstream analysis ##############
## set working directory
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lincRNA')
## load packages for data visulizations
library(ggplot2)
library(ggrepel)

### Summaries
#### The number of transcripts in intergenic and intronic regions is 16176
#### The number of transcripts in intergenic and intronic regions with exons >= 2 is 5600

#==============================================================================
#    Novel transcripts with >=2 exons in intergenic and intronic regions
#==============================================================================
## reads.csv detects the number of columns according to the first five rows.
## To prepare for the dataframe, first you need to know the max number of exons for the transcripts
num_exons_2plus <- read.csv('number_exons_2plus.csv', header = F)
max(num_exons_2plus) # the maximum number of exons is 35
ggplot(data = num_exons_2plus, aes(x = V1)) + geom_histogram(bins = 35)
summary(num_exons_2plus) 
### the mean exon numbers for new intergenic/intronic transcripts with >= 2 exons is 3.337, the median is 2


### Filter out those with transcript length < 200
length_exons_2plus <- read.table('exon_2plus.csv', header = FALSE, sep = ",", 
                                 col.names = paste0("V",seq_len(35)), fill = TRUE)
head(length_exons_2plus)
table(rowSums(length_exons_2plus, na.rm = T) < 200)
which(rowSums(length_exons_2plus, na.rm = T) < 200)
#### In transcripts with exons >= 2 (5600), transcripts with length < 200 = 7, >= 200 is 5593

### Check the length (mean, median, min and max) of the transcripts with exon >=2 and length >200nt
summary(rowSums(length_exons_2plus[rowSums(length_exons_2plus, na.rm = T) >= 200,], na.rm = T)) # length summary
summary(num_exons_2plus[which(rowSums(length_exons_2plus, na.rm = T) >= 200),]) # exon number summary


#==============================================================================
#    Novel transcripts with in intergenic and intronic regions (exon = 1 or >= 2)
#==============================================================================

num_exons_inter_intro <- read.csv('number_exons_assembly_intergenic_intronic.csv', header = F)
ggplot(data = num_exons_inter_intro, aes(x = V1)) + geom_histogram(bins = 35)
summary(num_exons_inter_intro)
### the mean exon numbers for intergenic/intronic transcripts is 1.809, the median is 1


### Filter out those with transcripts length < 200
length_exons_assembly_intergenic_intronic <- read.table('exon_assembly_intergenic_intronic.csv', header = FALSE, sep = ",", 
                                                        col.names = paste0("V",seq_len(35)), fill = TRUE)
table(rowSums(length_exons_assembly_intergenic_intronic, na.rm = T) < 200)
index_length_lessThan_200 <- which(rowSums(length_exons_assembly_intergenic_intronic, na.rm = T) < 200)
index_length_lessThan_200
### Check the length (mean, median, min and max) of the transcripts length >200nt
summary(rowSums(length_exons_assembly_intergenic_intronic[rowSums(length_exons_assembly_intergenic_intronic, na.rm = T) >= 200,], na.rm = T)) # length summary
summary(num_exons_inter_intro [which(rowSums(length_exons_assembly_intergenic_intronic, na.rm = T) >= 200),]) # exon number summary



#=====================================================================================
#    Clean your expression matrix with transcripts with exon >= 2 and length > 200nt
#=====================================================================================
## set working directory
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_data/WGCNA')


##################### Read in data (SampleGroup) ########################
## Load in sample data
sample <- read.csv('SampleGroup.csv', header = T, row.names = 1, colClasses = 'factor')

##################### Read in data (TPM) and data visualization ########################
## Load in count table based on TPM normalization (of all the 98528 transcripts)
transcripts <- read.table('merge_tpm.txt', header = T, row.names = 1)
intergenic_intronic <- read.table('transcript_intergenic_intronic_bedtools.txt')
intergenic_intronic <- as.vector(as.matrix(intergenic_intronic))
colnames(transcripts) <- paste("DCD002", c(492:511),"SQ", sep = "")
head(transcripts); dim(transcripts)

## Choose the one in the intergenic and intronic regions
intergenic_intronic_transcript <- transcripts[row.names(transcripts) %in% intergenic_intronic,]
dim(intergenic_intronic_transcript)

## Get rid of those with length less than 200
intergenic_intronic_transcript_moreThan_200 <- intergenic_intronic_transcript[-index_length_lessThan_200,]
head(intergenic_intronic_transcript_moreThan_200)

## Get rid of those with exons less than 2
exon_1 <- read.csv('exon_1_intergenic_intronic.csv', header = F)
intergenic_intronic_exon_highThan1_length_highThan_200 <- intergenic_intronic_transcript_moreThan_200[!row.names(intergenic_intronic_transcript_moreThan_200)  
                                                                                                      %in% exon_1$V1,]
dim(intergenic_intronic_exon_highThan1_length_highThan_200)


### Here we got the final expression matrix of TPM of intergenic and intronic regions of transcripts with >= 2 exons and length higher than 200bp
### the expression matrix is store in the object of 'intergenic_intronic_exon_highThan1_length_highThan_200'


############# STEP 1: Expression data and Phenotype data Preparation ##############
## Load package and basic setting
library(WGCNA)
options(stringsAsFactors = F)



####### read in the expression matrix ############
# counts/FPKM/RPKM/normalized counts are all can be used in WGCNA, the batch information is better to be taken into account
# Instead of using NumReads from Salmon, it is better to use TPM 
# All you need to change is the expression matrix
# data1  <- read.table('merge_tpm.txt', header = T, row.names = 1)
# data1 <- normalized_counts
data1 <- intergenic_intronic_exon_highThan1_length_highThan_200   # For TPM
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

### How do we pick up the threshold?
### The idea is we need to pick up the genes with certain amounts of expression level and are highly variable. The thing is for lowly-expressed genes, the variable
### (coeffiecient of variance) is commonly larger than those highly expressed genes because of technical noise. Therefore, using var or CV and set them at a certain
### value is not an optimal way to get the results. The good strategy which is also applied in different clever packages is to perform a regression model and select
### those genes above the regression line (more precisely, the prediction value). Here I do the regression model based on lowess/loess regression of log2_mean (x
### axis) and log2_CV (y axis). Log-transformation is used for scale the data and make it more plottable. But before plotting and perform the regression, I need to
### add two columns of log2_mean and log2_CV for each transcript.

n <- nrow(datExpr0); n
datExpr0[n+1,] <- apply(datExpr0[c(1:n),],2, function(x){log2(mean(x)+1)}); dim(datExpr0)
datExpr0[n+2,] <- apply(datExpr0[c(1:n),], 2, function(x){log2(sd(x)/mean(x)+1)}) ; dim(datExpr0)# 使用变异系数（coefficient of variance, CV)较正表达量高低对变异度的影响
datExpr1 <- as.data.frame(t(datExpr0))
names(datExpr1)[21] <- 'log2_mean';names(datExpr1)[22] <- 'log2_CV'; names(datExpr1)
head(datExpr1)[21]; head(datExpr1[22]); colnames(datExpr1)

### Use loess model to fit the curve in order to select the highly variable genes
#### Plotting the regression line and label the highly variable genes based on certain threshold (though subjective)
p <- ggplot(datExpr1, aes(x = log2_mean, y = log2_CV))+ geom_point() + 
  geom_smooth(span = 0.2, method = 'loess', na.rm = T) + 
  geom_smooth(method = lm, col = 'red', na.rm = T) + 
  ylim(c(0.4,2.6)) +
  geom_vline(xintercept = seq(0,2,0.2), col = 'darkgreen', lty = 2) +
  theme_classic() +
  geom_text_repel(data=subset(datExpr1, datExpr1$log2_mean > 2 & datExpr1$log2_CV> 1.25), 
                  aes(label=row.names(datExpr1[datExpr1$log2_mean > 2 & datExpr1$log2_CV> 1.25,])), 
                  col= 'black', nudge_x = 0.5); p

p <- ggplot(datExpr1, aes(x = log2_mean, y = log2_CV))+ geom_point() + 
  geom_smooth(span = 0.2, method = 'loess', na.rm = T) + 
  geom_smooth(method = lm, col = 'red', na.rm = T) + 
  ylim(c(0.4,2.6)) +
  geom_vline(xintercept = seq(0,2,0.2), col = 'darkgreen', lty = 2) +
  theme_classic() +     ## theme_classic() is from ggplot2 package and it is the only suitable theme for scientific publication.
  geom_label(data=subset(datExpr1, datExpr1$log2_mean > 4 & datExpr1$log2_CV> 1.25), 
                  aes(label=row.names(datExpr1[datExpr1$log2_mean > 4 & datExpr1$log2_CV> 1.25,])), 
                  col= 'black', 
                  nudge_x = 0.3,nudge_y = 0.05, fill = 'green'); p

#### Use loess regression to fit the model
model_xlog2mean_ylog2CV <- loess(datExpr1$log2_CV ~ datExpr1$log2_mean, span = 0.2, method = 'loess')
summary(model_xlog2mean_ylog2CV)
#### Use prediction to predict the y value (log2_CV) for each x (log2_mean)
prediction <- predict(object = model_xlog2mean_ylog2CV, data.frame(datExpr1$log2_mean), se = T)
head(prediction$fit)  ## get the y value (log2_CV) point prediction
#### Further filtering according to the predicted y value point prediction
datExpr0 <- datExpr1[datExpr1$log2_CV > prediction$fit & datExpr1$log2_mean > 2,1:20]; dim(datExpr0)  ## setting log2_mean > 2, I only get 109 condidate transcripts.

filtered_TPM_normalized_counts <- datExpr0
head(filtered_TPM_normalized_counts)
dim(filtered_TPM_normalized_counts)

### After selection of HVG using loess, let's have a look at the clustering and see whether the filtering highly variable genes are good enough to distinguish 
### different samples (supervised learning)
library('gplots')
library('pheatmap')
library('amap')
library('RColorBrewer')
pearson_cor <- as.matrix(cor(datExpr0, method = 'pearson'))
head(pearson_cor)
hc <- hcluster(t(datExpr0), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)

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


##################### Visualizing the gene network #######################
nSelect <- 109
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
# 比如对module = ‘blue’ 的模块进行分析
# 'blue' module gene
module <- 'blue'
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
blue_module_index <- which(moduleColors == 'blue') # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[blue_module_index])
length(rownames(filtered_TPM_normalized_counts)[blue_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
blue_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[blue_module_index]
length(blue_module_transcriptName); blue_module_transcriptName
# 'blue module 有12个基因



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
## Reset workding directory to export the node and edge data
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_data/WGCNA/WGCNA')
## Select each module
for (mod in 1:nrow(table(moduleColors)))
{
  
  modules = names(table(moduleColors))[mod]
  # Select module probes
  probes = colnames(t(datExpr0))
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



#======================================================================================
#
#                            Cytoscape network analysis 
#  (MCODE and cytohubba packages) for hub network and hub gene identification 
#
#======================================================================================



#======================================================================================
#
#                               Data frame preparation  
#
#======================================================================================
### So, Now we get the gene list, next we need to do is to get several information about the genes for downstream analysis.
### There are several things we need to do: (1), we need to pick up the promoter regions and use HOMER software (perl script)
### for TFBS prediction; (2), we need to pick up the flanking genes and perfrom enrichment analysis, and see whether the
### flanking genes are responsible for certain functions; (3), we need to get the expression matrix for the transcripts for
### each sample for downstream circos plot visualiation.

### So, what we need to do now is to get the coordinate of each transcripts, we have the coordinate and stored in assembly.bed
### file and let's first have a look at it in R; I have put the chr, start, end, transcript_id in a new file named 
### assembly_coordinate.txt in the working directory (with all 98528 total transcripts information).
coordiate_total_transcript <- read.table('assembly_coordiate.txt', header = F, sep = ' ')
colnames(coordiate_total_transcript) <- c('chr', 'start', 'end', 'transcript_ID')
head(coordiate_total_transcript); dim(coordiate_total_transcript);View(coordiate_total_transcript)


### Pick up the transcripts information of beta gene-module, the transcripts names are stored in 'blue_module_transcriptName' object
blue_beta_coordinate <- coordiate_total_transcript[coordiate_total_transcript$transcript_ID %in% as.vector(blue_module_transcriptName),]
blue_beta_coordinate
#### Here we see that there are several transcripts that are not matched to the autosomal chromosome, but now we don't filter them

### Next, we add the expression information, the TPM expression matrix is stored in 'data1' object
blue_beta_expr <- data1[row.names(data1) %in% blue_beta_coordinate$transcript_ID,]
blue_beta_expr <- cbind(blue_beta_coordinate, blue_beta_expr)
#### The coordination and expression matrix are stored in 'blue_beta_expr' object

################## Pick up the promoter regions ###################
### Normally we define the promoter regions as 500 bp upstream of TSS. Because this is unstranded library, we need to check
### the both ends. It is good to use 'dplyr' package here.
library('dplyr')
promoter_positive_strand <- mutate(.data = blue_beta_coordinate, Start = start -500, End = start)[, c(1,5,6,4)]
promoter_positive_strand
promoter_negative_strand <- mutate(.data = blue_beta_coordinate, Start = end, End = end + 500)[,c(1,5,6,4)]
promoter_negative_strand


### Write the table into txt form and the files would be used as input for HOMER
### Or we can use 'TFBSTools' and 'JASPAR2018' packages, together with 'Biostring'






################## Pick up flanking genes ###################
### Normally, the lincRNAs perform their jobs in cis, which means they have their functions locally. It is very important to 
### identify their flanking genes and the enrichment analysis will give some hints. We can use biomaRt to get the flanking 
### genes (entrezgene). But the first thing we need to do is to pick up the flanking regions. Normally, we check 10-20 kb regions
### So, we can divide our work into two steps:
### STEP1: get the flank coordiantes in the genome
### STEP2: use biomaRt to get the flank gene entrezgene id






