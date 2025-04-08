
### Deal arguements
args <- commandArgs(T)

indir   <- args[1]
sample  <- args[2]
add_lib <- args[3]
outdir  <- args[4]
rate    <- args[5]

if ( is.null( file ) | is.na( file ) ){
		warning( "\n  Usage : Seurat.R <parameter.yaml> (<outdir>)\n" )
		quit()
}



library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)
#plan("sequential")

source(add_lib, chdir = TRUE)

### Let's shake it
if ( ! is.na(outdir) ) setwd(outdir)

### Creat Seurat Object
message( "==>Reading 10x data<==" )
obj <- MakeSeuratObj(data_name = sample, data_dir = indir)


### Normalization Data
message( "==>Normalization Data<==" )
obj <- DoNormalization(obj, vars.regress = "none", is_SCTransform = FALSE, is.check = FALSE, scale.only.var.genes = TRUE)

### Reduce dimension
message( "==>Reduce dimension<==" )
obj <- DoDimReduc(obj, is.checkpca = FALSE, check_duplicates = FALSE)
dims <- seq(obj@reductions$pca)


### get pN
pN <- 0.25

### get pK
sweep.res.list <- paramSweep_v3(obj, PCs = dims, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)[1]]))

pdf(paste0("pK.", sample, ".pdf"))
plot(x = as.vector(bcmvn$pK), y = bcmvn$BCmetric, type = "b", xlab = "pK", ylab = "BCmetric", col = "blue", pch = 19)
abline(v = pK, lty = 2, col = "red")
dev.off()
WriteTable(bcmvn, file = paste0("pK.", sample, ".xls"))

### get nExp 
if ( is.na(rate) ) 
		rate <- 7.6 * 10^-6 * ncol(obj) + 5.27 * 10^-4
nExp_poi <- round(as.numeric(rate) * ncol(obj)) 

if( exists("seurat_clusters", obj@meta.data) ){
		annotations <- obj@meta.data$seurat_clusters
		homotypic.prop <- modelHomotypic(annotations)
		nExp_poi <- round(nExp_poi * (1 - homotypic.prop))
}

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
obj <- doubletFinder_v3(obj, PCs = dims, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(obj@meta.data)[grep('pANN', colnames(obj@meta.data))] <- "pANN"
colnames(obj@meta.data)[grep('DF.classifications', colnames(obj@meta.data))] <- "classifications"

embeddings <- cbind(obj@reductions[["tsne"]]@cell.embeddings, obj@reductions[["umap"]]@cell.embeddings)
data <- obj@meta.data[,c("pANN", "classifications")]
WriteTable(cbind(Cells = rownames(data), data), file = paste0("DF.classify.", sample, ".xls"))
WriteTable(cbind(Cells = rownames(data), data, embeddings), file = paste0("DF.classify.", sample, ".tmp"))

p1 <- DimPlot(obj, reduction = "umap", group.by = "classifications", cols = c("Singlet" = "black", "Doublet" = "red")) + dot_theme_default() + ggtitle(NULL)
ggsave(p1, file = paste0("DF.classify.UMAP.", sample, ".pdf"), width = 6, height = 5)
p2 <- DimPlot(obj, reduction = "tsne", group.by = "classifications", cols = c("Singlet" = "black", "Doublet" = "red")) + dot_theme_default() + ggtitle(NULL)
ggsave(p2, file = paste0("DF.classify.tSNE.", sample, ".pdf"), width = 6, height = 5)



PlotFeaturePlot(obj, "pANN", 'umap', outfile = paste0("pANN.UMAP.", sample, ".pdf"))
PlotFeaturePlot(obj, "pANN", 'tsne', outfile = paste0("pANN.tSNE.", sample, ".pdf"))

p <- ggplot(obj@meta.data, aes(x = pANN)) + geom_histogram() + bar_theme_default()
ggsave(p, file = paste0("pANN.hist.", sample, ".pdf"), width = 7, height = 7)

PlotBasicStat(obj, "BasicInfo", group.point.by = "classifications", group.point.color = c("Singlet" = "black", "Doublet" = "red"))

