
### Deal arguements
args <- commandArgs(T)

file    <- args[1]
outdir  <- args[2]
add_lib <- args[3]
if ( is.null( file ) | is.na( file ) ){
		warning( "\n  Usage : Seurat.R <parameter.yaml> (<outdir>)\n" )
		quit()
}


### Loading Library
handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parameter <- yaml::yaml.load_file( file, handlers = handlers)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)
#plan("sequential")

#source("/home/xushuyang/Pipeline/R/SCellWare/R/Seurat_lib.R", chdir = T)
if ( ! is.na(add_lib) ) source(add_lib, chdir = T)


### Let's shake it
if ( ! is.na(outdir) ) setwd(outdir)

### Creat Seurat Object
message( "==>Reading 10x data<==" )
obj <- MakeSeuratObj(parameter)
obj <- OverrideFeatures(obj, parameter)

### Add some check flag
message( "==>Adding MetaData<==" )
if ( ! is.null(parameter$marker$mito_list) ) obj <- StatFeatures(obj, parameter$marker$mito_list, col.name = "percent.mito", stat_pct = T, add_to_pdata = T)
if ( ! is.null(parameter$marker$plastid_list) ) obj <- StatFeatures(obj, parameter$marker$plastid_list, col.name = "percent.plastid", stat_pct = T, add_to_pdata = T)
if ( ! is.null(parameter$marker$expected) )  obj <- StatFeatures(obj, parameter$marker$expected,  col.name = "expected.marker")
if ( ! is.null(parameter$marker$excluded) )  obj <- StatFeatures(obj, parameter$marker$excluded,  col.name = "exclude.marker")
if ( ! is.null(parameter$marker$more) )      obj@misc[["more.marker"]] <- FindFeaturesID(obj, parameter$marker$more, unlist = FALSE)
WriteTable(tibble::rownames_to_column(obj@meta.data, var = "Cells"), file = "metadata.xls")

### Data Stat - before filter
message( "==>Stat before BasicInfo<==" )
if ( is.null(parameter$filter$filter.cells) ) {
		PlotBasicStat(obj, "BasicInfo", nRow = 1)
		if ( ! is.null(parameter$Groups) ) {
				PlotBasicStat(obj, "BasicInfo.groups", group.by = "Groups", nRow = 1)
		}
} else {
		obj[["DF"]] <- "Singlet"
		filter.cells <- readLines(parameter$filter$filter.cells)
		obj@meta.data[filter.cells, "DF"] <- "Doublet"
		PlotBasicStat(obj, "BasicInfo", group.point.by = "DF", group.point.color = c("Singlet" = "black", "Doublet" = "red"), nRow = 1)
		if ( ! is.null(parameter$Groups) ) {
				PlotBasicStat(obj, "BasicInfo.groups", group.by = "Groups", group.point.by = "DF", group.point.color = c("Singlet" = "black", "Doublet" = "red"), nRow = 1)
		}
}


### Filter
message( "==>Filter<==" )
obj <- FilterGenes(obj, parameter)
obj <- FilterCells(obj, parameter, do.stat = FALSE)
StatFilterCells(obj, group.by = "orig.ident", outfile = "Filter.stat.xls")
if ( ! is.null(parameter$Groups) ) {
		StatFilterCells(obj, group.by = "Groups", outfile = "Filter.stat.groups.xls")
}

### Data Stat - after filter
message( "==>Stat after BasicInfo<==" )
PlotBasicStat(obj, "AfterFilter.BasicInfo", nRow = 1)
if ( ! is.null(parameter$Groups) ) {
		PlotBasicStat(obj, "AfterFilter.BasicInfo.groups", group.by = "Groups", nRow = 1)
}

### Normalization Data
message( "==>Normalization Data<==" )
obj <- DoNormalization(obj, parameter, is_SCTransform = FALSE,
		scale.only.var.genes = TRUE,
		vfeature.must = obj@misc[['expected.marker']],
		vfeature.remove = obj@misc[['exclude.marker']])

### Reduce dimension
message( "==>Reduce dimension<==" )
obj <- DoDimReduc(obj)

### Find clusters
message( "==>Find clusters<==" )
obj <- DoFindClusters(obj, reduction = "pca", dims = NULL, resolution = parameter$cluster_resolution)

if ( ! is.null(parameter$integration$method) ) {
		if ( length(table(obj[["orig.ident"]])) > 1 ) {
				## check before integration visualization
				obj[["beforeInteg.cluster"]] <- Idents(object = obj)
				PlotCluster(obj, reduction = 'umap_RNA', outpref = "UMAP_before" )
				PlotCluster(obj, reduction = 'tsne_RNA', outpref = "tSNE_before" )
				if ( ! is.null(parameter$Groups) ) {
						PlotCluster(obj, reduction = 'umap_RNA', outpref = "UMAP_before.groups", split.by = "Groups", p1.group.by = "Groups" )
						PlotCluster(obj, reduction = 'tsne_RNA', outpref = "tSNE_before.groups", split.by = "Groups", p1.group.by = "Groups" )
				}

				### Integration
				message( "==> Do Integration <==" )
				if ( parameter$integration$method == "CCA" ) {
						obj <- DoIntegration(obj, split.by = "orig.ident")
						obj <- DoDimReduc(obj)
						obj <- DoFindClusters(obj, reduction = "pca", resolution = parameter$cluster_resolution)
				} else {
						obj <- RunHarmony(obj, group.by.vars = "orig.ident", project.dim = FALSE, assay.use = DefaultAssay(obj))
						obj <- DoDimReduc(obj, reduction = "harmony", reduction.surfix = "harmony")
						obj <- DoFindClusters(obj, reduction = "harmony", resolution = parameter$cluster_resolution)
				}
		}
}

## Draw t-SNE plot
message( "==>Draw t-SNE plot<==" )
PlotCluster(obj, reduction = 'umap', outpref = "UMAP" )
PlotCluster(obj, reduction = 'tsne', outpref = "tSNE" )
if ( ! is.null(parameter$Groups) ) {
		PlotCluster(obj, reduction = 'umap', outpref = "UMAP.groups", split.by = "Groups", p1.group.by = "Groups" )
		PlotCluster(obj, reduction = 'tsne', outpref = "tSNE.groups", split.by = "Groups", p1.group.by = "Groups" )
}

### Save data object 
message( "==>Output obj.Rda<==" )
DefaultAssay(obj) <- "RNA"
save(obj, file = "obj.Rda")

### stat table 
message( "==>Stat table<==" )
StatCluster(obj)
if ( ! is.null(parameter$Groups) ) {
		StatCluster(obj, "Groups")
}
CalAvgExp(obj)
CalAvgExp(obj, group.by = "orig.ident", outfile = "AllGene.avg_exp.Samples.xls")
CalPctExp(obj, outfile = "AllGene.avg_pct.xls")
CalPctExp(obj, group.by = "orig.ident", outfile = "AllGene.avg_pct.Samples.xls")
ListCellCluster(obj)
PlotPresetMarker(obj)


### Find maker genes
message( "==>Find maker genes<==" )
obj.markers <- DoFindAllMarkers(obj, parameter)

message( "==>Output markers.Rda<==" )
save( obj.markers, file = "markers.Rda" )
#obj.markers$gene <- ChangeOUTName(obj.markers$gene, object@misc$fdata)

## stat marker
message( "==>stat marker<==" )
StatMarker(obj.markers, color = obj@misc$color.cluster)
ListMarker(obj, obj.markers)

### Top marker
message( "==>display top markers<==" )
top <- FindTopMarker(obj.markers, top_num = parameter$heatmap$top, object = obj)
PlotAboutFeatures(obj, features = unique(top$gene), outpref = "Top")

dir.create("DensityPlot/", showWarnings = F, recursive = T)
unlink("DensityPlot/*", recursive = T)
PlotDensityPlot(obj, unique(top$gene), reduction = 'umap', outpref = "DensityPlot/")

dir.create("ExpPlot/", showWarnings = F, recursive = T)
unlink("ExpPlot/*", recursive = T)
PlotFeaturePlot(obj, unique(top$gene), reduction = 'umap', outpref = "ExpPlot/ExpPlot", is.combine = FALSE)

dir.create("ViolinPlot/", showWarnings = F, recursive = T)
unlink("ViolinPlot/*", recursive = T)
PlotVlnPlot(obj, unique(top$gene), outpref = "ViolinPlot/ViolinPlot")

### Hasta la vista, baby
message( "==>All Done!<==" )


