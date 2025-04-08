
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

source("Seurat_lib.R", chdir = T)

.PlotClusterStat = function(object, stat.what = "seurat_clusters", group.by = "orig.ident", color.st = NULL, color.gb = NULL, outpref = NULL, ...){
		if ( class(object) == "Seurat" ) {
				metadata <- object@meta.data
		} else {
				metadata <- object
		}
		if ( is.null(color.st) ) {
				if ( "misc" %in% slotNames(object) && exists(stat.what, object@misc) ) {
						color.st <- object@misc[[stat.what]]
				} else {
				color.st <- switch(stat.what,
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		if ( is.null(color.gb) ) {
				if ( "misc" %in% slotNames(object) && exists(group.by, object@misc) ) {
						color.gb <- object@misc[[group.by]]
				} else {
				color.gb <- switch(group.by, 
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		name.st <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		name.gb <- switch(group.by,  "seurat_clusters" = "Cluster", "orig.ident" = "Samples", group.by)
		stat.what <- as.name(stat.what)
		group.by  <- as.name(group.by)

		stat_sample <- metadata %>%
				group_by(!! name.gb := !! group.by, !! name.st := !! stat.what) %>%
				summarise("Number of cells" = n())

		p <- list()
		p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), y = ~ `Number of cells`, fill = as.name(name.st)))
		p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), y = ~ `Number of cells`, fill = as.name(name.gb)))
		if ( ! is.null(color.st) ) p[["by"]] <- p[["by"]] + scale_fill_manual(values = color.st)
		if ( ! is.null(color.gb) ) p[["in"]] <- p[["in"]] + scale_fill_manual(values = color.gb)
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )

		if ( is.null(outpref) ) {
				outpref <- paste0( name.st, ".stat")
		}
		for ( i in names(p) ) {
				p[[i]] <- p[[i]] + bar_theme_default() + theme(plot.margin = margin(l = 20))
				if (i == "by") w = length(unique(stat_sample[, name.gb])) * 1.9 + 3.5
				if (i == "in") w = length(unique(stat_sample[, name.st])) * 1.9 + 3
				ggsave( p[[i]] + geom_stack, file = paste0( outpref, ".", i, name.gb, ".pdf"), height = 6, width = w )
				ggsave( p[[i]] + geom_fill + ylab("Fraction of Cells"),  file = paste0( outpref, ".", i, name.gb, ".pct.pdf"), height = 6, width = w )
		}
}

StatCluster = function(object, group.by = "orig.ident", outpref = "Cluster.stat", stat.what = "seurat_clusters", assay = DefaultAssay(object), ...){
		.StatCluster(object, stat.what = stat.what, outpref = outpref, assay = assay)
		.StatCluster_by(object, stat.what = stat.what, group.by = group.by, outpref = outpref)
		.PlotClusterStat(object, stat.what = stat.what, group.by = group.by, outpref = outpref, ...)
}

obj = Load("obj_renamed.Rda")

# 1. umap
#obj$order_clusters = factor(as.character(obj$seurat_clusters), levels = names(sort(table(obj$seurat_clusters), decreasing = T)))
obj$order_clusters = obj$seurat_clusters
obj@misc$color.cluster2 = obj@misc$color.cluster[levels(obj$order_clusters)]

PlotCluster(obj, outpref = "UMAP", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, p2.label = F, outpref = "UMAP.nolabel", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, outpref = "UMAP.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, p2.label = F, outpref = "UMAP.nolabel.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", outpref = "TSNE", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", p2.label = F, outpref = "TSNE.nolabel", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", outpref = "TSNE.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", p2.label = F, outpref = "TSNE.nolabel.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters", p2.color = obj@misc$color.cluster2)

# 2. Stat
StatCluster(obj, stat.what = "seurat_clusters", color.st = obj@misc$color.cluster2)
StatCluster(obj, "Groups", stat.what = "seurat_clusters", color.st = obj@misc$color.cluster2)

cell = data.frame(
	Cells = Cells(obj),
	Cluster = obj$seurat_clusters,
	Groups = obj$Groups
)
color = obj@misc$color.cluster2

library(ggrepel)
for (i in levels(cell$Groups)) {
  df = cell[cell$Groups == i, ]
  pdata = df %>% group_by(Cluster) %>% summarise(Freq = length(Cells), .groups = 'drop')
 
  pdata$fraction = pdata$Freq / sum(pdata$Freq)
  pdata$ymax = cumsum(pdata$fraction)
  pdata$ymin = c(0, pdata$ymax[-nrow(pdata)])
  pdata$Cluster2 = paste0(round(pdata$fraction * 100, digits = 2), '%')
#  names(color) = pdata$Cluster2
  p = ggplot(pdata, aes(fill = Cluster, ymax = ymax, ymin = ymin, xmax = 3.5, xmin = 2, col = Cluster)) +
    geom_rect() + coord_polar('y') + xlim(c(0, 4)) +
    scale_fill_manual(values = color[pdata$Cluster]) + 
    scale_color_manual(values = color[pdata$Cluster]) +
    geom_label_repel(aes(label = Cluster2, x = 3.5, y = (ymax + ymin)/2), show_guide  = FALSE, color = 'black', size = 4,  max.overlaps = Inf) +
    labs(fill = 'Cluster', color = 'Cluster') +
    theme_void() +
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(paste0(i, '.cluster_pie.pdf'), p, w = 7, h = 5)
}

