
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

library(reshape2)
library(ComplexHeatmap)

source("Seurat_lib.R", chdir = T)

genes = readLines("genes.list")

df = "filter.annot.xls.tmp" %>%
    data.table::fread(sep = "\t", header = TRUE, stringsAsFactors = F, data.table = F, fill = TRUE) %>%
    filter(Index %in% genes) %>%
    tibble::column_to_rownames("Index")

avg = df[, 1:2]
rna = ScaleData(avg)
rownames(rna) = paste0("m/z: ", df$mz)

cols = circlize::colorRamp2(c(-2, 0, 2), c("#3976B0", "white", "#BC323B"))
hmp = Heatmap(
  rna, col = cols,
  name = " ", 
  show_row_names = T, 
  cluster_columns = F,
  cluster_rows = T,
  column_title_side = "bottom",
  column_title_rot = 0, column_names_rot = 90
)
pdf("Heatmap.pdf", width = 4, height = 5)
draw(hmp, padding = unit(c(25, 2, 2, 2), "mm"))
dev.off()


