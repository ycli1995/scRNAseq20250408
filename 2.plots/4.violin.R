
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

library(reshape2)
library(ggpubr)

source("Seurat_lib.R", chdir = T)

obj = Load("obj_renamed.Rda")

genes = readLines("genes.list")

obj = subset(obj, seurat_clusters %in% c("Mesenchymal cells",  "Epithelial cells", "Perivascular cells"))
obj@meta.data = droplevels(obj@meta.data)

features = FindFeaturesID(obj, genes)
df = FetchData(obj, c("Groups", "seurat_clusters", features)) %>%
  rename(FindFeaturesName(obj, features, "name") %>% setNames(names(.), .)) %>%
  reshape2::melt(c("Groups", "seurat_clusters")) %>%
  rename(GeneName = variable, Expression = value)

for (i in levels(obj$seurat_clusters)) {
  i2 = gsub("\\s|\\/", "_", i)
  pos = position_dodge(width = 1)
  for (j in unique(df$GeneName)) {
    df2 = df %>%
      filter(GeneName %in% j & seurat_clusters %in% i)
    p = ggplot(df2, aes(Groups, Expression, fill = Groups)) +
      geom_violin(position = pos, width = 1) +
      geom_boxplot(outlier.shape = NA, position = pos, width=0.1, fill = "white") +
      scale_fill_manual(values = obj@misc$color.group) +
      labs(x = "") +
      bar_theme_default() +
      RotatedAxis()
    h = 4.5
    w = 4
    ggsave(p, filename = paste0(j, ".", i2, ".Violin.pdf"), height = h, width = w)

    p = ggplot(df2, aes(Groups, Expression, fill = Groups)) +
      geom_violin(position = pos, width = 1) +
      #geom_boxplot(outlier.shape = NA, position = pos, width=0.2) +
      scale_fill_manual(values = obj@misc$color.group) +
      labs(x = "") +
      bar_theme_default() +
      RotatedAxis()
    h = 4.5
    w = 4
    ggsave(p, filename = paste0(j, ".", i2, ".Violin.no_box.pdf"), height = h, width = w)
  }
}

