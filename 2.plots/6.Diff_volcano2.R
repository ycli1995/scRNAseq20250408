
me = sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
lib.dir <- dirname(dirname(normalizePath(dirname(me))))

library(Seurat)
library(dplyr)
library(rlang)
library(ggplot2)
library(data.table)

source(file.path(lib.dir, "seuobj_lib.R"), chdir = TRUE)

args <- commandArgs(T)

file <- args[1]
in_file <- args[2]
outfile <- args[3]

if (is.null(file) | is.na(file)) {
  warning("\n  Usage: Rscript ", me, " <parameter.yaml>\n")
  quit()
}

file <- tools::file_path_as_absolute(file)
in_file <- tools::file_path_as_absolute(in_file)

handlers <- YAML_HANDLERS
param <- yaml::yaml.load_file(file, handlers = handlers)

param$cluster_column <- param$cluster_column %||% "Cluster"
param$remove.infinite <- param$remove.infinite %||% FALSE

markers = fread(in_file, header = T, sep = "\t", stringsAsFactors = F, data.table = F)

if (nrow(markers) == 0) {
  warning("\n  No differential genes!")
  quit()
}

markers[[param$cluster_column]] = factor(
  markers[[param$cluster_column]], 
  levels = stringr::str_sort(unique(markers[[param$cluster_column]]))
)

## Rename clusters
if (!is.null(param$cluster_rename)) {
  markers[[param$cluster_column]] = ReplaceEntries(markers[[param$cluster_column]], param$cluster_rename)
}

## Set colors
all_color = SetColor(markers[[param$cluster_column]], "tsne", "set1")
if (!is.null(param$cluster_colors)) {
  param$cluster_colors = unlist(param$cluster_colors)
  all_color[names(param$cluster_colors)] = param$cluster_colors
}

## Subset clusters
all_cl = param$clusters
if (!is.null(all_cl)) {
  all_cl = levels(markers[[param$cluster_column]])
}
markers$cluster = markers[[param$cluster_column]]
markers = droplevels(markers)
all_color = all_color[levels(markers$cluster)]

## Format df
markers$gene = markers[[param$gene_column]]
markers$log2FC = markers[[param$log2FC_column]]
markers$p_val = markers[[param$pval_column]]
markers$p_val_adj = markers[[param$qval_column]]

df = markers[, c("gene", "log2FC", "p_val", "p_val_adj", "cluster")]

## Deal with infinite values
if (param$remove.infinite) {
  df = df[!is.infinite(df$log2FC), , drop = FALSE]
  df = droplevels(df)
} else {
  df$log2FC[is.infinite(df$log2FC) & df$log2FC > 0] = max(df$log2FC[!is.infinite(df$log2FC) & df$log2FC > 0])
  df$log2FC[is.infinite(df$log2FC) & df$log2FC < 0] = min(df$log2FC[!is.infinite(df$log2FC) & df$log2FC < 0])
}

## Plot
## https://github.com/junjunlab/scRNAtoolVis/blob/master/R/jjVolcano.R
jjVolcano <- function(
  diffData = NULL,
  type = c("updown", "sig"),
  
  log2FC.cutoff = 0.25,
  sig.by = c("pvalue", "qvalue", "both"),
  pvalue.cutoff = 0.05,
  qvalue.cutoff = 0.05,
  
  order.by = "log2FC", # c("log2FC","p_val")
  topGeneN = 5,
  color.div = c('#0099CC','#CC3333'),
  
  tile.col = NULL,
  cluster.order = NULL,
  
  back.col = 'grey93',
  pSize = 0.5,
  
  base_size = 14,
  
  polar = FALSE,
  expand = c(-1, 1),
  flip = FALSE,
  ...
) {
  type <- match.arg(type)
  sig.by <- match.arg(sig.by)
  col.type = type

  filter.p <- switch(
    EXPR = sig.by,
    "pvalue" = TRUE,
    "qvalue" = FALSE,
    "both" = TRUE
  )
  filter.q <- switch(
    EXPR = sig.by,
    "pvalue" = FALSE,
    "qvalue" = TRUE,
    "both" = FALSE
  )
  
  # filter data
  diff.marker <- diffData %>%
    dplyr::filter(abs(log2FC) >= log2FC.cutoff)
  if (filter.p) {
    diff.marker <- diff.marker %>%
      dplyr::filter(p_val < pvalue.cutoff)
  }
  if (filter.q) {
    diff.marker <- diff.marker %>%
      dplyr::filter(p_val_adj < qvalue.cutoff)
  }
  
  # assign type
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(log2FC >= log2FC.cutoff, "sigUp", "sigDown"))
  diff.marker <- switch(
    EXPR = sig.by,
    "pvalue" = diff.marker %>% 
      dplyr::mutate(
        type2 = ifelse(p_val < pvalue.cutoff,
                       paste("P < ", pvalue.cutoff, sep = ''),
                       paste("P >= ", pvalue.cutoff, sep = ''))
        ),
    "qvalue" = diff.marker %>% 
      dplyr::mutate(
        type2 = ifelse(p_val_adj < qvalue.cutoff,
                       paste("P_adj < ", qvalue.cutoff, sep = ''),
                       paste("P_adj >= ", qvalue.cutoff, sep = ''))
      ),
    "both" = diff.marker %>% 
      dplyr::mutate(
        type2 = ifelse(p_val < pvalue.cutoff & p_val_adj < qvalue.cutoff,
                       paste("P < ", pvalue.cutoff, " & P_adj < ", qvalue.cutoff, sep = ''),
                       paste("P >= ", pvalue.cutoff, " & P_adj >= ", sep = ''))
      )
  )
  
  # cluster orders
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(as.character(diff.marker$cluster), levels = cluster.order)
  }
  
  # get background cols
  back.data <- purrr::map_df(unique(diff.marker$cluster), function(x){
    tmp <- diff.marker %>% dplyr::filter(cluster == x)
    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$log2FC) - 0.2,
      max = max(tmp$log2FC) + 0.2
    )
    return(new.tmp)
  })
  
  # get top gene
  top.marker.tmp <- diff.marker %>%
    dplyr::group_by(cluster)
  
  # order
  # if(length(order.by) == 1){
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::group_by(cluster) %>%
  #     dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  #
  # }else{
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_head(n = topGeneN)
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_tail(n = topGeneN)
  # }
  
  top.marker.max <- top.marker.tmp %>%
    dplyr::slice_max(n = topGeneN, order_by = get(order.by), with_ties = F)
  top.marker.min <- top.marker.tmp %>%
    dplyr::slice_min(n = topGeneN, order_by = get(order.by), with_ties = F)
  
  # combine
  top.marker <- rbind(top.marker.max, top.marker.min)
  
  if (is.null(x = tile.col)) {
    tile.col <- scales::hue_pal()(length(unique(diff.marker$cluster)))
  }
  
  # ====================================================================
  # plot
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, y = log2FC)) +
    # add back cols
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = min), fill = back.col) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = max), fill = back.col)
  
  # ap1 <- paste("adjust Pvalue >= ",adjustP.cutoff,sep = '')
  # ap2 <- paste("adjust Pvalue < ",adjustP.cutoff,sep = '')
  
  # color type
  if (col.type == "updown") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize) +
      ggplot2::scale_color_manual(values = c("sigDown" = color.div[1], "sigUp" = color.div[2]))
  } else if (col.type == "adjustP") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type2), size = pSize) +
      ggplot2::scale_color_manual(values = c(color.div[2], color.div[1]))
  }
  
  # theme details
  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank()) +
    ggplot2::xlab('Clusters') + ggplot2::ylab('log2FC') + 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5, fill = NA)))
  
  # add tile
  p4 <- p3 +
    ggplot2::geom_tile(data = diff.marker, ggplot2::aes(x = cluster,y = 0,fill = cluster),
                       color = 'black',
                       height = log2FC.cutoff * 2,
                       #alpha = 0.3,
                       show.legend = T) +
    ggplot2::scale_fill_manual(values = tile.col) +
	#ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) +
    # add gene label
    ggrepel::geom_text_repel(data = top.marker,
                             ggplot2::aes(x = cluster,y = log2FC, label = gene),
                             max.overlaps = 50,
                             ...)
  
  # whether coord_plolar
  if (polar == TRUE) {
    p5 <- p4 +
      #geomtextpath::geom_textpath(ggplot2::aes(x = cluster,y = 0,label = cluster)) +
      ggplot2::scale_y_continuous(n.breaks = 6,
                                  expand = ggplot2::expansion(mult = expand)) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
                     legend.title = ggplot2::element_blank()) +
      ggplot2::coord_polar(clip = 'off',theta = 'x')
  } else {
    # whether flip plot
    if (flip == TRUE) {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = cluster,y = 0,label = cluster)) +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::coord_flip()
    }else{
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
      #  ggplot2::geom_text(ggplot2::aes(x = cluster,y = 0,label = cluster)) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
  }
  return(p5)
}

p = jjVolcano(df, type = param$type, log2FC.cutoff = param$log2FC_thres, pvalue.cutoff = param$pval_thres, qvalue.cutoff = param$qval_thres, topGeneN = param$top_n, color.div = param$color.div, tile.col = all_color, sig.by = param$sig.by)

#w = max(7, 2 + length(unique(df$cluster)) * 0.08 * max(nchar(as.character(df$cluster))))
w = max(7, 1.5 + length(unique(df$cluster)) * 0.75 + 0.08 * max(nchar(as.character(df$cluster))))
ggsave(p, file = outfile, height = 8, width = w, limitsize = FALSE)

