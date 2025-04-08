Q3Normalize <- function(count) {
  qt <- sapply(count, quantile, probs = 0.75, names = F)
  mean_Q3 <- exp(mean(log(qt)))
  data <- as.data.frame(as.matrix(sapply(names(qt), function(i) mean_Q3 / qt[[i]] * count[[i]])))
  rownames(data) <- rownames(count)
  return(data)
}

DotPlot <- function(meta, var = "SegmentDisplayName", outfile = "QC_Aligned_Rate.pdf"){
		dt <- meta %>% select(SegmentDisplayName, RawReads, AlignedReads) %>%
				mutate(rate = AlignedReads/RawReads)
		p <- ggplot(dt, aes(x = RawReads, y = rate)) +
				geom_point(color = "royalblue") +
				geom_hline(yintercept = 0.8, linetype = "dotted", color = "brown") +
				geom_text(aes(x = mean(range(RawReads)), y = 0.82, label = "Aligned Rate = 80%"), vjust = 0) +
				scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
				labs(x = "Raw Reads", y = "Aligned Rate(%)") +
				theme_minimal()
		ggsave(p, file = outfile, width = 5, height = 5)
}

BarPlot <- function(count, meta, feature = "Negative Probe",
                        var = "SegmentDisplayName", var.color = "royalblue",
                        group.by = "SegmentLabel", group.cols = NULL,
                        out.id = "Neg", outpref = paste0("QC_", out.id, "_barplot"),
                        scale_y = NULL, yintercept = NULL) {
  if ( ! all(colnames(count) %in% meta[[var]]) ) {
    stop()
  }
  input <- count[feature,, drop = FALSE]
  input <- exp(colMeans(log(input)))
  dt <- tibble::rownames_to_column(as.data.frame(input), var = var) %>% left_join(x = meta)
  dt[[var]] <- factor(dt[[var]], levels = dt[[var]][order(dt[["input"]], decreasing = FALSE)])
  
  y    <- "input"
  ylab <- if(length(feature) == 1) paste0(feature, " Counts") else paste0(out.id, "_geomean")
  if (! is.null(scale_y)) {
    y    <- paste0(scale_y, "(", y,    ")")
    ylab <- paste0(scale_y, "(", ylab, ")")
  }
  
  p <- ggplot(dt, aes_string(x = var, y = y)) +
    labs(x = "AOIs", y = ylab) +
    theme_minimal() +
    coord_flip()
  
  if (! is.null(yintercept)) {
    p <- p + geom_hline(yintercept = yintercept, linetype = "dotted")
  }
  
  p1 <- p + geom_bar(fill = var.color, stat = "identity")
  ggsave(p1, file = paste0(outpref, ".pdf"), width = 7, height = 7) ## TODO
  
  for ( i in seq(group.by) ) {
	name <- group.by[[i]]
    color <- if ( is.null(group.cols) || ! exists(name, group.cols) ) {
      SetColor(dt[[name]], tag = i)
    } else {
      group.cols[[name]]
    }
    p2 <- p + geom_bar(aes_string(fill = name), stat = "identity") +
      scale_fill_manual(values = color)
    ggsave(p2, file = paste0(outpref, ".", name,".pdf"), width = 7, height = 7) ## TODO
  }
}

BarPlot.HK <- function(count, meta, HK.gene, var = "SegmentDisplayName",
                           group.by = "SegmentLabel", group.cols = NULL,
                           outpref = "QC_HK_barplot") {
  if ( ! all(colnames(count) %in% meta[[var]]) ) {
    stop()
  }
  input <- count[HK.gene,,drop = FALSE]
  if ( length(HK.gene) == 1 ) {
    input <- exp(colMeans(log(input))) ## same as t(input)
    ylab <- paste0(HK.gene, " Counts")
  }else {
    input <- exp(colMeans(log(input)))
    ylab <- paste0(HK.gene, "_geomean")
  }
  dt <- tibble::rownames_to_column(as.data.frame(input), var = var) %>% left_join(x = meta)
  dt[[var]] <- factor(dt[[var]], levels = dt[[var]][order(dt[["input"]], decreasing = FALSE)])
  
  p <- ggplot(dt, aes_string(x = var, y = "log2(input)")) +
    labs(x = "AOIs", y = "log2(HK_geomean)") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          legend.direction = "horizontal", 
          legend.position = c(0, 1),
          legend.justification = c(0,1))
  
  for ( i in seq(group.by) ) {
	name <- group.by[[i]]
    color <- if ( is.null(group.cols) || ! exists(name, group.cols) )
      SetColor(dt[[name]], tag = i)
    else 
      group.cols[[name]]
    p1 <- p + geom_bar(aes_string(fill = name), stat = "identity") +
      scale_fill_manual(values = color)
    ggsave(p1, file = paste0(outpref, ".", name,".pdf"), width = 7, height = 7)
  }
}

BarPlot.Q3 <- function(count, color = "#3E70C0", outfile = "QC_Q3_barplot.pdf") {
  qt <- sapply(count, quantile, probs = 0.75, names = F)
  p <- ggplot(as.data.frame(qt), aes(x = log2(qt))) +
    geom_density(fill = color, alpha = 0.5) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.5,
                   fill = color, color = "black") +
    geom_density(alpha = 0.5) +
    xlim(0, 15) +
    xlab("log2 Q3(Counts)") +
    theme_minimal()
  ggsave(p, file = outfile, width = 7, height = 7)
}

HeatmapPlot <- function(data, meta, var = "SegmentDisplayName", 
                        group.by = "SegmentLabel", group.cols = NULL,
                        is.scale = TRUE, max.value = 2, 
                        color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                        show_rownames = FALSE, outpref = "heatmap", plot.each.labels = TRUE){
  color <- colorRampPalette(color)(100)
  annot <- tibble::column_to_rownames(meta, var = var)[, group.by, drop = FALSE]
  annot_color <- sapply(seq(group.by), function(i) {
	name <- group.by[[i]]
    if ( is.null(group.cols) || ! exists(name, group.cols) ) {
      SetColor(annot[[name]], tag = i)
    } else {
      group.cols[[name]]
    }}, simplify = FALSE)
  names(annot_color) <- group.by
  filename <- ifelse(is.null(outpref), NA, paste0(outpref, "_overall.pdf"))
  
  if ( is.scale ) {
    mat <- t(scale(t(data)))
  }
  if ( ! is.null(max.value) ) {
    mat[mat > max.value] <- max.value
    mat[mat < -max.value] <- -max.value
  }
  h = 7
  w = max(4.5, 1.3 + ncol(mat) * 0.3)
  pheatmap::pheatmap(mat, show_rownames = show_rownames, color = color,
                     annotation_col = annot,
                     annotation_colors = annot_color,
                     filename = filename, width = w, height = h)
  if ( plot.each.labels ) {
    for ( i in group.by) {
      mat2 = mat[,order(annot[,i])]
      h = 7
      w = max(4.5, 1.3 + ncol(mat2) * 0.3)
      filename <- ifelse(is.null(outpref), NA, paste0(outpref, "_label.", i, ".pdf"))
      pheatmap::pheatmap(mat2, cluster_cols = FALSE,
                         show_rownames = show_rownames, color = color, 
                         annotation_col = annot[,i,drop = FALSE],
                         annotation_colors = annot_color[i],
                         filename = filename)
    }
  }
}

CorrelationPlot <- function(data, color = rev(RColorBrewer::brewer.pal(11, "RdBu")), outfile = "CorrelationPlot.pdf") {
  cor <- cor(data, data)
  dt <- reshape2::melt(cor)
  
  ft <- ggplot2::.pt*72.27/96
  size <- ncol(data) * 1.3 + 6
  p <- ggplot(dt, aes(x = Var1, y = Var2)) +
    geom_point(aes(color = value, size = abs(value))) +
    scale_color_gradientn(NULL, colours = rev(color), limits = c(-1,1))+
    scale_size(range = c(0, size), guide = NULL) +
    geom_tile(fill = NA, color = "black")+
    scale_y_discrete(NULL, limits = rev)+
    scale_x_discrete(NULL, position = "top")+
    coord_equal()+
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "red"),
          axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
          legend.key.height = grid::unit((size+1)/ggplot2::.pt * ncol(cor), "bigpts"))
  
  source("ggplot_save.R")
  SavePlot(p, plot.size = (size + 1) * ft * ncol(cor), units = "bigpts", file = outfile, is.legend.overlap = F)
}

PCAPlot <- function(data, meta,
                    var = "SegmentDisplayName", var.color = "royalblue",
                    group.by = "SegmentLabel", group.cols = NULL,
                    nfeature = 15, outpref = "PCA") {
  annot <- tibble::column_to_rownames(meta, var = var)[, group.by, drop = F]
  npcs <- 10
  pca.results <- stats::prcomp(x = t(data), rank. = npcs)
  feature.loadings <- as.data.frame(pca.results$rotation)
  feature.loadings$feature <- rownames(feature.loadings)
  cell.embeddings <- as.data.frame(pca.results$x)
  cell.embeddings <- cbind(cell.embeddings, annot)
  
  p <- ggplot(head(feature.loadings, nfeature), aes(x = PC1, y = PC2)) +
    geom_segment(aes(group = feature, xend = 0, yend = 0), arrow = arrow(ends = "first")) +
    ggrepel::geom_text_repel(aes(label = feature), color = var.color, min.segment.length = 9999) +
    theme_bw()
  ggsave(p, file = paste0(outpref, "_arrow_", nfeature, "_feature.pdf"), width = 7, height = 7)
  
  for ( i in seq(group.by)) {
	name <- group.by[[i]]
    color <- if ( is.null(group.cols) || ! exists(name, group.cols) ) {
      SetColor(annot[[name]], tag = i)
    } else {
      group.cols[[name]]
    }    
    p1 <- ggplot(cell.embeddings, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(color = name))+
      scale_color_manual(values = color) +
      geom_hline(yintercept = 0,linetype = "dashed") +
      geom_vline(xintercept = 0,linetype = "dashed") +
      theme_bw()
    ggsave(p1, file = paste0(outpref, ".", name, ".pdf"), width = 7, height = 7)
  }
}

ComparePlot <- function(dt, genes, compare.by = "SegmentLabel", test.method = NULL,
                        color = NULL, group.by = "Classification", 
                        ncol = NULL, nrow = NULL, plot.size = 2.25,
                        outpref = "Signature.test") {
  if ( is.null(test.method) ) test.method <- if (length(unique(dt2[[compare.by]])) > 2) "kruskal.test" else "wilcox.test"
  result <- compare_means(formula = as.formula(paste0("value ~ ", compare.by)),
                          data = dt2, method = test.method,
                          group.by = group.by)
  result <- result[, c(group.by, "p", "p.adj", "p.signif", "method")]
  write.table(result, file = paste0(outpref, ".xls"), quote = F, row.names = F, sep = "\t")
  
  num <- length(unique(dt2[[group.by]]))
  if ( ! is.null(ncol) ) {
    nrow <- ceiling(num / ncol)
  } else if ( ! is.null(nrow) ) {
    ncol <- ceiling(num / nrow)
  } else {
    ncol <- ceiling(sqrt(num))
    nrow <- ceiling(num / ncol)
  }
  dt2[[group.by]] <- stringr::str_wrap(dt2[[group.by]], 20)
  p <- ggboxplot(dt2, x = compare.by, y = "value", fill = compare.by, palette = color) + 
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method = test.method) +
    facet_wrap(group.by, scales = "free", strip.position = "left", nrow = nrow, ncol = ncol) + 
    labs(x = NULL, y = NULL) + 
    coord_cartesian(clip = "off") + 
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
	  axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )
  width <- max(7, plot.size * ncol * max(1, log(length(unique(dt2[[compare.by]])))/log(3)) )
  height <- max(7, plot.size * nrow)
  ggsave(p, file = paste0(outpref, ".pdf"), width = width, height = height, limitsize = FALSE)
}

SetColor <- function(x, type = "tsne", tag = 1L, ...) {
  if ( file.exists("Colors.R") ) {
    source("Colors.R", chdir = TRUE)
  }
 
  if ( ! is.factor(x) ) x <- as.factor(x)
  if ( nlevels(x) == 2 && tag == 1 ) {
    color <- c("royalblue4", "red3")
  } else {
	if ( is.numeric(tag) ) {
	  tag <- max(tag, length(color.list[[type]]))
	  tag <- names(color.list[[type]])[tag]
    }
    color <- if ( exists("fetch_color") ) fetch_color(n = nlevels(x), type = type, tag = tag, ...)
    else rainbow(n = nlevels(x))
  }
  names(color) <- levels(x)
  return(color)
}
