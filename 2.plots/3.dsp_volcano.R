
args = commandArgs(TRUE)

infile = args[1]

library(dplyr)
library(ggplot2)
library(data.table)

Pvalue = 0.05
Log2FC = 1

xlabs = "log2FoldChange"
ylabs = "-log10(p.adj)"
colorlabs = xlabs
sizelabs = ylabs

df = infile %>%
  fread(sep = "\t", header = TRUE, data.table = FALSE, stringsAsFactors = FALSE) %>%
  mutate(group = ifelse(`log2(fc)` > 0, yes = "up", no = "down") %>% factor()) %>%
  mutate(log2fc = `log2(fc)`) %>%
  mutate(logFDR = -log10(FDR))

top10 = df %>%
  filter(FDR < Pvalue) %>%
  group_by(group) %>%
  top_n(10, wt = abs(log2fc))

df$label = NA
df$label[df$id %in% top10$id] = df$id[df$id %in% top10$id]

max.logFDR = ceiling(max(df$logFDR[is.finite(df$logFDR)]) / 5) * 5
max.logFDR = min(ceiling(median(top10$logFDR)), max.logFDR)
b = ceiling((max.logFDR / 5) / 10) * 10

df$logFDR[df$logFDR > max.logFDR] = max.logFDR
df$log2fc[df$log2fc > 10] = 10
df$log2fc[df$log2fc < -10] = -10

colors = RColorBrewer::brewer.pal(name = "Spectral", n = 9) %>% rev()

p = ggplot(df, aes(x = log2fc, y = logFDR)) +
  geom_point(aes(color = log2fc, size = logFDR)) +
  scale_color_gradientn(colors = colors) +
  scale_size_continuous(limits = c(0, max.logFDR), breaks = seq(0, max.logFDR, b)) +
  geom_hline(yintercept = -log10(Pvalue), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1) * Log2FC, linetype = "dashed", color = "grey40") +
  xlim(-10, 10) + ylim(0, max.logFDR) +
  labs(x = xlabs, y = ylabs, color = colorlabs, size = sizelabs) +
  theme_light() +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = 20)

w = 7
h = 5.5
outfile = gsub("\\.xls", ".Volcano.pdf", infile)
ggsave(p, filename = outfile, width = w, height = h)

