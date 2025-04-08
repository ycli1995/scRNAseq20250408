
library(dplyr)
library(data.table)
library(ggplot2)

top.n = 20

GO = "GO.bar_Gradient.xls" %>%
  fread(sep = "\t", stringsAsFactors = F, data.table = F, header = T) %>%
  mutate(type = class)

KO = "KO.bar_Gradient.xls" %>%
  fread(sep = "\t", stringsAsFactors = F, data.table = F, header = T) %>%
  mutate(type = "KEGG")

df = rbind(GO, KO) %>%
  mutate(logFDR = -log10(Pvalue))

FDR.colors = RColorBrewer::brewer.pal(name = "YlOrRd", n = 9) %>% head(8)
type.colors = RColorBrewer::brewer.pal(name = "Set1", n = 4) %>% setNames(unique(df$type))

# KEGG
i = "KEGG"
i2 = gsub("\\s|\\/", "_", i)
df2 = df %>%
  filter(type == i) %>%
  top_n(top.n, wt = logFDR) %>%
  mutate(Descrption = factor(Descrption, rev(unique(Descrption))))
p1 = ggplot(df2, aes(y = ratio, x = Descrption)) +
  geom_bar(aes(fill = logFDR), stat = "identity") +
  scale_fill_gradientn(colors = FDR.colors) +
  geom_text(aes(y = 0, label = Descrption), hjust = 0) +
  coord_flip() +
  labs(y = "Gene ratio", x = "Description", fill = "-log10(Pvalue)", title = paste(i, "enrichment")) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))
h = 0.5 + 0.26 * nrow(df2)
w = 2 + 0.08 * max(nchar(levels(df2$Descrption)))
outfile = paste0("Top", top.n, ".BarPlot.", i2, ".pdf")
ggsave(p1, filename = outfile, height = h, width = w)

# GO
top.n = 5
i2 = "GO"
df2 = df %>%
  filter(!type == i) %>%
  group_by(type) %>%
  top_n(top.n, wt = logFDR) %>%
  mutate(Descrption = factor(Descrption, rev(unique(Descrption))))

max.logFDR = ceiling(max(df2$logFDR))
step.logFDR = round(max.logFDR / 5, digits = 2)
break.logFDR = seq(0, max.logFDR, step.logFDR)

max.count = ceiling(max(df2$ratio * 100))
step.count = round(max.count / 5, digits = 2)
break.count = seq(0, max.count, step.count)
fold = max.count / max.logFDR

p1 = ggplot(df2) +
  ggforce::geom_link(
    aes(x = 0, y = Descrption, xend = logFDR, yend = Descrption, color = type, linewidth = after_stat(index)),
    n = 500,
    show.legend = F
  ) +
  geom_point(aes(x = logFDR, y = Descrption), color = "black", fill = "white", size = 6, shape = 21) +
  geom_line(aes(x = ratio * 100 / fold, y = Descrption, group = 1), orientation = "y", linewidth = 1, color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./fold, labels = scales::label_percent(), name = "Percent of geneRatio")) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black", size = 12), plot.title = element_text(hjust = 0.5)) +
  labs(x = "-log10(Pvalue)", y = "", title = paste(i2, "enrichment")) +
  scale_color_manual(values = type.colors)
h = 1 + 0.3 * nrow(df2)
w = 3 + 0.07 * max(nchar(levels(df2$Descrption)))
outfile = paste0("Top", top.n, ".MeteorPlot.", i2, ".pdf")
ggsave(p1, filename = outfile, height = h, width = w)

