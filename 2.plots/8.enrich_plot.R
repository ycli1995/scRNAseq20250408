
library(dplyr)
library(data.table)
library(ggplot2)

top.n = 20

KO = "KO.bar_Gradient.xls" %>%
  fread(sep = "\t", stringsAsFactors = F, data.table = F, header = T) %>%
  mutate(type = "KEGG")

df = KO %>%
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
top.n = nrow(df2)
p1 = ggplot(df2, aes(y = ratio, x = Descrption)) +
  geom_bar(aes(fill = logFDR), stat = "identity") +
  scale_fill_gradientn(colors = FDR.colors) +
  geom_text(aes(y = 0, label = Descrption), hjust = 0) +
  coord_flip() +
  labs(y = "Gene ratio", x = "Description", fill = "-log10(Pvalue)", title = paste(i, "enrichment")) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))
h = max(0.5 + 0.26 * nrow(df2), 3)
w = 2 + 0.08 * max(nchar(levels(df2$Descrption)))
outfile = paste0("Top", top.n, ".BarPlot.", i2, ".pdf")
ggsave(p1, filename = outfile, height = h, width = w)

