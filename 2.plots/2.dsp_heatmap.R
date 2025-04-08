
getopt <- function(){
		library(optparse)
		bin <- dirname(normalizePath(getopt::get_Rscript_filename()))
		option_list <- list(
				make_option(c("-i", "--infile"), type = "character",
						help = "pristine xlsx data."),
				make_option(c("-c", "--count"), type = "character",
						help = "TargetCountMatrix. overwrite which in --infile." ),
				make_option(c("-m", "--metadata"), type = "character", 
						help = "SegmentProperties. overwrite which in --infile." ),
				make_option(c("-f", "--fdata"), type = "character", 
						help = "TargetProperties. overwrite which in --infile." ),
				make_option(c("-g", "--groups"), type = "character",
						help = "Addtional groups info."), 
				make_option(c("-o", "--outdir"), type = "character", default = "./",
						help = "outdir. default '%default'."),
				make_option(c("-a", "--add_lib"), type = "character", default = paste0(bin, "/DSP_lib.R"),
						help = "additional lib script. default '%default'." )
#				make_option(c("-s", "--signiture"), type = "character", default = paste0(bin, "/../database/signiture.annot.yaml"),
#						help = "signiture annot info. default '%default'.")
			)
		opts.obj <- OptionParser(option_list = option_list, description = "e.g.
    Rscipt %prog -i <xlsx> (-c <count> -m <metadata> -f <fdata>)
        OR
    Rscipt %prog -c <count> -m <metadata> -f <fdata>")
		opts <- parse_args(opts.obj)

		if ( is.null(opts$infile) && ( is.null(opts$count) || is.null(opts$metadata) || is.null(opts$fdata) ) ) {
				print_help(opts.obj)
				q()
		}
		return(opts)
}

opts <- getopt()


##################
### start here ###
library(dplyr)
library(ggplot2)

if ( ! is.null(opts$add_lib) ) {
		source(opts$add_lib, chdir = T)
}

if ( ! is.null(opts$infile) ) {
		library(readxl)
		infile <- opts$infile
		count  <- read_xlsx(infile, sheet = "TargetCountMatrix") %>% tibble::column_to_rownames(var = "TargetName")
		meta   <- read_xlsx(infile, sheet = "SegmentProperties") 
		fdata  <- read_xlsx(infile, sheet = "TargetProperties")		
}

if ( ! is.null(opts$count) ) {
		count <- read.table(opts$count, sep = "\t", quote = "", header = T, row.names = 1, check.names = F)
}
if ( ! is.null(opts$metadata) ) {
		meta  <- read.table(opts$metadata, sep = "\t", quote = "", header = T, check.names = F)
}
if ( ! is.null(opts$fdata) ) {
		fdata <- read.table(opts$fdata, sep = "\t", quote = "", header = T, check.names = F)
}

meta$SegmentDisplayName <- gsub("[ |,]+", "_", meta$SegmentDisplayName)
colnames(count) <- gsub("[ |,]+", "_", colnames(count))

group_label <- "SegmentLabel"
#group_label <- c("SegmentLabel", "SlideName")
if ( ! is.null(opts$groups) ) {
		group_info <- read.table(opts$groups, sep = "\t", quote = "", header = T, check.names = F, colClasses = 'character')
		colnames(group_info)[1] <- "SegmentDisplayName"
		group_label <- colnames(group_info)[-1]
		meta <- left_join(group_info, meta, by = "SegmentDisplayName", suffix = c('', '.y'))
}

if ( ! is.null(opts$outdir) ) {
		setwd(opts$outdir)
}

### 1.Quality_Control
#DotPlot(meta)

## Neg
Neg.gene <- if ( "Negative Probe" %in% rownames(count) ) "Negative Probe" else "NegProbe-WTX"
#BarPlot(count, meta, feature = Neg.gene, yintercept = 10, group.by = group_label)

## HK
HK.gene <- subset(fdata, CodeClass == "Control")[["TargetName"]]
#if ( length(HK.gene) > 0 ) {
#	BarPlot(count, meta, feature = HK.gene, out.id = "HK", scale_y = "log2", group.by = group_label)
#}

## Q3
#BarPlot.Q3(count)


### 2.Overall_Analysis
data  <- Q3Normalize(count)
## heatmap
HeatmapPlot(data, meta, group.by = group_label)

## cor
CorrelationPlot(data)

## pca
#PCAPlot(data, meta, group.by = group_label)


