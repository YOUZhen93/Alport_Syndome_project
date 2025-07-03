# Alport Syndrome Figures

# shell command to generate input file:
#
# bioawk -c fastx '{print length($seq)"\t"meanqual($qual)}' ${fastq}
library(optparse)

Ranno_list <- list(
    make_option(c("-i", "--input"), type="character", default=FALSE,
        help="input stats file columned by length and average qs, delimited by tab"),
    make_option(c("-o", "--outfile"), type="character", default="Nanopore_readsLen_QS",
        help = "output file prefix")
    )

opt <- parse_args(OptionParser(option_list = Ranno_list))

if (is.na(c(opt$input))) {
  stop("Error: QC files are required! ")
}

options(stringsAsFactors=F)

suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggExtra))

input = opt$input
output = opt$outfile
out1 = paste0(output, "_readsLen_QS_2d.pdf")
out2 = paste0(output, "_QC_stats_metrics.txt")
gradient_set = c('#68ACF0', '#6DA7EC', '#71A2E7', '#769CE3', '#7A97DF', '#7F92DB', '#848DD6', '#8888D2', '#8D83CE', '#927DCA', '#9678C5', '#9B73C1', '#9F6EBD', '#A469B8', '#A964B4', '#AD5EB0', '#B259AC', '#B754A7', '#BB4FA3', '#C04A9F', '#C4459A', '#C93F96', '#CE3A92', '#D2358E', '#D73089', '#DC2B85', '#E02681', '#E5207D', '#E91B78', '#EE1674')

N50_cal <- function(num) {
l1 = rev(sort(num))
l2 = cumsum(l1) <= sum(l1)/2
print(paste("number needs to reach N50:", sum(l2)))
N50 = l1[sum(l2)]
return(N50)
}


x = fread(input, header=F)
x = x[which(x$V1 > 100 & x$V2 > 7), ]
tot_obs = nrow(x)
tot_bases = sum(x$V1)
yield = round(tot_bases/1000/1000/1000, 3)
max_len = max(x$V1)
mean_len = round(mean(x$V1),3)
median_len = median(x$V1)
N50 = N50_cal(as.numeric(x$V1))
mean_q = round(mean(x$V2),3)
median_q = median(x$V2)

stats1 = data.frame(V1 = c("sample", "total_reads", "total_bases", "yield (GB)", "max_length", "mean_length", "median_length", "N50", "mean_quality", "median_quality"), V2 = c(output, tot_obs, tot_bases, yield, max_len, mean_len, median_len, N50, mean_q, median_q))
write.table(stats1, out2, sep="\t", quote=F, col.names=F, row.names=F)

N50=format(N50, nsmall=0, big.mark=",")
x$V1 = log10(x$V1)
hists = hist(x$V1, breaks = 40)
max_his = max(hists$density)
coef = max(x$V2) / max_his
pdf(out1, w = 10, h = 7)
ggplot(x, aes(x=V1, y = V2)) +
geom_bin2d(aes(x = V1, , y = V2), bins = 200) +
geom_histogram(aes(x = V1, y=after_stat(coef*density)), alpha=0, color="black", position="identity", fill = NA, bins = 40) +
scale_fill_gradientn(colours = gradient_set, values = scales::rescale(c(-0.1, 0, 0.05, 0.1, 0.5, 1))) +
scale_y_continuous(sec.axis = sec_axis(~. / coef, name = "Density")) +
scale_x_continuous(limits = c(2,6), breaks = c(2,3,4,5,6), labels=c("0.1K","1K","10K","100K", "1,000K")) + theme_bw() +
theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm")) +
geom_vline(aes(xintercept = 2), colour="black", size=0.7) + geom_hline(aes(yintercept = 7), colour="black", size=0.7) +
annotate('text', x=2, label="\nLen=100", y=40, colour="black", angle=90, size=8) +
annotate('text', x=5, label=paste0("N50 = ",N50, "\n"), y=40, colour="black", size=8) +
annotate('text', x=4, label="QS=7\n", y=7, colour="black", size=8) + xlab("Log10(length)") +
ylab("Average QS") + theme(axis.title.y = element_text(size=28,color="black", margin = margin(t = 0, r = 2, b = 0, l = 0, unit="mm")),
    axis.title.x = element_text(size=28,color="black", margin = margin(t = 2, r = 0, b = 0, l = 0, unit="mm")),
    axis.title.y.right = element_text(size=28,color="black", margin = margin(t = 0, r = 0, b = 0, l = 2, unit="mm")),
    axis.text = element_text(size=20,color="black"),  axis.ticks = element_line(size=1)) +
theme(axis.ticks.length=unit(0.2, "cm")) + theme(panel.background = element_rect(fill = "white", colour = "white", size = 2),
    panel.border = element_rect(size = 1, color="black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(legend.key.size = unit(0.7, 'cm'), legend.title = element_text(size=22), legend.text = element_text(size=18))
dev.off()

