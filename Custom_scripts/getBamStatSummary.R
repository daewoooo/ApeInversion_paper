## Load required libraries ##
#############################
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(Biostrings) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38.masked) )

## Get genome size ##
#####################
#example.bam <- "/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/GG09x02PE20396.sort.mdup.bam"
#file.header <- Rsamtools::scanBamHeader(example.bam)[[1]]
#chroms <- file.header$targets
genome <- BSgenome.Hsapiens.UCSC.hg38.masked
genomeSize <- 0
for (chr in paste0('chr', c(1:22, 'X'))) {
  chr.len <- length(genome[[chr]])
  gap.len <- maskedwidth(masks(genome[[chr]])[1])
  genomeSize <- genomeSize + (chr.len - gap.len)
}

## Parameters ##
################
read.len <- 80
outputfolder <- "/home/porubsky/WORK/Great_apes/Data_stats"
if (!dir.exists(outputfolder)) {
  dir.create(outputfolder)
}

## Load required data ##
########################
## Chimpanzee
chimpanzee.bams <- list.files("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/", pattern = "\\.bam$", full.names = TRUE)
message("Calculating BAM stats for chimpanzee ...")
chimpanzee.bams.stat <- list()
for (i in seq_along(chimpanzee.bams)) {
  bam <- chimpanzee.bams[i]
  bam.stat <- bam2stat(bamfile = bam, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
  chimpanzee.bams.stat[[i]] <- bam.stat
}
chimpanzee.bams.stat <- do.call(rbind, chimpanzee.bams.stat)
chimpanzee.bams.stat$depth <- (chimpanzee.bams.stat$total.reads * read.len) / genomeSize
chimpanzee.bams.stat$covered.pos.perc <- (chimpanzee.bams.stat$covered.pos / genomeSize) * 100

#Load merged bam
merged.bam <- "/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/Dorien_GRCh38_merged/Dorien_GRCh38_selected.bam"
bam.stat <- bam2stat(bamfile = merged.bam, chunkSize = 10000, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
bam.stat$depth <- (bam.stat$total.reads * read.len) / genomeSize
bam.stat$covered.pos.perc <- (bam.stat$covered.pos / genomeSize) * 100
bam.stat$filename <- 'total'
#Export results
final.df <- rbind(chimpanzee.bams.stat, bam.stat)
final.df$ID <- 'chimpanzee'
destination <- file.path(outputfolder, 'chimpanzee_selected_bamStat.txt')
write.table(final.df, file = destination, quote = FALSE, row.names = FALSE)


## Bonobo
bonobo.bams <- list.files("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/", pattern = "\\.bam$", full.names = TRUE)
message("Calculating BAM stats for bonobo ...")
bonobo.bams.stat <- list()
for (i in seq_along(bonobo.bams)) {
  bam <- bonobo.bams[i]
  bam.stat <- bam2stat(bamfile = bam, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
  bonobo.bams.stat[[i]] <- bam.stat
}
bonobo.bams.stat <- do.call(rbind, bonobo.bams.stat)
bonobo.bams.stat$depth <- (bonobo.bams.stat$total.reads * read.len) / genomeSize
bonobo.bams.stat$covered.pos.perc <- (bonobo.bams.stat$covered.pos / genomeSize) * 100

#Load merged bam
merged.bam <- "/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/Bonobo_GRCh38_merged/bonobo_GRCh38_selected_merged.bam"
bam.stat <- bam2stat(bamfile = merged.bam, chunkSize = 10000, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
bam.stat$depth <- (bam.stat$total.reads * read.len) / genomeSize
bam.stat$covered.pos.perc <- (bam.stat$covered.pos / genomeSize) * 100
bam.stat$filename <- 'total'
#Export results
final.df <- rbind(bonobo.bams.stat, bam.stat)
final.df$ID <- 'bonobo'
destination <- file.path(outputfolder, 'bonobo_selected_bamStat.txt')
write.table(final.df, file = destination, quote = FALSE, row.names = FALSE)

## Gorilla
gorilla.bams <- list.files("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/", pattern = "\\.bam$", full.names = TRUE)
message("Calculating BAM stats for gorilla ...")
gorilla.bams.stat <- list()
for (i in seq_along(gorilla.bams)) {
  bam <- gorilla.bams[i]
  bam.stat <- bam2stat(bamfile = bam, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
  gorilla.bams.stat[[i]] <- bam.stat
}
gorilla.bams.stat <- do.call(rbind, gorilla.bams.stat)
gorilla.bams.stat$depth <- (gorilla.bams.stat$total.reads * read.len) / genomeSize
gorilla.bams.stat$covered.pos.perc <- (gorilla.bams.stat$covered.pos / genomeSize) * 100

#Load merged bam
merged.bam <- "/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/Gorilla_GRCh38_merged/gorilla_GRCh38_selected_merged.bam"
bam.stat <- bam2stat(bamfile = merged.bam, chunkSize = 10000, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
bam.stat$depth <- (bam.stat$total.reads * read.len) / genomeSize
bam.stat$covered.pos.perc <- (bam.stat$covered.pos / genomeSize) * 100
bam.stat$filename <- 'total'
#Export results
final.df <- rbind(gorilla.bams.stat, bam.stat)
final.df$ID <- 'gorilla'
destination <- file.path(outputfolder, 'gorilla_selected_bamStat.txt')
write.table(final.df, file = destination, quote = FALSE, row.names = FALSE)

## Orangutan
orangutan.bams <- list.files("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/", pattern = "\\.bam$", full.names = TRUE)
message("Calculating BAM stats for orangutan ...")
orangutan.bams.stat <- list()
for (i in seq_along(orangutan.bams)) {
  bam <- orangutan.bams[i]
  bam.stat <- bam2stat(bamfile = bam, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
  orangutan.bams.stat[[i]] <- bam.stat
}
orangutan.bams.stat <- do.call(rbind, orangutan.bams.stat)
orangutan.bams.stat$depth <- (orangutan.bams.stat$total.reads * read.len) / genomeSize
orangutan.bams.stat$covered.pos.perc <- (orangutan.bams.stat$covered.pos / genomeSize) * 100

#Load merged bam
merged.bam <- "/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/Orangutan_GRCh38_merged/orangutan_GRCh38_selected_merged.bam"
bam.stat <- bam2stat(bamfile = merged.bam, chunkSize = 10000, min.mapq = 10, chromosomes = paste0('chr', c(1:22, 'X')), filt.alt = FALSE, filt.flag = 3328)
bam.stat$depth <- (bam.stat$total.reads * read.len) / genomeSize
bam.stat$covered.pos.perc <- (bam.stat$covered.pos / genomeSize) * 100
bam.stat$filename <- 'total'
#Export results
final.df <- rbind(orangutan.bams.stat, bam.stat)
final.df$ID <- 'orangutan'
destination <- file.path(outputfolder, 'orangutan_selected_bamStat.txt')
write.table(final.df, file = destination, quote = FALSE, row.names = FALSE)


## Load BAM stat data ##
chimpanzee.stat <- read.table(file.path(outputfolder, "chimpanzee_selected_bamStat.txt"), header=TRUE)
bonobo.stat <- read.table(file.path(outputfolder, "bonobo_selected_bamStat.txt"), header=TRUE)
gorilla.stat <- read.table(file.path(outputfolder, "gorilla_selected_bamStat.txt"), header=TRUE)
orangutan.stat <- read.table(file.path(outputfolder, "orangutan_selected_bamStat.txt"), header=TRUE)

total.counts <- rbind(chimpanzee.stat[nrow(chimpanzee.stat),],
                      bonobo.stat[nrow(bonobo.stat),],
                      gorilla.stat[nrow(gorilla.stat),],
                      orangutan.stat[nrow(orangutan.stat),])

single.cell.counts <- rbind(chimpanzee.stat[-nrow(chimpanzee.stat),],
                            bonobo.stat[-nrow(bonobo.stat),],
                            gorilla.stat[-nrow(gorilla.stat),],
                            orangutan.stat[-nrow(orangutan.stat),])

library.counts.perIndivid <- table(single.cell.counts$ID)

## Plot BAM stats
message("Preparing plots ...")
suppressPackageStartupMessages( library(reshape2) )
plt.df <- total.counts[c('total.reads', 'depth', 'covered.pos.perc', 'ID')]
colnames(plt.df) <- c('Total # of reads', 'Depth of coverage', '% of genome covered', 'ID')
plt.df$ID <- paste0(plt.df$ID, "\n(n=", library.counts.perIndivid, ")")

plt.df <- melt(plt.df, measure.vars = c('Total # of reads', 'Depth of coverage', '% of genome covered'), id.vars = 'ID')
plt <- ggplot(plt.df) +
  geom_col(aes(x=ID, y=value, fill=ID)) +
  scale_fill_manual(values =  c('#3182bd','#31a354','#8856a7','#e6550d'), name="") +
  scale_y_continuous(labels = comma) +
  facet_grid(variable ~ ., scales = 'free') +
  xlab("") + ylab("") +
  theme(legend.position = "none")
## Save plot
destination <- file.path(outputfolder, "total_stat.pdf")
ggsave(plt, filename = destination, width = 5, height = 6)  

## Plot total coverage per library
my_theme <- theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())

p1 <- ggplot(single.cell.counts) +
  geom_boxplot(aes(x=ID, y=total.reads, fill=ID)) +
  scale_fill_manual(values =  c('#3182bd','#31a354','#8856a7','#e6550d'), guide="none") +
  ylab("Total # of reads") + xlab("")
p2 <- ggplot(single.cell.counts) +
  geom_boxplot(aes(x=ID, y=depth, fill=ID)) +
  scale_fill_manual(values =  c('#3182bd','#31a354','#8856a7','#e6550d'), guide="none") +
  ylab("Depth of coverage") + xlab("")
p3 <- ggplot(single.cell.counts) +
  geom_boxplot(aes(x=ID, y=covered.pos.perc, fill=ID)) +
  scale_fill_manual(values =  c('#3182bd','#31a354','#8856a7','#e6550d'), guide="none") +
  ylab("% of genome covered") + xlab("")

plt <- plot_grid(p1, p2, p3, nrow = 1)
## Save plot
destination <- file.path(outputfolder, "perCell_stat.pdf")
ggsave(plt, filename = destination, width = 10, height = 5)  

message("DONE!!!")