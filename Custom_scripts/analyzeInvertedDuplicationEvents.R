## Load required libraries
library(primatR)
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(UpSetR)
library(VennDiagram)
library(BSgenome.Hsapiens.UCSC.hg38)

## Set output directory
outputdirectory <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/"

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load centromere Track
cent <- read.table("/home/porubsky/WORK/Great_apes/Annotations/centromeres_GRCh38.bed.gz")
cent.gr <- GRanges(seqnames=cent$V2, ranges=IRanges(start=cent$V3, end=cent$V4))

## Blacklisted regions ##
#########################
## +/- 500kb from the centromeres
#blacklist1 <- GRanges(seqnames=cent$V2, ranges=IRanges(start=cent$V3-500000, end=cent$V4+500000))
## 5Mb from each chromosome end
#chrom.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, 'X'))]
#blacklist2 <- GRanges(seqnames=names(chrom.len), ranges=IRanges(start = 1, end = 5000000))
#blacklist3 <- GRanges(seqnames=names(chrom.len), ranges=IRanges(start = chrom.len - 5000000, end = chrom.len))
## Merge all blacklists together
#blacklist <- reduce(c(blacklist1, blacklist2, blacklist3))
## Blacklist regions that are likely interchromosomal links
links.folder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/"
## Get gorilla blacklisted regions
gorilla.links <- get(load(file.path(links.folder, 'gorilla/gorilla_GRCh38align_links.RData'))) 
gorilla.links <- gorilla.links$inter.links
ranges1 <- gorilla.links[,0]
ranges2 <- gorilla.links$to.gr
blacklist.gorilla <- reduce(c(ranges1, ranges2))
## Get orangutan blacklisted regions
orangutan.links <- get(load(file.path(links.folder, 'orangutan/orangutan_GRCh38align_links.RData'))) 
orangutan.links <- orangutan.links$inter.links
ranges1 <- orangutan.links[,0]
ranges2 <- orangutan.links$to.gr
blacklist.orangutan <- reduce(c(ranges1, ranges2))
## Get chimpanzee blacklisted regions
chimpanzee.links <- get(load(file.path(links.folder, 'chimpanzee/chimpanzee_GRCh38align_links.RData'))) 
chimpanzee.links <- chimpanzee.links$inter.links
ranges1 <- chimpanzee.links[,0]
ranges2 <- chimpanzee.links$to.gr
blacklist.chimpanzee <- reduce(c(ranges1, ranges2))
## Get bonobo blacklisted regions
bonobo.links <- get(load(file.path(links.folder, 'bonobo/bonobo_GRCh38align_links.RData'))) 
bonobo.links <- bonobo.links$inter.links
ranges1 <- bonobo.links[,0]
ranges2 <- bonobo.links$to.gr
blacklist.bonobo <- reduce(c(ranges1, ranges2))

## Read in all called inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Filter only inverted duplications invDups
all.invDups.gr <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Load manually selected duplicated regions with direct orientation
directDups <- read.table(file.path(outputdirectory, "directDups_list.txt"), stringsAsFactors = FALSE, header = FALSE)
directDups.gr <- GRanges(seqnames=directDups$V1, ranges=IRanges(start=directDups$V2, end=directDups$V3), ID=directDups$V4)
## Keep only events >=10kb
directDups.gr <- directDups.gr[width(directDups.gr) >= 10000]

## Plot invDup counts and size distribution for inverted dups ##
################################################################
all.invDups.df <- as.data.frame(all.invDups.gr)
colors = c('#3182bd','#31a354','#8856a7','#e6550d')
plt1 <- all.invDups.df %>% group_by(ID) %>% summarise(count=n()) %>% 
        ggplot() + geom_col(aes(x=ID, y=count, fill=ID)) +
        scale_fill_manual(values = colors, guide='none') +
        xlab("") + 
        geom_text(aes(x=ID, y=count, label=count), vjust=0)+
        scale_x_discrete(labels = c('B', 'C', 'G', 'O'))
plt2 <- all.invDups.df %>% group_by(ID) %>% summarise(invBases=sum(width)) %>% 
        ggplot() + geom_col(aes(x=ID, y=invBases, fill=ID)) +
        scale_fill_manual(values = colors, guide='none') +
        xlab("") +
        ylab("Total inverted bases (bp)") +
        scale_y_continuous(labels = comma) +
        scale_x_discrete(labels = c('B', 'C', 'G', 'O'))
  
final.plt <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(0.8, 1))
destination <- file.path(outputdirectory, "invDupCountsPerIndivid.pdf")
ggsave(final.plt, filename = destination, width = 5, height = 5, useDingbats=FALSE)
## Plot invDup size distribution
#destination <- file.path(outputdirectory, "invDupSizeDist.pdf")
#invDup.sizeDist <- rangesSizeDistribution(all.invDups.gr)
#ggsave(invDup.sizeDist, filename = destination, device = 'pdf', width = 10, height = 10)

## Plot proportion of direct vs inverted duplications ##
########################################################
## Remove calls from the blacklisted regions
#all.invDups.gr <- subsetByOverlaps(all.invDups.gr, blacklist, invert = TRUE)
#directDups.gr <- subsetByOverlaps(directDups.gr, blacklist, invert = TRUE)

## Remove previously obtained blacklisted region for each individual based on inter-chromosomal link analysis
gorilla.invDups.gr <- subsetByOverlaps(all.invDups.gr[all.invDups.gr$ID == 'gorilla'], blacklist.gorilla, invert = TRUE)
orangutan.invDups.gr <- subsetByOverlaps(all.invDups.gr[all.invDups.gr$ID == 'orangutan'], blacklist.orangutan, invert = TRUE)
chimpanzee.invDups.gr <- subsetByOverlaps(all.invDups.gr[all.invDups.gr$ID == 'chimpanzee'], blacklist.chimpanzee, invert = TRUE)
bonobo.invDups.gr <- subsetByOverlaps(all.invDups.gr[all.invDups.gr$ID == 'bonobo'], blacklist.bonobo, invert = TRUE)
all.invDups.filt.gr <- sort(c(gorilla.invDups.gr, orangutan.invDups.gr, chimpanzee.invDups.gr, bonobo.invDups.gr))

invDups.counts <- as.data.frame(table(all.invDups.filt.gr$ID))
dirDups.counts <- as.data.frame(table(directDups.gr$ID))
## Calculate percentage of both direct and inverted dups
perc.direct <- (dirDups.counts$Freq / (invDups.counts$Freq + dirDups.counts$Freq))
perc.inverted <- (invDups.counts$Freq / (invDups.counts$Freq + dirDups.counts$Freq))
## Prepare data.frame for plotting
dirDups.counts$perc <- perc.direct
invDups.counts$perc <- perc.inverted
dirDups.counts$Dir <- 'direct'
invDups.counts$Dir <- 'inverted'
plt.df <- rbind(dirDups.counts, invDups.counts)

## Calculate significance of a difference between count of direct vs inverted invDups
signif.results.all <- list()
for (i in unique(plt.df$Var1)) {
  sub.tab <- plt.df[plt.df$Var1 == i,]
  data.matrix <- matrix(c(sub.tab$Freq[sub.tab$Dir == 'inverted'],
                          sum(sub.tab$Freq)/2,
                          sub.tab$Freq[sub.tab$Dir == 'direct'],
                          sum(sub.tab$Freq)/2),
                        nrow = 2
  )
  ## Calculate Chi-square
  chisq.pval <- chisq.test(data.matrix)$p.value
  suppressWarnings( fish.oddsr <- fisher.test(data.matrix)$estimate )
  signif.results <- data.frame(ID=i, chisq.pval=chisq.pval, fish.oddsr=unname(fish.oddsr))
  signif.results.all[[length(signif.results.all) + 1]] <- signif.results
}
signif.results.all <- do.call(rbind, signif.results.all)
## Do multiple testing correction (Bonferroni)
signif.results.all$chisq.pval.BonfCorr <- p.adjust(p = signif.results.all$chisq.pval, method = 'bonferroni', n = length(signif.results.all$chisq.pval))
signif.results.all$pVal <- format(signif.results.all$chisq.pval.BonfCorr, scientific = TRUE, digits = 2)

plt3 <- ggplot(plt.df, aes(x=Var1, y=Freq, group=Var1)) +
  geom_col(aes(fill=Dir)) +
  geom_text(data=signif.results.all, aes(x=ID, y=Inf, label=pVal), inherit.aes = FALSE, vjust=1) +
  geom_text(aes(group=Var1, label=Freq), position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("Counts") +
  scale_fill_manual(values = c('darkslategray4', 'indianred4'), name="Direction") +
  scale_x_discrete(labels = c('B', 'C', 'G', 'O'))

plt4 <- ggplot(plt.df) +
  geom_col(aes(x=Var1, y=perc, fill=Dir)) +
  xlab("") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c('darkslategray4', 'indianred4')) +
  scale_x_discrete(labels = c('B', 'C', 'G', 'O'))

final.plt <- plot_grid(plt3, plt4, nrow = 1)
destination <- file.path(outputdirectory, "invervedVSdirect_duplications.pdf")
ggsave(final.plt, filename = destination, width = 10, height = 5, useDingbats=FALSE)
destination <- file.path(outputdirectory, "invervedVSdirect_duplications.txt")
write.table(plt.df, file = destination, append = FALSE, quote = FALSE, row.names = FALSE)

## Make plot for paper
#final.plt <- plot_grid(plt1, plt3, ncol = 1, align = 'v', axis = "rl")
#destination <- file.path(outputdirectory, "invDups_summary.pdf")
#ggsave(final.plt, filename = destination, device = 'pdf', width = 5, height = 8)

## Export merged inverted and direct duplication as bed formated file
directDups.gr$ID <- paste0('dirDup_', directDups.gr$ID)
all.invDups.gr$ID <- paste0('invDup_', all.invDups.gr$ID)
all.dups.gr <- c(directDups.gr[,'ID'], all.invDups.gr[,'ID'])
seqlevels(all.dups.gr) <- paste0('chr', c(1:22,'X'))
all.dups.gr <- sort(all.dups.gr)
all.dups.df <- as.data.frame(all.dups.gr)[,c('seqnames','start','end','ID')]
destination <- file.path(outputdirectory, "dirDupANDinvDup.regions.all.bed")
write.table(x = all.dups.df, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE)

## Find shared vs unique inverted duplications ## 
#################################################
apes.invDups.overlaps <- getDisjointOverlapsWeighted(all.invDups.gr, percTh = 50)
overlaps <- split(apes.invDups.overlaps$sub.group, apes.invDups.overlaps$ID)

overlaps.df <- fromList(overlaps)
overlaps.df$n <- sample(1:nrow(overlaps.df))

destination <- file.path(outputDirectory, "upsetR_invDups.pdf")
pdf(destination, width = 8, height = 4, useDingbats=FALSE)
upset(overlaps.df,
      order.by = 'freq', 
      nsets = length(overlaps), 
      sets.bar.color = c("#e6550d","#8856a7","#31a354","#3182bd"),
      sets = c("bonobo", "chimpanzee", "gorilla", "orangutan"),
      queries = list(
        list(query = intersects, params = list("bonobo"), color="#3182bd", active = T),
        list(query = intersects, params = list("chimpanzee"), color="#31a354", active = T),
        list(query = intersects, params = list("gorilla"), color="#8856a7", active = T),
        list(query = intersects, params = list("orangutan"), color="#e6550d", active = T)
      )
    )
dev.off()
## Prepare evolutionary tree based on shared invDups ##
#######################################################
# overlaps.m <- t(overlaps.df)
# ## Get counts of shared events
# #vennPartitions <- get.venn.partitions(overlaps)
# node1 <- length(intersect(intersect(intersect(overlaps[['orangutan']], overlaps[['gorilla']]), overlaps[['bonobo']]), overlaps[['chimpanzee']]))  
# node2 <- length(intersect(intersect(overlaps[['gorilla']], overlaps[['bonobo']]), overlaps[['chimpanzee']]))
# node3 <- length(intersect(overlaps[['bonobo']], overlaps[['chimpanzee']]))
# ## Get counts of unique events
# tip1 <- 6
# tip2 <- 13
# tip3 <- 96
# tip4 <- 120
# ## Contruct the tree
# hc.tree <- hclust(dist(overlaps.m))
# hc.tree <- as.phylo(hc.tree)
# plt <- ggtree(hc.tree) 
# ## Prepare extra annotation data.frames
# annot.nodes <- data.frame(node=unique(hc.tree$edge[,1]), text=c(node1, node2, node3), stringsAsFactors = FALSE)
# annot.tips <- data.frame(node=1:4, good=c(tip1, tip2, tip3, tip4), stringsAsFactors = FALSE)
# ## Add extra annotation to the tree
# plt <- plt %<+% annot.nodes + geom_nodelab(aes(node=y, label=text), geom = "label")
# plt <- plt %<+% annot.tips + geom_tiplab(aes(label=good, subset=isTip), geom = "label", hjust = 1.5) + geom_tiplab()
# plt <- plt + xlim_tree(xlim = c(1:max(hc.tree$edge.length)+2))

## Number of invDups given total size of SDs per chromosomes ##
###############################################################
plt.df <- as.data.frame(all.invDups.gr)
sd.df <- as.data.frame(seg.dup.gr)
sd.size.per.chr <- sd.df %>% group_by(seqnames) %>% summarise(sd.size=sum(width))
data.tab <- plt.df %>% group_by(.dots='seqnames') %>% summarize(counts = n())
data.tab$sd.size <- sd.size.per.chr$sd.size[match(data.tab$seqnames, sd.size.per.chr$seqnames)]

l.mod = lm(counts ~ sd.size, data=data.tab)
conf.level <- predict(l.mod, interval="prediction", level = 0.95)
data.tab <- data.tab %>% mutate(resid=abs(resid(l.mod)), fitted=fitted(l.mod))
data.tab <- cbind(data.tab, conf.level)

r.sq <- paste0('R^2=', round(summary(l.mod)$r.squared, 3))

plt <- data.tab %>% ggplot() + geom_line(aes_string(x='sd.size', y='fitted')) +
  geom_line(aes_string(x='sd.size', y='lwr'), color = "red", linetype = "dashed") +
  geom_line(aes_string(x='sd.size', y='upr'), color = "red", linetype = "dashed") +
  geom_point(aes_string(x='sd.size', y='counts', color='resid')) +
  geom_text(aes_string(x='sd.size', y='counts', label='seqnames'), vjust=-0.5, hjust=-0.1) +
  geom_text(aes(x=-Inf, y=Inf, label=r.sq), inherit.aes = F, hjust=-0.5, vjust=1) +
  scale_colour_gradient(low="blue", high="red") +
  scale_x_continuous(labels = comma) +
  labs(x="Total human SD bases per Chromosome (bp)", y='Number of inverted duplications', colour="Residuals")

destination <- file.path(outputdirectory, "invDups_vs_SDsize.pdf")
ggsave(plt, filename = destination, device = 'pdf', width = 7, height = 5, useDingbats=FALSE)

#plt <- data.tab %>% ggplot() + geom_point(aes(x=sd.size, y=counts)) + 
#  geom_text(data=data.tab, aes(x=sd.size, y=counts, label=seqnames), vjust=-0.5, hjust=0.1) +
#  scale_x_continuous(labels = comma) +  
#  xlab("Total human SD bases per Chromosome (bp)") +
#  ylab("# of invDups per Chromosome") +
#  theme_bw()

## Plot genome-wide distribution of inverted and direct duplications ##
#######################################################################
## NOTE: load helper function at the end of this script !!!
## Put all duplications together
invertDups.gr <- all.invDups.gr[,'ID']
invertDups.gr$Direction <- 'inverted'
directDups.gr$Direction <- 'direct'
all.dups.gr <- c(invertDups.gr, directDups.gr)
## Recode NHP IDs to shorter names
all.dups.gr$ID <- recode(all.dups.gr$ID, 'chimpanzee' = 'C', 'bonobo' = 'B', 'gorilla' = 'G', 'orangutan' = 'O')

## Select only certain chromosomes
select.chr <- c('chr5','chr7','chr10','chr16','chr17')
dups.gr <- keepSeqlevels(all.dups.gr, value = select.chr, pruning.mode = 'coarse')
sub.segDup.gr <- keepSeqlevels(seg.dup.gr, value = select.chr, pruning.mode = 'coarse')
sub.segDup.gr <- reduce(sub.segDup.gr)
sub.segDup.gr <- sub.segDup.gr[width(sub.segDup.gr) >= 10000] ## Keep only SDs 10b and longer!!!
dups.gr <- sort(dups.gr)
sub.segDup.gr <- sort(sub.segDup.gr)
## Add human SD annotation
sub.segDup.gr$ID <- 'H'
sub.segDup.gr$Direction <- 'SD'
dups.gr <- sort(c(dups.gr, sub.segDup.gr))

## Transfer chromosome coords into genome-wide coords
#seqlevels(dups.gr) <- paste0("chr", c(1:22, "X"))
seqlengths(dups.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(dups.gr)]
#seqlengths(sub.segDup.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(sub.segDup.gr)]
dups.gr <- transCoord(dups.gr)
#sub.segDup.gr <- transCoord(sub.segDup.gr)
## Get positions of ends of each chromosome to plot lones between the chromosomes
cum.seqlengths <- cumsum(as.numeric(seqlengths(dups.gr)))
if (length(cum.seqlengths) > 1) {
  chr.lines <- cum.seqlengths[-length(cum.seqlengths)]
  names(chr.lines) <- seqlevels(dups.gr)[-length(seqlevels(dups.gr))]
} else {
  chr.lines <- 0 
}
## Get positions of each chromosomes names
cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(dups.gr) ) )
#names(chr.label.pos) <- gsub("chr", "", names(chr.label.pos))

## Prepare final plot
my_theme <- theme( strip.text.y = element_text(angle = 360),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   legend.position = "none")

## Plot human SDs
segDup.df <- as.data.frame(dups.gr[dups.gr$ID == 'H'])
annot <- ggplot() +
  geom_rect(data=segDup.df, aes(xmin=start.genome, xmax=end.genome, ymin=0, ymax=1), fill='orange') +
  geom_vline(xintercept = chr.lines) +
  scale_x_continuous(expand = c(0,0)) + my_theme

## Plot NHP dir+inv duplications
plt.df <- as.data.frame(dups.gr[dups.gr$ID != 'H'])
plt <- ggplot() + 
  geom_vline(xintercept = chr.lines, color="gray43") +
  geom_point(data=plt.df, aes(x=start.genome, y=1, color=Direction, shape=Direction), position = position_jitter(height = 0.5), stroke=5) +
  scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) +
  scale_color_manual(values = c('darkslategray4', 'indianred4', 'orange')) +
  scale_shape_manual(values=c(62, 60, 108)) +
  xlab("Chromosome") +
  facet_grid(ID ~ .) +
  my_theme
## Arrange plots into a sigle column
suppressPackageStartupMessages( library(ggbio) )
final.plt <- tracks(annot, plt, heights = c(0.2, 1))
## Export final plot
destination <- file.path(outputdirectory, "selectedChrDistribution_DirVSInvDups.pdf")
ggsave(final.plt, filename = destination, width = 12, height = 4, useDingbats=FALSE)

## HELPER FUNCTION ##
transCoord <- function(gr) {
  cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
  cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
  gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  return(gr)
}
