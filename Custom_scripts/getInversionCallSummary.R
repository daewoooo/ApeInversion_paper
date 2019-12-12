## Load required libraries ##
#############################
library(VennDiagram)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyr)
library(reshape2)
library(primatR)
library(ggtree)
library(UpSetR)
library(gtools)
library(ape)

outputDirectory <- "/home/porubsky/WORK/Great_apes/Summary_plots/"

## Load required data ##
########################
## segDup track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load repeatMasker track
rmsk <- read.table("/home/porubsky/WORK/Great_apes/Annotations/repeatMasker_GRCh38.bed.gz", header=FALSE, stringsAsFactors = FALSE)
rmsk.gr <- GRanges(seqnames=rmsk$V6, ranges=IRanges(start=rmsk$V7, end=rmsk$V8), strand=rmsk$V10, repName=rmsk$V11, repClass=rmsk$V12, repFamily=rmsk$V13)
rmsk.THE1.gr <- rmsk.gr[grep(rmsk.gr$repName, pattern = 'THE1')]
ranges2UCSC(rmsk.THE1.gr, outputDirectory = "/home/porubsky/WORK/Great_apes/Annotations/", index = 'THE1elems', id.field = 'repName')

## NCBI refSeqGene Track
#genes <- read.table(gzfile("/home/porubsky/WORK/Great_apes/Annotations/NCBI_refSeqGenes_GRCh38_curated.gz"))
#genes.gr <- GRanges(seqnames=genes$V3, strand=genes$V4, ranges=IRanges(start=genes$V5, end=genes$V6), geneName=genes$V13)
#genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = 'coarse')
#genes.gr.reduced <- split(genes.gr, genes.gr$geneName)
#genes.gr.reduced <- endoapply(genes.gr.reduced, collapseOverlaps)

## Load gene lists
vega68 <- read.table(gzfile("/home/porubsky/WORK/Great_apes/Annotations/vega68_GRCh38_geneList_mart_export.txt.gz"), header=TRUE, sep=",", stringsAsFactors = FALSE)
gencode29 <- read.table("/home/porubsky/WORK/Great_apes/Annotations/gencode.v29.annotation..processed.txt", header = TRUE)
vega68.gr <- GRanges(seqnames=vega68$Chromosome.scaffold.name, ranges=IRanges(start=vega68$Gene.start..bp., end=vega68$Gene.end..bp.), gene.name=vega68$Gene.name, gene.type=vega68$Gene.type, gene.id=vega68$Gene.stable.ID)
gencode29.gr <- GRanges(seqnames=gencode29$chr, ranges=IRanges(start=gencode29$start, end=gencode29$end), gene.name=gencode29$gene_name, gene.type=gencode29$gene_type, gene.id=gencode29$gene_id)
vega68PLUSgencode29.gr <- suppressWarnings( c(vega68.gr, gencode29.gr) )
vega68PLUSgencode29.gr <- keepStandardChromosomes(vega68PLUSgencode29.gr, pruning.mode = 'coarse')
genes.gr <- reduce(vega68PLUSgencode29.gr)

## Load geneHancer track
geneHancer <- read.table("/home/porubsky/WORK/Great_apes/Annotations/geneHancerRegElements_GRCh38.bed.gz", header=FALSE, stringsAsFactors = FALSE)
geneHancer.gr <-GRanges(seqnames=geneHancer$V1, ranges=IRanges(start=geneHancer$V2, end=geneHancer$V3), name=geneHancer$V4, type=geneHancer$V11, categ=geneHancer$V12)
# Keep only enhancer regions
geneHancer.gr <- geneHancer.gr[grep(geneHancer.gr$type, pattern = 'Enhancer')]

## Load complete Strand-seq dataset ##
######################################
## Load complete dataset including putative human specific inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))

## Data plotting before filtering ##
####################################
plotColumnCounts(all.inversion.calls.filt, colName = 'gen', facetID = 'ID', colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
plotColumnCounts(all.inversion.calls.filt, colName = 'SVclass', facetID = 'ID', colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
#rangesSizeDistribution(all.inversion.calls.filt, plotUniqueBases=FALSE)
#rangesSizeDistribution(all.inversion.calls.filt, plotUniqueBases=TRUE)

## Compile Strand-seq callset for subsequent analysis ##
########################################################
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
## Load simple inversion with annotated potential genome misorients
#all.SimpleInversion.calls.filt.annot <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.SimpleInversion.calls.filt.annot.RData"))
## Retain only inverted duplications
all.invertedDups.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Color set ##
## All great apes
# c('#3182bd','#31a354','#8856a7','#a6611a','#e6550d')
## Non-human primates only
# c('#3182bd','#31a354','#8856a7','#e6550d') 

## Data plotting after filtering ##
###################################
## Make sure GRanges obejct is properly ordered
seqlevels(all.SimpleInversion.calls.filt) <- gtools::mixedsort(seqlevels(all.SimpleInversion.calls.filt))
seqlevels(all.invertedDups.calls.filt) <- gtools::mixedsort(seqlevels(all.invertedDups.calls.filt))
#plotColumnCounts(all.SimpleInversion.calls.filt, colName = 'gen')
## Plot simple inversion size distribution and genotype distribution
plt1 <- plotColumnCounts(all.SimpleInversion.calls.filt, colName = 'gen', facetID = 'ID', colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
plt2 <- rangesSizeDistribution(all.SimpleInversion.calls.filt, plotUniqueBases=FALSE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
final.plt <- plot_grid(plt2, plt1, nrow = 1, rel_widths = c(4,1))
destination <- file.path(outputDirectory, "simpleInversionSizeDistAndGenotypes.pdf")
ggsave(filename = destination, plot = final.plt, device = 'pdf', width = 9, height = 5, useDingbats=FALSE)
rangesSizeDistribution(all.SimpleInversion.calls.filt, plotUniqueBases=TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
## Plot genotype frequency and total inverted bases per chromosome
plt1 <- plotColumnCountsPerChr(all.SimpleInversion.calls.filt, colName = 'gen', facetID = 'ID', normChrSize = TRUE)
plt2 <- basesPerGenotypePerChr(all.SimpleInversion.calls.filt, normChrSize = TRUE)
final.plt <- plot_grid(plt1, plt2, ncol = 1, align = 'v')
destination <- file.path(outputDirectory, "genotypeFreqAndTotalSimpleINVbasesPerChromosome.pdf")
ggsave(filename = destination, plot = final.plt, device = 'pdf', width = 10, height = 7, useDingbats=FALSE)
## Plot distribution of simple inversion and inverted duplication given the chromosome size
plt1 <- eventsPerChrSizeScatter(gr = all.SimpleInversion.calls.filt, bsgenome = BSgenome.Hsapiens.UCSC.hg38, lm = TRUE)
plt2 <- eventsPerChrSizeScatter(gr = all.SimpleInversion.calls.filt, bsgenome = BSgenome.Hsapiens.UCSC.hg38, colBy = 'gen', facetID = 'gen')
#destination <- file.path(outputDirectory, "simpleINV_vs_CHRsize.pdf")
#ggsave(filename = destination, plot = plt1, width = 7, height = 5)
destination <- file.path(outputDirectory, "simpleINV_vs_CHRsize_perGEN.pdf")
ggsave(filename = destination, plot = plt2, width = 8, height = 6, useDingbats=FALSE)

plt3 <- eventsPerChrSizeScatter(all.invertedDups.calls.filt, bsgenome = BSgenome.Hsapiens.UCSC.hg38, lm = TRUE)
destination <- file.path(outputDirectory, "invDups_vs_CHRsize.pdf")
ggsave(filename = destination, plot = plt3, width = 8, height = 6, useDingbats=FALSE)

## Number of invDups given total size of SDs per chromosomes ##
###############################################################
plt.df <- as.data.frame(all.invertedDups.calls.filt)
sd.df <- as.data.frame(seg.dup.gr)
sd.size.per.chr <- sd.df %>% group_by(seqnames) %>% summarise(sd.size=sum(width))
data.tab <- plt.df %>% group_by(.dots='seqnames') %>% summarize(counts = n())
data.tab$sd.size <- sd.size.per.chr$sd.size[match(data.tab$seqnames, sd.size.per.chr$seqnames)]

l.mod = lm(counts ~ sd.size, data=data.tab)
conf.level <- predict(l.mod, interval="prediction", level = 0.95)
data.tab <- data.tab %>% mutate(resid=abs(resid(l.mod)), fitted=fitted(l.mod))
data.tab <- cbind(data.tab, conf.level)

r.sq <- paste0('R^2=', round(summary(l.mod)$r.squared, 3))
r.sq.df <- as.data.frame(r.sq)

plt4 <- data.tab %>% ggplot() + geom_line(aes_string(x='sd.size', y='fitted')) +
  geom_line(aes_string(x='sd.size', y='lwr'), color = "red", linetype = "dashed") +
  geom_line(aes_string(x='sd.size', y='upr'), color = "red", linetype = "dashed") +
  geom_point(aes_string(x='sd.size', y='counts', color='resid')) +
  geom_text(aes_string(x='sd.size', y='counts', label='seqnames'), vjust=-0.5, hjust=-0.1) +
  geom_text(data = r.sq.df, aes(x=-Inf, y=Inf, label=r.sq), inherit.aes = F, hjust=-0.5, vjust=1) +
  scale_colour_gradient(low="blue", high="red") +
  scale_x_continuous(labels = comma) +
  labs(x="Total human SD bases per Chromosome (bp)", y='Number of inverted duplications', colour="Residuals")
## Export plot
final.plt <- plot_grid(plt1, plt4, nrow = 1)
destination <- file.path(outputDirectory, "simpleINV_vs_CHRsize_invDup_vs_SDcont.pdf")
ggsave(filename = destination, device = 'pdf', plot = final.plt, width = 12, height = 4, useDingbats=FALSE)


## Plot overlap of inversions with genes
#Reduce multiple isoforms of the same gene!!!
#plotOverlapWithRanges(query.gr=all.SimpleInversion.calls.filt.noHS, subject.gr=genes.gr)

## Plot inversion size distribution per chromosome ##
#####################################################
size.dist.df <- data.frame(ID=seqnames(all.SimpleInversion.calls.filt), size=width(all.SimpleInversion.calls.filt))
size.dist.df <- size.dist.df %>% group_by(ID) %>% mutate(mean.size = mean(size))
size.dist.df$ID <- factor(size.dist.df$ID, levels = paste0('chr' ,c(1:22,'X')))
plt <- ggplot(size.dist.df) + 
  geom_point(aes(x=ID, y=size), position = position_jitter(width = 0.2)) +
  geom_point(aes(x=ID, y=mean.size), shape=95, color='red', size=10) +
  scale_y_continuous(trans = 'log10', breaks=c(1000, 10000, 100000, 1000000, 10000000, 1000000000), labels=c('1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb')) +
  xlab("") +
  ylab("Inversion size") +
  geom_hline(yintercept = c(1000, 100000), linetype='dashed', color='blue')
destination <- file.path(outputDirectory, "invSizeDistPerChrom.pdf")
ggsave(filename = destination, plot = plt, device='pdf', width = 12, height = 4, useDingbats=FALSE)

## Get inversion summary per chromosome given chromosome the size and SD content ##
###################################################################################
df <- as.data.frame(all.inversion.calls.filt)
df$seqnames <- factor(df$seqnames, levels = paste0('chr', c(1:22,'X')))
df <- df %>% group_by(seqnames, SVclass) %>% summarise(count=n())
df <- tidyr::spread(df, SVclass, count)
df2 <- as.data.frame(seg.dup.gr[seqnames(seg.dup.gr) %in% paste0('chr', c(1:22,'X'))])
df2$seqnames <- factor(df2$seqnames, levels = paste0('chr', c(1:22,'X')))
df2 <- df2 %>% group_by(seqnames) %>% summarise(SDsize=sum(width))
df$SDsize <- df2$SDsize
df$ChrSize <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[df$seqnames]
## Add mean inversion and invDup size
df.tmp <- as.data.frame(all.inversion.calls.filt)
df.tmp$seqnames <- factor(df.tmp$seqnames, levels = paste0('chr', c(1:22,'X')))
df.tmp <- df.tmp %>% group_by(seqnames, SVclass) %>% summarise(Mean.size=mean(width))
df.tmp <- tidyr::spread(df.tmp, SVclass, Mean.size)
colnames(df.tmp) <- paste0(colnames(df.tmp),'.meanSize')
df <- cbind(df, df.tmp)
## Calculate variable relationship using linear model (lm)
#gl.mod <- glm(INV ~ SDsize + ChrSize + INV.meanSize, data = df)
#l.mod <- lm(INV ~ SDsize + ChrSize + INV.meanSize, data = df)
l.mod <- lm(INV ~ SDsize + ChrSize, data = df)
suppressWarnings( conf.level <- predict(l.mod, interval="prediction", level = 0.95) )
df$resid <- abs(resid(l.mod))
df$fitted <- fitted(l.mod)

r.sq <- paste0('R^2=', round(summary(l.mod)$r.squared, 3))
r.sq.df <- as.data.frame(r.sq)

plt.df <- cbind(as.data.frame(df), conf.level)
plt <- ggplot(plt.df) + geom_line(aes_string(x='ChrSize', y='fitted')) +
  geom_line(aes_string(x='ChrSize', y='lwr'), color = "red", linetype = "dashed") +
  geom_line(aes_string(x='ChrSize', y='upr'), color = "red", linetype = "dashed") +
  geom_point(aes_string(x='ChrSize', y='INV', color='resid')) +
  geom_text(aes_string(x='ChrSize', y='INV', label='seqnames'), vjust=-0.5, hjust=-0.1) +
  geom_text(data = r.sq.df, aes(x=-Inf, y=Inf, label=r.sq), inherit.aes = F, hjust=-0.5, vjust=1) +
  scale_colour_gradient(low="blue", high="red") +
  scale_x_continuous(labels = comma) +
  labs(x="Chromosome size (bp)", y='counts', colour="Residuals")

## Plot genome-wide ideogram ##
###############################
plt <- genomewideRangesIdeo(all.SimpleInversion.calls.filt, userTrack = seg.dup.gr, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), bsgenome = BSgenome.Hsapiens.UCSC.hg38)
destination <- file.path(outputDirectory, "simpleINV_genomewide.pdf")
ggsave(filename = destination, plot = plt, device = 'pdf', width = 8, height = 4)
plt <- genomewideRangesIdeo(all.invertedDups.calls.filt, userTrack = seg.dup.gr, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), bsgenome = BSgenome.Hsapiens.UCSC.hg38)
destination <- file.path(outputDirectory, "invDups_genomewide.pdf")
ggsave(filename = destination, plot = plt, device = 'pdf', width = 8, height = 4, useDingbats=FALSE)

## Measure the distance of simple inversions to human SDs   ##
## Plot size distribution of inversion calls per individual ##
##############################################################
all.SimpleInversion.calls.filt.dists <- range2rangeDistance(gr=all.SimpleInversion.calls.filt, userTrack=seg.dup.gr, allow.overlap = TRUE)
SDflank <- which(all.SimpleInversion.calls.filt.dists$leftDist < 5000 & all.SimpleInversion.calls.filt.dists$rightDist < 5000 & all.SimpleInversion.calls.filt.dists$leftDist >= 0 & all.SimpleInversion.calls.filt.dists$rightDist >= 0)
## Plot size distribution of all simple inversions
plt1 <- rangesSizeDistribution(gr = all.SimpleInversion.calls.filt, violin = TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), title = "All simple inversions")
## Plot size distribution of simple inversions 10kb and longer in size
plt2 <- rangesSizeDistribution(gr = all.SimpleInversion.calls.filt[width(all.SimpleInversion.calls.filt) >= 10000], violin = TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), title = "Simple inversions 10kb and longer")
## Plot size distribution of inversions flanked by SDs
plt3 <- rangesSizeDistribution(gr = all.SimpleInversion.calls.filt[SDflank], violin = TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), title = "Simple inversion flanked by SDs")
## Plot size distribution of inversions NOT flanked by SDs
plt4 <- rangesSizeDistribution(gr = all.SimpleInversion.calls.filt[-SDflank], violin = TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), title = "Simple inversion not flanked by SDs")
## Plot size distribution of invDups
plt5 <- rangesSizeDistribution(gr = all.invertedDups.calls.filt, violin = TRUE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), title = "Inverted duplications")
final.plt <- plot_grid(plt1, plt2, plt3, plt4, plt5, ncol = 1)
destination <- file.path(outputDirectory, "simpleINV_sizeDist_violins.pdf")
ggsave(filename = destination, plot = final.plt, device = 'pdf', width = 9, height = 15, useDingbats=FALSE)

# SimpleInversion.SDflank <- all.SimpleInversion.calls.filt[SDflank]
# SimpleInversion.noSDflank <- all.SimpleInversion.calls.filt[-SDflank]
# SimpleInversion.SDflank.GENcounts <- table(SimpleInversion.SDflank$gen)
# SimpleInversion.noSDflank.GENcounts <- table(SimpleInversion.noSDflank$gen)
# (SimpleInversion.SDflank.GENcounts[2] / sum(SimpleInversion.SDflank.GENcounts)) * 100
# (SimpleInversion.noSDflank.GENcounts[2] / sum(SimpleInversion.noSDflank.GENcounts)) * 100

## Plot size distribution of simple inversions (flanked & not flanked by SD) and inverted duplications ##
#########################################################################################################
## Subset the data
simpleINV.SDflank <- all.SimpleInversion.calls.filt[SDflank]
simpleINV.SDflank$SVclass <- 'SDflankINV'
simpleINV.noSDflank <- all.SimpleInversion.calls.filt[-SDflank]
simpleINV.noSDflank$SVclass <- 'noSDflankINV' 
plt.gr <- c(simpleINV.SDflank, simpleINV.noSDflank, all.invertedDups.calls.filt)
## Construct plotting table
plt.df <- as.data.frame(plt.gr)
plt.df$SVclass <- factor(plt.df$SVclass, levels=c('SDflankINV','noSDflankINV','invDup'))
plt.df <- plt.df[order(plt.df$width),]
xaxis.labels <- paste0(levels(plt.df$SVclass), "\n(n=", table(plt.df$SVclass), ")")
## Plot the data
plt <- ggplot(plt.df, aes(x=SVclass, y=width, fill=SVclass)) + geom_violin(trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  #geom_dotplot(aes(x=SVclass, y=width), binaxis='y', stackdir='center', dotsize=0.05, binwidth = 1) +
  scale_y_continuous(breaks=c(1000,10000,100000,1000000), labels = comma, trans = 'log10') +
  scale_x_discrete(labels=xaxis.labels) +
  geom_hline(yintercept = c(1000, 10000, 1000000), linetype="dashed") +
  ylab("Inversion size (log10)") +
  xlab("") + 
  scale_fill_brewer(palette="Dark2") + 
  theme_minimal() +
  stat_summary(fun.data='mean_sdl', size=0.5, geom="pointrange", color="white")
## Test if the size distribution of SD flanked vs not-flanked is significant
wilcox.p.val <- wilcox.test(width(simpleINV.SDflank), width(simpleINV.noSDflank), paired = FALSE, alternative = 'greater')
#kruskal.p.val <- kruskal.test(list(width(simpleINV.noSDflank), width(simpleINV.SDflank)))
destination <- file.path(outputDirectory, "SVclasses_sizeDist.pdf")
ggsave(filename = destination, plot = plt, device = 'pdf', width = 6, height = 4, useDingbats=FALSE)

## Measure the distance of human specific inversions to gene enhancers ... ##
#############################################################################
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
## Run permutation test
HSinv.vs.geneHancer.enrich <- permDistanceToFeature(query.gr = HSinv.gr, feature.gr = geneHancer.gr, nperm = 1000, randomize = 'query', mask.gr = seg.dup.gr, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
## Plot results
dataset <- HSinv.vs.geneHancer.enrich
ttest.p.val <- format(round(dataset$ttest.p.val, digits = 6), scientific = FALSE)
wilcox.p.val <- format(round(dataset$wilcox.p.val, digits = 6), scientific = FALSE)

plt.df <- data.frame(dist = c(dataset$observed.dist, dataset$rand.dist), 
                     ID = rep(c('observed', 'permuted'), c(length(dataset$observed.dist), length(dataset$rand.dist))))

plt1 <- ggplot(plt.df, aes(x=ID, y=dist, fill=ID)) + 
  geom_violin(trim = FALSE, scale = TRUE) +
  ylab("Closest gene distance (bp)") +
  xlab("") +
  ggtitle("Test", subtitle = paste0("t-test: ", ttest.p.val, "\n", "wilcox-test: ", wilcox.p.val))
plt2 <- ggplot(plt.df, aes(x=dist, group=ID, fill=ID)) + 
  geom_histogram(bins = 100) + 
  scale_y_continuous(trans='log10') +
  ggtitle("Test", subtitle = paste0("t-test: ", ttest.p.val, "\n", "wilcox-test: ", wilcox.p.val))


## Plot distance of HS and polymorphic inversions to enhancers  ##
##################################################################
outputDirectory <- "/home/porubsky/WORK/Great_apes/Ape_paper/Plots/"
## Load predicted HS invertsions  
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
## Load putative polymorphic inversions
polymorphic.inversions <- get(load("/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites/shared_HGSVC&NHP_sites.RData"))
## Load geneHancer track
geneHancer <- read.table("/home/porubsky/WORK/Great_apes/Annotations/geneHancerRegElements_GRCh38.bed.gz", header=FALSE, stringsAsFactors = FALSE)
geneHancer.gr <-GRanges(seqnames=geneHancer$V1, ranges=IRanges(start=geneHancer$V2, end=geneHancer$V3), name=geneHancer$V4, type=geneHancer$V11, categ=geneHancer$V12)
## Keep only enhancer regions
geneHancer.gr <- geneHancer.gr[grep(geneHancer.gr$type, pattern = 'Enhancer')]
geneHancer.gr <- geneHancer.gr[,0]
geneHancer.gr$ID <- 'enhancer'
## Get distance of HS and polymorphic inv to enhancers
HSinv.gr <- HSinv.gr[,0]
HSinv.gr$ID <- 'HSinv'
polymorphic.inversions <- polymorphic.inversions[,0]
polymorphic.inversions$ID <- 'polymINV'

roi.inv <- sort(c(HSinv.gr, polymorphic.inversions))
hits.enhancer <- findOverlaps(roi.inv, geneHancer.gr)
dist.to.enhancer <- getMinDist(gr = roi.inv, userTrack = geneHancer.gr)
## Export HS inversion with distance to above mentioned features
roi.inv$dist.to.enhancer <- dist.to.enhancer
## Set regions that overlap directly with a feature to zero
roi.inv$dist.to.enhancer[unique(queryHits(hits.enhancer))] <- 0
## Export results
destination <- file.path(outputDirectory, "HSinvAndPolyInv_dist2enhancers.RData")
save(roi.inv, file = destination)
## Sort by chromosome
seqlevels(roi.inv) <- mixedsort(seqlevels(roi.inv))
roi.inv <- sort(roi.inv)
## Prepare data for plotting
plt.df <- as.data.frame(roi.inv)
plt.df$roi <- as.character(roi.inv[,0])

plt.df <- reshape2::melt(plt.df, measure.vars='dist.to.enhancer')
plt <- plt.df %>% mutate(enhancer.overlap = value == 0) %>% group_by(ID, enhancer.overlap) %>% summarise(count=n()) %>% 
  ggplot() + geom_col(aes(x=ID, y=count, fill=enhancer.overlap)) +
  scale_fill_manual(values=brewer.pal(n = 4, name = "Dark2")) +
  geom_text(aes(x=ID, y=count, label=count), vjust=1, position = position_stack(vjust = .5)) +
  theme_bw() +
  xlab("")
## Save the plot
destination <- file.path(outputDirectory, "HSinvAndPolyInv_dist2enhancers.pdf")
ggsave(plt, filename = destination, width = 4, height = 6, useDingbats=FALSE)

# plt <- ggplot(plt.df, aes(x=roi, y=value, color=ID, group=ID)) + 
#   geom_point(position=position_dodge(width=0.5), size=3) +
#   scale_color_manual(values=brewer.pal(n = 4, name = "Dark2")) +
#   scale_y_continuous(breaks=c(1000, 10000, 100000, 1000000), labels = comma, trans = 'log10') +
#   geom_hline(yintercept = c(1000, 10000, 100000, 1000000), linetype="dashed") +
#   theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
#   xlab("") +
#   ylab("Log10 Distance to enhancer (bp)")
#   #coord_flip()
# ## Save the plot
# destination <- file.path(outputDirectory, "HSinvAndPolyInv_dist2enhancers.pdf")
# ggsave(plt, filename = destination, width = 12, height = 4, useDingbats=FALSE)


## Inversions and distribution of THE1 element ##
#################################################
# SimpleInversion <- all.SimpleInversion.calls.filt
# SimpleInversion <- SimpleInversion[width(SimpleInversion) < 1000000]
# SimpleInversion$THE1.hits <- countOverlaps(SimpleInversion, rmsk.THE1.gr, maxgap = 10000)
# SimpleInversion$THE1.hits.norm <- (SimpleInversion$THE1.hits / width(SimpleInversion))*1000000
# 
# invertedDups <- all.invertedDups.calls.filt
# invertedDups$THE1.hits <- countOverlaps(invertedDups, rmsk.THE1.gr)
# invertedDups$THE1.hits.norm <- (invertedDups$THE1.hits / width(invertedDups))*1000000
# 
# all.ranges <- c(SimpleInversion, invertedDups)
# 
# plt.df <- as.data.frame(all.ranges)
# plt.df <- plt.df[plt.df$THE1.hits > 0,]
# plt.df <- plt.df[order(plt.df$THE1.hits.norm, decreasing = TRUE),]
# plt.df$x <- 1:nrow(plt.df)
# 
# ggplot(plt.df) + 
#   geom_col(aes(x=x, y=THE1.hits.norm, fill=SVclass)) +
#   scale_fill_manual(values = c('darkkhaki', 'royalblue3')) +
#   ylab("THE1 elem norm hits") +
#   xlab("Ordered INV & invDups")




            