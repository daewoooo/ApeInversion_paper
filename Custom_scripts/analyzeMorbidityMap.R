## Load required libraries
library(primatR)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

outputDirectory <- "/home/porubsky/WORK/Great_apes/Annotations/Morbidity_map/"

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Morbidity map data [gold standard]
morbMap.gold <- read.table("/home/porubsky/WORK/Great_apes/Annotations/Morbidity_map/Sig37disorders_hg19_lifted2hg38.bed", header = FALSE)
morbMap.gold.gr <- GRanges(seqnames=morbMap.gold$V1, ranges=IRanges(start=morbMap.gold$V2, end=morbMap.gold$V3), name=morbMap.gold$V4)
## Export Morbidity map region as UCSC browser formatted file
ranges2UCSC(gr = morbMap.gold.gr, outputDirectory = outputDirectory, index = "morbidityMapCoe_gold_grch38", colorRGB = "234,56,87")

## Morbidity map data [all]
#morbMap.all <- read.table("/home/porubsky/WORK/Great_apes/Annotations/Morbidity_map/Signature_29085SamplesRareUniqueCalls_lifted_to_hg38.bed")
#morbMap.all.gr <- GRanges(seqnames=morbMap.all$V1, ranges=IRanges(start=morbMap.all$V2, end=morbMap.all$V3), ID=morbMap.all$V4, CN=morbMap.all$V5)

## Load complete dataset including putative human specific inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
all.invertedDups.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Load putative human specific inversions
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
## Load inversion hotspots
hotspots.gr <- get(load("/home/porubsky/WORK/Great_apes/Hotspots_breakpoints/hotspots.RData"))
## Load putative polymorphic sites
polymorph.gr <- get(load("/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites/shared_HGSVC&NHP_sites.RData"))

## Find overlaps with Morbidity map
## Get reciprocal overlap of simple INV with gold standard morbidty map (50% threshold)
#all.SimpleInversion.calls.filt <- getReciprocalOverlaps(query = all.SimpleInversion.calls.filt, subject = morbMap.gold.gr, report = 'query')
#SimpleInversion.calls.morbMap50 <- all.SimpleInversion.calls.filt[all.SimpleInversion.calls.filt$perc.overlap >= 50]
## Get reciprocal overlap of invDup with gold standard morbidty map
#all.invertedDups.calls.filt <- getReciprocalOverlaps(query = all.invertedDups.calls.filt, subject = morbMap.gold.gr, report = 'query')
#invertedDups.calls.morbMap50 <- all.invertedDups.calls.filt[all.invertedDups.calls.filt$perc.overlap >= 50]
## Get reciprocal overlap all INV calls with all morbidty map regions
#all.inversion.calls.filt <- getReciprocalOverlaps(query = all.inversion.calls.filt, subject = morbMap.all.gr, report = 'query')
#inversion.calls.filt.allmorbMap50 <- all.inversion.calls.filt[all.inversion.calls.filt$perc.overlap >= 50]
## Get reciprocal overlap of CNV regions against simple inversions
morbMap.gold.gr.overlap <- getReciprocalOverlaps(query = morbMap.gold.gr, subject = all.SimpleInversion.calls.filt, report = 'query')
## Export overlap of morbid CNVs with simple INVs
morbMap.gold.gr.overlap.df <- as.data.frame(morbMap.gold.gr.overlap)
destination <- file.path(outputDirectory, 'CNVregions_simpleINV_overlap.csv')
write.table(x = morbMap.gold.gr.overlap.df, file = destination, quote = FALSE, row.names = FALSE, sep = ',')

## Annotate simple inversions overlapping with morbidity map
morbMap.gold.gr.overlap$annot <- 'simpleINV'
hits <- findOverlaps(morbMap.gold.gr.overlap, HSinv.gr)
morbMap.gold.gr.overlap$annot[unique(queryHits(hits))] <- 'HSinv'
hits <- findOverlaps(morbMap.gold.gr.overlap, polymorph.gr)
morbMap.gold.gr.overlap$annot[unique(queryHits(hits))] <- 'PolymInv'
morbMap.gold.gr.overlap$annot[morbMap.gold.gr.overlap$perc.overlap == 0] <- NA

## Regenptype morbidity map regions using NHP composite files
## Load ape's composite files
chimp.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
chimp.data$ID <- 'chimpanzee'
bonobo.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
bonobo.data$ID <- 'bonobo'
gorilla.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
gorilla.data$ID <- 'gorilla'
orangutan.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
orangutan.data$ID <- 'orangutan'

## Re-genotype morbid CNVs ##
#############################
## Get CNV regions with at least 50% reciprocal overlap
morbMap.gold.gr.overlap.filt <- morbMap.gold.gr.overlap[morbMap.gold.gr.overlap$perc.overlap >= 50]
## Remove flanking SDs from morbid CNV regions
morbMap.gold.gr.noSD <- primatR::subtractRegions(gr = morbMap.gold.gr.overlap.filt, remove.gr = reduce(seg.dup.gr), mode = 'flanks')
## Set ranges to be re-genotyped
regions.to.genotype <- morbMap.gold.gr.noSD
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

#genotyped.regions[which(nhp.alleles.dir[,4] == 2)]

## Recode genotypes to number in order to calculate allele frequencies
nhp.genotypes <- mcols(genotyped.regions)[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]
nhp.alleles.inv <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
nhp.alleles.dir <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
## Calculate allele freaquencies for NHP genotypes per site
#nhp.total.inv.alleles <- rowSums(nhp.alleles.inv)
#nhp.total.dir.alleles <- rowSums(nhp.alleles.dir)
#nhp.inv.allele.freq <- nhp.total.inv.alleles / (nhp.total.inv.alleles + nhp.total.dir.alleles)
## Calculate allele frequencies for NHP genotypes per individual
nhp.total.inv.alleles <- colSums(nhp.alleles.inv)
nhp.total.dir.alleles <- colSums(nhp.alleles.dir)
nhp.inv.allele.freq <- nhp.total.inv.alleles / (nhp.total.inv.alleles + nhp.total.dir.alleles)
nhp.inv.allele.freq.df <- as.data.frame(nhp.inv.allele.freq)
destination <- file.path(outputDirectory, 'overlapped_CNVregions_InvAlleleFreq.txt')
write.table(x = nhp.inv.allele.freq.df, file = destination, quote = FALSE, row.names = FALSE)

## Re-genotype morbid CNVs in HGSVC ##
######################################
HGSVC.samples <- c('HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240')
genotyped.regions.hgsvc <- morbMap.gold.gr.noSD
for (sample in HGSVC.samples) {
  sample.data <- get(load(paste0("/home/porubsky/WORK/Great_apes/ChromosomeX/Phased_HGSVC/",sample, "_compositeFile/data/compositeFile.RData")))
  sample.data$ID <- sample

  ## Genotype directional data
  directional.reads <- sample.data$fragments
  ## Directionality of HGSVC composite files is flipped!!!
  strand(directional.reads) <- ifelse(strand(directional.reads) == '+', '-', '+')
  genotyped.regions.hgsvc <- genotypeRegions(regions = genotyped.regions.hgsvc, directional.reads = directional.reads, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = sample)
}
hgsvc.genotypes <- mcols(genotyped.regions.hgsvc)[,c('genoT_HG00512','genoT_HG00513','genoT_HG00731','genoT_HG00732','genoT_NA19238','genoT_NA19239')]
hgsvc.alleles.inv <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
hgsvc.alleles.dir <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
## Calculate allele frequencies for NHP genotypes per individual
hgsvc.total.inv.alleles <- colSums(hgsvc.alleles.inv)
hgsvc.total.dir.alleles <- colSums(hgsvc.alleles.dir)
hgsvc.inv.allele.freq <- hgsvc.total.inv.alleles / (hgsvc.total.inv.alleles + hgsvc.total.dir.alleles)
hgsvc.inv.allele.freq.df <- as.data.frame(hgsvc.inv.allele.freq) 
  
plt.df <- data.frame(ID=as.character(morbMap.gold.gr.noSD),
           HGSVC=-rowSums(hgsvc.alleles.inv),
           NHP=rowSums(nhp.alleles.inv))
plt.df <- reshape2::melt(plt.df)
ggplot(plt.df) + 
  geom_col(aes(x=ID, y=value, fill=variable)) +
  scale_fill_manual(values = brewer.pal(n=3, name = 'Dark2')) +
  coord_flip() +
  theme_bw() +
  ylab("Inverted loci count") +
  xlab("Pathogenic CNVs")

  
## Get list of inverted regions in respect to CNV morbid list
inv.orang <- genotyped.regions[which(nhp.alleles.inv[,4] == 2)]
inv.orang.df <- as.data.frame(inv.orang)
destination <- file.path(outputDirectory, 'CNVregions_invertedInOrang.bed')
write.table(x = inv.orang.df, file = destination, quote = FALSE, row.names = FALSE)

## Plot hotspots genome-wide
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
seq.len <- seqlengths(bsgenome)[paste0('chr', c(1:22, 'X'))]
ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
ideo.df$seqnames <- factor(ideo.df$seqnames, levels=paste0('chr', c(1:22, 'X')))

plt.df <- as.data.frame(hotspots.gr)
plt1 <- ggplot() + geom_rect(data = ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill="white", color="black") +
  facet_grid(seqnames ~ ., switch = 'y') +
  geom_rect(data=plt.df , aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=num.events), fill='gray') +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text.y = element_text(angle = 180))

annot.df <- as.data.frame(morbMap.gold.gr.overlap.filt)
#annot.df <- as.data.frame(inversion.calls.filt.allmorbMap50)
annot.df$midpoint <- annot.df$start + ((annot.df$end - annot.df$start)/2)
annot.df$midpoint.toGR <- annot.df$toGR.start + ((annot.df$toGR.end - annot.df$toGR.start)/2)
plt1 <- plt1 + geom_point(data=annot.df, aes(x=midpoint, y=0.5, color=annot), inherit.aes = FALSE) +
  scale_color_manual(values = c('olivedrab1', 'darkorange1', 'deepskyblue3','darkorchid1')) +
  theme(legend.position = "bottom")
plt1 <- plt1 + geom_point(data=annot.df, aes(x=midpoint.toGR, y=1.5), shape=25, color='red', inherit.aes = FALSE)

## Make a summary plot of represented SV classes
plt2 <- annot.df %>% group_by(annot) %>% summarise(count=n()) %>%
  ggplot(aes(x=annot, y=count, fill=annot)) + 
  geom_col() +
  #scale_fill_manual(values = c('olivedrab1', 'darkorange1', 'deepskyblue3','darkorchid1')) +
  scale_fill_manual(values = gray.colors(3), guide='none') +
  geom_text(aes(label=count), vjust=0) +
  xlab("") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# plt3 <- as.data.frame(SimpleInversion.calls.morbMap50) %>% group_by(ID) %>% summarise(count=n()) %>%
#   ggplot(aes(x=ID, y=count, fill=ID)) + 
#   geom_col() +
#   scale_fill_manual(values = c('#3182bd','#31a354','#8856a7','#e6550d'), guide='none') +
#   scale_x_discrete(labels=c('B', 'C', 'G', 'O')) +
#   xlab("")

nhp.inv.allele.freq.df$ID <- rownames(nhp.inv.allele.freq.df)
plt4 <- ggplot(data=nhp.inv.allele.freq.df, aes(x=ID, y=nhp.inv.allele.freq, fill=ID)) + 
  geom_col() +
  scale_fill_manual(values = c('#3182bd','#31a354','#8856a7','#e6550d'), guide='none') +
  scale_x_discrete(labels=c('B', 'C', 'G', 'O')) +
  xlab("") +
  ylab("Inverted loci frequency over CN (n=16)")

#final.plt <- plot_grid(plt1, plt2, plt3, plt4, nrow = 1, rel_widths = c(4, 1.5, 1, 1))
#plt2 <- plt2 + coord_flip()
#plt4 <- plt4 + coord_flip()
final.plt <- plot_grid(plt1, plt2, plt4, nrow = 1, rel_widths = c(4,1,1), align = 'h', axis = 'b') 
## Export final plots
ggsave(filename = file.path(outputDirectory, 'morbMapOverlap_simpINV.pdf'), plot = final.plt, width = 10, height = 6, useDingbats=FALSE)
#ggsave(filename = file.path(outputDirectory, 'morbMapOverlap_simpINV_genomewide.pdf'), plot = plt1, width = 10, height = 6)

## Calculate enrichment of simple inversions overlapping with morbid CNVs ##
############################################################################
## Get reciprocal overlap of CNV regions against simple inversions
morbMap.gold.gr.overlap <- getReciprocalOverlaps(query = morbMap.gold.gr, subject = all.SimpleInversion.calls.filt, report = 'query')
observed.count <- length(morbMap.gold.gr.overlap[morbMap.gold.gr.overlap$perc.overlap >= 50])

## Create GRCh38 mask from centromeric regions and gaps
## Add centromere information
cent.df <- read.table("/home/porubsky/WORK/Saarclust_project/centromeres_GRCh38.bed")
cent.df <- cent.df[,c(2,3,4)]
colnames(cent.df) <- c('seqnames', 'start', 'end')
cent.gr <- makeGRangesFromDataFrame(cent.df)
## Add gap information
gaps.df <- read.table("/home/porubsky/WORK/Saarclust_project/gaps_GRCh38")
gaps.df <- gaps.df[,c(2,3,4)]
colnames(gaps.df) <- c('seqnames', 'start', 'end')
gaps.gr <- makeGRangesFromDataFrame(gaps.df)
hg38.mask <- sort(c(cent.gr, gaps.gr))

## Keep only simple inversion no longer than longest morbid CNV
## Not needed as we require 50 reciprocal overlap
#max.len.morbid <- max(width(morbMap.gold.gr))
#simpleInversion.calls.maxMorbid <- all.SimpleInversion.calls.filt[width(all.SimpleInversion.calls.filt) <= max.len.morbid]

## Check if simple inversions overlap more often with morbid CNVs
randomized.count <- list()
for (i in 1:100) {
  message("Working on iteration: ", i)
  ## Randomize morbid CNV regions
  morbMap.gold.gr.rand <- randomizeRanges(gr = morbMap.gold.gr, bsgenome = BSgenome.Hsapiens.UCSC.hg38, mask.gr = hg38.mask)
  ## Calculate overlaps of size selected inversions with randomized morbid CNVs
  overlap <- getReciprocalOverlaps(query = morbMap.gold.gr.rand, subject = all.SimpleInversion.calls.filt, report = 'query')
  random.count <- length(overlap[overlap$perc.overlap >= 50])
  randomized.count[[i]] <- random.count
}
## Calculate statistical significance
random.counts <- unlist(randomized.count, use.names = FALSE)
zscore <- (observed.count - mean(random.counts)) / sd(random.counts)
pvalue2sided <- 2 * pnorm(-abs(zscore))
## Calculate enrichment factor
enrich.fact <- observed.count / mean(random.counts)
