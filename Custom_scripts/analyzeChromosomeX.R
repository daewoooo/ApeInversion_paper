## Chromosome X polymorphic events ##
#####################################
## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(UpSetR) )
suppressPackageStartupMessages( library(reshape2) )
suppressPackageStartupMessages( library(biovizBase) )
suppressPackageStartupMessages( library(ggbio) )
suppressPackageStartupMessages( library(gtools) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38) )

message("Analyzing Chromosome X inversions ...")

outputDirectory <- "/home/porubsky/WORK/Great_apes/ChromosomeX/"

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load Strand-seq simple INV calls ##
######################################
all.SimpleInversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.SimpleInversion.calls.filt.annot.RData"))
## Remove predicted misorients !!!
all.SimpleInversion.calls.filt <- all.SimpleInversion.calls.filt[all.SimpleInversion.calls.filt$misorient == FALSE]

## Load HGSCV simple inversions ##
##################################
hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
## Construct GRanges object
hgsvc.simple.inv.gr <- GRanges(seqnames=hgsvc.simple.inv$seqnames, 
                               ranges=IRanges(start = hgsvc.simple.inv$innerBP_start, end = hgsvc.simple.inv$innerBP_end)
)
## Add genotype information
mcols(hgsvc.simple.inv.gr) <- hgsvc.simple.inv[,c(10:ncol(hgsvc.simple.inv ))]
## Filter inversion with at least 1000kb of unique sequence
hgsvc.simple.inv.gr <- getSegDupOverlaps(query.gr = hgsvc.simple.inv.gr, subject.gr = seg.dup.gr)
hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[hgsvc.simple.inv.gr$TotalUniqueBases >= 1000]

## Load ape's composite files ##
################################
chimp.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
chimp.data$ID <- 'chimpanzee'
bonobo.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
bonobo.data$ID <- 'bonobo'
gorilla.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
gorilla.data$ID <- 'gorilla'
orangutan.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
orangutan.data$ID <- 'orangutan'

## Load predicted hotspots
hotspots <- get(load("/home/porubsky/WORK/Great_apes/Hotspots_breakpoints/hotspots.RData"))
hotspots.chrX <- hotspots[seqnames(hotspots) == 'chrX']

## Remove putative misorients from HGSVC calls !!!
misorients.annot <- "/home/porubsky/WORK/Great_apes/Ashley_inversions/Misorients/"
HGSVC.misorients <- read.table(file.path(misorients.annot, "misAssem2remove.NHP.csv"), sep=',', header = TRUE, stringsAsFactors = FALSE)
HGSVC.misorients.gr <- GRanges(seqnames = HGSVC.misorients$seqnames, ranges=IRanges(start = HGSVC.misorients$start, end = HGSVC.misorients$end), toRemove = HGSVC.misorients$toRemove)
HGSVC.misorients.gr <- HGSVC.misorients.gr[HGSVC.misorients.gr$toRemove == TRUE]
hits <- suppressWarnings( findOverlaps(hgsvc.simple.inv.gr, HGSVC.misorients.gr) )
hgsvc.simple.inv.gr$misorient <- FALSE
hgsvc.simple.inv.gr$misorient[queryHits(hits)] <- TRUE
hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[hgsvc.simple.inv.gr$misorient == FALSE]

## Check for different inversion hotspots ##
############################################
## Keep only chromosome X
chrX.nhp.simpleInv <- all.SimpleInversion.calls.filt[seqnames(all.SimpleInversion.calls.filt) == 'chrX']
chrX.nhp.simpleInv.regions <- reduce(chrX.nhp.simpleInv)
chrX.hgsvc.simpleInv <- hgsvc.simple.inv.gr[seqnames(hgsvc.simple.inv.gr) == 'chrX']
hgsvc.genotypes <- mcols(chrX.hgsvc.simpleInv[,1:9])

## Genotype inverted regions on ChrX in NHP ##
##############################################
## Set ranges to be re-genotyped
regions.to.genotype <- chrX.nhp.simpleInv.regions[,0]
regions.to.genotype <- subtractRegions(regions.to.genotype, remove.gr = reduce(seg.dup.gr), mode = 'flanks')
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

nhp.genotypes <- mcols(genotyped.regions)[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]

## Recode genotypes to number in order to calculate allele frequencies
hgsvc.alleles.inv <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
hgsvc.alleles.dir <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
nhp.alleles.inv <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
nhp.alleles.dir <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
## Calculate allele freaquencies for HGSVC and NHP genotypes
hgsvc.total.inv.alleles <- rowSums(hgsvc.alleles.inv)
hgsvc.total.dir.alleles <- rowSums(hgsvc.alleles.dir)
nhp.total.inv.alleles <- rowSums(nhp.alleles.inv)
nhp.total.dir.alleles <- rowSums(nhp.alleles.dir)
hgsvc.inv.allele.freq <- hgsvc.total.inv.alleles / (hgsvc.total.inv.alleles + hgsvc.total.dir.alleles)
nhp.inv.allele.freq <- nhp.total.inv.alleles / (nhp.total.inv.alleles + nhp.total.dir.alleles)

## Plot allele frequencies
hgsvc.df <- as.data.frame(chrX.hgsvc.simpleInv[,0])
hgsvc.df$ID <- 'HGSVC'
nhp.df <- as.data.frame(chrX.nhp.simpleInv.regions[,0])
nhp.df$ID <- 'NHP'
hgsvc.df$allele.freq <- hgsvc.inv.allele.freq
nhp.df$allele.freq <- nhp.inv.allele.freq * -1
plt.df <- rbind(hgsvc.df, nhp.df)
plt.df$midpoint <- plt.df$start + ((plt.df$end - plt.df$start) / 2)
plt.df$ymin <- rep(c(0.1, -0.1), table(plt.df$ID))
plt.df$value <- plt.df$allele.freq + plt.df$ymin

# chrX.hgsvc.simple.inv.gr$count <- apply(hgsvc.gens, 1, function(x) length(x[unlist(x) != 'REF']))
# chrX.hgsvc.simple.inv.gr$count.norm <- chrX.hgsvc.simple.inv.gr$count / 9
# chrX.hgsvc.simple.inv.gr$ID <- 'HGSVC'
# ## Count overlapping inversion in NHP
# disjoin.ranges <- disjoin(chrX.nhp.simpleInv)
# disjoin.ranges$count <- countOverlaps(disjoin.ranges, chrX.nhp.simpleInv)
# disjoin.ranges$count.norm <- disjoin.ranges$count / 4
# disjoin.ranges$count.norm <- disjoin.ranges$count.norm * -1
# disjoin.ranges$ID <- 'NHP'
# ## Prepare data for plotting
# plt.gr <- c(chrX.hgsvc.simple.inv.gr[,c('count.norm', 'ID')], disjoin.ranges[,c('count.norm', 'ID')])
# plt.df <- as.data.frame(plt.gr)
# plt.df$midpoint <- plt.df$start + ((plt.df$end - plt.df$start) / 2)
# plt.df$ymin <- rep(c(0.1, -0.1), table(plt.df$ID))

## Make banded ideogram ##
##########################
suppressMessages( chrX.hg38IdeogramCyto <- getIdeogram("hg38", cytobands = TRUE, subchr = 'chrX') )
ideo.df <- as.data.frame(chrX.hg38IdeogramCyto)
ideo <- ggplot() + 
  geom_rect(ideo.df, aes(xmin=start, xmax=end, ymin=-0.1, ymax=0.1, fill=gieStain), color='black', show.legend = FALSE) + 
  scale_fill_giemsa() +
  scale_x_continuous(labels = comma) +
  xlab("Genomic Position (bp)") +
  ylab("Inverted loci frequency") +
  ggtitle("Chromosome X inversion distribution")
## Add inversion counts bars
#plt <- ideo + geom_linerange(data=plt.df, aes(x=midpoint, ymin=ymin, ymax=count.norm, color=ID), inherit.aes = FALSE) +
plt <- ideo + geom_linerange(data=plt.df, aes(x=midpoint, ymin=ymin, ymax=value, color=ID), inherit.aes = FALSE) +
  scale_color_manual(values = c('sienna1', 'slateblue3'))
## Add position of predicted hotspots
hotspots.chrX.df <- as.data.frame(hotspots.chrX)
plt <- plt + geom_rect(data=hotspots.chrX.df, aes(xmin=start, xmax=end, ymin=1, ymax=1.1), fill='gray')
## Save final plot
ggsave(filename = file.path(outputDirectory, 'chromosomeX_invDistribution.pdf'), plot = plt, width = 12, height = 4, useDingbats=FALSE)

## Add ChromosomeX strata
#chrX.strata <- read.table("/home/porubsky/WORK/Great_apes/ChromosomeX/pandey_etal_2013_chrXstrata.csv", header = TRUE, sep = ',')
#plt2 <- ggplot() + geom_rect(chrX.strata, aes(xmin=start, xmax=end, ymin=1, ymax=1.1, fill=factor(stratum))) +
#  scale_fill_manual(values = brewer.pal(n = 12, name = 'Set3'), name="Evol_Stratum")

## Load Recomb rates for chromosome X
#chrX.recomb.hg19 <- read.table("/home/porubsky/WORK/Great_apes/Annotations/chrX_recombRate_hg19.bed", header=TRUE)
#chrX.recomb.hg38pos <- read.table("/home/porubsky/WORK/Great_apes/Annotations/chrX_recombRate_hg19tohg38_posOnly.bed")
#chrX.recomb.gr <- GRanges(seqnames=chrX.recomb.hg38pos$V1, ranges=IRanges(start=chrX.recomb.hg38pos$V2, end=chrX.recomb.hg38pos$V3))
#chrX.recomb.gr$decodeFemale <- chrX.recomb.hg19$decodeAvg
#chrX.recomb.gr$fill <- findInterval(chrX.recomb.gr$decodeFemale, c(1,2,3,4))

#chrX.recomb.df <- as.data.frame(chrX.recomb.gr)
#plt3 <- ggplot() + geom_rect(chrX.recomb.df, aes(xmin=start, xmax=end, ymin=1, ymax=1.1, fill=fill)) +
#  scale_fill_gradient(low = 'white', high = 'black')

#plot_grid(plt3, plt, plt2, ncol = 1, align = 'v', rel_heights = c(0.5, 3, 0.5))


## Regenotype human inversions in non-human primates
## Genotype regions from NA19240 ##
###################################
## Get common set of non-overlapping regions
chrX.regions <- reduce(c(chrX.hgsvc.simpleInv[,0], chrX.nhp.simpleInv[,0]))
## Set ranges to be re-genotyped
regions.to.genotype <- chrX.regions[,0]
regions.to.genotype <- subtractRegions(regions.to.genotype, remove.gr = reduce(seg.dup.gr), mode = 'flanks')
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

## Compile genotypes for all individuals ##
chrX.genotypes <- genotyped.regions
mcols(chrX.genotypes) <- mcols(chrX.genotypes )[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]
names(mcols(chrX.genotypes)) <- gsub(names(mcols(chrX.genotypes)), pattern = 'genoT_', replacement = '')
## Add HGSVC genotypes
hits <- findOverlaps(chrX.genotypes, chrX.hgsvc.simpleInv)
hgsvc.gens <- matrix(rep('REF', length(chrX.genotypes) * 9, ncol = 9), nrow = length(chrX.genotypes))
colnames(hgsvc.gens) <- names(mcols(chrX.hgsvc.simpleInv)[1:9])
hgsvc.gens.df <- data.frame(hgsvc.gens, stringsAsFactors = FALSE)
hgsvc.gens.df[queryHits(hits),] <- as.data.frame(mcols(chrX.hgsvc.simpleInv)[1:9])
mcols(chrX.genotypes) <- cbind(mcols(chrX.genotypes), hgsvc.gens.df)

## Plot genotypes
plt.df <- as.data.frame(mcols(chrX.genotypes))
plt.df$ID <- as.character(chrX.genotypes[,0])
plt.df$midpoint <- start(chrX.genotypes) + ((end(chrX.genotypes) - start(chrX.genotypes))/2)
plt.df <- melt(plt.df, id.vars = c('ID','midpoint'))
plt.df$ID <- factor(plt.df$ID, as.character(chrX.genotypes[,0]))

plt2 <- ggplot(plt.df) + 
  geom_tile(aes(x=ID, y=variable, fill=value)) +
  scale_fill_manual(values = c("cornflowerblue","coral3", "white", "gray"), name="Genotype") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  xlab("") +
  ylab("")
## Save final plot
ggsave(filename = file.path(outputDirectory, 'genotChrXInvertedRegions.pdf'), plot = plt2, width = 10, height = 5, useDingbats=FALSE)


## Get fold enrichment of ChrX inversions in comparison to autosomes ##
#######################################################################
## Load inversion calls
data <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
simple.inv <- data[data$SVclass == 'INV']

## Get chromosome lengths
chr.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, 'X'))]
## Get inversion counts
inv.counts <- data.frame(chr = as.character(runValue(seqnames(simple.inv))), count = runLength(seqnames(simple.inv)), stringsAsFactors = FALSE)
inv.counts <- inv.counts[mixedorder(inv.counts$chr),]
## Add chromosome length
inv.counts$chr.len <- chr.len[names(chr.len) %in% inv.counts$chr]
## normalize inv.counts (INV.count per chr.length per Mb)
inv.counts$count.norm <- (inv.counts$count / inv.counts$chr.len) * 1000000
## Separate counts for autosomes and chrX
inv.counts.autosomes <- inv.counts[inv.counts$chr %in% paste0('chr', c(1:22)),]
inv.counts.chrX <- inv.counts[inv.counts$chr == 'chrX',]

## Get fold enrichment
mean.enrich <- inv.counts.chrX$count / mean(inv.counts.autosomes$count)
mean.enrich.norm <- inv.counts.chrX$count.norm / mean(inv.counts.autosomes$count.norm)

## create null expectation using only autosome counts to convert to normalised counts (i.e., z scores)
zscores.aut.count <- ( (inv.counts.autosomes$count - mean(inv.counts.autosomes$count)) / sd(inv.counts.autosomes$count) )
zscores.aut.count.norm <- ( (inv.counts.autosomes$count.norm - mean(inv.counts.autosomes$count.norm)) / sd(inv.counts.autosomes$count.norm) )

## z-score for chrX with respect to autosomes
zscore.chrX <- ( (inv.counts.chrX$count - mean(inv.counts.autosomes$count)) / sd(inv.counts.autosomes$count) )
zscore.chrX.norm <- ( (inv.counts.chrX$count.norm - mean(inv.counts.autosomes$count.norm)) / sd(inv.counts.autosomes$count.norm) )

## Convert to p-value (single-sided)
pvalue = pnorm(-abs(zscore.chrX))
pvalue.norm = pnorm(-abs(zscore.chrX.norm))

## Final answer
print(paste("    Enrichment factor is", mean.enrich.norm,"and P-value =", pvalue.norm))


## Number of HETs on chromosome X in comparison to autosomes ##
###############################################################
## Keep only HET inversions
simple.inv.het <- simple.inv[simple.inv$gen == 'HET']
## Get chromosome lengths
chr.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, 'X'))]
## Get inversion counts
inv.counts <- data.frame(chr = as.character(runValue(seqnames(simple.inv.het))), count = runLength(seqnames(simple.inv.het)), stringsAsFactors = FALSE)
inv.counts <- inv.counts[mixedorder(inv.counts$chr),]
## Add chromosome length
inv.counts$chr.len <- chr.len[names(chr.len) %in% inv.counts$chr]
## normalize inv.counts (INV.count per chr.length per Mb)
inv.counts$count.norm <- (inv.counts$count / inv.counts$chr.len) * 1000000
## Separate counts for autosomes and chrX
inv.counts.autosomes <- inv.counts[inv.counts$chr %in% paste0('chr', c(1:22)),]
inv.counts.chrX <- inv.counts[inv.counts$chr == 'chrX',]

## Get fold enrichment
mean.enrich <- inv.counts.chrX$count / mean(inv.counts.autosomes$count)
mean.enrich.norm <- inv.counts.chrX$count.norm / mean(inv.counts.autosomes$count.norm)

## create null expectation using only autosome counts to convert to normalised counts (i.e., z scores)
zscores.aut.count <- ( (inv.counts.autosomes$count - mean(inv.counts.autosomes$count)) / sd(inv.counts.autosomes$count) )
zscores.aut.count.norm <- ( (inv.counts.autosomes$count.norm - mean(inv.counts.autosomes$count.norm)) / sd(inv.counts.autosomes$count.norm) )

## z score for chrX with respect to autosomes
zscore.chrX <- ( (inv.counts.chrX$count - mean(inv.counts.autosomes$count)) / sd(inv.counts.autosomes$count) )
zscore.chrX.norm <- ( (inv.counts.chrX$count.norm - mean(inv.counts.autosomes$count.norm)) / sd(inv.counts.autosomes$count.norm) )

## Convert to p-value (single-sided)
pvalue = pnorm(-abs(zscore.chrX))
pvalue.norm = pnorm(-abs(zscore.chrX.norm))

## Final answer
#print(paste("    Enrichment factor is", mean.enrich.norm,"and P-value =", pvalue.norm))

message("DONE!!!")