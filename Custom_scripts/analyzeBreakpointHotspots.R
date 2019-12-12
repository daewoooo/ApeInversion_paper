## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38) )

message("Analyzing breakpoint hotspots ...")

outputDirectory <- "/home/porubsky/WORK/Great_apes/Hotspots_breakpoints/"

# ## Load all called inversions
# all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
# ## Retain only simple inversions
# all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
# ## Retain only inverted duplications
# all.invertedDups.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']
# 
# ## Load HGSCV simple inversions
# hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
# ## Construct GRanges object
# hgsvc.simple.inv.gr <- GRanges(seqnames=hgsvc.simple.inv$seqnames, 
#                                ranges=IRanges(start = hgsvc.simple.inv$innerBP_start, end = hgsvc.simple.inv$innerBP_end)
# )
# ## Add genotype information
# mcols(hgsvc.simple.inv.gr) <- hgsvc.simple.inv[,c(10:ncol(hgsvc.simple.inv ))]
# ## Filter inversion with at least 1000kb of unique sequence
# hgsvc.simple.inv.gr <- getSegDupOverlaps(query.gr = hgsvc.simple.inv.gr, subject.gr = seg.dup.gr)
# hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[hgsvc.simple.inv.gr$TotalUniqueBases >= 1000]

## Load putative polymorphic inversions
polymorphic.inversions <- get(load("/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites/shared_HGSVC&NHP_sites.RData"))

# ## Merge NHP and HGSVC simple inversion calls
# all.SimpleInversion.calls.filt <- all.SimpleInversion.calls.filt[,0]
# all.SimpleInversion.calls.filt$ID <- 'NHP'
# hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[,0]
# hgsvc.simple.inv.gr$ID <- 'HGSVC'
# simpleINV.greatApes <- c(all.SimpleInversion.calls.filt, hgsvc.simple.inv.gr)
# 
# ## Load known misorients detected based on HGSVC data
# misorients.annot <- "/home/porubsky/WORK/Great_apes/Ashley_inversions/Misorients/"
# HGSVC.misorients <- read.table(file.path(misorients.annot, "misAssem2remove.NHP.csv"), sep=',', header = TRUE, stringsAsFactors = FALSE)
# HGSVC.misorients.gr <- GRanges(seqnames = HGSVC.misorients$seqnames, ranges=IRanges(start = HGSVC.misorients$start, end = HGSVC.misorients$end), toRemove = HGSVC.misorients$toRemove)
# HGSVC.misorients.gr <- HGSVC.misorients.gr[HGSVC.misorients.gr$toRemove == TRUE]

## Load predicted human specific inversions
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))

# ## Collapse ranges that overlaps with predicted human specific inversions
# overlaps2HS <- getReciprocalOverlaps(query = simpleINV.greatApes, subject = HSinv.gr, report = 'query')
# simpleINV.greatApes.filt <- overlaps2HS[!overlaps2HS$perc.overlap >= 50, 1]
# simpleINV.greatApes.filt <- c(simpleINV.greatApes.filt, HSinv.gr[,0])

# ## Remove known misorients/misassemblies
# overlaps2misO <- getReciprocalOverlaps(query = simpleINV.greatApes.filt, subject = HGSVC.misorients.gr, report = 'query')
# simpleINV.greatApes.filt <- overlaps2misO[!overlaps2misO$perc.overlap >= 50, 1]

# ## Assigne chromosome lengths
# seqlengths(simpleINV.greatApes.filt) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(simpleINV.greatApes.filt)]

## Test hotspots based on simple inversion breakpoints in non-redundant dataset
hgsvc.simple.inv.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/HGSVC.nonred.filt.RData"))
hgsvc.simple.inv.gr$ID <- "HGSVC"
NHP.nonred.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/NHP.nonred.filt.RData"))
all.simple.INV.gr <- sort(c(hgsvc.simple.inv.gr[,'ID'], NHP.nonred.gr[,'ID']))
# simpleINV.greatApes <- getDisjointOverlapsWeighted(simpleINV.greatApes, percTh = 50)
# simpleINV.greatApes.collapsed <- collapseBins(simpleINV.greatApes, id.field = 6)
simpleINV.breakpoints <- getRegionBoundaries(all.simple.INV.gr)
seqlengths(simpleINV.breakpoints) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(simpleINV.breakpoints)]
simpleINV.breakpoints <- sort(simpleINV.breakpoints)
hotspots <- primatR::hotspotter(gr = simpleINV.breakpoints[,0], bw = 2000000, pval = 5e-10)
#hotspots <- primatR::hotspotter(gr = sort(test.regions[,0]), bw = 1000000, pval = 5e-10)

## Test hotspots in redundant dataset
#hotspots <- primatR::hotspotter(gr = simpleINV.greatApes.filt[,0], bw = 1000000, pval = 5e-10)
#hotspots$ID <- 'break.hotspots'

## Plot hotspots genome-wide
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
seq.len <- seqlengths(bsgenome)[paste0('chr', c(1:22, 'X'))]
ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
ideo.df$seqnames <- factor(ideo.df$seqnames, levels=paste0('chr', c(1:22, 'X')))

plt.df <- as.data.frame(hotspots)
plt <- ggplot() + geom_rect(data = ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill="white", color="black") +
  facet_grid(seqnames ~ ., switch = 'y') +
  geom_rect(data=plt.df , aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=num.events)) +
  scale_fill_gradient(low = "dodgerblue1", high = "red") +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text.y = element_text(angle = 180))

## Export results
destination <- file.path(outputDirectory, "hotspots.RData")
save(hotspots, file = destination)
destination <- file.path(outputDirectory, "hotspots_genomeWidePlot.pdf")
ggsave(plt, filename = destination, width = 10, height = 6, useDingbats=FALSE)

## Add position of likely polymorphic inversions
polymorphic.inversions.df <- as.data.frame(polymorphic.inversions)
polymorphic.inversions.df$midpoint <- polymorphic.inversions.df$start + ((polymorphic.inversions.df$end - polymorphic.inversions.df$start)/2)
polymorphic.inversions.df$ID <- 'sharedINV'
annot.df <- polymorphic.inversions.df

## Get number of hotspots overlapping with polymorphic inversions
counts <- countOverlaps(hotspots, polymorphic.inversions)
counts <- counts[counts > 0]

plt <- plt + geom_point(data=annot.df, aes(x=midpoint, y=0.5, color=ID), inherit.aes = FALSE) +
  scale_color_manual(values = c('chartreuse3', 'darkgoldenrod3'))

destination <- file.path(outputDirectory, "hotspots_genomeWidePlot_annot.pdf")
ggsave(plt, filename = destination, width = 10, height = 6, useDingbats=FALSE)


message("DONE!!!")
