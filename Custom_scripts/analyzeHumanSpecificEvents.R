## Load required libraries ##
#############################
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(gtools) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38) )

outputDirectory <- "/home/porubsky/WORK/Great_apes/Human_specific_events/"

message("Searching for human specific inversions ...")

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Read in all called inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))

## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']

## Load manually selected human specific inversion [outdated!!!]
#HSinv <- read.table(file.path(outputDirectory, "manuallySelected_humanSpecific.bed"))
#HSinv.gr <- GRanges(seqnames=HSinv$V1, ranges=IRanges(start=HSinv$V2, end=HSinv$V3), SVclass=HSinv$V4, Note=HSinv$V5)
#HSinv.gr <- HSinv.gr[HSinv.gr$Note == 'HS' | HSinv.gr$Note == 'HS?']
#HSinv.gr <- HSinv.gr[HSinv.gr$Note == 'varGen'| HSinv.gr$Note == 'near']
#HSinv.gr <- keepSeqlevels(HSinv.gr, value = unique(seqnames(HSinv.gr)), pruning.mode = 'coarse')

## Get human specific inversions ##
###################################
simpleInversion.comp <- getDisjointOverlapsWeighted(gr = all.SimpleInversion.calls.filt, percTh = 50) #required 50% reciprocal overlap
simpleInversion.comp.grl <- split(simpleInversion.comp, simpleInversion.comp$sub.group)
mask <- which(lengths(simpleInversion.comp.grl) >= 4) #Select regions that appear in all 4 individuals
human.specific.grl <- simpleInversion.comp.grl[mask]
## Filter out loci with various genotypes (Likely polymorphic)
## Keep only loci where all NHP are HOMs
mask <- sapply(human.specific.grl, function(x) all(x$gen == 'HOM'))
human.specific.grl <- human.specific.grl[mask]
human.specific.gr <- unlist(human.specific.grl, use.names = FALSE)
HSinv.plt <- genomewideRangesIdeo(gr = human.specific.gr, userTrackGeom = 'point', bsgenome = BSgenome.Hsapiens.UCSC.hg38)
destination <- file.path(outputDirectory, "putative_humanSpecificEvents_genomewide.pdf")
ggsave(HSinv.plt, filename = destination, width = 12, height = 7, limitsize = FALSE, useDingbats=FALSE)
## Reduce likely human specific inversions into a reduced regions
human.specific.regions <- reduce(human.specific.gr)

## Check putative human specific inversions against known misorients ##
#######################################################################
misorients.annot <- "/home/porubsky/WORK/Great_apes/Ashley_inversions/Misorients/"
HGSVC.misorients <- read.table(file.path(misorients.annot, "misAssem2remove.NHP.csv"), sep=',', header = TRUE, stringsAsFactors = FALSE)
HGSVC.misorients.gr <- GRanges(seqnames = HGSVC.misorients$seqnames, ranges=IRanges(start = HGSVC.misorients$start, end = HGSVC.misorients$end), toRemove = HGSVC.misorients$toRemove)
HGSVC.misorients.gr <- HGSVC.misorients.gr[HGSVC.misorients.gr$toRemove == TRUE]

## Plot putative human specific INV against known genome misorients
#HSinvPlusMisorients.plt <- genomewideRangesIdeo(gr = human.specific.gr, userTrack = misorients.gr, userTrackGeom = 'point', bsgenome = BSgenome.Hsapiens.UCSC.hg38)
#destination <- file.path(outputDirectory, "putative_humanSpecificEvents_vs_knownMisorients_genomewide.pdf")
#ggsave(HSinvPlusMisorients.plt, filename = destination, width = 12, height = 7, limitsize = FALSE)

## Check putative human specific inversions against HGSVC data ##
#################################################################
## Load HGSCV data
hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
## Construct GRanges object
hgsvc.simple.inv.gr <- GRanges(seqnames=hgsvc.simple.inv$seqnames, 
                                       ranges=IRanges(start = hgsvc.simple.inv$innerBP_start, end = hgsvc.simple.inv$innerBP_end))
## Get overlaps between putative human specific inversion and HGSCV polymorphic inversions
HSinv.vs.HGSVC <- getReciprocalOverlaps(query = human.specific.regions, subject = hgsvc.simple.inv.gr, thresh = 50, report = 'query', index = 'HGSVC')

## Annotate putative human specific inversion as known misorients as been published on ##
#########################################################################################
## Get info on alread published regions from Stuart's table
stuart.annot <- read.table(file.path(outputDirectory, "human_spec_inversions_annotTable.csv"), sep = ',', header = TRUE, stringsAsFactors = FALSE)
stuart.annot.gr <- GRanges(seqnames = stuart.annot$chr, ranges=IRanges(start = stuart.annot$start, end = stuart.annot$end), note = stuart.annot$Reported.Previously, comments = stuart.annot$Comments)
published.gr <- stuart.annot.gr[stuart.annot.gr$comments != 'reference error' & nchar(stuart.annot.gr$note) > 0]
## Annotate predicted human specific regions
hits <- suppressWarnings( findOverlaps(human.specific.regions, HGSVC.misorients.gr) )
human.specific.regions$misorient <- FALSE
human.specific.regions$misorient[queryHits(hits)] <- TRUE
hits <- suppressWarnings( findOverlaps(human.specific.regions, published.gr) )
human.specific.regions$published <- FALSE
human.specific.regions$published[queryHits(hits)] <- TRUE
human.specific.regions$HGSVCpolymorph <- FALSE
human.specific.regions$HGSVCpolymorph[HSinv.vs.HGSVC$perc.overlap_HGSVC > 50] <- TRUE #Mark HS events that are polymorphic in humans

## Plot putative human specific mark published and misoriented regions
known.or.misorients.gr <- human.specific.regions[human.specific.regions$misorient == TRUE | human.specific.regions$published == TRUE]
unknown.HSinv.plt <- genomewideRangesIdeo(gr = human.specific.gr, userTrack = known.or.misorients.gr, userTrackGeom = 'point', bsgenome = BSgenome.Hsapiens.UCSC.hg38, colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
destination <- file.path(outputDirectory, "putative_humanSpecificEvents_vs_knownEventsAndMisorients_genomewide.pdf")
ggsave(unknown.HSinv.plt, filename = destination, width = 12, height = 7, limitsize = FALSE, useDingbats=FALSE)

## Export annotated human specific regions
human.specific.regions.df <- as.data.frame(human.specific.regions)
destination <- file.path(outputDirectory, "putative_humanSpecific_annot.txt")
write.table(as.data.frame(human.specific.regions), file = destination, quote = FALSE, row.names = FALSE)
destination <- file.path(outputDirectory, "putative_humanSpecific_annot.RData")
save(human.specific.regions, file = destination)
destination <- file.path(outputDirectory, "putative_humanSpecific_annot.bed")
write.table(as.data.frame(human.specific.regions)[,c(1,2,3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
## Export true human specific inversions
human.specific.regions.trueSet <- human.specific.regions[human.specific.regions$misorient == FALSE & human.specific.regions$HGSVCpolymorph == FALSE]
destination <- file.path(outputDirectory, "human.specific.regions.trueSet.RData")
save(human.specific.regions.trueSet, file = destination)
ranges2UCSC(gr = human.specific.regions.trueSet, outputDirectory = outputDirectory, index = "predicted_humanSpecificInversions", colorRGB = "80,34,25")

## Plot venn of predicted human specific inversions ##
## This part has to be run manually!!!
## Load Vennerable package first!!!
venn.df <- mcols(human.specific.regions)
## Prepare list of shared rows
overlaps.list <- list()
overlaps.list[['misorients']] <- which(venn.df$misorient == TRUE)
overlaps.list[['published']] <- which(venn.df$published == TRUE)
overlaps.list[['HGSVCpolymorph']] <- which(venn.df$HGSVCpolymorph == TRUE)
overlaps.list[['HumanSpecific']] <- 1:nrow(venn.df)
## Make Venn using Vennerable package [run Manually]
suppressPackageStartupMessages( library(Vennerable) )
venn.data <- Venn(overlaps.list)
destination <- file.path(outputDirectory, "HSinv_venn.pdf")
pdf(destination, useDingbats=FALSE)
plot(venn.data, type='squares')
dev.off()

## Mark misorients in the final Ape inversion callset ##
########################################################
ranges2mark <- subsetByOverlaps(human.specific.gr, human.specific.regions[human.specific.regions$misorient == TRUE])
all.SimpleInversion.calls.filt.annot <- all.SimpleInversion.calls.filt
all.SimpleInversion.calls.filt.annot$misorient <- FALSE
all.SimpleInversion.calls.filt.annot$misorient[ranges2mark$idx] <- TRUE
destination <- "/home/porubsky/WORK/Great_apes/Final_INV_calls/all.SimpleInversion.calls.filt.annot.RData"
save(all.SimpleInversion.calls.filt.annot, file = destination)

## Extract and plot reads for each putative human specific event ##
###################################################################
## Resize manually selected HS inversions to twice of their size
HSinv.regions <- resizeRanges(human.specific.regions, times = 2, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
## Import reads from manually selected HS inversions
composite.folder <- "/home/porubsky/WORK/Great_apes/Composite_files"
chimp.region.reads <- importReadsFromComposite(compositeFile = file.path(composite.folder, "syncReads_chimpanzee_final.RData"), regions = HSinv.regions, ID = 'chimpanzee')
bonobo.region.reads <- importReadsFromComposite(compositeFile = file.path(composite.folder, "syncReads_bonobo_final.RData"), regions = HSinv.regions, ID = 'bonobo')
gorilla.region.reads <- importReadsFromComposite(compositeFile = file.path(composite.folder, "syncReads_gorilla_final.RData"), regions = HSinv.regions, ID = 'gorilla')
orangutan.region.reads <- importReadsFromComposite(compositeFile = file.path(composite.folder, "syncReads_orangutan_final.RData"), regions = HSinv.regions, ID = 'orangutan')
human.region.reads <- importReadsFromComposite(compositeFile = file.path(composite.folder, "syncFrags_NA19240.RData"), regions = HSinv.regions, ID = 'NA19240')
## Export disjoint position of every single read for plotting
chimp.region.reads.grl <- coveragePerRegion(grl = split(chimp.region.reads, chimp.region.reads$region.ID), ID = 'chimpanzee', coverage = FALSE)
bonobo.region.reads.grl <- coveragePerRegion(grl = split(bonobo.region.reads, bonobo.region.reads$region.ID), ID = 'bonobo', coverage = FALSE)
gorilla.region.reads.grl <- coveragePerRegion(grl = split(gorilla.region.reads, gorilla.region.reads$region.ID), ID = 'gorilla', coverage = FALSE)
orangutan.region.reads.grl <- coveragePerRegion(grl = split(orangutan.region.reads, orangutan.region.reads$region.ID), ID = 'orangutan', coverage = FALSE)
human.region.reads.grl <- coveragePerRegion(grl = split(human.region.reads, human.region.reads$region.ID), ID = 'NA19240', coverage = FALSE)
chimp.plt.gr <- unlist(chimp.region.reads.grl)
bonobo.plt.gr <- unlist(bonobo.region.reads.grl)
gorilla.plt.gr <- unlist(gorilla.region.reads.grl)
orangutan.plt.gr <- unlist(orangutan.region.reads.grl)
human.plt.gr <- unlist(human.region.reads.grl)
## Prepare data.frame for plotting
plt.gr <- c(chimp.plt.gr, bonobo.plt.gr, gorilla.plt.gr, orangutan.plt.gr, human.plt.gr)
plt.df <- as(plt.gr, 'data.frame') 
## Construct a plot
plt <- ggplot(plt.df) + 
          geom_linerange(aes(ymin=start, ymax=end, x=level, color=strand), size=3) + 
          coord_flip() + 
          facet_wrap(ID ~ region.ID, scales = 'free', ncol = length(HSinv.regions)) +
          scale_color_manual(values = c("#678B8B","#F3A561"))
## Save the plot
destination <- file.path(outputDirectory, "humanSpecificEvents_readCoveragePlot.pdf")
ggsave(plt, filename = destination, width = 80, height = 10, limitsize = FALSE, useDingbats=FALSE)

message("DONE!!!")

## Plot distance of HS inversions to enhancers and differentialy expressed genes ##
###################################################################################
if (FALSE) {
  ## Load predicted HS invertsions  
  HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
  ## Load geneHancer track
  geneHancer <- read.table("/home/porubsky/WORK/Great_apes/Annotations/geneHancerRegElements_GRCh38.bed.gz", header=FALSE, stringsAsFactors = FALSE)
  geneHancer.gr <-GRanges(seqnames=geneHancer$V1, ranges=IRanges(start=geneHancer$V2, end=geneHancer$V3), name=geneHancer$V4, type=geneHancer$V11, categ=geneHancer$V12)
  ## Keep only enhancer regions
  geneHancer.gr <- geneHancer.gr[grep(geneHancer.gr$type, pattern = 'Enhancer')]
  ## Keep only ELITE enhancers
  #geneHancer.gr <- geneHancer.gr[geneHancer.gr$categ == 'Elite']
  
  ## Load Alex's gene list
  enriched.gene.list <- get(load("/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/upregul.gene.list.RData"))
  to.collapse <- GRangesList()
  for (i in seq_along(enriched.gene.list)) {
    ID <- names(enriched.gene.list[i])
    ## Skip over CONCERTED gene list
    if (ID == 'CONCERTED') {next}
    gr <- unlist(enriched.gene.list[[i]], use.names = FALSE)
    gr$ID <- paste0(ID, "_", gr$upregul.in)
    to.collapse[[length(to.collapse) + 1]] <- gr
  }
  enriched.genes.gr <- unlist(to.collapse, use.names = FALSE)
  ## Get overlap with Alex's gene list and gene enhancers +/- 1kb
  hits.upHuman <- findOverlaps(HSinv.gr, enriched.genes.gr[enriched.genes.gr$upregul.in == 'human'])
  hits.upChimp <- findOverlaps(HSinv.gr, enriched.genes.gr[enriched.genes.gr$upregul.in == 'human'])
  hits.enhancer <- findOverlaps(HSinv.gr, geneHancer.gr)
  dist.to.upHuman <- getMinDist(gr = HSinv.gr, userTrack = enriched.genes.gr[enriched.genes.gr$upregul.in == 'human'])
  dist.to.upChimp <- getMinDist(gr = HSinv.gr, userTrack = enriched.genes.gr[enriched.genes.gr$upregul.in == 'chimp'])
  dist.to.enhancer <- getMinDist(gr = HSinv.gr, userTrack = geneHancer.gr)
  ## Export HS inversion with distance to above mentioned features
  HSinvs.dist2feature.gr <- HSinv.gr[,0]
  ## Do not use Pollen's data
  #HSinvs.dist2feature.gr$dist.to.upHuman <- dist.to.upHuman
  #HSinvs.dist2feature.gr$dist.to.upChimp <- dist.to.upChimp
  HSinvs.dist2feature.gr$dist.to.enhancer <- dist.to.enhancer
  ## Set regions that overlap directly with a feature to zero
  #HSinvs.dist2feature.gr$dist.to.upHuman[unique(queryHits(hits.upHuman))] <- 0
  #HSinvs.dist2feature.gr$dist.to.upChimp[unique(queryHits(hits.upChimp))] <- 0
  HSinvs.dist2feature.gr$dist.to.enhancer[unique(queryHits(hits.enhancer))] <- 0
  ## Export results
  destination <- file.path(outputDirectory, "human.specific.regions.trueSet.dist2feature.RData")
  save(HSinvs.dist2feature.gr, file = destination)
  ## Sort by chromosome
  seqlevels(HSinvs.dist2feature.gr) <- mixedsort(seqlevels(HSinvs.dist2feature.gr))
  HSinvs.dist2feature.gr <- sort(HSinvs.dist2feature.gr)
  ## Prepare data for plotting
  plt.df <- as.data.frame(HSinvs.dist2feature.gr)
  ## Sort by chromosome
  plt.df$ID <- factor(as.character(HSinvs.dist2feature.gr[,0]), levels = as.character(HSinvs.dist2feature.gr[,0]))
  #plt.df <- reshape2::melt(plt.df, measure.vars=c('dist.to.upHuman','dist.to.upChimp','dist.to.enhancer'))
  plt.df <- reshape2::melt(plt.df, measure.vars='dist.to.enhancer')
  plt3 <- ggplot(plt.df, aes(x=ID, y=value, color=variable, group=variable)) + 
    geom_point(position=position_dodge(width=0.5), size=3) +
    scale_color_manual(values=brewer.pal(n = 4, name = "Dark2"), guide="none") +
    scale_y_continuous(breaks=c(1000, 10000, 100000, 1000000), labels = comma, trans = 'log10') +
    geom_hline(yintercept = c(1000, 10000, 100000, 1000000), linetype="dashed") +
    #theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab("") +
    ylab("Log10 Distance to enhancer (bp)") +
    coord_flip()
  ## Save the plot
  destination <- file.path(outputDirectory, "humanSpecificEvents_dist2features.pdf")
  ggsave(plt3, filename = destination, width = 8, height = 5, useDingbats=FALSE)
}