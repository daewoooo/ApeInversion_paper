## Export lineage specific inversions for great apes ##
#######################################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Params
thresh <- 50 #required reciprocal overlap (%)
## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/MasterTables"
if (!dir.exists(outputfolder)) {
  dir.create(outputfolder)
}

## Load Strand-seq calls after filtering and merging with validated uncertain calls ##
######################################################################################
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
# retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
# retain only inverted duplications
all.invertedDups.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Get comparison all simple inversions against each other
all.SimpleInversion.calls.filt.comp <- getDisjointOverlapsWeighted(gr = all.SimpleInversion.calls.filt, percTh = 50)
inversion.calls.comp.perRegion <- split(all.SimpleInversion.calls.filt.comp, all.SimpleInversion.calls.filt.comp$sub.group)
## Select inversions that do not overlap with any other inversion
inversion.calls.unique <- inversion.calls.comp.perRegion[lengths(inversion.calls.comp.perRegion) == 1]
inversion.calls.unique <- unlist(inversion.calls.unique, use.names = FALSE)

## Plot lineage specific inversions
seqlevels(inversion.calls.unique) <- gtools::mixedsort(seqlevels(inversion.calls.unique))
genomewideRangesIdeo(inversion.calls.unique, userTrack = seg.dup.gr, colors = c('#3182bd','#31a354','#8856a7','#e6550d'), bsgenome = BSgenome.Hsapiens.UCSC.hg38)
plt1 <- plotColumnCounts(inversion.calls.unique, colName = 'gen', facetID = 'ID', colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
plt2 <- rangesSizeDistribution(inversion.calls.unique, plotUniqueBases=FALSE, colors = c('#3182bd','#31a354','#8856a7','#e6550d'))
plot_grid(plt2, plt1, nrow = 1, rel_widths = c(4,1))
plt1 <- plotColumnCountsPerChr(inversion.calls.unique, colName = 'gen', facetID = 'ID', normChrSize = TRUE)
plt2 <- basesPerGenotypePerChr(inversion.calls.unique, normChrSize = TRUE)
plot_grid(plt1, plt2, ncol = 1, align = 'v')

## Chimpanzee ##
################
stranS.callsChimp.gr <- all.inversion.calls.filt.grl[['chimpanzee']]
other.apes.grl <- all.inversion.calls.filt.grl[names(all.inversion.calls.filt.grl) != 'chimpanzee']

message("Comparing chimpanzee data ...")
query <- stranS.callsChimp.gr
for (i in seq_along(other.apes.grl)) {
  ID <- names(other.apes.grl[i])
  message("    Processing ", ID, " data ...")
  subject <- other.apes.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = ID)
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
