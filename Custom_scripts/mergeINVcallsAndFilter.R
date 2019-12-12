## Load required libraries ##
#############################
suppressPackageStartupMessages( library(primatR) )

outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls"

## segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load Strand-seq calls ##
###########################
message("Preparing final inversion callset ...")
all.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/allGreatApesFinalINVcalls.RData"))
validUncertain.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/validUncertain_calls.RData"))

## Filter Strand-seq inversion calls ##
#######################################
all.gr <- all.gr[grep("segDup|CEN", all.gr$SVclass, ignore.case = TRUE, invert = TRUE)] #remove all call other then inversions
all.gr <- all.gr[grep("\\?", all.gr$SVclass, invert = TRUE)] #remove all uncertain calls

## Merge confident and valid uncertain callsets ##
##################################################
all.inversion.calls <- sort(c(all.gr, validUncertain.gr))

## Get overlaps with Human segmental duplication track ##
#########################################################
all.inversion.calls <- getSegDupOverlaps(query.gr = all.inversion.calls, subject.gr = seg.dup.gr)

## Filter data ##
#################
all.inversion.calls.filt <- all.inversion.calls
## Keep inversion calls (INV or invDup) that have at last 1kb of unique sequence
all.inversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$TotalUniqueBases >= 1000]
## Keep only inverted duplications 10kb longer
all.inversion.calls.filt <- all.inversion.calls.filt[!(all.inversion.calls.filt$SVclass == 'invDup' & width(all.inversion.calls.filt) < 10000)]

## Export filtered INV calls ##
###############################
destination <- file.path(outputfolder, 'all.inversion.calls.filt.RData')
save(all.inversion.calls.filt, file = destination)

## Export UCSC bed formated files ##
####################################
index <- "Final_INVcalls_filtered"
inv.calls.ucsc <- as.data.frame(all.inversion.calls.filt)
inv.calls.ucsc$ucsc.id <- paste0(inv.calls.ucsc$SVclass, "_", inv.calls.ucsc$ID)
inv.calls.ucsc$score <- 0
inv.calls.ucsc$gen <- dplyr::recode(inv.calls.ucsc$gen, 'HET' = '+', 'HOM' = '-')
## Select columns to export
inv.calls.ucsc <- inv.calls.ucsc[,c('seqnames', 'start', 'end', 'ucsc.id', 'score', 'gen')]
## Write to a file
ucsc.header <- paste0('track name=', index,' ROIs description=', index, '_Bed_of_inverted_regions visibility=dense colorByStrand=\"28,144,153 221,28,119\"')
destination <- file.path(outputfolder, 'all.inversion.calls.filt.ucsc.bed')
write.table(ucsc.header, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
write.table(inv.calls.ucsc, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')

message("DONE!!!")