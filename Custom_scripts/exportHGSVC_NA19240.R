## Compile human specific inversion callset [NA19240] ##
########################################################

library(primatR)

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls"

message("Loading INV calls ...")
exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Composite_files/NA19240_CompositeFile_breakpoints/syncReads_Multireads_NA19240_breakPoints_checked.bed", 
               outputfolder = outputfolder,
               index = "NA19240"
)

message("Exporting INV calls ...")
stranS.callsHuman <- read.table(file.path(outputfolder, "NA19240_INVcalls.bed"), skip = 1)
stranS.callsHuman.gr <- GRanges(seqnames=stranS.callsHuman$V1, ranges=IRanges(start=stranS.callsHuman$V2, end=stranS.callsHuman$V3))
stranS.callsHuman.gr$gen <- stranS.callsHuman$V6
stranS.callsHuman.gr$gen <- dplyr::recode(stranS.callsHuman.gr$gen, '+' = 'HET', '-' = 'HOM')
stranS.callsHuman.gr$ID <- 'NA19240'
stranS.callsHuman.gr$SVclass <- stranS.callsHuman$V4
## Filter NA19240 inversion list
stranS.callsHuman.gr <- stranS.callsHuman.gr[grep("segDup|CEN", stranS.callsHuman.gr$SVclass, ignore.case = TRUE, invert = TRUE)] #remove all call other then inversions
stranS.callsHuman.gr <- stranS.callsHuman.gr[grep("\\?", stranS.callsHuman.gr$SVclass, invert = TRUE)] #remove all uncertain calls
## Calculate overlap with segmental duplications
stranS.callsHuman.gr <- getSegDupOverlaps(query.gr = stranS.callsHuman.gr, subject.gr = seg.dup.gr)
## Keep inversion calls (INV or invDup) that have at last 1kb of unique sequence
stranS.callsHuman.gr <- stranS.callsHuman.gr[stranS.callsHuman.gr$TotalUniqueBases >= 1000]
## Keep only inverted duplications 10kb longer
stranS.callsHuman.gr <- stranS.callsHuman.gr[!(stranS.callsHuman.gr$SVclass == 'invDup' & width(stranS.callsHuman.gr) < 10000)]
## Mark misorients in the final human inversion callset
misorients.annot <- "/home/porubsky/WORK/Great_apes/Ashley_inversions/Misorients/"
HGSVC.misorients <- read.table(file.path(misorients.annot, "misAssem2remove.NHP.csv"), sep=',', header = TRUE, stringsAsFactors = FALSE)
HGSVC.misorients.gr <- GRanges(seqnames = HGSVC.misorients$seqnames, ranges=IRanges(start = HGSVC.misorients$start, end = HGSVC.misorients$end), toRemove = HGSVC.misorients$toRemove)
HGSVC.misorients.gr <- HGSVC.misorients.gr[HGSVC.misorients.gr$toRemove == TRUE]
hits <- suppressWarnings( findOverlaps(stranS.callsHuman.gr, HGSVC.misorients.gr) )
stranS.callsHuman.gr$misorient <- FALSE
stranS.callsHuman.gr$misorient[queryHits(hits)] <- TRUE
## Export human inversion callset
save(stranS.callsHuman.gr, file = file.path(outputfolder, "human.NA19240.filt.RData"))

#=====================================================
# Old analysis
#=====================================================
# ## Load composite files
# chimp.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
# chimp.data$ID <- 'chimpanzee'
# bonobo.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
# bonobo.data$ID <- 'bonobo'
# gorilla.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
# gorilla.data$ID <- 'gorilla'
# orangutan.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
# orangutan.data$ID <- 'orangutan'
# na19240.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncFrags_NA19240.RData"))
# na19240.data$ID <- 'NA19240'
# 
# ## Load HGSCV data [NA19240]
# hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
# #hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions_manualCheck.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
# hgsvc.simple.inv.NA19240 <- hgsvc.simple.inv[hgsvc.simple.inv$NA19240 == 'HOM' | hgsvc.simple.inv$NA19240 == 'HET',]
# ## Filter out manually excluded regions
# #hgsvc.simple.inv.NA19240 <- hgsvc.simple.inv.NA19240[hgsvc.simple.inv.NA19240$FILTER != TRUE,]
# ## Keep only inversion supported by StrandS
# hgsvc.simple.inv.NA19240 <- hgsvc.simple.inv.NA19240[grepl(hgsvc.simple.inv.NA19240$platforms, pattern = 'Ss'),]
# #hgsvc.simple.inv.NA19240 <- hgsvc.simple.inv.NA19240[grepl(hgsvc.simple.inv.NA19240$platforms, pattern = 'Ss') | grepl(hgsvc.simple.inv.NA19240$support, pattern = 'Ss'),]
# ## Construct GRanges object
# hgsvc.simple.inv.NA19240.gr <- GRanges(seqnames=hgsvc.simple.inv.NA19240$seqnames, 
#                                        ranges=IRanges(start = hgsvc.simple.inv.NA19240$innerBP_start, end = hgsvc.simple.inv.NA19240$innerBP_end),
#                                        gen=hgsvc.simple.inv.NA19240$NA19240,
#                                        ID='NA19240',
#                                        SVclass='INV')
# ## Get overlaps with Human segmental duplication track
# hgsvc.simple.inv.NA19240.gr <- getSegDupOverlaps(query.gr = hgsvc.simple.inv.NA19240.gr, subject.gr = seg.dup.gr)
# ## Keep inversion calls that have at last 1kb of unique sequence and do not overlap with SD more than 90%
# hgsvc.simple.inv.NA19240.gr <- hgsvc.simple.inv.NA19240.gr[hgsvc.simple.inv.NA19240.gr$TotalUniqueBases >= 1000 & hgsvc.simple.inv.NA19240.gr$SDTrackPercOverlap <= 90]
# 
# ## Genotype regions from NA19240 ##
# ###################################
# ## Set ranges to be re-genotyped
# regions.to.genotype <- hgsvc.simple.inv.NA19240.gr[,0]
# ## Genotype chimpanzee data
# genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
# ## Genotype bonobo data
# genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
# ## Genotype gorilla data
# genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
# ## Genotype orangutan data
# genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')
# 
# ## Extract genotypes for all apes
# genoT.df <- mcols(genotyped.regions)[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]
# ## Filter out regions that couldn't not been genotyped or are all in reference orientation
# mask <- apply(genoT.df, 1, function(x) all(x == 'lowReads') | all(x == 'REF'))
# ## Filter NA19240 callset
# genotyped.regions <- genotyped.regions[!mask]
# hgsvc.simple.inv.NA19240.gr <- hgsvc.simple.inv.NA19240.gr[!mask]
# 
# ## Exclude SD region from called regions [OPTIONAL]
# #hgsvc.simple.inv.NA19240.gr <- removeSDflanks(gr = hgsvc.simple.inv.NA19240.gr, sd.track = reduce(seg.dup.gr))
# 
# ## Load complete dataset including putative human specific inversions
# all.inversion.calls.filt.noHS <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.noHS.RData"))
# 
# ## Find missed calls from NA19240 that could have been genotyped in at least one ape
# missed.NA19240 <- subsetByOverlaps(hgsvc.simple.inv.NA19240.gr, all.inversion.calls.filt.noHS, invert = TRUE)
# 
# ## Add human specific inversion not present in the HGSVC callset
# ## Load human specific inversions
# human.specific.regions.trueSet <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
# human.specific.regions.trueSet.toadd <- human.specific.regions.trueSet[,0]
# human.specific.regions.trueSet.toadd$gen <- 'HOM'
# human.specific.regions.trueSet.toadd$ID <- 'NA19240'
# human.specific.regions.trueSet.toadd$SVclass <- 'INV'
# human.specific.regions.trueSet.toadd <- getSegDupOverlaps(query.gr = human.specific.regions.trueSet.toadd, subject.gr = seg.dup.gr)
# hgsvc.simple.inv.NA19240.gr <- c(human.specific.regions.trueSet.toadd, hgsvc.simple.inv.NA19240.gr) #add unique human specific inversion to NA19240 callset
# seqlevels(hgsvc.simple.inv.NA19240.gr) <- paste0('chr', c(1:22, 'X'))
# hgsvc.simple.inv.NA19240.gr <- sort(hgsvc.simple.inv.NA19240.gr)
# 
# ## Export human inversion callset
# save(hgsvc.simple.inv.NA19240.gr, file = "/home/porubsky/WORK/Great_apes/Final_INV_calls/human.NA19240.plusHS.RData")
# 
