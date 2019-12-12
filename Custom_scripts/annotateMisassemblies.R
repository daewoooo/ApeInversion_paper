## Remove calls overlapping with predicted reference misassemblies (misorients) ##
##################################################################################
## Load required libraries
suppressPackageStartupMessages( library(primatR) )

message("Annotating predicted human misassemblies ...")

## Load complete dataset including putative human specific inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']

## Load annotation of human misorients
misorients.annot <- "/home/porubsky/WORK/Great_apes/Ashley_inversions/Misorients/"
#GRashley.misorients <- read.table(file.path(misorients.annot, "Supplemental_Tables_S3_potential_misorients.csv"), sep = ',', header = TRUE, stringsAsFactors = FALSE)
#GRashley.misorients.gr <- GRanges(seqnames = GRashley.misorients$Chr, ranges=IRanges(start = GRashley.misorients$Start, end = GRashley.misorients$End))
HGSVC.misorients <- read.table(file.path(misorients.annot, "Misorients_HOMinAll_Feb202017_HGSVC.txt"), header = TRUE, stringsAsFactors = FALSE)
HGSVC.misorients.gr <- GRanges(seqnames = HGSVC.misorients$chr, ranges=IRanges(start = HGSVC.misorients$start, end = HGSVC.misorients$end))
HGSVC.valid.misorients <- read.table(file.path(misorients.annot, "suppTable11_validatedMisassemblies.txt"), header = TRUE, stringsAsFactors = FALSE)
HGSVC.valid.misorients.gr <- GRanges(seqnames = HGSVC.valid.misorients$chr, ranges=IRanges(start = HGSVC.valid.misorients$start, end = HGSVC.valid.misorients$end))
## Take HGSVC validated misorients not present in comprehensive HGSVC list of misorients
HGSVC.valid.misorients.gr <- subsetByOverlaps(HGSVC.valid.misorients.gr, HGSVC.misorients.gr, invert = TRUE)
misorients.gr <- suppressWarnings( c(HGSVC.misorients.gr, HGSVC.valid.misorients.gr) )
misorients.gr <- reduce(misorients.gr)

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

## Genotype likely misassembled regions in NHP ##
#################################################
## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Find ranges having 50% reciprocal overlap with predicted misorients
misorients.gr.inNHP <- getReciprocalOverlaps(query = misorients.gr, subject = all.SimpleInversion.calls.filt, report = 'query', thresh = 50)
misorients.gr.inNHP <- misorients.gr.inNHP[misorients.gr.inNHP$perc.overlap >= 50]

## Set of ranges to be re-genotyped
regions.to.genotype <- misorients.gr.inNHP[,0]
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

## Remove regions being HOM in all NHP ##
#########################################
col.idx <- grep(names(mcols(genotyped.regions)), pattern = 'genoT')
misAssem2remove <- as.data.frame(genotyped.regions[,col.idx])
write.table(misAssem2remove, file = file.path(misorients.annot, 'misAssem2remove.NHP.csv'), quote = FALSE, sep = ',', row.names = FALSE)
