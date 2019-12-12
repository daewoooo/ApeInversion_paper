## BioNano calls for Chimp, Bonobo, Orangutan & Gorilla ##
##########################################################

##Load required libraries
suppressPackageStartupMessages( library(GenomicRanges) )

## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/BioNano_calls"

message("Reading in BioNano inversion calls ...")

## BioNano data for chimpanzee ##
#################################
BioNano.callsChimpanzee <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/chimp_twoEnzymes_sop_cluster_molecule_variant_summary.txt", skip = 1, header = TRUE, comment.char = '&', stringsAsFactors = FALSE)
## Keep SV types: inversion, duplication_inverted  & translocation_intrachr
BioNano.callsChimpanzee <- BioNano.callsChimpanzee[BioNano.callsChimpanzee$type %in% c('inversion', 'duplication_inverted', 'translocation_intrachr'),]
BioNano.callsChimpanzee.gr <- GRanges(seqnames=paste0('chr', BioNano.callsChimpanzee$chrom1), ranges=IRanges(start=BioNano.callsChimpanzee$position1, end=BioNano.callsChimpanzee$position2))
BioNano.callsChimpanzee.gr$gen <- BioNano.callsChimpanzee$type
BioNano.callsChimpanzee.gr$ID <- 'chimpanzee'
BioNano.callsChimpanzee.gr$study <- 'BioNano'

## BioNano data for bonobo ##
#################################
BioNano.callsBonobo <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/bonobo_twoEnzymes_sop_cluster_molecule_variant_summary.txt", skip = 1, header = TRUE, comment.char = '&', stringsAsFactors = FALSE)
## Keep SV types: inversion, duplication_inverted  & translocation_intrachr
BioNano.callsBonobo <- BioNano.callsBonobo[BioNano.callsBonobo$type %in% c('inversion', 'duplication_inverted', 'translocation_intrachr'),]
BioNano.callsBonobo.gr <- GRanges(seqnames=paste0('chr', BioNano.callsBonobo$chrom1), ranges=IRanges(start=BioNano.callsBonobo$position1, end=BioNano.callsBonobo$position2))
BioNano.callsBonobo.gr$gen <- BioNano.callsBonobo$type
BioNano.callsBonobo.gr$ID <- 'bonobo'
BioNano.callsBonobo.gr$study <- 'BioNano'

## BioNano data for gorilla ##
#################################
BioNano.callsGorilla <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/gorilla_dle1_sop_cluster_molecule_variant_summary.txt", skip = 1, header = TRUE, comment.char = '&', stringsAsFactors = FALSE)
## Keep SV types: inversion, duplication_inverted  & translocation_intrachr
BioNano.callsGorilla <- BioNano.callsGorilla[BioNano.callsGorilla$type %in% c('inversion', 'duplication_inverted', 'translocation_intrachr'),]
BioNano.callsGorilla.gr <- GRanges(seqnames=paste0('chr', BioNano.callsGorilla$chrom1), ranges=IRanges(start=BioNano.callsGorilla$position1, end=BioNano.callsGorilla$position2))
BioNano.callsGorilla.gr$gen <- BioNano.callsGorilla$type
BioNano.callsGorilla.gr$ID <- 'gorilla'
BioNano.callsGorilla.gr$study <- 'BioNano'

## BioNano data for orangutan ##
#################################
BioNano.callsOrangutan <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/orangutan_twoEnzymes_sop_cluster_molecule_variant_summary.txt", skip = 1, header = TRUE, comment.char = '&', stringsAsFactors = FALSE)
## Keep SV types: inversion, duplication_inverted  & translocation_intrachr
BioNano.callsOrangutan <- BioNano.callsOrangutan[BioNano.callsOrangutan$type %in% c('inversion', 'duplication_inverted', 'translocation_intrachr'),]
BioNano.callsOrangutan.gr <- GRanges(seqnames=paste0('chr', BioNano.callsOrangutan$chrom1), ranges=IRanges(start=BioNano.callsOrangutan$position1, end=BioNano.callsOrangutan$position2))
BioNano.callsOrangutan.gr$gen <- BioNano.callsOrangutan$type
BioNano.callsOrangutan.gr$ID <- 'orangutan'
BioNano.callsOrangutan.gr$study <- 'BioNano'

## Export all BioNano data into RData object
BioNano.invCalls <- suppressWarnings( c(BioNano.callsChimpanzee.gr, BioNano.callsBonobo.gr, BioNano.callsGorilla.gr, BioNano.callsOrangutan.gr) )
destination <- file.path(outputfolder, "BioNano.invCalls.RData")
save(BioNano.invCalls, file = destination)

message("DONE!!!")

# ####################################################################################################################
# ## Below is reading old BioNano data ## 
# ####################################################################################################################
# library(dplyr)
# 
# ## BioNano data for chimpanzee and orangutan (Zev_2018) ##
# ##########################################################
# BioNano.callsZev <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/aar6343_TableS12_3_BioNano.csv", header = TRUE, sep = ",", skip = 1)
# BioNano.callsZev.gr <- GRanges(seqnames=BioNano.callsZev$chr, ranges=IRanges(start=BioNano.callsZev$GRCh38.start, end=BioNano.callsZev$GRCh38.end))
# BioNano.callsZev.gr$gen <-BioNano.callsZev$Genotype
# BioNano.callsZev.gr$ID <- BioNano.callsZev$Species
# BioNano.callsZev.gr$study <- 'BioNano'
# BioNano.callsZev.gr$gen <- recode(BioNano.callsZev.gr$gen, heterozygous = "HET", homozygous = "HOM")
# 
# ## BioNano data for bonobo ##
# #############################
# BioNano.callsBonobo <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/BioNano_Bonobo/Bonobo_SVMerge_hg38_mergedSV_invOnly.txt", comment.char = '#', header = TRUE)
# ## Mask negative ranges and ranges of size 1
# mask <- which(BioNano.callsBonobo$RefEndPos < BioNano.callsBonobo$RefStartPos | BioNano.callsBonobo$RefEndPos == BioNano.callsBonobo$RefStartPos)
# BioNano.callsBonobo <- BioNano.callsBonobo[-mask,]
# ## Create GRanges object
# BioNano.callsBonobo.gr <- GRanges(seqnames=paste0('chr', BioNano.callsBonobo$RefcontigID1), ranges=IRanges(start=BioNano.callsBonobo$RefStartPos, end=BioNano.callsBonobo$RefEndPos))
# BioNano.callsBonobo.gr$gen <- BioNano.callsBonobo$Zygosity
# BioNano.callsBonobo.gr$gen <- recode(BioNano.callsBonobo.gr$gen, heterozygous = "HET", homozygous = "HOM")
# BioNano.callsBonobo.gr$ID <- 'bonobo'
# BioNano.callsBonobo.gr$study <- 'BioNano'
# 
# ## BioNano data for gorilla ##
# ##############################
# BioNano.callsGorilla <- read.table("/home/porubsky/WORK/Great_apes/BioNano_calls/BioNano_Gorilla_Kamilah/kamillah_dle1_sop_cluster_inv_minSize_1_rSize_50_percent_summary.txt", header=TRUE)
# ## Insert 'chr' if missing in chromosome name
# mask <- which(!grepl('chr', BioNano.callsGorilla$chrom))
# BioNano.callsGorilla$chrom <- sub(pattern='^', replacement='chr', BioNano.callsGorilla$chrom[mask])
# ## Set proper chromosome name for X and Y
# BioNano.callsGorilla$chrom[BioNano.callsGorilla$chrom == 'chr23'] <- 'chrX'
# BioNano.callsGorilla$chrom[BioNano.callsGorilla$chrom == 'chr24'] <- 'chrY'
# ## Create GRanges object
# BioNano.callsGorilla.gr <- GRanges(seqnames=BioNano.callsGorilla$chrom, ranges=IRanges(start=BioNano.callsGorilla$start, end=BioNano.callsGorilla$end))
# BioNano.callsGorilla.gr$gen <- ""
# BioNano.callsGorilla.gr$ID <- 'gorilla'
# BioNano.callsGorilla.gr$study <- 'BioNano'
# 
# ## Export all BioNano data into RData object
# BioNano.invCalls <- suppressWarnings( c(BioNano.callsZev.gr, BioNano.callsBonobo.gr, BioNano.callsGorilla.gr) )
# destination <- file.path(outputfolder, "BioNano.invCalls.RData")
# save(BioNano.invCalls, file = destination)
