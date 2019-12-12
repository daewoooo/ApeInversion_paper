## Export manually processed inversion calls ##
###############################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )

outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls"

message("Loading INV calls ...")
exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_chimp_breakPoints_checked.bed", 
               outputfolder = outputfolder,
               index = "chimpanzee"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_bonobo_breakPoints_checked.bed", 
               outputfolder = outputfolder,
               index = "bonobo"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_gorilla_breakPoints_checked.bed", 
               outputfolder = outputfolder,
               index = "gorilla"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_orangutan_breakPoints_checked.bed", 
               outputfolder = outputfolder,
               index = "orangutan"
)

## Compile all inversions into a single RData file
message("Exporting INV calls ...")
stranS.callsChimp <- read.table(file.path(outputfolder, "chimpanzee_INVcalls.bed"), skip = 1)
stranS.callsChimp.gr <- GRanges(seqnames=stranS.callsChimp$V1, ranges=IRanges(start=stranS.callsChimp$V2, end=stranS.callsChimp$V3))
stranS.callsChimp.gr$gen <- stranS.callsChimp$V6
stranS.callsChimp.gr$gen <- dplyr::recode(stranS.callsChimp.gr$gen, '+' = 'HET', '-' = 'HOM')
stranS.callsChimp.gr$ID <- 'chimpanzee'
stranS.callsChimp.gr$SVclass <- stranS.callsChimp$V4

stranS.callsBonobo <- read.table(file.path(outputfolder, "bonobo_INVcalls.bed"), skip = 1)
stranS.callsBonobo.gr <- GRanges(seqnames=stranS.callsBonobo$V1, ranges=IRanges(start=stranS.callsBonobo$V2, end=stranS.callsBonobo$V3))
stranS.callsBonobo.gr$gen <- stranS.callsBonobo$V6
stranS.callsBonobo.gr$gen <- dplyr::recode(stranS.callsBonobo.gr$gen, '+' = 'HET', '-' = 'HOM')
stranS.callsBonobo.gr$ID <- 'bonobo'
stranS.callsBonobo.gr$SVclass <- stranS.callsBonobo$V4

stranS.callsOrangutan <- read.table(file.path(outputfolder, "orangutan_INVcalls.bed"), skip = 1)
stranS.callsOrangutan.gr <- GRanges(seqnames=stranS.callsOrangutan$V1, ranges=IRanges(start=stranS.callsOrangutan$V2, end=stranS.callsOrangutan$V3))
stranS.callsOrangutan.gr$gen <- stranS.callsOrangutan$V6
stranS.callsOrangutan.gr$gen <- dplyr::recode(stranS.callsOrangutan.gr$gen, '+' = 'HET', '-' = 'HOM')
stranS.callsOrangutan.gr$ID <- 'orangutan'
stranS.callsOrangutan.gr$SVclass <- stranS.callsOrangutan$V4

stranS.callsGorilla <- read.table(file.path(outputfolder, "gorilla_INVcalls.bed"), skip = 1)
stranS.callsGorilla.gr <- GRanges(seqnames=stranS.callsGorilla$V1, ranges=IRanges(start=stranS.callsGorilla$V2, end=stranS.callsGorilla$V3))
stranS.callsGorilla.gr$gen <- stranS.callsGorilla$V6
stranS.callsGorilla.gr$gen <- dplyr::recode(stranS.callsGorilla.gr$gen, '+' = 'HET', '-' = 'HOM')
stranS.callsGorilla.gr$ID <- 'gorilla'
stranS.callsGorilla.gr$SVclass <- stranS.callsGorilla$V4

all.inversion.calls <- c(stranS.callsChimp.gr, stranS.callsBonobo.gr, stranS.callsOrangutan.gr, stranS.callsGorilla.gr)
destination <- file.path(outputfolder, 'allGreatApesFinalINVcalls.RData')
save(all.inversion.calls, file = destination)

message("DONE!!!")