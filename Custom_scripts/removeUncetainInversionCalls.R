## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Load Strand-seq calls ##
###########################
all.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/allGreatApesFinalINVcalls.RData"))
## Get uncertain calls only
uncertainCalls.gr <- all.gr[grep("\\?", all.gr$SVclass)]

## Load calls from orthogonal technologies ##
#############################################
pbsv.calls <- get(load("/home/porubsky/WORK/Great_apes/PacBio_calls/PBSV_calls_greatApes.RData"))
delly.calls <- get(load("/home/porubsky/WORK/Great_apes/Delly_calls/DELLY_calls_greatApes.RData"))
sniffles.calls <- get(load("/home/porubsky/WORK/Great_apes/PacBio_calls/SNIFFLES_calls_greatApes.RData"))
orthogonal.calls <- c(pbsv.calls, delly.calls, sniffles.calls)

## Split calls by individual
uncertainCalls.grl <- split(uncertainCalls.gr, uncertainCalls.gr$ID)
orthogonalCalls.grl <- split(orthogonal.calls, orthogonal.calls$ID)

## Validate uncertain calls ##
##############################
## Get reciprocal overlap
message("Validating uncertain inversion calls ...")
thresh <- 50
chimpanzee.uncertain <- getReciprocalOverlaps(query = uncertainCalls.grl[[2]], subject = orthogonalCalls.grl[[2]], thresh = thresh, report = 'query', index = 'uncertain.calls')
bonobo.uncertain <- getReciprocalOverlaps(query = uncertainCalls.grl[[1]], subject = orthogonalCalls.grl[[1]], thresh = thresh, report = 'query', index = 'uncertain.calls')
gorilla.uncertain <- getReciprocalOverlaps(query = uncertainCalls.grl[[3]], subject = orthogonalCalls.grl[[3]], thresh = thresh, report = 'query', index = 'uncertain.calls')
orangutan.uncertain <- getReciprocalOverlaps(query = uncertainCalls.grl[[4]], subject = orthogonalCalls.grl[[4]], thresh = thresh, report = 'query', index = 'uncertain.calls')

## Only inversion with 50% reciprocal overlap are considered valid
chimpanzee.uncertain.valid <- chimpanzee.uncertain[chimpanzee.uncertain$perc.overlap_uncertain.calls >= 50,]
bonobo.uncertain.valid <- bonobo.uncertain[bonobo.uncertain$perc.overlap_uncertain.calls >= 50,]
gorilla.uncertain.valid <- gorilla.uncertain[gorilla.uncertain$perc.overlap_uncertain.calls >= 50,]
orangutan.uncertain.valid <- orangutan.uncertain[orangutan.uncertain$perc.overlap_uncertain.calls >= 50,]

## Export validated uncertain inversions
outputfolder <- "/home/porubsky/WORK/Great_apes/Final_INV_calls/"
valid.uncertainCalls.gr <- c(chimpanzee.uncertain.valid, bonobo.uncertain.valid, gorilla.uncertain.valid, orangutan.uncertain.valid)
destination <- file.path(outputfolder, "validUncertain_calls_reciprocalOverlap.RData")
save(valid.uncertainCalls.gr, file = destination)

valid.uncertainCalls.gr$SVclass <- gsub(valid.uncertainCalls.gr$SVclass, pattern = "\\?", replacement = "")
valid.uncertainCalls.gr <- valid.uncertainCalls.gr[,1:3]
destination <- file.path(outputfolder, "validUncertain_calls.RData")
save(valid.uncertainCalls.gr , file = destination)

## Get proportion of validated uncertain calls
perc.valid <- round( (length(valid.uncertainCalls.gr)  / length(uncertainCalls.gr)) * 100, digits = 2 )
message("Validated ", perc.valid, "% of all uncertain calls")
message("DONE!!!")
