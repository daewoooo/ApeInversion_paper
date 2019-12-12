## 10XGenomics calls for Chimp, Bonobo, Orangutan & Gorilla ##
##########################################################

##Load required libraries
library(dplyr)
library(gsubfn)

## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/TenX_calls/"

## 10X data for chimpanzee (Clint) ##
#####################################
tenX.callsChimp <- read.table("/home/porubsky/WORK/Great_apes/TenX_calls/INV_large_SVs_filtered.tab", stringsAsFactors = FALSE)
## Filter only inversions
tenX.callsChimp <- tenX.callsChimp[tenX.callsChimp$V5 == '<INV>',]
## Get inversion coordinates
end.coord <- strapplyc(tenX.callsChimp$V8, "^END=(\\d+)", simplify = T)
end.coord <- as.numeric(end.coord)
## Construct GRanges object
tenX.callsChimp.gr <- GRanges(seqnames=tenX.callsChimp$V1, ranges=IRanges(start=tenX.callsChimp$V2, end=end.coord))
tenX.callsChimp.gr$gen <- 'BioNano.callsZev$Genotype'
tenX.callsChimp.gr$ID <- 'chimpanzee'
tenX.callsChimp.gr$study <- 'tenX'

## Export all 10X calls into RData object
tenX.invCalls <- suppressWarnings( c(tenX.callsChimp.gr) )
destination <- file.path(outputfolder, "tenX.invCalls.RData")
save(tenX.invCalls, file = destination)

## 10X VALOR for chimpanzee (Clint) ##
######################################
valor.callsChimp <- read.table("/home/porubsky/WORK/Great_apes/TenX_calls/clint_valor_predicted_svs.bedpe", stringsAsFactors = FALSE)
valor.callsChimp.gr <- GRanges(seqnames=valor.callsChimp$V1, ranges=IRanges(start=valor.callsChimp$V2, end=valor.callsChimp$V6))
valor.callsChimp.gr <- sort(valor.callsChimp.gr) 
#ranges2UCSC(gr = valor.callsChimp.gr, outputDirectory = "/home/porubsky/WORK/Great_apes/TenX_calls/", index = "Valor_clint_calls", colorRGB = "178,102,255")  
