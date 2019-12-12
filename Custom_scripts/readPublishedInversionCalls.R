## This script reads in previously published inversion calls ##
###############################################################

## Load required libraries
suppressPackageStartupMessages( library(GenomicRanges) )

## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/Published_data"

message("Reading in published inversion calls ...")

## Feuk_2015 (available only for chimp) ##
##########################################
feuk.callsChimp <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Feuk_2005/panTroINVcalls_Feuk_TableS2.txt")
feuk.callsChimp.gr <- GRanges(seqnames=feuk.callsChimp$V1, ranges=IRanges(start=feuk.callsChimp$V2, end=feuk.callsChimp$V3))
feuk.callsChimp.gr$gen <- ""
feuk.callsChimp.gr$ID <- 'chimpanzee'
feuk.callsChimp.gr$study <- 'feuk_2005'

## Catacchio_2018 ##
####################
catacchio.calls <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Catacchio_2018/Supplemental_Tables_S3.csv", sep = ",", header=TRUE, skip = 1, stringsAsFactors=FALSE)
gen.position <- strsplit(catacchio.calls$Mapping..GRCh38.hg38., split="[\\:|\\-]+")
gen.position <- do.call(rbind, gen.position)
## Load Chimpanzee data
mask <- catacchio.calls$panTro5.calls != '' | grepl(catacchio.calls$Chimpanzee.Validation, pattern = 'inv')
gen.position.panTro <- gen.position[mask,]
gen.position.panTro[,1] <- gsub(gen.position.panTro[,1], pattern = 'Chr', replacement = 'chr')
catacchio.callsChimp.gr <- GRanges(seqnames=gen.position.panTro[,1], ranges=IRanges(start=as.numeric(gen.position.panTro[,2]), end=as.numeric(gen.position.panTro[,3])))
catacchio.callsChimp.gr$gen <- ""
catacchio.callsChimp.gr$ID <- 'chimpanzee'
catacchio.callsChimp.gr$study <- "Catacchio_2018"
## Load Gorilla data
mask <- catacchio.calls$gorGor4.calls != '' | grepl(catacchio.calls$Gorilla.Validation, pattern = 'inv')
gen.position.gorGor <- gen.position[mask,]
gen.position.gorGor[,1] <- gsub(gen.position.gorGor[,1], pattern = 'Chr', replacement = 'chr')
catacchio.callsGor.gr <- GRanges(seqnames=gen.position.gorGor[,1], ranges=IRanges(start=as.numeric(gen.position.gorGor[,2]), end=as.numeric(gen.position.gorGor[,3])))
catacchio.callsGor.gr$gen <- ""
catacchio.callsGor.gr$ID <- 'gorilla'
catacchio.callsGor.gr$study <- "Catacchio_2018"
## Load Orangutan data
mask <- catacchio.calls$ponAbe2.calls != '' | grepl(catacchio.calls$Orangutan.Validation, pattern = 'inv')
gen.position.ppy <- gen.position[mask,]
gen.position.ppy[,1] <- gsub(gen.position.ppy[,1], pattern = 'Chr', replacement = 'chr')
catacchio.callsOrang.gr <- GRanges(seqnames=gen.position.ppy[,1], ranges=IRanges(start=as.numeric(gen.position.ppy[,2]), end=as.numeric(gen.position.ppy[,3])))
catacchio.callsOrang.gr$gen <- ""
catacchio.callsOrang.gr$ID <- 'orangutan'
catacchio.callsOrang.gr$study <- "Catacchio_2018"

## Zev_2018 (large validated calls) ##
######################################
zev.callsValid <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Zev_2018/aar6343_TableS12_1.csv", sep=",", header=TRUE, skip=1, stringsAsFactors=FALSE)
zev.callsValid.gr <- GRanges(seqnames=paste0('chr', zev.callsValid$GRCH38_chr), ranges=IRanges(start=zev.callsValid$start, end=zev.callsValid$end), ID=zev.callsValid$species)
## Load Chimpanzee data
zev.callsChimp.gr <- zev.callsValid.gr[zev.callsValid.gr$ID == 'chimpanzee']
zev.callsChimp.gr <- zev.callsChimp.gr[,0]
zev.callsChimp.gr$gen <- ""
zev.callsChimp.gr$ID <- "chimpanzee"
zev.callsChimp.gr$study <- "Zev_2018"
## Load Orangutan data
zev.callsOrang.gr <- zev.callsValid.gr[zev.callsValid.gr$ID == 'orangutan']
zev.callsOrang.gr <- zev.callsOrang.gr[,0]
zev.callsOrang.gr$gen <- ""
zev.callsOrang.gr$ID <- "orangutan"
zev.callsOrang.gr$study <- "Zev_2018"

## Prakash_1982 (large pericentromeric inv) ##
##############################################
prakash.calls <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Prakash_1982/pericenINV_coord_shwethaM_prakashStudy.csv", sep = ",")
prakash.calls.gr <- GRanges(seqnames=prakash.calls$V1, ranges=IRanges(start=prakash.calls$V2, end=prakash.calls$V3))
prakash.calls.gr$gen <- ""
prakash.calls.gr$ID <- "PericenINV"
prakash.calls.gr$study <- "Prakash_1982"

## Szamalek_2018 ##
###################
szamalek.calls <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Szamalek_2006/szamalek_2006_pericenINV.bed", header = TRUE)
szamalek.calls.gr <- GRanges(seqnames=szamalek.calls$chr.1, ranges=IRanges(start=szamalek.calls$start.1, end=szamalek.calls$end.2))
szamalek.calls.gr$gen <- ""
szamalek.calls.gr$ID <- "PericenINV"
szamalek.calls.gr$study <- "Szamalek_2006"


## Export all published data into RData object
published.invCalls <- suppressWarnings( c(feuk.callsChimp.gr, 
                                          catacchio.callsChimp.gr,
                                          catacchio.callsGor.gr,
                                          catacchio.callsOrang.gr, 
                                          zev.callsChimp.gr, 
                                          zev.callsOrang.gr, 
                                          prakash.calls.gr, 
                                          szamalek.calls.gr) )
destination <- file.path(outputfolder, "published.invCalls.RData")
save(published.invCalls, file = destination)

message("DONE!!!")
