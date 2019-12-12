## Load calls from orthogonal technologies ##
#############################################

## Load required libraries
suppressPackageStartupMessages( library(gsubfn) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(GenomicRanges) )

message("Reading in PB&IL inversion calls ...")

## Load PBSV calls
#chimapnzee
pbsv.calls.chimpanzee <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/chimpanzee/chimpanzee_pbsv.vcf", stringsAsFactors = FALSE)
end.pos <- as.numeric( strapplyc(as.character(pbsv.calls.chimpanzee$V8), "END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(pbsv.calls.chimpanzee$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
pbsv.calls.chimpanzee.gr <- GRanges(seqnames=pbsv.calls.chimpanzee$V1, ranges=IRanges(start=pbsv.calls.chimpanzee$V2, end=end.pos), gen=gen, ID='chimpanzee', study='PBSV')

#bonobo
pbsv.calls.bonobo <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/bonobo/bonobo_pbsv.vcf", stringsAsFactors = FALSE)
end.pos <- as.numeric( strapplyc(as.character(pbsv.calls.bonobo$V8), "END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(pbsv.calls.bonobo$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
pbsv.calls.bonobo.gr <- GRanges(seqnames=pbsv.calls.bonobo$V1, ranges=IRanges(start=pbsv.calls.bonobo$V2, end=end.pos), gen=gen, ID='bonobo', study='PBSV')

#gorilla
pbsv.calls.gorilla <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/gorilla/gorilla_pbsv.vcf", stringsAsFactors = FALSE)
end.pos <- as.numeric( strapplyc(as.character(pbsv.calls.gorilla$V8), "END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(pbsv.calls.gorilla$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
pbsv.calls.gorilla.gr <- GRanges(seqnames=pbsv.calls.gorilla$V1, ranges=IRanges(start=pbsv.calls.gorilla$V2, end=end.pos), gen=gen, ID='gorilla', study='PBSV')

#orangutan
pbsv.calls.orangutan <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/orangutan/orangutan_pbsv.vcf", stringsAsFactors = FALSE)
end.pos <- as.numeric( strapplyc(as.character(pbsv.calls.orangutan$V8), "END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(pbsv.calls.orangutan$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
pbsv.calls.orangutan.gr <- GRanges(seqnames=pbsv.calls.orangutan$V1, ranges=IRanges(start=pbsv.calls.orangutan$V2, end=end.pos), gen=gen, ID='orangutan', study='PBSV')

#Export all PBSV calls
outputfolder <- "/home/porubsky/WORK/Great_apes/PacBio_calls"
all.pbsv.calls <- c(pbsv.calls.chimpanzee.gr, pbsv.calls.bonobo.gr, pbsv.calls.gorilla.gr, pbsv.calls.orangutan.gr)
destination <- file.path(outputfolder, 'PBSV_calls_greatApes.RData')
save(all.pbsv.calls, file = destination)


## Load SNIFFLES calls
#chimapnzee
sniffles.calls.chimpanzee <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/chimpanzee/chimpanzee_sniffles.vcf", stringsAsFactors = FALSE)
sniffles.calls.chimpanzee <- filter(sniffles.calls.chimpanzee, V5 == '<INV>', V7 == 'PASS')
precision <- strapplyc(as.character(sniffles.calls.chimpanzee$V8), "^(\\w+);", simplify = T)
sniffles.calls.chimpanzee <- sniffles.calls.chimpanzee[precision == 'PRECISE',]
end.pos <- as.numeric( strapplyc(as.character(sniffles.calls.chimpanzee$V8), "END=(\\d+)", simplify = T) )
sniffles.calls.chimpanzee.gr <- GRanges(seqnames=sniffles.calls.chimpanzee$V1, ranges=IRanges(start=sniffles.calls.chimpanzee$V2, end=end.pos))
sniffles.calls.chimpanzee.gr$gen <- ''
sniffles.calls.chimpanzee.gr$ID <- 'chimpanzee'
sniffles.calls.chimpanzee.gr$study <- 'SNIFFLES'

#bonobo
sniffles.calls.bonobo <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/bonobo/bonobo_sniffles.vcf", stringsAsFactors = FALSE)
sniffles.calls.bonobo <- filter(sniffles.calls.bonobo, V5 == '<INV>', V7 == 'PASS')
precision <- strapplyc(as.character(sniffles.calls.bonobo$V8), "^(\\w+);", simplify = T)
sniffles.calls.bonobo <- sniffles.calls.bonobo[precision == 'PRECISE',]
end.pos <- as.numeric( strapplyc(as.character(sniffles.calls.bonobo$V8), "END=(\\d+)", simplify = T) )
sniffles.calls.bonobo.gr <- GRanges(seqnames=sniffles.calls.bonobo$V1, ranges=IRanges(start=sniffles.calls.bonobo$V2, end=end.pos))
sniffles.calls.bonobo.gr$gen <- ''
sniffles.calls.bonobo.gr$ID <- 'bonobo'
sniffles.calls.bonobo.gr$study <- 'SNIFFLES'

#gorilla
sniffles.calls.gorilla <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/gorilla/gorilla_sniffles.vcf", stringsAsFactors = FALSE)
sniffles.calls.gorilla <- filter(sniffles.calls.gorilla, V5 == '<INV>', V7 == 'PASS')
precision <- strapplyc(as.character(sniffles.calls.gorilla$V8), "^(\\w+);", simplify = T)
sniffles.calls.gorilla <- sniffles.calls.gorilla[precision == 'PRECISE',]
end.pos <- as.numeric( strapplyc(as.character(sniffles.calls.gorilla$V8), "END=(\\d+)", simplify = T) )
sniffles.calls.gorilla.gr <- GRanges(seqnames=sniffles.calls.gorilla$V1, ranges=IRanges(start=sniffles.calls.gorilla$V2, end=end.pos))
sniffles.calls.gorilla.gr$gen <- ''
sniffles.calls.gorilla.gr$ID <- 'gorilla'
sniffles.calls.gorilla.gr$study <- 'SNIFFLES'

#orangutan
sniffles.calls.orangutan <- read.table("/home/porubsky/WORK/Great_apes/PacBio_calls/orangutan/orangutan_sniffles.vcf", stringsAsFactors = FALSE)
sniffles.calls.orangutan <- filter(sniffles.calls.orangutan, V5 == '<INV>', V7 == 'PASS')
precision <- strapplyc(as.character(sniffles.calls.orangutan$V8), "^(\\w+);", simplify = T)
sniffles.calls.orangutan <- sniffles.calls.orangutan[precision == 'PRECISE',]
end.pos <- as.numeric( strapplyc(as.character(sniffles.calls.orangutan$V8), "END=(\\d+)", simplify = T) )
sniffles.calls.orangutan.gr <- GRanges(seqnames=sniffles.calls.orangutan$V1, ranges=IRanges(start=sniffles.calls.orangutan$V2, end=end.pos))
sniffles.calls.orangutan.gr$gen <- ''
sniffles.calls.orangutan.gr$ID <- 'orangutan'
sniffles.calls.orangutan.gr$study <- 'SNIFFLES'

#Export all SNIFFLES calls
outputfolder <- "/home/porubsky/WORK/Great_apes/PacBio_calls/"
all.sniffles.calls <- c(sniffles.calls.chimpanzee.gr, sniffles.calls.bonobo.gr, sniffles.calls.gorilla.gr, sniffles.calls.orangutan.gr)
all.sniffles.calls <- keepStandardChromosomes(all.sniffles.calls, pruning.mode = 'coarse')
destination <- file.path(outputfolder, 'SNIFFLES_calls_greatApes.RData')
save(all.sniffles.calls, file = destination)


## Load Delly calls
#chimpanzee
delly.calls.chimpanzee <- read.table("/home/porubsky/WORK/Great_apes/Delly_calls/chimpanzee_clint_delly.vcf")
delly.calls.chimpanzee <- delly.calls.chimpanzee[delly.calls.chimpanzee$V7 == 'PASS',]
end.pos <- as.numeric( strapplyc(as.character(delly.calls.chimpanzee$V8), ";END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(delly.calls.chimpanzee$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
delly.calls.chimpanzee.gr <- GRanges(seqnames=delly.calls.chimpanzee$V1, ranges=IRanges(start=delly.calls.chimpanzee$V2, end=end.pos), gen=gen, ID='chimpanzee', study='DELLY')
delly.calls.chimpanzee.gr <- keepSeqlevels(delly.calls.chimpanzee.gr, paste0('chr', c(1:22, 'X')), pruning.mode = 'coarse')

#bonobo
delly.calls.bonobo <- read.table("/home/porubsky/WORK/Great_apes/Delly_calls/bonobo_delly.vcf")
delly.calls.bonobo <- delly.calls.bonobo[delly.calls.bonobo$V7 == 'PASS',]
end.pos <- as.numeric( strapplyc(as.character(delly.calls.bonobo$V8), ";END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(delly.calls.bonobo$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
delly.calls.bonobo.gr <- GRanges(seqnames=delly.calls.bonobo$V1, ranges=IRanges(start=delly.calls.bonobo$V2, end=end.pos), gen=gen, ID='bonobo', study='DELLY')
delly.calls.bonobo.gr <- keepSeqlevels(delly.calls.bonobo.gr, paste0('chr', c(1:22, 'X')), pruning.mode = 'coarse')

#gorilla
delly.calls.gorilla <- read.table("/home/porubsky/WORK/Great_apes/Delly_calls/gorilla_kamilah_delly.vcf")
delly.calls.gorilla <- delly.calls.gorilla[delly.calls.gorilla$V7 == 'PASS',]
end.pos <- as.numeric( strapplyc(as.character(delly.calls.gorilla$V8), ";END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(delly.calls.gorilla$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
delly.calls.gorilla.gr <- GRanges(seqnames=delly.calls.gorilla$V1, ranges=IRanges(start=delly.calls.gorilla$V2, end=end.pos), gen=gen, ID='gorilla', study='DELLY')
delly.calls.gorilla.gr <- keepSeqlevels(delly.calls.gorilla.gr, paste0('chr', c(1:22, 'X')), pruning.mode = 'coarse')

#orangutan
delly.calls.orangutan <- read.table("/home/porubsky/WORK/Great_apes/Delly_calls/orangutan_susie_delly.vcf")
delly.calls.orangutan <- delly.calls.orangutan[delly.calls.orangutan$V7 == 'PASS',]
end.pos <- as.numeric( strapplyc(as.character(delly.calls.orangutan$V8), ";END=(\\d+)", simplify = T) )
gen <- strapplyc(as.character(delly.calls.orangutan$V10), "^(.+?):", simplify = T)
gen <- do.call(rbind, strsplit(gen, "/"))
gen <- ifelse(gen[,1] == gen[,2], 'HOM', 'HET')
delly.calls.orangutan.gr <- GRanges(seqnames=delly.calls.orangutan$V1, ranges=IRanges(start=delly.calls.orangutan$V2, end=end.pos), gen=gen, ID='orangutan', study='DELLY')
delly.calls.orangutan.gr <- keepSeqlevels(delly.calls.orangutan.gr, paste0('chr', c(1:22, 'X')), pruning.mode = 'coarse')

#Export all PBSV calls
outputfolder <- "/home/porubsky/WORK/Great_apes/Delly_calls/"
all.delly.calls <- c(delly.calls.chimpanzee.gr, delly.calls.bonobo.gr, delly.calls.gorilla.gr, delly.calls.orangutan.gr)
destination <- file.path(outputfolder, 'DELLY_calls_greatApes.RData')
save(all.delly.calls, file = destination)

message("DONE!!!")