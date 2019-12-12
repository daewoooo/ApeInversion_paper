## Prepare dot plots from nucmer outputs ##
###########################################

## Load required libraries
library(primatR)

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)
seg.dup.gr <- reduce(seg.dup.gr)

## Process through all nucmer outputfiles
## chimpanzee ##
message("Preparing chimpanzee dotplots ...")
nuc.files <- list.files("/home/porubsky/WORK/Great_apes/Nucmer_dotplots/chimpanzee/nucmer_out/", pattern = "\\.coords", full.names = TRUE)
dotplots <- list()
for (file in nuc.files) {
  filename <- basename(file)
  message("  Processing ", filename, " ...")
  
  ## Get region to highlight
  region <- as.numeric( strsplit(filename, "_|\\.")[[1]][2:3] )
  ## Plot the data
  plt <- plotNucmerCoords(nucmer.coords = file, genome.coord = TRUE, highlight.pos = region, title = filename, sd.track = seg.dup.gr)
  dotplots[[length(dotplots) + 1]] <- plt
}
## Export final plots in PDF
message("  Printing to PDF ...")
filename = "/home/porubsky/WORK/Great_apes/Nucmer_dotplots/chimpanzee/chimpanzee_dotplots.pdf"
grDevices::pdf(filename, width=5, height=5)
bquiet = lapply(dotplots, print)
d <- grDevices::dev.off()


## bonobo ##
message("Preparing bonobo dotplots ...")
nuc.files <- list.files("/home/porubsky/WORK/Great_apes/Nucmer_dotplots/bonobo/nucmer_out/", pattern = "\\.coords", full.names = TRUE)
dotplots <- list()
for (file in nuc.files) {
  filename <- basename(file)
  message("  Processing ", filename, " ...")
  
  ## Get region to highlight
  region <- as.numeric( strsplit(filename, "_|\\.")[[1]][2:3] )
  ## Plot the data
  plt <- plotNucmerCoords(nucmer.coords = file, genome.coord = TRUE, highlight.pos = region, title = filename, sd.track = seg.dup.gr)
  dotplots[[length(dotplots) + 1]] <- plt
}
## Export final plots in PDF
message("  Printing to PDF ...")
filename = "/home/porubsky/WORK/Great_apes/Nucmer_dotplots/bonobo/bonobo_dotplots.pdf"
grDevices::pdf(filename, width=5, height=5)
bquiet = lapply(dotplots, print)
d <- grDevices::dev.off()


## gorilla ##
message("Preparing gorilla dotplots ...")
nuc.files <- list.files("/home/porubsky/WORK/Great_apes/Nucmer_dotplots/gorilla/nucmer_out/", pattern = "\\.coords", full.names = TRUE)
dotplots <- list()
for (file in nuc.files) {
  filename <- basename(file)
  message("  Processing ", filename, " ...")
  
  ## Get region to highlight
  region <- as.numeric( strsplit(filename, "_|\\.")[[1]][2:3] )
  ## Plot the data
  plt <- plotNucmerCoords(nucmer.coords = file, genome.coord = TRUE, highlight.pos = region, title = filename, sd.track = seg.dup.gr)
  dotplots[[length(dotplots) + 1]] <- plt
}
## Export final plots in PDF
message("  Printing to PDF ...")
filename = "/home/porubsky/WORK/Great_apes/Nucmer_dotplots/gorilla/gorilla_dotplots.pdf"
grDevices::pdf(filename, width=5, height=5)
bquiet = lapply(dotplots, print)
d <- grDevices::dev.off()


## orangutan ##
message("Preparing orangutan dotplots ...")
nuc.files <- list.files("/home/porubsky/WORK/Great_apes/Nucmer_dotplots/orangutan/nucmer_out/", pattern = "\\.coords", full.names = TRUE)
dotplots <- list()
for (file in nuc.files) {
  filename <- basename(file)
  message("  Processing ", filename, " ...")
  
  ## Get region to highlight
  region <- as.numeric( strsplit(filename, "_|\\.")[[1]][2:3] )
  ## Plot the data
  plt <- plotNucmerCoords(nucmer.coords = file, genome.coord = TRUE, highlight.pos = region, title = filename, sd.track = seg.dup.gr)
  dotplots[[length(dotplots) + 1]] <- plt
}
## Export final plots in PDF
message("  Printing to PDF ...")
filename = "/home/porubsky/WORK/Great_apes/Nucmer_dotplots/orangutan/orangutan_dotplots.pdf"
grDevices::pdf(filename, width=5, height=5)
bquiet = lapply(dotplots, print)
d <- grDevices::dev.off()
