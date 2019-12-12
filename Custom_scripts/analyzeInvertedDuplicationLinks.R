## This script detects links between donor and accptor sites of inverted duplications ##
########################################################################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38) )

message("Analysing inter-chromosomal links between invDups ...")

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Analyze links for gorilla ##
###############################
inputfolder = "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/gorilla/bams_GRCh38align/"
linkedReads.gorilla <- readPairsAsLinks(inputfolder = inputfolder
                                        , min.mapq = 60
                                        , filt.flag = 3328
                                        , min.reads = 10
                                        , chromosomes = paste0('chr', c(1:22, 'X'))
                                        , bsgenome = BSgenome.Hsapiens.UCSC.hg38
                                        , blacklist = seg.dup.gr)

signif.links.gorilla <- processReadLinks(gr.links = linkedReads.gorilla
                                         , min.reads = 10
                                         , chromosomes = paste0('chr', c(1:22, 'X'))
                                         , bsgenome = BSgenome.Hsapiens.UCSC.hg38)

plt.gorilla <- plotLinks(links = signif.links.gorilla, 
                         chromosomes = paste0('chr', c(1:22, 'X')),
                         index = 'Gorilla',
                         bsgenome = BSgenome.Hsapiens.UCSC.hg38)

destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/gorilla/"
ggsave(plt.gorilla, file = file.path(destination, "gorilla_GRCh38align_links.pdf"))
destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/gorilla/"
save(signif.links.gorilla, file = file.path(destination, "gorilla_GRCh38align_links.RData"))
write.table(as.data.frame(signif.links.gorilla$intra.links), file = file.path(destination, "gorilla_GRCh38align_intraLinks.txt"), quote = FALSE, row.names = FALSE)
write.table(as.data.frame(signif.links.gorilla$inter.links), file = file.path(destination, "gorilla_GRCh38align_interLinks.txt"), quote = FALSE, row.names = FALSE)

## gorilla export bedgraphs
bams <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
# Set outputdirectory
outputdirectory <- file.path(inputfolder, "bams_GRCh38align_bedgraphs/")
if (!dir.exists(outputdirectory)) {
  dir.create(outputdirectory)
}

for (bam in bams) {
  bam.name <- basename(bam)
  message("Exporting bedgraph for ", bam.name)
  ## Export split read maps to a bedGraph
  exportBedGraph(bamfile = bam
                 , outputdirectory = outputdirectory
                 , mapq = 60
                 , filt.flag = 3328
                 , min.read.len = 2000
                 , blacklist = seg.dup.gr
  )
}


## Analyze links for orangutan ##
#################################
inputfolder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/orangutan/bams_GRCh38align/"
linkedReads.orangutan <- readPairsAsLinks(inputfolder = inputfolder 
                                        , min.mapq = 60
                                        , filt.flag = 3328
                                        , min.reads = 10
                                        , chromosomes = paste0('chr', c(1:22, 'X'))
                                        , bsgenome = BSgenome.Hsapiens.UCSC.hg38
                                        , blacklist = seg.dup.gr)

signif.links.orangutan <- processReadLinks(gr.links = linkedReads.orangutan
                                         , min.reads = 10
                                         , chromosomes = paste0('chr', c(1:22, 'X'))
                                         , bsgenome = BSgenome.Hsapiens.UCSC.hg38)

plt.orangutan <- plotLinks(links = signif.links.orangutan,
                           chromosomes = paste0('chr', c(1:22, 'X')),
                           index = 'Orangutan',
                           bsgenome = BSgenome.Hsapiens.UCSC.hg38)

destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/orangutan/"
ggsave(plt.orangutan , file = file.path(destination, "orangutan_GRCh38align_links.pdf"))
destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/orangutan"
save(signif.links.orangutan, file = file.path(destination, "orangutan_GRCh38align_links.RData"))
write.table(signif.links.orangutan$intra.links, file = file.path(destination, "orangutan_GRCh38align_intraLinks.txt"), quote = FALSE, row.names = FALSE)
write.table(signif.links.orangutan$inter.links, file = file.path(destination, "orangutan_GRCh38align_interLinks.txt"), quote = FALSE, row.names = FALSE)

## orangutan export bedgraphs
bams <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
# Set outputdirectory
outputdirectory <- file.path(inputfolder, "bams_GRCh38align_bedgraphs/")
if (!dir.exists(outputdirectory)) {
  dir.create(outputdirectory)
}

for (bam in bams) {
  bam.name <- basename(bam)
  message("Exporting bedgraph for ", bam.name)
  ## Export split read maps to a bedGraph
  exportBedGraph(bamfile = bam
                 , outputdirectory = outputdirectory
                 , mapq = 60
                 , filt.flag = 3328
                 , min.read.len = 2000
                 , blacklist = seg.dup.gr
  )
}


## Analyze links for chimpanzee ##
##################################
inputfolder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/chimpanzee/bams_GRCh38align/"
test <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/chimpanzee/bams_GRCh38align/TEST/"
linkedReads.chimpanzee <- readPairsAsLinks(inputfolder = inputfolder
                                       , min.mapq = 60
                                       , filt.flag = 3328
                                       , min.reads = 10
                                       , chromosomes = paste0('chr', c(1:22, 'X'))
                                       , bsgenome = BSgenome.Hsapiens.UCSC.hg38
                                       , blacklist = seg.dup.gr)

signif.links.chimpanzee <- processReadLinks(gr.links = linkedReads.chimpanzee
                                        , min.reads = 10
                                        , chromosomes = paste0('chr', c(1:22, 'X'))
                                        , bsgenome = BSgenome.Hsapiens.UCSC.hg38)

plt.chimpanzee <- plotLinks(links = signif.links.chimpanzee,
                            chromosomes = paste0('chr', c(1:22, 'X')),
                            index = 'Chimpanzee',
                            bsgenome = BSgenome.Hsapiens.UCSC.hg38)

destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/chimpanzee/"
ggsave(plt.chimpanzee, file = file.path(destination, "chimpanzee_GRCh38align_links.pdf"))
destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/chimpanzee/"
save(signif.links.chimpanzee, file = file.path(destination, "chimpanzee_GRCh38align_links.RData"))
write.table(signif.links.chimpanzee$intra.links, file = file.path(destination, "chimpanzee_GRCh38align_intraLinks.txt"), quote = FALSE, row.names = FALSE)
write.table(signif.links.chimpanzee$inter.links, file = file.path(destination, "chimpanzee_GRCh38align_interLinks.txt"), quote = FALSE, row.names = FALSE)

## chimpanzee export bedgraphs
bams <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
# Set outputdirectory
outputdirectory <- file.path(inputfolder, "bams_GRCh38align_bedgraphs/")
if (!dir.exists(outputdirectory)) {
  dir.create(outputdirectory)
}

for (bam in bams) {
  bam.name <- basename(bam)
  message("Exporting bedgraph for ", bam.name)
  ## Export split read maps to a bedGraph
  exportBedGraph(bamfile = bam
                 , outputdirectory = outputdirectory
                 , mapq = 60
                 , filt.flag = 3328
                 , min.read.len = 2000
                 , blacklist = seg.dup.gr
  )
}

## Analyze links for bonobo ##
##############################
inputfolder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/bonobo/bams_GRCh38align/"
linkedReads.bonobo <- readPairsAsLinks(inputfolder = inputfolder
                                        , min.mapq = 60
                                        , filt.flag = 3328
                                        , min.reads = 10
                                        , chromosomes = paste0('chr', c(1:22, 'X'))
                                        , bsgenome = BSgenome.Hsapiens.UCSC.hg38
                                        , blacklist = seg.dup.gr)

signif.links.bonobo <- processReadLinks(gr.links = linkedReads.bonobo
                                         , min.reads = 10
                                         , chromosomes = paste0('chr', c(1:22, 'X'))
                                         , bsgenome = BSgenome.Hsapiens.UCSC.hg38)

plt.bonobo <- plotLinks(links = signif.links.bonobo,
                        chromosomes = paste0('chr', c(1:22, 'X')),
                        index = 'Bonobo',
                        bsgenome = BSgenome.Hsapiens.UCSC.hg38)

destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/bonobo/"
ggsave(plt.bonobo, file = file.path(destination, "bonobo_GRCh38align_links.pdf"))
destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/bonobo/"
save(signif.links.bonobo, file = file.path(destination, "bonobo_GRCh38align_links.RData"))
write.table(signif.links.bonobo$intra.links, file = file.path(destination, "bonobo_GRCh38align_intraLinks.txt"), quote = FALSE, row.names = FALSE)
write.table(signif.links.bonobo$inter.links, file = file.path(destination, "bonobo_GRCh38align_interLinks.txt"), quote = FALSE, row.names = FALSE)

## bonobo export bedgraphs
bams <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
# Set outputdirectory
outputdirectory <- file.path(inputfolder, "bams_GRCh38align_bedgraphs/")
if (!dir.exists(outputdirectory)) {
  dir.create(outputdirectory)
}

for (bam in bams) {
  bam.name <- basename(bam)
  message("Exporting bedgraph for ", bam.name)
  ## Export split read maps to a bedGraph
  exportBedGraph(bamfile = bam
                 , outputdirectory = outputdirectory
                 , mapq = 60
                 , filt.flag = 3328
                 , min.read.len = 2000
                 , blacklist = seg.dup.gr
  )
}

## Compile all plots together
final.plot <- plot_grid(plt.chimpanzee, plt.bonobo, plt.gorilla, plt.orangutan, nrow = 2)
destination <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites"
ggsave(final.plot, file = file.path(destination, "invDupLinks_allApes.pdf"), device = 'pdf', width = 20, height = 12, useDingbats=FALSE)

message("DONE!!!")

