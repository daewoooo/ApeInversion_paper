## Load required libraries
library(primatR)
library(BSgenome.Hsapiens.UCSC.hg38)

## Load putative human specific inversions
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/putative_humanSpecific_annot.RData"))
# Remove known misorients
HSinv.gr <- HSinv.gr[HSinv.gr$misorient != TRUE]
HSinv.gr <- HSinv.gr[,0]

## Load putative polymorphic inversions
polymorphic.inversions <- get(load("/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites/putative_polymorphic_regions.RData"))

## Load chimpanzee inversion in a range from 100kb to 4MB
#all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
#all.inversion.calls.filt.grl <- split(all.inversion.calls.filt, all.inversion.calls.filt$ID)
#stranS.callsChimp.gr <- all.inversion.calls.filt.grl[['chimpanzee']]
#stranS.callsChimp.gr <- stranS.callsChimp.gr[width(stranS.callsChimp.gr) >= 100000 & width(stranS.callsChimp.gr) <= 4000000]

## Genomic regions to plot contact matrices in
#regions.gr <- stranS.callsChimp.gr
regions.gr <- c(HSinv.gr, polymorphic.inversions)

## Human HIC data NA12878
bamfile.human <- "/home/porubsky/WORK/Great_apes/HiC_analysis/GM12878_hic_mdup_filt_srt.bam"

## Chimpanzee
bamfile.chimp <- "/home/porubsky/WORK/Great_apes/HiC_analysis/chimpanzee_hic_mdup_filt_srt.bam"

region.plots <- list()
for (i in seq_along(regions.gr)) {
  inv.region <- regions.gr[i]
  inv.region.id <- as.character(inv.region)
  message("Processing region: ", inv.region.id)
  region.size <- width(inv.region)
  
  ## Extend inverted region left and right 1,5,10 or 20 times the size of inversion
  if (region.size >= 1000000) {
    region <- resizeRanges(inv.region, times = 1, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  } else if (region.size >= 100000 & region.size < 1000000) {
    region <- resizeRanges(inv.region, times = 5, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  } else if (region.size >= 50000 & region.size < 100000) {
    region <- resizeRanges(inv.region, times = 10, bsgenome = BSgenome.Hsapiens.UCSC.hg38)  
  } else if (region.size >= 10000 & region.size < 50000) {
    region <- resizeRanges(inv.region, times = 20, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  } else if (region.size < 10000) {
    region <- resizeRanges(inv.region, times = 40, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  }  
  ## Plot HIC contacts for human and non-human primate
  plot.l <- plotHICregional(bamfile = bamfile.human, region = region, min.mapq = 10, resolution = c(50000, 25000, 10000), highlight.pos = c(start(inv.region), end(inv.region)), bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  plt.human <- plot_grid(plotlist = plot.l, nrow = 1)
  plot.l <- plotHICregional(bamfile = bamfile.chimp, region = region, min.mapq = 10, resolution = c(50000, 25000, 10000), highlight.pos = c(start(inv.region), end(inv.region)), bsgenome = BSgenome.Hsapiens.UCSC.hg38)
  plt.chimp <- plot_grid(plotlist = plot.l, nrow = 1)
  ## Prepare final plot
  final.plt <- plot_grid(plt.human, plt.chimp, ncol = 1)
  region.plots[[i]] <- final.plt
}

message("Printing to PDF ...")
filename = "/home/porubsky/WORK/Great_apes/HiC_analysis/chimpanzee_vs_human_HICcontacts.pdf"
grDevices::pdf(filename, width=15, height=5)
bquiet = lapply(region.plots, print)
d <- grDevices::dev.off()
