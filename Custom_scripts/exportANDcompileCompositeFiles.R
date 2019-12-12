## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/Composite_files"
if (!dir.exists(outputfolder)) {
  dir.create(outputfolder)
}

## Merge sync reads Chimpanzee ##
#################################
message("Creating composite file for chimpanzee ...")
## Read in original sync reads and split them by chromosome
sync.all <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/syncFrags_1MbBin.RData"))
sync.all.grl <- split(sync.all, seqnames(sync.all))
## Read in sync reads prepared for individual chromosomes in order to avoid noise caused by complex rearrangements such as SCEs etc.
sync.chr1 <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr1_filtAlt/syncFrags_1MbBin_chr1.RData"))
sync.chr2 <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr2_filtAlt/syncFrags_1MbBin_chr2.RData"))
sync.chr5 <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr5_filtAlt/syncFrags_1MbBin_chr5.RData"))
sync.chr6 <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr6_filtAlt/syncFrags_1MbBin_chr6.RData"))
sync.chr17 <- get(load("/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr17_filtAlt/syncFrags_1MbBin_chr17.RData"))
## Replace original sync reads (for certain chromosomes) by new sync reads based on manual data curation. See selected libraries for each chromsome in analysis README 
sync.all.grl[['chr1']] <- sync.chr1
sync.all.grl[['chr2']] <- sync.chr2
sync.all.grl[['chr5']] <- sync.chr5
sync.all.grl[['chr6']] <- sync.chr6
sync.all.grl[['chr17']] <- sync.chr17
sync.all.new <- unlist(sync.all.grl, use.names = FALSE)
## Save compiled composite file
destination <- file.path(outputfolder, "syncReads_chimpanzee_final.RData")
save(sync.all.new, file = destination)
## Export composite file for UCSC visualisation
fragments2UCSC(index = 'ChimpanzeeCompositeFile', outputDirectory = outputfolder, fragments = sync.all.new)
## Plot final composite file
chimpanzee.ideo <- plotCompositeIdeo(sync.all.new, colors = c('#deebf7','#3182bd'))

## Merge sync reads Bonobo ##
#############################
message("Creating composite file for bonobo ...")
## Read in original sync reads and split them by chromosome
sync.all <- get(load("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/syncFrags_1MbBin_bonobo.RData"))
sync.all.grl <- split(sync.all, seqnames(sync.all))
## Read in sync reads prepared for individual chromosomes in order to avoid noise caused by complex rearrangements such as SCEs etc.
sync.chr2 <- get(load("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr2_filtAlt/syncFrags_1MbBin_chr2.RData"))
sync.chr5 <- get(load("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr5_filtAlt/syncFrags_1MbBin_chr5.RData"))
sync.chr6 <- get(load("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr6_filtAlt/syncFrags_1MbBin_chr6.RData"))
sync.chr17 <- get(load("/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_chr17_filtAlt/syncFrags_1MbBin_chr17.RData"))
## Replace original sync reads (for certain chromosomes) by new sync reads based on manual data curation. See selected libraries for each chromsome in analysis README 
sync.all.grl[['chr2']] <- sync.chr2
sync.all.grl[['chr5']] <- sync.chr5
sync.all.grl[['chr6']] <- sync.chr6
sync.all.grl[['chr17']] <- sync.chr17
sync.all.new <- unlist(sync.all.grl, use.names = FALSE)
## Save compiled composite file
destination <- file.path(outputfolder, "syncReads_bonobo_final.RData")
save(sync.all.new, file = destination)
## Export composite file for UCSC visualisation
fragments2UCSC(index = 'BonoboCompositeFile', outputDirectory = outputfolder, fragments = sync.all.new)
## Plot final composite file
bonobo.ideo <- plotCompositeIdeo(sync.all.new, colors = c('#e5f5e0','#31a354'))


## Merge sync reads Gorilla ##
##############################
message("Creating composite file for gorilla ...")
## Read in original sync reads and split them by chromosome
sync.all <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncFrags_1MbBin.RData"))
sync.all.grl <- split(sync.all, seqnames(sync.all))
## Read in sync reads prepared for individual chromosomes in order to avoid noise caused by complex rearrangements such as SCEs etc.
sync.chr1 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr1_filtAlt/syncFrags_1MbBin_chr1.RData"))
sync.chr2 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr2_filtAlt/syncFrags_1MbBin_chr2.RData"))
sync.chr3 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_500KbBin_chr3_filtAlt/syncFrags_500kbbBin_chr3.RData"))
sync.chr5 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr5_filtAlt/syncFrags_1MbBin_chr5.RData"))
sync.chr10 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr10_filtAlt/syncFrags_1MbBin_chr10.RData"))
sync.chr12 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr12_filtAlt/syncFrags_1MbBin_chr12.RData"))
sync.chr14 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr14_filtAlt/syncFrags_1MbBin_chr14.RData"))
sync.chr17 <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chr17_filtAlt/syncFrags_1MbBin_chr17.RData"))
sync.chrX <- get(load("/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_chrX_filtAlt/syncFrags_1MbBin_chrX.RData"))
## Replace original sync reads (for certain chromosomes) by new sync reads based on manual data curation. See selected libraries for each chromsome in analysis README 
sync.all.grl[['chr1']] <- sync.chr1
sync.all.grl[['chr2']] <- sync.chr2
sync.all.grl[['chr3']] <- sync.chr3
sync.all.grl[['chr5']] <- sync.chr5
sync.all.grl[['chr10']] <- sync.chr10
sync.all.grl[['chr12']] <- sync.chr12
sync.all.grl[['chr14']] <- sync.chr14
sync.all.grl[['chr17']] <- sync.chr17
sync.all.grl[['chrX']] <- sync.chrX
sync.all.new <- unlist(sync.all.grl, use.names = FALSE)
## Save compiled composite file
destination <- file.path(outputfolder, "syncReads_gorilla_final.RData")
save(sync.all.new, file = destination)
## Export composite file for UCSC visualisation
fragments2UCSC(index = 'GorillaCompositeFile', outputDirectory = outputfolder, fragments = sync.all.new)
## Plot final composite file
gorilla.ideo <- plotCompositeIdeo(sync.all.new, colors = c('#e0ecf4','#8856a7'))


## Merge sync reads Orangutan ##
################################
message("Creating composite file for orangutan ...")
## Read in original sync reads and split them by chromosome
sync.all <- get(load("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncFrags_1MbBin.RData"))
sync.all.grl <- split(sync.all, seqnames(sync.all))
## Read in sync reads prepared for individual chromosomes in order to avoid noise caused by complex rearrangements such as SCEs etc.
sync.chr2 <- get(load("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_chr2_filtAlt/syncFrags_1MbBin_chr2.RData"))
sync.chr7 <- get(load("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_chr7_filtAlt/syncFrags_1MbBin_chr7.RData"))
sync.chr11 <- get(load("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_chr11_filtAlt/syncFrags_1MbBin_chr11.RData"))
sync.chrX <- get(load("/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_chrX_filtAlt/syncFrags_1MbBin_chrX.RData"))
## Replace original sync reads (for certain chromosomes) by new sync reads based on manual data curation. See selected libraries for each chromsome in analysis README 
sync.all.grl[['chr2']] <- sync.chr2
sync.all.grl[['chr7']] <- sync.chr7
sync.all.grl[['chr11']] <- sync.chr11
sync.all.grl[['chrX']] <- sync.chrX
sync.all.new <- unlist(sync.all.grl, use.names = FALSE)
## Save compiled composite file
destination <- file.path(outputfolder, "syncReads_orangutan_final.RData")
save(sync.all.new, file = destination)
## Export composite file for UCSC visualisation
fragments2UCSC(index = 'OrangutanCompositeFile', outputDirectory = outputfolder, fragments = sync.all.new)
## Plot final composite file
orangutan.ideo <- plotCompositeIdeo(sync.all.new, colors = c('#fee6ce','#e6550d'))

## Export and save all plots ##
###############################
message("Plotting all composite files ...")
plt.list <- list(chimpanzee.ideo, bonobo.ideo, gorilla.ideo, orangutan.ideo)
## Save plots in RData file
destination <- file.path(outputfolder, "compositeFilesPlots.RData")
save(plt.list, file = destination)
## Save final plot in pdf
final.plot <- cowplot::plot_grid(plotlist = plt.list, nrow = 2)
destination <- file.path(outputfolder, "compositeFilesPlots.pdf")
ggsave(final.plot, file = destination, width = 20, height = 15)

message("DONE!!!")