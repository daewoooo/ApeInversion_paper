## Export breakpoints of inverted duplication calls ##
######################################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/"

message("Exporting inverted duplication calls ...")

## Read in all called inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))

## Filter only inverted duplications invDups
all.invDups.gr <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Export bed file for all great apes
all.invDups.df <- as(all.invDups.gr, 'data.frame')
all.invDups.df <- all.invDups.df[,c('seqnames', 'start','end', 'ID')]
destination <- file.path(outputfolder, 'invDup.regions.all.bed')
write.table(all.invDups.df, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE)

## Load extra manually selected breakpoints
extra.breaks <- read.table("/home/porubsky/WORK/Great_apes/InvertedDups_analysis/extra_invDups_breakpoints.txt", header = TRUE, stringsAsFactors = FALSE)
extra.breaks.l <- split(extra.breaks[,c(1:4)], extra.breaks$id)

## Export bed file per individual ##
####################################
## chimpanzee ##
chimpanzee.invDup <- all.invDups.gr[all.invDups.gr$ID == 'chimpanzee']
chimpanzee.invDup.df <- as.data.frame(chimpanzee.invDup)
chimpanzee.invDup.df$index <- paste0('roi', c(1:nrow(chimpanzee.invDup.df)))
chimpanzee.invDup.df <- dplyr::mutate(chimpanzee.invDup.df, name = paste(seqnames, SVclass, ID, index, sep = '_'))
chimpanzee.invDup.df.starts <- chimpanzee.invDup.df
chimpanzee.invDup.df.ends <- chimpanzee.invDup.df
chimpanzee.invDup.df.starts$end <- chimpanzee.invDup.df.starts$start + 1
chimpanzee.invDup.df.ends$start <- chimpanzee.invDup.df.ends$end - 1
chimpanzee.invDup.df.out <- rbind(chimpanzee.invDup.df.starts, chimpanzee.invDup.df.ends)
chimpanzee.invDup.df.out <- dplyr::arrange(chimpanzee.invDup.df.out , name)
chimpanzee.invDup.df.out <- dplyr::select(chimpanzee.invDup.df.out, seqnames, start, end, name)

## Add extra manually selected breakpoints
extra.breaks.chimp <- extra.breaks.l[['chimpanzee']]
chimpanzee.invDup.df.out <- rbind(chimpanzee.invDup.df.out, extra.breaks.chimp)

destination <- file.path(outputfolder, 'chimpanzee_invDup_breakPointList.bed')
write.table(chimpanzee.invDup.df.out, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## bonobo ##
bonobo.invDup <- all.invDups.gr[all.invDups.gr$ID == 'bonobo']
bonobo.invDup.df <- as.data.frame(bonobo.invDup)
bonobo.invDup.df$index <- paste0('roi', c(1:nrow(bonobo.invDup.df)))
bonobo.invDup.df <- dplyr::mutate(bonobo.invDup.df, name = paste(seqnames, SVclass, ID, index, sep = '_'))
bonobo.invDup.df.starts <- bonobo.invDup.df
bonobo.invDup.df.ends <- bonobo.invDup.df
bonobo.invDup.df.starts$end <- bonobo.invDup.df.starts$start + 1
bonobo.invDup.df.ends$start <- bonobo.invDup.df.ends$end - 1
bonobo.invDup.df.out <- rbind(bonobo.invDup.df.starts, bonobo.invDup.df.ends)
bonobo.invDup.df.out <- dplyr::arrange(bonobo.invDup.df.out , name)
bonobo.invDup.df.out <- dplyr::select(bonobo.invDup.df.out, seqnames, start, end, name)

## Add extra manually selected breakpoints
extra.breaks.bonobo <- extra.breaks.l[['bonobo']]
bonobo.invDup.df.out <- rbind(bonobo.invDup.df.out, extra.breaks.bonobo)

destination <- file.path(outputfolder, 'bonobo_invDup_breakPointList.bed')
write.table(bonobo.invDup.df.out, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## gorilla ##
gorilla.invDup <- all.invDups.gr[all.invDups.gr$ID == 'gorilla']
gorilla.invDup.df <- as.data.frame(gorilla.invDup)
gorilla.invDup.df$index <- paste0('roi', c(1:nrow(gorilla.invDup.df)))
gorilla.invDup.df <- dplyr::mutate(gorilla.invDup.df, name = paste(seqnames, SVclass, ID, index, sep = '_'))
gorilla.invDup.df.starts <- gorilla.invDup.df
gorilla.invDup.df.ends <- gorilla.invDup.df
gorilla.invDup.df.starts$end <- gorilla.invDup.df.starts$start + 1
gorilla.invDup.df.ends$start <- gorilla.invDup.df.ends$end - 1
gorilla.invDup.df.out <- rbind(gorilla.invDup.df.starts, gorilla.invDup.df.ends)
gorilla.invDup.df.out <- dplyr::arrange(gorilla.invDup.df.out , name)
gorilla.invDup.df.out <- dplyr::select(gorilla.invDup.df.out, seqnames, start, end, name)

## Add extra manually selected breakpoints
extra.breaks.gorilla <- extra.breaks.l[['gorilla']]
gorilla.invDup.df.out <- rbind(gorilla.invDup.df.out, extra.breaks.gorilla)

destination <- file.path(outputfolder, 'gorilla_invDup_breakPointList.bed')
write.table(gorilla.invDup.df.out, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## orangutan ##
orangutan.invDup <- all.invDups.gr[all.invDups.gr$ID == 'orangutan']
orangutan.invDup.df <- as.data.frame(orangutan.invDup)
orangutan.invDup.df$index <- paste0('roi', c(1:nrow(orangutan.invDup.df)))
orangutan.invDup.df <- dplyr::mutate(orangutan.invDup.df, name = paste(seqnames, SVclass, ID, index, sep = '_'))
orangutan.invDup.df.starts <- orangutan.invDup.df
orangutan.invDup.df.ends <- orangutan.invDup.df
orangutan.invDup.df.starts$end <- orangutan.invDup.df.starts$start + 1
orangutan.invDup.df.ends$start <- orangutan.invDup.df.ends$end - 1
orangutan.invDup.df.out <- rbind(orangutan.invDup.df.starts, orangutan.invDup.df.ends)
orangutan.invDup.df.out <- dplyr::arrange(orangutan.invDup.df.out , name)
orangutan.invDup.df.out <- dplyr::select(orangutan.invDup.df.out, seqnames, start, end, name)

## Add extra manually selected breakpoints
extra.breaks.orangutan <- extra.breaks.l[['orangutan']]
orangutan.invDup.df.out <- rbind(orangutan.invDup.df.out, extra.breaks.orangutan)

destination <- file.path(outputfolder, 'orangutan_invDup_breakPointList.bed')
write.table(orangutan.invDup.df.out, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

message("DONE!!!")

