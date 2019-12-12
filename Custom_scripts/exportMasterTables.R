## Export master table for each great ape ##
############################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Params
thresh <- 50 #required reciprocal overlap (%)
## Set output directory
outputfolder <- "/home/porubsky/WORK/Great_apes/MasterTables"
if (!dir.exists(outputfolder)) {
  dir.create(outputfolder)
}

## Read-in published data
published.invCalls <- get(load("/home/porubsky/WORK/Great_apes/Published_data/published.invCalls.RData"))
## Read-in BioNano data
bioNano.invCalls <- get(load("/home/porubsky/WORK/Great_apes/BioNano_calls/BioNano.invCalls.RData"))
## Read-in SmartieSV data
#smartieSV.invCalls <- get(load("/home/porubsky/WORK/Great_apes/SmartieSV_calls/smartieSV.invCalls.RData"))
## Read-in PBSV data
pbsv.invCalls <- get(load("/home/porubsky/WORK/Great_apes/PacBio_calls/PBSV_calls_greatApes.RData"))
## Read-in SNIFFLES data
sniffles.invCalls <- get(load("/home/porubsky/WORK/Great_apes/PacBio_calls/SNIFFLES_calls_greatApes.RData"))
## Read-in DELLY data
delly.invCalls <- get(load("/home/porubsky/WORK/Great_apes/Delly_calls/DELLY_calls_greatApes.RData"))

## Merge all available inversion calls from different sources than StrandS
suppressWarnings( all.avail.invCalls <- c(published.invCalls, bioNano.invCalls, pbsv.invCalls, sniffles.invCalls, delly.invCalls) )

## Load Strand-seq calls after filtering and merging with validated uncertain calls ##
######################################################################################
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
all.inversion.calls.filt.grl <- split(all.inversion.calls.filt, all.inversion.calls.filt$ID)

## Chimpanzee ##
################
stranS.callsChimp.gr <- all.inversion.calls.filt.grl[['chimpanzee']]
## Subset available inv calls for chimpanzee
avail.invCalls.chimp <- all.avail.invCalls[grepl(all.avail.invCalls$ID, pattern = "chimpanzee|PericenINV")]
avail.invCalls.chimp.grl <- split(avail.invCalls.chimp, avail.invCalls.chimp$study)

message("Comparing chimpanzee data ...")
query <- stranS.callsChimp.gr
for (i in seq_along(avail.invCalls.chimp.grl)) {
  ID <- unique(avail.invCalls.chimp.grl[[i]]$study)
  message("    Processing ", ID, " data ...")
  subject <- avail.invCalls.chimp.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = unique(subject$study))
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
## Add column saying which inversion have been found in published datasets
remove.sets <- which(grepl(names(perc.overlaps), pattern = 'BioNano|DELLY|PBSV|SNIFFLES'))
perc.overlaps.published <- perc.overlaps[,-remove.sets]
query$published <- apply(perc.overlaps.published, 1, function(x) any(x >= thresh))
## Add column saying which inversion have been validated
query$valid50 <- apply(perc.overlaps, 1, function(x) any(x >= 50))
query$valid90 <- apply(perc.overlaps, 1, function(x) any(x >= 90))
## Get WC counts for each inverted regions
message("    Counting W & C reads ...")
chimpanzee.composite <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
read.counts <- list()
for (i in seq_along(stranS.callsChimp.gr)) {
  inv.region <- stranS.callsChimp.gr[i]
  reads.subset <- subsetByOverlaps(chimpanzee.composite, inv.region)
  region.read.counts <- table(strand(reads.subset))[c('+','-')]
  read.counts[[i]] <- region.read.counts
}
read.counts.df <- do.call(rbind, read.counts)
colnames(read.counts.df) <- dplyr::recode(colnames(read.counts.df), '+' = 'Crick', '-' = 'Watson')
## Export chimpanzee comparison 
chimpanzee.df.export <- as(query, 'data.frame')
chimpanzee.df.export <- cbind(chimpanzee.df.export, read.counts.df)
destination <- file.path(outputfolder, 'chimpanzee_master_table.csv')
write.table(chimpanzee.df.export , file = destination, quote = FALSE, sep = ",", row.names = FALSE)


## Bonobo ##
############
stranS.callsBonobo.gr <- all.inversion.calls.filt.grl[['bonobo']]
## Subset available inv calls for bonobo
avail.invCalls.bonobo <- all.avail.invCalls[grepl(all.avail.invCalls$ID, pattern = "bonobo|PericenINV")]
avail.invCalls.bonobo.grl <- split(avail.invCalls.bonobo, avail.invCalls.bonobo$study)

message("Comparing bonobo data ...")
query <- stranS.callsBonobo.gr
for (i in seq_along(avail.invCalls.bonobo.grl)) {
  ID <- unique(avail.invCalls.bonobo.grl[[i]]$study)
  message("    Processing ", ID, " data ...")
  subject <- avail.invCalls.bonobo.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = unique(subject$study))
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
## Add column saying which inversion have been found in published datasets
remove.sets <- which(grepl(names(perc.overlaps), pattern = 'BioNano|DELLY|PBSV|SNIFFLES'))
perc.overlaps.published <- perc.overlaps[,-remove.sets]
query$published <- apply(perc.overlaps.published, 1, function(x) any(x >= thresh))
## Add column saying which inversion have been validated
query$valid50 <- apply(perc.overlaps, 1, function(x) any(x >= 50))
query$valid90 <- apply(perc.overlaps, 1, function(x) any(x >= 90))
## Get WC counts for each inverted regions
message("    Counting W & C reads ...")
bonobo.composite <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
read.counts <- list()
for (i in seq_along(stranS.callsBonobo.gr)) {
  inv.region <- stranS.callsBonobo.gr[i]
  reads.subset <- subsetByOverlaps(bonobo.composite, inv.region)
  region.read.counts <- table(strand(reads.subset))[c('+','-')]
  read.counts[[i]] <- region.read.counts
}
read.counts.df <- do.call(rbind, read.counts)
colnames(read.counts.df) <- dplyr::recode(colnames(read.counts.df), '+' = 'Crick', '-' = 'Watson')
## Export bonobo comparison 
bonobo.df.export <- as(query, 'data.frame')
bonobo.df.export <- cbind(bonobo.df.export, read.counts.df)
destination <- file.path(outputfolder, 'bonobo_master_table.csv')
write.table(bonobo.df.export , file = destination, quote = FALSE, sep = ",", row.names = FALSE)


## Orangutan ##
###############
stranS.callsOrangutan.gr <- all.inversion.calls.filt.grl[['orangutan']]
## Subset available inv calls for orangutan
avail.invCalls.orangutan <- all.avail.invCalls[grepl(all.avail.invCalls$ID, pattern = "orangutan|PericenINV")]
avail.invCalls.orangutan.grl <- split(avail.invCalls.orangutan, avail.invCalls.orangutan$study)

message("Comparing orangutan data ...")
query <- stranS.callsOrangutan.gr
for (i in seq_along(avail.invCalls.orangutan.grl)) {
  ID <- unique(avail.invCalls.orangutan.grl[[i]]$study)
  message("    Processing ", ID, " data ...")
  subject <- avail.invCalls.orangutan.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = unique(subject$study))
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
## Add column saying which inversion have been found in published datasets
remove.sets <- which(grepl(names(perc.overlaps), pattern = 'BioNano|DELLY|PBSV|SNIFFLES'))
perc.overlaps.published <- perc.overlaps[,-remove.sets]
query$published <- apply(perc.overlaps.published, 1, function(x) any(x >= thresh))
## Add column saying which inversion have been validated
query$valid50 <- apply(perc.overlaps, 1, function(x) any(x >= 50))
query$valid90 <- apply(perc.overlaps, 1, function(x) any(x >= 90))
## Get WC counts for each inverted regions
message("    Counting W & C reads ...")
orangutan.composite <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
read.counts <- list()
for (i in seq_along(stranS.callsOrangutan.gr)) {
  inv.region <- stranS.callsOrangutan.gr[i]
  reads.subset <- subsetByOverlaps(orangutan.composite, inv.region)
  region.read.counts <- table(strand(reads.subset))[c('+','-')]
  read.counts[[i]] <- region.read.counts
}
read.counts.df <- do.call(rbind, read.counts)
colnames(read.counts.df) <- dplyr::recode(colnames(read.counts.df), '+' = 'Crick', '-' = 'Watson')
## Export orangutan comparison 
orangutan.df.export <- as(query, 'data.frame')
orangutan.df.export <- cbind(orangutan.df.export, read.counts.df)
destination <- file.path(outputfolder, 'orangutan_master_table.csv')
write.table(orangutan.df.export , file = destination, quote = FALSE, sep = ",", row.names = FALSE)


## Gorilla ##
#############
stranS.callsGorilla.gr <- all.inversion.calls.filt.grl[['gorilla']]
## Subset available inv calls for gorilla
avail.invCalls.gorilla <- all.avail.invCalls[grepl(all.avail.invCalls$ID, pattern = "gorilla|PericenINV")]
avail.invCalls.gorilla.grl <- split(avail.invCalls.gorilla, avail.invCalls.gorilla$study)

message("Comparing gorilla data ...")
query <- stranS.callsGorilla.gr
for (i in seq_along(avail.invCalls.gorilla.grl)) {
  ID <- unique(avail.invCalls.gorilla.grl[[i]]$study)
  message("    Processing ", ID, " data ...")
  subject <- avail.invCalls.gorilla.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = unique(subject$study))
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
## Add column saying which inversion have been found in published datasets
remove.sets <- which(grepl(names(perc.overlaps), pattern = 'BioNano|DELLY|PBSV|SNIFFLES'))
perc.overlaps.published <- perc.overlaps[,-remove.sets]
query$published <- apply(perc.overlaps.published, 1, function(x) any(x >= thresh))
## Add column saying which inversion have been validated
query$valid50 <- apply(perc.overlaps, 1, function(x) any(x >= 50))
query$valid90 <- apply(perc.overlaps, 1, function(x) any(x >= 90))
## Get WC counts for each inverted regions
message("    Counting W & C reads ...")
gorilla.composite <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
read.counts <- list()
for (i in seq_along(stranS.callsGorilla.gr)) {
  inv.region <- stranS.callsGorilla.gr[i]
  reads.subset <- subsetByOverlaps(gorilla.composite, inv.region)
  region.read.counts <- table(strand(reads.subset))[c('+','-')]
  read.counts[[i]] <- region.read.counts
}
read.counts.df <- do.call(rbind, read.counts)
colnames(read.counts.df) <- dplyr::recode(colnames(read.counts.df), '+' = 'Crick', '-' = 'Watson')
## Export orangutan comparison 
gorilla.df.export <- as(query, 'data.frame')
gorilla.df.export <- cbind(gorilla.df.export, read.counts.df)
destination <- file.path(outputfolder, 'gorilla_master_table.csv')
write.table(gorilla.df.export , file = destination, quote = FALSE, sep = ",", row.names = FALSE)


## Export unvalidated simple inversions and inverted duplications ##
####################################################################
message("Reporting unvalidated calls for dotplot analysis ...")
chimpanzee.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/chimpanzee_master_table.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
bonobo.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/bonobo_master_table.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
gorilla.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/gorilla_master_table.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
orangutan.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/orangutan_master_table.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

## Chimpanzee
# Recode T/F into valid vs unvalid
chimpanzee.master.table$valid <- dplyr::recode(as.character(chimpanzee.master.table$valid50), "TRUE"="valid", "FALSE"="unvalid")
chimpanzee.master.table.unvalid <- chimpanzee.master.table[chimpanzee.master.table$valid50 == 'unvalid' & chimpanzee.master.table$SVclass == 'INV',]
destination <- file.path(outputfolder, 'chimpanzee_unvalid_simpleINV.bed')
write.table(chimpanzee.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
chimpanzee.master.table.unvalid <- chimpanzee.master.table[chimpanzee.master.table$valid == 'unvalid' & chimpanzee.master.table$SVclass == 'invDup',]
destination <- file.path(outputfolder, 'chimpanzee_unvalid_invDup.bed')
write.table(chimpanzee.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Bonobo
# Recode T/F into valid vs unvalid
bonobo.master.table$valid <- dplyr::recode(as.character(bonobo.master.table$valid50), "TRUE"="valid", "FALSE"="unvalid")
bonobo.master.table.unvalid <- bonobo.master.table[bonobo.master.table$valid50 == 'unvalid' & bonobo.master.table$SVclass == 'INV',]
destination <- file.path(outputfolder, 'bonobo_unvalid_simpleINV.bed')
write.table(bonobo.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
bonobo.master.table.unvalid <- bonobo.master.table[bonobo.master.table$valid == 'unvalid' & bonobo.master.table$SVclass == 'invDup',]
destination <- file.path(outputfolder, 'bonobo_unvalid_invDup.bed')
write.table(bonobo.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Gorilla
# Recode T/F into valid vs unvalid
gorilla.master.table$valid <- dplyr::recode(as.character(gorilla.master.table$valid50), "TRUE"="valid", "FALSE"="unvalid")
gorilla.master.table.unvalid <- gorilla.master.table[gorilla.master.table$valid50 == 'unvalid' & gorilla.master.table$SVclass == 'INV',]
destination <- file.path(outputfolder, 'gorilla_unvalid_simpleINV.bed')
write.table(gorilla.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
gorilla.master.table.unvalid <- gorilla.master.table[gorilla.master.table$valid == 'unvalid' & gorilla.master.table$SVclass == 'invDup',]
destination <- file.path(outputfolder, 'gorilla_unvalid_invDup.bed')
write.table(gorilla.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Orangutan
# Recode T/F into valid vs unvalid
orangutan.master.table$valid <- dplyr::recode(as.character(orangutan.master.table$valid50), "TRUE"="valid", "FALSE"="unvalid")
orangutan.master.table.unvalid <- orangutan.master.table[orangutan.master.table$valid50 == 'unvalid' & orangutan.master.table$SVclass == 'INV',]
destination <- file.path(outputfolder, 'orangutan_unvalid_simpleINV.bed')
write.table(orangutan.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
orangutan.master.table.unvalid <- orangutan.master.table[orangutan.master.table$valid == 'unvalid' & orangutan.master.table$SVclass == 'invDup',]
destination <- file.path(outputfolder, 'orangutan_unvalid_invDup.bed')
write.table(orangutan.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

message("DONE!!!")
