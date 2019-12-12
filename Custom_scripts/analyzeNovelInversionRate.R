## Evaluate novel inversion rate based on Strand-seq inv calls ##
#################################################################
## Load required libraries
suppressPackageStartupMessages( library(primatR) )

## Load published data
published.invCalls <- get(load("/home/porubsky/WORK/Great_apes/Published_data/published.invCalls.RData"))
published.invCalls.df <- as.data.frame(published.invCalls)
published.summary <-
  published.invCalls.df %>%
  group_by(study, ID) %>%
  summarise(count=n())

## Load complete dataset including putative human specific inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
all.inversion.calls.filt.grl <- split(all.inversion.calls.filt, all.inversion.calls.filt$ID)
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
## Load set of validated simple inversions detected by Strand-seq
strandseq.valid.simpleINV <- get(load("/home/porubsky/WORK/Great_apes/MasterTables/strandseq.valid.simpleINV.RData"))

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

## Validate raw callset from published inversions ##
####################################################
## Load Zev_2018 smartieSV calls
zev.smartieSV.gr <- get(load("/home/porubsky/WORK/Great_apes/SmartieSV_calls/smartieSV.invCalls.RData"))
zev.smartieSV.gr$study <- 'zev_2018_all'

## Load Feuk_2005 all calls
feuk.callsChimp <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Feuk_2015/TableS_ All Putative_1576_Inversion_Regions.csv", sep = ",", header = TRUE)
feuk.callsChimp.gr <- GRanges(seqnames=feuk.callsChimp$Human.Chr, ranges=IRanges(start=feuk.callsChimp$Human.Start, end=feuk.callsChimp$Human.Stop))
feuk.callsChimp.gr$gen <- ""
feuk.callsChimp.gr$ID <- 'chimpanzee'
feuk.callsChimp.gr$study <- 'feuk_2005_all'

## Merge unfiltered callset from Zev_2018 and Feuk_2005
suppressWarnings( published.unfiltered.calls <- c(zev.smartieSV.gr, feuk.callsChimp.gr) )

## Merge all available inversion calls from different sources than StrandS
suppressWarnings( all.avail.invCalls <- c(bioNano.invCalls, pbsv.invCalls, sniffles.invCalls, delly.invCalls) )

## Chimpanzee ##
################
unfilt.callsChimp.gr <- published.unfiltered.calls[published.unfiltered.calls$ID == "chimpanzee"]
## Subset available inv calls for chimpanzee
avail.invCalls.chimp <- all.avail.invCalls[grepl(all.avail.invCalls$ID, pattern = "chimpanzee")]
avail.invCalls.chimp.grl <- split(avail.invCalls.chimp, avail.invCalls.chimp$study)

message("Comparing chimpanzee data ...")
query <- unfilt.callsChimp.gr
for (i in seq_along(avail.invCalls.chimp.grl)) {
  ID <- unique(avail.invCalls.chimp.grl[[i]]$study)
  message("    Processing ", ID, " data ...")
  subject <- avail.invCalls.chimp.grl[[i]]
  query <-  getReciprocalOverlaps(query = query, subject = subject, thresh = thresh, report = 'query', index = unique(subject$study))
}
perc.overlaps <- mcols( query[,which(grepl(names(mcols(query)), pattern = 'perc.overlap'))] )
query$valid <- apply(perc.overlaps, 1, function(x) any(x >= thresh))
chimpanzee.df.export <- as(query, 'data.frame')
