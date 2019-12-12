## Load required libraries ##
#############################
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(UpSetR) )
#Upset plots have to be prepared manually

message("Analyzing validated inversion calls ...")

## segDup track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load exported master tables ##
#################################
chimpanzee.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/chimpanzee_master_table_checked.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
bonobo.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/bonobo_master_table_checked.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
gorilla.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/gorilla_master_table_checked.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
orangutan.master.table <- read.table("/home/porubsky/WORK/Great_apes/MasterTables/orangutan_master_table_checked.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

#chimpanzee.master.table.long <- melt(chimpanzee.master.table,
#                                     id.vars = c('seqnames','start','end','width','gen','ID','SVclass'),
#                                     measure.vars = c('perc.overlap_BioNano','perc.overlap_Catacchio_2018','perc.overlap_DELLY','perc.overlap_PBSV','perc.overlap_Prakash_1982','perc.overlap_SNIFFLES','perc.overlap_Szamalek_2006','perc.overlap_Zev_2018','perc.overlap_feuk_2015')
#)

## Set parameters
#dataset.ids <- c("BioNano","Catacchio_2018","DELLY","PBSV","Prakash_1982","SNIFFLES","Szamalek_2006","Zev_2018","feuk_2015", "Strand-seq")
thresh <- 50
export.unvalids <- "/home/porubsky/WORK/Great_apes/MasterTables/"

## Inversion validation summary ##
##################################
all.plots <- list()
## Chimpanzee ##
valid.chimpanzee <- vector()
for (i in 1:nrow(chimpanzee.master.table)) {
  row <- chimpanzee.master.table[i,]
  perc.overlap.idx <- grep(names(row), pattern = 'perc.overlap')
  perc.overlaps <- row[perc.overlap.idx]
  IDs <- gsub(names(perc.overlaps), pattern = 'perc.overlap_', replacement = '')
  valid.IDs <- IDs[which(perc.overlaps >= thresh)]
  if (length(valid.IDs) > 0) {
    to.add <- rep(i, length(valid.IDs))
    names(to.add) <- valid.IDs
    valid.chimpanzee <- c(valid.chimpanzee, to.add)
    names(to.add) <- "Strand-seq"
    valid.chimpanzee <- c(valid.chimpanzee, to.add)
  } else {
    to.add <- i
    names(to.add) <- "Strand-seq"
    valid.chimpanzee <- c(valid.chimpanzee, to.add)
  }  
}
valid.chimpanzee.list <- split(valid.chimpanzee, names(valid.chimpanzee))
# Plot the overlaps
#upset(fromList(valid.chimpanzee.list), order.by = 'freq', nsets = length(valid.chimpanzee.list))
# Get total number of validated inversions including dotlots and manuall BN calls
chimpanzee.master.table$valid.total <- apply(chimpanzee.master.table[c('published','valid50', 'valid.dotplot', 'valid.BN.manual')], 1, any)
chimpanzee.master.table$valid.total[is.na(chimpanzee.master.table$valid.total)] <- FALSE
# Recode T/F into valid vs unvalid
chimpanzee.master.table$valid.total <- recode(as.character(chimpanzee.master.table$valid.total), "TRUE"="valid", "FALSE"="unvalid")
# Plot size distribution of valid vs unvalid
plt1 <- chimpanzee.master.table %>% 
  ggplot() + geom_boxplot(aes(x=valid.total, y=width, fill=valid.total)) +
  scale_y_continuous(breaks=c(10000,100000,1000000), labels = comma) +
  coord_trans(y='log10') +
  ylab("Log10 Inversion size (bp)") +
  geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed") + 
  scale_fill_manual(values = c('gold3', 'lightcyan4'), name="") +
  xlab("") +
  theme_bw()
# Plot SVclass of valid vs invalid
plt2 <- chimpanzee.master.table %>%
  group_by(SVclass, valid.total) %>% 
  summarise(count=n()) %>% 
  mutate(ID = paste0(SVclass,"_",valid.total)) %>%
  ggplot(aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
  geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightblue3', 'lightblue4', 'gold3', 'gold4'), name="") +
  xlab("")
merge.plt <- plot_grid(plt1, plt2, nrow = 1)
title <- ggdraw() + draw_label("Chimpanzee validation", fontface='bold')
final.plt <- plot_grid(title, merge.plt, ncol=1, rel_heights=c(0.1, 1))
all.plots[[length(all.plots) + 1]] <- final.plt
# Export unvalidated simple inversions as a bed file
chimpanzee.master.table.unvalid <- chimpanzee.master.table[chimpanzee.master.table$valid.total == 'unvalid' & chimpanzee.master.table$SVclass == 'INV',]
destination <- file.path(export.unvalids, 'chimpanzee_remainUnvalid_simpleINV.bed')
write.table(chimpanzee.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Bonobo ##
valid.bonobo <- vector()
for (i in 1:nrow(bonobo.master.table)) {
  row <- bonobo.master.table[i,]
  perc.overlap.idx <- grep(names(row), pattern = 'perc.overlap')
  perc.overlaps <- row[perc.overlap.idx]
  IDs <- gsub(names(perc.overlaps), pattern = 'perc.overlap_', replacement = '')
  valid.IDs <- IDs[which(perc.overlaps >= thresh)]
  if (length(valid.IDs) > 0) {
    to.add <- rep(i, length(valid.IDs))
    names(to.add) <- valid.IDs
    valid.bonobo <- c(valid.bonobo, to.add)
    names(to.add) <- "Strand-seq"
    valid.bonobo <- c(valid.bonobo, to.add)
  } else {
    to.add <- i
    names(to.add) <- "Strand-seq"
    valid.bonobo <- c(valid.bonobo, to.add)
  }  
}
valid.bonobo.list <- split(valid.bonobo, names(valid.bonobo))
# Plot the overlaps
#upset(fromList(valid.bonobo.list), order.by = 'freq', nsets = length(valid.bonobo.list))
# Get total number of validated inversions including dotlots and manuall BN calls
bonobo.master.table$valid.total <- apply(bonobo.master.table[c('published','valid50', 'valid.dotplot', 'valid.BN.manual')], 1, any)
bonobo.master.table$valid.total[is.na(bonobo.master.table$valid.total)] <- FALSE
# Recode T/F into valid vs unvalid
bonobo.master.table$valid.total <- recode(as.character(bonobo.master.table$valid.total), "TRUE"="valid", "FALSE"="unvalid")
# Plot size distribution of valid vs invalid
plt1 <- bonobo.master.table %>% 
  ggplot() + geom_boxplot(aes(x=valid.total, y=width, fill=valid.total)) +
  scale_y_continuous(breaks=c(10000,100000,1000000), labels = comma) +
  coord_trans(y='log10') +
  ylab("Log10 Inversion size (bp)") +
  geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed") + 
  scale_fill_manual(values = c('gold3', 'lightcyan4'), name="") +
  xlab("") +
  theme_bw()
# Plot SVclass of valid vs invalid
plt2 <- bonobo.master.table %>%
  group_by(SVclass, valid.total) %>% 
  summarise(count=n()) %>% 
  mutate(ID = paste0(SVclass,"_",valid.total)) %>%
  ggplot(aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
  geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightblue3', 'lightblue4', 'gold3', 'gold4'), name="") +
  xlab("")
merge.plt <- plot_grid(plt1, plt2, nrow = 1)
title <- ggdraw() + draw_label("Bonobo validation", fontface='bold')
final.plt <- plot_grid(title, merge.plt, ncol=1, rel_heights=c(0.1, 1))
all.plots[[length(all.plots) + 1]] <- final.plt
# Export unvalidated simple inversions as a bed file
bonobo.master.table.unvalid <- bonobo.master.table[bonobo.master.table$valid.total == 'unvalid' & bonobo.master.table$SVclass == 'INV',]
destination <- file.path(export.unvalids, 'bonobo_remainUnvalid_simpleINV.bed')
write.table(bonobo.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Gorilla ##
valid.gorilla <- vector()
for (i in 1:nrow(gorilla.master.table)) {
  row <- gorilla.master.table[i,]
  perc.overlap.idx <- grep(names(row), pattern = 'perc.overlap')
  perc.overlaps <- row[perc.overlap.idx]
  IDs <- gsub(names(perc.overlaps), pattern = 'perc.overlap_', replacement = '')
  valid.IDs <- IDs[which(perc.overlaps >= thresh)]
  if (length(valid.IDs) > 0) {
    to.add <- rep(i, length(valid.IDs))
    names(to.add) <- valid.IDs
    valid.gorilla <- c(valid.gorilla, to.add)
    names(to.add) <- "Strand-seq"
    valid.gorilla <- c(valid.gorilla, to.add)
  } else {
    to.add <- i
    names(to.add) <- "Strand-seq"
    valid.gorilla <- c(valid.gorilla, to.add)
  }  
}
valid.gorilla.list <- split(valid.gorilla, names(valid.gorilla))
# Plot the overlaps
#upset(fromList(valid.gorilla.list), order.by = 'freq', nsets = length(valid.gorilla.list))
# Get total number of validated inversions including dotlots and manuall BN calls
gorilla.master.table$valid.total <- apply(gorilla.master.table[c('published','valid50', 'valid.dotplot', 'valid.BN.manual')], 1, any)
gorilla.master.table$valid.total[is.na(gorilla.master.table$valid.total)] <- FALSE
# Recode T/F into valid vs unvalid
gorilla.master.table$valid.total <- recode(as.character(gorilla.master.table$valid.total), "TRUE"="valid", "FALSE"="unvalid")
# Plot size distribution of valid vs invalid
plt1 <- gorilla.master.table %>% 
  ggplot() + geom_boxplot(aes(x=valid.total, y=width, fill=valid.total)) +
  scale_y_continuous(breaks=c(10000,100000,1000000), labels = comma) +
  coord_trans(y='log10') +
  ylab("Log10 Inversion size (bp)") +
  geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed") + 
  scale_fill_manual(values = c('gold3', 'lightcyan4'), name="") +
  xlab("") +
  theme_bw()
# Plot SVclass of valid vs invalid
plt2 <- gorilla.master.table %>%
  group_by(SVclass, valid.total) %>% 
  summarise(count=n()) %>% 
  mutate(ID = paste0(SVclass,"_",valid.total)) %>%
  ggplot(aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
  geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightblue3', 'lightblue4', 'gold3', 'gold4'), name="") +
  xlab("")
merge.plt <- plot_grid(plt1, plt2, nrow = 1)
title <- ggdraw() + draw_label("Gorilla validation", fontface='bold')
final.plt <- plot_grid(title, merge.plt, ncol=1, rel_heights=c(0.1, 1))
all.plots[[length(all.plots) + 1]] <- final.plt
# Export unvalidated simple inversions as a bed file
gorilla.master.table.unvalid <- gorilla.master.table[gorilla.master.table$valid.total == 'unvalid' & gorilla.master.table$SVclass == 'INV',]
destination <- file.path(export.unvalids, 'gorilla_remainUnvalid_simpleINV.bed')
write.table(gorilla.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Orangutan ##
valid.orangutan <- vector()
for (i in 1:nrow(orangutan.master.table)) {
  row <- orangutan.master.table[i,]
  perc.overlap.idx <- grep(names(row), pattern = 'perc.overlap')
  perc.overlaps <- row[perc.overlap.idx]
  IDs <- gsub(names(perc.overlaps), pattern = 'perc.overlap_', replacement = '')
  valid.IDs <- IDs[which(perc.overlaps >= thresh)]
  if (length(valid.IDs) > 0) {
    to.add <- rep(i, length(valid.IDs))
    names(to.add) <- valid.IDs
    valid.orangutan <- c(valid.orangutan, to.add)
    names(to.add) <- "Strand-seq"
    valid.orangutan <- c(valid.orangutan, to.add)
  } else {
    to.add <- i
    names(to.add) <- "Strand-seq"
    valid.orangutan <- c(valid.orangutan, to.add)
  }  
}
valid.orangutan.list <- split(valid.orangutan, names(valid.orangutan))
# Plot the overlaps
#upset(fromList(valid.orangutan.list), order.by = 'freq', nsets = length(valid.orangutan.list))
# Get total number of validated inversions including dotplots and manuall BN calls
orangutan.master.table$valid.total <- apply(orangutan.master.table[c('published','valid50', 'valid.dotplot', 'valid.BN.manual')], 1, any)
orangutan.master.table$valid.total[is.na(orangutan.master.table$valid.total)] <- FALSE
# Recode T/F into valid vs unvalid
orangutan.master.table$valid.total <- recode(as.character(orangutan.master.table$valid.total), "TRUE"="valid", "FALSE"="unvalid")
# Plot size distribution of valid vs invalid
plt1 <- orangutan.master.table %>% 
  ggplot() + geom_boxplot(aes(x=valid.total, y=width, fill=valid.total)) +
  scale_y_continuous(breaks=c(10000,100000,1000000), labels = comma) +
  coord_trans(y='log10') +
  ylab("Log10 Inversion size (bp)") +
  geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed") + 
  scale_fill_manual(values = c('gold3', 'lightcyan4'), name="") +
  xlab("") +
  theme_bw()
# Plot SVclass of valid vs invalid
plt2 <- orangutan.master.table %>%
  group_by(SVclass, valid.total) %>% 
  summarise(count=n()) %>% 
  mutate(ID = paste0(SVclass,"_",valid.total)) %>%
  ggplot(aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
  geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightblue3', 'lightblue4', 'gold3', 'gold4'), name="") +
  xlab("")
merge.plt <- plot_grid(plt1, plt2, nrow = 1)
title <- ggdraw() + draw_label("Orangutan validation", fontface='bold')
final.plt <- plot_grid(title, merge.plt, ncol=1, rel_heights=c(0.1, 1))
all.plots[[length(all.plots) + 1]] <- final.plt
# Export unvalidated simple inversions as a bed file
orangutan.master.table.unvalid <- orangutan.master.table[orangutan.master.table$valid.total == 'unvalid' & orangutan.master.table$SVclass == 'INV',]
destination <- file.path(export.unvalids, 'orangutan_remainUnvalid_simpleINV.bed')
write.table(orangutan.master.table.unvalid[,c(1:3)], file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Prepare final plot and save
final.plt <- plot_grid(plotlist = all.plots, ncol = 2)
destination <- file.path(export.unvalids, 'INVvalidation_plots.RData')
save(final.plt, file = destination)
destination <- file.path(export.unvalids, 'INVvalidation_plots.pdf')
ggsave(final.plt, filename = destination, width = 15, height = 10)

## Make summary stat for all calls ##
#####################################
allCalls.master.table <- rbind(chimpanzee.master.table[c('seqnames','start','end','width','gen','ID','SVclass','published','valid.total')], 
                                bonobo.master.table[c('seqnames','start','end','width','gen','ID','SVclass','published','valid.total')], 
                                gorilla.master.table[c('seqnames','start','end','width','gen','ID','SVclass','published','valid.total')], 
                                orangutan.master.table[c('seqnames','start','end','width','gen','ID','SVclass','published','valid.total')])

plt1 <- allCalls.master.table %>% 
  ggplot() + geom_boxplot(aes(x=valid.total, y=width, fill=valid.total)) +
  scale_y_continuous(breaks=c(10000,100000,1000000), labels = comma) +
  coord_trans(y='log10') +
  ylab("Log10 Inversion size (bp)") +
  geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed") + 
  scale_fill_manual(values = c('gold3', 'lightcyan4'), name="") +
  xlab("") +
  theme_bw()
# Plot SVclass of valid vs invalid
plt2 <- allCalls.master.table %>%
  group_by(SVclass, valid.total) %>% 
  summarise(count=n()) %>% 
  mutate(ID = paste0(SVclass,"_",valid.total)) %>%
  ggplot(aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
  geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightblue3', 'lightblue4', 'gold3', 'gold4'), name="") +
  xlab("")
summary.plt <- plot_grid(plt1, plt2, nrow = 1)
## Prepare final plot and save
destination <- file.path(export.unvalids, 'summary_INVvalidation.pdf')
ggsave(summary.plt, filename = destination, width = 8, height = 4)
## Export summary table of all validated inversion calls
destination <- file.path(export.unvalids, 'allCalls.master.table.validatio.csv')
write.table(allCalls.master.table, file = destination, quote = FALSE, row.names = FALSE, sep = ",")

## Extra calculations [OPTIONAL] ##
###################################
## Get percentage of HET inversions among unvalidated inversions
# unvalidSimpleINV <- allCalls.master.table[allCalls.master.table$valid.total == 'unvalid' & allCalls.master.table$SVclass == 'INV',]
# gens <- table(unvalidSimpleINV$gen)
# het.perc <- (gens[1] / sum(gens)) * 100
## Get proportions on SD flanked inversion among unvalidated calls
# unvalidSimpleINV.gr <- GRanges(seqnames=unvalidSimpleINV$seqnames, ranges=IRanges(start=unvalidSimpleINV$start, end=unvalidSimpleINV$end), gen=unvalidSimpleINV$gen)
# unvalidSimpleINV.dists <- range2rangeDistance(gr=unvalidSimpleINV.gr, userTrack=seg.dup.gr, allow.overlap = TRUE)
# SDflank <- which(unvalidSimpleINV.dists$leftDist < 5000 & unvalidSimpleINV.dists$rightDist < 5000 & unvalidSimpleINV.dists$leftDist >= 0 & unvalidSimpleINV.dists$rightDist >= 0)
# unvalidSimpleINV.gr$SDflank <- FALSE
# unvalidSimpleINV.gr[SDflank]$SDflank <- TRUE
# (length(unvalidSimpleINV.gr$SDflank[unvalidSimpleINV.gr$SDflank]) / length(unvalidSimpleINV.gr)) * 100
# unvalidSimpleINV.gr[unvalidSimpleINV.gr$gen == 'HET' | unvalidSimpleINV.gr$SDflank == TRUE]

message("DONE!!!")

## Run this part manually!!!
if (FALSE) {
  ## New inversions vs published data ##
  ######################################
  ## Load published data
  published.invCalls <- get(load("/home/porubsky/WORK/Great_apes/Published_data/published.invCalls.RData"))
  published.invCalls.df <- as.data.frame(published.invCalls)
  ## Summarize Strand-seq validated inversions
  strandseq.valid.simpleINV <- allCalls.master.table[allCalls.master.table$SVclass == 'INV' & allCalls.master.table$valid.total == 'valid',]
  
  published.summary <-
    published.invCalls.df %>%
    group_by(study, ID) %>%
    summarise(count=n())
  
  strandseq.summary <-
    strandseq.valid.simpleINV %>%
    filter(SVclass == 'INV') %>%  
    group_by(ID, SVclass) %>%
    summarise(count=n()) %>%
    mutate(study='Strand-seq') %>%
    select(study, ID, count)
  
  ## Export counts of validated Strand-seq invertsions
  destination <- file.path(export.unvalids, 'strandseq.valid.simpleINV.counts.txt')
  write.table(strandseq.summary, file = destination, quote = FALSE, row.names = FALSE)
  
  inversion.summary <- rbind(published.summary, strandseq.summary)
  inversion.summary <- inversion.summary[inversion.summary$ID != 'PericenINV',]
  
  ## Plot summary plot of validated inversion
  plt <- ggplot(inversion.summary, aes(x=ID, y=count, fill=study)) + 
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_text(aes(label=count), position = position_dodge2(width = 0.9, preserve = "single"), vjust=1) +
    scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
    xlab("") +
    ylab("Validated inversions")
  ## Save the plot
  destination <- file.path(export.unvalids, 'validStrandS_vs_publishedData.pdf')
  ggsave(plt, filename = destination, width = 7, height = 4)
  
  ## Get summary of likely novel simple inversions
  strandseq.summary.novel <-
    strandseq.valid.simpleINV %>%
    filter(SVclass == 'INV') %>%  
    group_by(ID, published) %>%
    summarise(count=n())
  
  strandseq.summary.valid <- allCalls.master.table %>%
    filter(SVclass == 'INV') %>% 
    group_by(SVclass, valid.total) %>% 
    summarise(count=n()) %>% 
    mutate(ID = paste0(SVclass,"_",valid.total))
  
  ## Plot summary of validated simple inversions
  plt1 <- ggplot(strandseq.summary.valid, aes(x=SVclass, y=count)) + geom_col(aes(fill=ID)) +
    geom_text(aes(group=ID, label=count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c('lightblue3', 'lightblue4'), name="") +
    xlab("") + 
    coord_flip() +
    theme_bw()
  
  ## Plot summary plot of likely novel inversions
  plt2 <- ggplot(strandseq.summary.novel, aes(x=ID, y=count, fill=published)) + 
    geom_col(position = position_stack()) +
    geom_text(aes(label=count), position = position_stack(), vjust=0.5, hjust=1) +
    scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("") +
    coord_flip() +
    theme_bw()
  
  final.plt <- plot_grid(plt1, plt2, ncol = 1, align = 'v', rel_heights = c(1,2))
  
  ## Save the plot
  destination <- file.path(export.unvalids, 'validInversions_publishedVSnovel.pdf')
  ggsave(final.plt, filename = destination, width = 8, height = 4)
  
  ## Plot summary table of validated inversion by Strand-seq and published data
  library(gridExtra)
  inversion.summary.tab <- reshape2::dcast(inversion.summary , study ~ ID)
  # Set theme to allow for plotmath expressions
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  tbl <- tableGrob(inversion.summary.tab, rows=NULL, theme=tt)
  # Plot chart and table into one object
  grid.arrange(tbl,
               nrow=1,
               as.table=TRUE)
  
  ## Plot proportion of novel inversions
  upset(fromList(valid.chimpanzee.list), order.by = 'freq', nsets = length(valid.chimpanzee.list))
}  
