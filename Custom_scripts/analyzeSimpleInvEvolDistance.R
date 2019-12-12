## Compare all pairs of individuals and estimate evolutionary distances  among Great Apes ##
############################################################################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(UpSetR) )
suppressPackageStartupMessages( library(VennDiagram) )
suppressPackageStartupMessages( library(ggtree) )
suppressPackageStartupMessages( library(ape) )
suppressPackageStartupMessages( library(dplyr) )

outputDirectory <- "/home/porubsky/WORK/Great_apes/Evolutionary_distance/"

message("Estimating evolutionary distances among Great Apes ...")

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load simple inversion without previously established human specific events and potential misorients
all.SimpleInversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.SimpleInversion.calls.filt.annot.RData"))
## Load simple inversion callset for human (NA19240)
hgsvc.NA19240.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/human.NA19240.filt.RData"))
hgsvc.simpleINV.NA19240.gr <- hgsvc.NA19240.gr[hgsvc.NA19240.gr$SVclass != 'invDup']
## Remove marked misorients from NHP and HGSVC callset
all.SimpleInversion.calls.filt <- all.SimpleInversion.calls.filt[all.SimpleInversion.calls.filt$misorient == FALSE]
hgsvc.simpleINV.NA19240.gr <- hgsvc.simpleINV.NA19240.gr[hgsvc.simpleINV.NA19240.gr$misorient == FALSE]
hgsvc.simpleINV.NA19240.gr <- hgsvc.simpleINV.NA19240.gr[, 1:3]

## Remove HS inversions from NHP set and add them to NA19240
HSinv.gr <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
overlaps2HS <- getReciprocalOverlaps(query = all.SimpleInversion.calls.filt, subject = HSinv.gr, report = 'query')
all.SimpleInversion.calls.noHS <- overlaps2HS[!overlaps2HS$perc.overlap >= 50, 1:3]
HSinv2add.gr <- HSinv.gr[,0]
HSinv2add.gr$gen <- 'HOM'
HSinv2add.gr$ID <- 'NA19240'
HSinv2add.gr$SVclass <- 'INV'
hgsvc.simpleINV.NA19240.plusHS <- sort(c(hgsvc.simpleINV.NA19240.gr, HSinv2add.gr))

## Merge calls and split by ID into GRangesList
all.SimpleInversion.allGreatApes <- sort(c(all.SimpleInversion.calls.noHS[,1:3], hgsvc.simpleINV.NA19240.plusHS))
all.SimpleInversion.allGreatApes.grl <- split(all.SimpleInversion.allGreatApes, all.SimpleInversion.allGreatApes$ID)

## Compare all possible pairs of great ape individuals ## [OPTIONAL]
#########################################################
message("    Comparing all pairs of individuals ...")
all.comparisons <- list()
all.comparisons.widths <- list()
## bonobo.vs.chimp ##
chimp.vs.bonobo <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['chimpanzee']], all.SimpleInversion.allGreatApes.grl[['bonobo']]), percTh = 50)
intersect <- split(chimp.vs.bonobo$sub.group, chimp.vs.bonobo$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['bonobo']]) - shared)
second <- abs(length(intersect[['chimpanzee']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['bonobo.vs.chimp']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(chimp.vs.bonobo), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[chimp.vs.bonobo$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['bonobo.vs.chimp']] <- inv.size.df 

## chimp.vs.gorilla ##
chimp.vs.gorilla <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['chimpanzee']], all.SimpleInversion.allGreatApes.grl[['gorilla']]), percTh = 50)
intersect <- split(chimp.vs.gorilla$sub.group, chimp.vs.gorilla$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['chimpanzee']]) - shared)
second <- abs(length(intersect[['gorilla']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['chimp.vs.gorilla']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(chimp.vs.gorilla), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[chimp.vs.gorilla$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['chimp.vs.gorilla']] <- inv.size.df 

## chimp.vs.orangutan ##
chimp.vs.orangutan <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['chimpanzee']], all.SimpleInversion.allGreatApes.grl[['orangutan']]), percTh = 50)
intersect <- split(chimp.vs.orangutan$sub.group, chimp.vs.orangutan$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['chimpanzee']]) - shared)
second <- abs(length(intersect[['orangutan']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['chimp.vs.orangutan']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(chimp.vs.orangutan), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[chimp.vs.orangutan$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['chimp.vs.orangutan']] <- inv.size.df

## bonobo.vs.gorilla ##
bonobo.vs.gorilla <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['bonobo']], all.SimpleInversion.allGreatApes.grl[['gorilla']]), percTh = 50)
intersect <- split(bonobo.vs.gorilla$sub.group, bonobo.vs.gorilla$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['bonobo']]) - shared)
second <- abs(length(intersect[['gorilla']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['bonobo.vs.gorilla']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(bonobo.vs.gorilla), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[bonobo.vs.gorilla$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['bonobo.vs.gorilla']] <- inv.size.df

## bonobo.vs.orangutan ##
bonobo.vs.orangutan <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['bonobo']], all.SimpleInversion.allGreatApes.grl[['orangutan']]), percTh = 50)
intersect <- split(bonobo.vs.orangutan$sub.group, bonobo.vs.orangutan$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['bonobo']]) - shared)
second <- abs(length(intersect[['orangutan']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['bonobo.vs.orangutan']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(bonobo.vs.orangutan), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[bonobo.vs.orangutan$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['bonobo.vs.orangutan']] <- inv.size.df

## gorilla.vs.orangutan ##
gorilla.vs.orangutan <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['gorilla']], all.SimpleInversion.allGreatApes.grl[['orangutan']]), percTh = 50)
intersect <- split(gorilla.vs.orangutan$sub.group, gorilla.vs.orangutan$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['gorilla']]) - shared)
second <- abs(length(intersect[['orangutan']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['gorilla.vs.orangutan']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(gorilla.vs.orangutan), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[gorilla.vs.orangutan$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['gorilla.vs.orangutan']] <- inv.size.df

## gorilla.vs.na19240 ##
gorilla.vs.na19240 <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['gorilla']], all.SimpleInversion.allGreatApes.grl[['NA19240']]), percTh = 50)
intersect <- split(gorilla.vs.na19240$sub.group, gorilla.vs.na19240$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['gorilla']]) - shared)
second <- abs(length(intersect[['NA19240']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['gorilla.vs.na19240']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(gorilla.vs.na19240), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[gorilla.vs.na19240$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['gorilla.vs.na19240']] <- inv.size.df

## orangutan.vs.na19240 ##
orangutan.vs.na19240 <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['orangutan']], all.SimpleInversion.allGreatApes.grl[['NA19240']]), percTh = 50)
intersect <- split(orangutan.vs.na19240$sub.group, orangutan.vs.na19240$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['orangutan']]) - shared)
second <- abs(length(intersect[['NA19240']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['orangutan.vs.na19240']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(orangutan.vs.na19240), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[orangutan.vs.na19240$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['orangutan.vs.na19240']] <- inv.size.df

## bonobo.vs.na19240 ##
bonobo.vs.na19240 <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['bonobo']], all.SimpleInversion.allGreatApes.grl[['NA19240']]), percTh = 50)
intersect <- split(bonobo.vs.na19240$sub.group, bonobo.vs.na19240$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['bonobo']]) - shared)
second <- abs(length(intersect[['NA19240']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['bonobo.vs.na19240']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(bonobo.vs.na19240), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[bonobo.vs.na19240$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['bonobo.vs.na19240']] <- inv.size.df

## chimpanzee.vs.na19240 ##
chimpanzee.vs.na19240 <- getDisjointOverlapsWeighted(c(all.SimpleInversion.allGreatApes.grl[['chimpanzee']], all.SimpleInversion.allGreatApes.grl[['NA19240']]), percTh = 50)
intersect <- split(chimpanzee.vs.na19240$sub.group, chimpanzee.vs.na19240$ID)
shared <- length(intersect(intersect[[1]],intersect[[2]]))
first <- abs(length(intersect[['chimpanzee']]) - shared)
second <- abs(length(intersect[['NA19240']]) - shared)
overlap <- data.frame(first=first, second=second, shared=shared)
all.comparisons[['chimpanzee.vs.na19240']] <- overlap

overlap.object <- VennDiagram::get.venn.partitions(intersect)
shared.idx <- overlap.object$..values..[[1]]
inv.size.df <- data.frame(size=width(chimpanzee.vs.na19240), ID='unique', stringsAsFactors = FALSE)
inv.size.df$ID[chimpanzee.vs.na19240$sub.group %in% shared.idx] <- 'shared'
all.comparisons.widths[['chimpanzee.vs.na19240']] <- inv.size.df

## Compile all comparisons into a table
all.comparisons.m <- do.call(rbind, all.comparisons)
all.comparisons.df <- as.data.frame(all.comparisons.m)
colnames(all.comparisons.df) <- c('first','second','shared')
# Plot the overlaps
all.comparisons.df.perc <- round((all.comparisons.df / rowSums(all.comparisons.df)) * 100, digits = 1)
all.comparisons.df.perc <- all.comparisons.df.perc[order(all.comparisons.df.perc$shared),]
all.comparisons.df.perc$ID <- factor(rownames(all.comparisons.df.perc), levels = rownames(all.comparisons.df.perc))
plt.df <- reshape2::melt(all.comparisons.df.perc, measure.vars = c('first','second','shared'))
plt.df$variable <- factor(plt.df$variable, levels = c('second','shared','first'))
plt1 <- ggplot(plt.df, aes(x=ID, y=value, group=variable)) + 
  geom_col(aes(fill=variable)) +
  geom_text(aes(group=variable, label=value), position = position_stack(vjust = 0.5)) +
  ylab("% overlaping simple inversions") +
  xlab("Comparison") +
  scale_fill_manual(values = c('lightblue1','lightblue4', 'lightblue3'), name='') +
  coord_flip() + 
  theme(legend.position="top")

## Compile size distribution into a table
# all.comparisons.widths.m <- do.call(rbind, all.comparisons.widths)
# all.comparisons.widths.df <- as.data.frame(all.comparisons.widths.m)
# all.comparisons.widths.df$Comparison <- gsub(rownames(all.comparisons.widths.df), pattern = "\\.\\d+", replacement = "")
# plt.df <- all.comparisons.widths.df %>% group_by(Comparison, ID) %>% summarise(mean=mean(size), median=median(size))
# plt.df$x <- recode(plt.df$ID, 'shared'=1, 'unique'=2) 
# plt.df$y <- rep(plt.df$median[plt.df$ID == 'shared'], each=2)
# plt.df$yend <- rep(plt.df$median[plt.df$ID == 'unique'], each=2)
# plt2 <- ggplot() +
#   geom_segment(data=plt.df, aes(x = 1, y = y, xend = 2, yend = yend, group=Comparison), color='gray') +
#   geom_point(data=plt.df, aes(x=x, y=median, color=ID, size=mean), shape=15, inherit.aes = FALSE) +
#   facet_grid(Comparison ~ ., scales = 'free') +
#   scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2), labels = c('shared', 'unique')) +
#   scale_color_manual(values = c('lightblue4', 'lightblue1')) +
#   xlab("") +
#   ylab("Median inversion size") +
#   theme(strip.text.y = element_text(angle = 1))
# 
# final.plt <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(2,1))

#p <- venn.diagram(intersect, fill=1:length(intersect), alpha = 0.25, filename = NULL, main = "Chimpanzee vs Bonobo", print.mode="percent")
#grid.newpage()
#grid.draw(p)

## Plot contingency table of shared inversions ##
#################################################
## Match all possible pairs of great apes
all.pairs <- as.matrix( expand.grid(x=unique(all.SimpleInversion.allGreatApes$ID), y=unique(all.SimpleInversion.allGreatApes$ID)) )
## Prepare data matrix
IDs <- c('NA19240','chimpanzee','bonobo','gorilla','orangutan')
comparison.m <- matrix(data = rep(0, nrow(all.pairs)), nrow = length(IDs), ncol = length(IDs))
rownames(comparison.m) <- IDs
colnames(comparison.m) <- IDs

for (i in 1:nrow(all.pairs)) {
  pair <- all.pairs[i,]
  if (pair[1] == pair[2]) {
    shared <- 0
  } else {
    message("    Comparing ", pair[1], " vs ",pair[2])
    gr1 <- all.SimpleInversion.allGreatApes[all.SimpleInversion.allGreatApes$ID == pair[1]]
    gr2 <- all.SimpleInversion.allGreatApes[all.SimpleInversion.allGreatApes$ID == pair[2]]

    overlaps <- getDisjointOverlapsWeighted(c(gr1, gr2), percTh = 50)
    
    intersect <- split(overlaps$sub.group, overlaps$ID)
    shared <- length(intersect(intersect[[1]],intersect[[2]]))
  }
  comparison.m[pair[1], pair[2]] <- shared
}

## Prepare plot
plt.df <- reshape2::melt(comparison.m)
plt3 <- ggplot(plt.df) + 
  geom_tile(aes(x=Var1, y=Var2, fill=value)) +
  geom_text(aes(x=Var1, y=Var2, label=value), color="white") +
  scale_fill_gradient(name = "Shared INVs", trans = "log", low = "gray34", high = "orange2") +
  xlab("") +
  ylab("")

suppressWarnings( final.plt <- plot_grid(plt1, plt3, nrow = 1, align = 'h', axis = 't') )
## Export final plot
destination <- file.path(outputDirectory, "pairwiseOverlaps_greatApes.pdf")
ggsave(filename = destination, plot = final.plt, width = 15, height = 5, useDingbats=FALSE)

## Compare all great apes and construct phylogeny ##
####################################################
## For upsetR plot keep HSinversions and misorients in their original callset within NHP and NA19240!!!
## Load simple inversion without previously established human specific events and potential misorients
#all.SimpleInversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.SimpleInversion.calls.filt.annot.RData"))
## Load simple inversion callset for human (NA19240)
#hgsvc.NA19240.gr <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/human.NA19240.filt.RData"))
#hgsvc.simpleINV.NA19240.gr <- hgsvc.NA19240.gr[hgsvc.NA19240.gr$SVclass != 'invDup']
## Calculate overlaps for a set of genomic ranges
greatApes.inv.overlaps <- getDisjointOverlapsWeighted(c(all.SimpleInversion.calls.filt, hgsvc.simpleINV.NA19240.gr), percTh = 50)
## Split unique overlaps by individual
overlaps <- split(greatApes.inv.overlaps$sub.group, greatApes.inv.overlaps$ID)
## Prepare UpsetR plot
overlaps.df <- fromList(overlaps)
overlaps.df$n <- sample(1:nrow(overlaps.df))

destination <- file.path(outputDirectory, "upsetR_greatApesOverlaps.pdf")
pdf(destination, width = 8, height = 5, useDingbats=FALSE)
upset(overlaps.df,
      order.by = 'freq', 
      nsets = length(overlaps), 
      sets.bar.color = c("#e6550d","#8856a7","#31a354","#3182bd","#a6611a"),
      sets = c("bonobo", "chimpanzee", "gorilla", "NA19240", "orangutan"),
      queries = list(
        list(query = intersects, params = list("NA19240"), color='#a6611a', active = T),
        list(query = intersects, params = list("bonobo"), color="#3182bd", active = T),
        list(query = intersects, params = list("chimpanzee"), color="#31a354", active = T),
        list(query = intersects, params = list("gorilla"), color="#8856a7", active = T),
        list(query = intersects, params = list("orangutan"), color="#e6550d", active = T)
      )
)
dev.off()

## Get likely recurrent inversions ##
#####################################
vennPartitions <- get.venn.partitions(overlaps)
vennPartitions$..set..
#c(13,4,15,6,5,2,8,7,14)
select.idx <- unlist(vennPartitions$..values..[c(13,4,15,6,5,2,8,7,14)])
recurr.inv <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% select.idx]
## Calculate distance to known SDs
recurr.inv.SDdist <- range2rangeDistance(gr=reduce(recurr.inv), userTrack=seg.dup.gr, allow.overlap = TRUE)
## Export likely recurrent loci
destination <- file.path(outputDirectory, "putative_recurrent_sites.RData")
save(recurr.inv.SDdist, file = destination)

## Construct matrix of shared inverted regions and report genotype for each shared inversion ##
###############################################################################################
## Here same filtering applied as for upsetR plot (see above) !!!
## Calculate overlaps for a set of genomic ranges
#ranges2compare <- c(all.SimpleInversion.calls.filt, hgsvc.simpleINV.NA19240.gr)
## Here human specific inversions are assigned to human and removed from NHPs
ranges2compare <- all.SimpleInversion.allGreatApes
## Keep only autosomal chromosomes
ranges2compare <- ranges2compare[seqnames(ranges2compare) != 'chrX']
greatApes.inv.overlaps <- getDisjointOverlapsWeighted(gr = ranges2compare, percTh = 50)
## Collapse overlapping regions into a nonredundant dataset
NHP.nonred.gr <- collapseBins(greatApes.inv.overlaps[order(greatApes.inv.overlaps$sub.group)], id.field = ncol(mcols(greatApes.inv.overlaps)))
NHP.nonred.gr <- NHP.nonred.gr[,0]
NHP.nonred.gr <- getSegDupOverlaps(query.gr = NHP.nonred.gr, subject.gr = seg.dup.gr)
## Add lineage specific information
greatApes.inv.overlaps.grl <- split(greatApes.inv.overlaps, greatApes.inv.overlaps$sub.group)
lineages <- lapply(greatApes.inv.overlaps.grl, function(x) paste(x$ID, collapse = "-"))
NHP.nonred.gr$lineages <- lineages
## Export results
destination <- file.path(outputDirectory, "shared.invRegions.autosomes.bed")
NHP.nonred.df <- as.data.frame(NHP.nonred.gr)
write.table(NHP.nonred.df, file = destination, quote = FALSE, row.names = FALSE, sep = "\t")

## Split unique overlaps by individual
overlaps <- split(greatApes.inv.overlaps$sub.group, greatApes.inv.overlaps$ID)
## Prepare UpsetR plot
overlaps.df <- fromList(overlaps)
comp <- overlaps.df
## Split ranges by individual
overlaps.grl <- split(greatApes.inv.overlaps, greatApes.inv.overlaps$ID)

comp$all.inv.ids <- unique(unlist(overlaps, use.names = FALSE))
comp$chr <- as.character(seqnames(greatApes.inv.overlaps[match(comp$all.inv.ids, greatApes.inv.overlaps$sub.group)]))
comp$bonobo[comp$all.inv.ids %in% overlaps[['bonobo']]] <- as.character(overlaps.grl[['bonobo']]$gen[order(match(overlaps[['bonobo']], comp$all.inv.ids))])
comp$chimpanzee[comp$all.inv.ids %in% overlaps[['chimpanzee']]] <- as.character(overlaps.grl[['chimpanzee']]$gen[order(match(overlaps[['chimpanzee']], comp$all.inv.ids))])
comp$gorilla[comp$all.inv.ids %in% overlaps[['gorilla']]] <- as.character(overlaps.grl[['gorilla']]$gen[order(match(overlaps[['gorilla']], comp$all.inv.ids))])
comp$NA19240[comp$all.inv.ids %in% overlaps[['NA19240']]] <- as.character(overlaps.grl[['NA19240']]$gen[order(match(overlaps[['NA19240']], comp$all.inv.ids))])
comp$orangutan[comp$all.inv.ids %in% overlaps[['orangutan']]] <- as.character(overlaps.grl[['orangutan']]$gen[order(match(overlaps[['orangutan']], comp$all.inv.ids))])
comp.m <- apply(comp, 2, function(x) dplyr::recode(x, 'HET' = '1', 'HOM' = '2'))
## Export results
destination <- file.path(outputDirectory, "shared.inv.perGen.autosomes.tsv")
write.table(comp.m, file = destination, quote = FALSE, row.names = FALSE, sep = "\t")

## Prepare evolutionary tree based on shared simple inversions ##
#################################################################
greatApes.inv.overlaps <- getDisjointOverlapsWeighted(all.SimpleInversion.allGreatApes, percTh = 50)
## Split unique overlaps by individual
overlaps <- split(greatApes.inv.overlaps$sub.group, greatApes.inv.overlaps$ID)

overlaps.df <- fromList(overlaps)
overlaps.m <- t(overlaps.df)
## Get counts of shared events
vennPartitions <- get.venn.partitions(overlaps)
overlap <- data.frame(set = vennPartitions$..set.., count = vennPartitions$..count..)
node.counts <- overlap$count[c(1,17,21,29)]
tip.counts <- overlap$count[c(31,30,28,24,16)]
## Contruct the tree
hc.tree <- hclust(dist(overlaps.m))
hc.tree <- as.phylo(hc.tree)
## Prepare extra annotation data.frames
annot.nodes <- data.frame(node=unique(hc.tree$edge[,1]), text=node.counts, stringsAsFactors = FALSE)
annot.tips <- data.frame(node=1:length(tip.counts), good=tip.counts, stringsAsFactors = FALSE)
## Add genotype information on unique inversions
chimp.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[30])]
chimp.only.gen <- table(chimp.only$gen)
bonobo.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[31])]
bonobo.only.gen <- table(bonobo.only$gen)
gorilla.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[28])]
gorilla.only.gen <- table(gorilla.only$gen)
orangutan.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[16])]
orangutan.only.gen <- table(orangutan.only$gen)
human.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[24])]
human.only.gen <- table(human.only$gen)

#annot.tip.genoT.hom <- data.frame(node=rep(1:length(tip.counts)),
#                                  HOM.t=c(bonobo.only.gen['HOM'], chimp.only.gen['HOM'], gorilla.only.gen['HOM'], human.only.gen['HOM'], orangutan.only.gen['HOM']))
#annot.tip.genoT.het <- data.frame(node=rep(1:length(tip.counts)),
#                                  HET.t=c(bonobo.only.gen['HET'], chimp.only.gen['HET'], gorilla.only.gen['HET'], human.only.gen['HET'], orangutan.only.gen['HET']))

## Add genotype information on shared inversions
node.plots <- list()
CBHGO.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[1])]
#CBHGO.only.gen <- table(CBHGO.only$gen)
CBHGO.only.gen <- as.data.frame(CBHGO.only) %>% group_by(ID, gen) %>% summarise(count=n())
CBHGO.only.gen$ID <- dplyr::recode(CBHGO.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
node.plots[['6']] <- ggplot(CBHGO.only.gen) + geom_col(aes(y=count, x=ID, fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none')
CBHG.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[17])]
#CBHG.only.gen <- table(CBHG.only$gen)
CBHG.only.gen <- as.data.frame(CBHG.only) %>% group_by(ID, gen) %>% summarise(count=n())
CBHG.only.gen$ID <-  dplyr::recode(CBHG.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
node.plots[['7']] <- ggplot(CBHG.only.gen) + geom_col(aes(y=count, x=ID, fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none')
CBH.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[21])]
#CBH.only.gen <- table(CBH.only$gen)
CBH.only.gen <- as.data.frame(CBH.only) %>% group_by(ID, gen) %>% summarise(count=n())
CBH.only.gen$ID <-  dplyr::recode(CBH.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
node.plots[['8']] <- ggplot(CBH.only.gen) + geom_col(aes(y=count, x=ID, fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none')
CB.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[29])]
#CB.only.gen <- table(CB.only$gen)
CB.only.gen <- as.data.frame(CB.only) %>% group_by(ID, gen) %>% summarise(count=n())
CB.only.gen$ID <-  dplyr::recode(CB.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
node.plots[['9']] <- ggplot(CB.only.gen) + geom_col(aes(y=count, x=ID, fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none')

#annot.node.genoT.hom <- data.frame(node=unique(hc.tree$edge[,1]),
#                                  HOM.n=c(CBHGO.only.gen['HOM'], CBHG.only.gen['HOM'], CBH.only.gen['HOM'], CB.only.gen['HOM']))
#annot.node.genoT.het <- data.frame(node=unique(hc.tree$edge[,1]),
#                                 HET.n=c(CBHGO.only.gen['HET'], CBHG.only.gen['HET'], CBH.only.gen['HET'], CB.only.gen['HET']))

## Add genotype information on unique inversions
tip.plots <- list()
C.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[30])]
C.only.gen <- as.data.frame(C.only) %>% group_by(ID, gen) %>% summarise(count=n())
C.only.gen$ID <-  dplyr::recode(C.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
tip.plots[['2']] <- ggplot(C.only.gen, aes(y=count, x=ID, group=gen)) + geom_col(aes(fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none') + geom_text(aes(group=gen, label = count), color='white', position = position_stack(vjust = 0.5))
B.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[31])]
B.only.gen <- as.data.frame(B.only) %>% group_by(ID, gen) %>% summarise(count=n())
B.only.gen$ID <-  dplyr::recode(B.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
tip.plots[['1']] <- ggplot(B.only.gen, aes(y=count, x=ID, group=gen)) + geom_col(aes(fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none') + geom_text(aes(group=gen, label = count), color='white', position = position_stack(vjust = 0.5))
G.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[28])]
G.only.gen <- as.data.frame(G.only) %>% group_by(ID, gen) %>% summarise(count=n())
G.only.gen$ID <-  dplyr::recode(G.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
tip.plots[['3']] <- ggplot(G.only.gen, aes(y=count, x=ID, group=gen)) + geom_col(aes(fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none') + geom_text(aes(group=gen, label = count), color='white', position = position_stack(vjust = 0.5))
O.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[16])]
O.only.gen <- as.data.frame(O.only) %>% group_by(ID, gen) %>% summarise(count=n())
O.only.gen$ID <-  dplyr::recode(O.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
tip.plots[['5']] <- ggplot(O.only.gen, aes(y=count, x=ID, group=gen)) + geom_col(aes(fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none') + geom_text(aes(group=gen, label = count), color='white', position = position_stack(vjust = 0.5))
H.only <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[24])]
H.only.gen <- as.data.frame(H.only) %>% group_by(ID, gen) %>% summarise(count=n())
H.only.gen$ID <-  dplyr::recode(H.only.gen$ID, 'bonobo' = 'B', 'chimpanzee' = 'C', 'gorilla' = 'G', 'NA19240' = 'H', 'orangutan' = 'O')
tip.plots[['4']] <- ggplot(H.only.gen, aes(y=count, x=ID, group=gen)) + geom_col(aes(fill=gen)) + xlab("") + ylab("") + scale_fill_manual(values=c("red", "blue"), guide='none') + geom_text(aes(group=gen, label = count), color='white', position = position_stack(vjust = 0.5))

## Add extra annotation to the tree
plt <- ggtree(hc.tree)
plt <- plt %<+% annot.nodes + geom_nodelab(aes(label=text), geom = "label")
#plt <- plt %<+% annot.node.genoT.hom + geom_nodelab(aes(label=HOM.n), color="red", geom = "label", hjust=-0.5, vjust=0)
#plt <- plt %<+% annot.node.genoT.het + geom_nodelab(aes(label=HET.n), color="blue", geom = "label", hjust=-0.5, vjust=1)
plt <- plt %<+% annot.tips + geom_tiplab(aes(label=good, subset=isTip), geom = "label", hjust = 2) 
#plt <- plt %<+% annot.tip.genoT.hom + geom_tiplab(aes(label=HOM.t, subset=isTip), color="red", geom = "label", hjust = 1, vjust=0)
#plt <- plt %<+% annot.tip.genoT.het + geom_tiplab(aes(label=HET.t, subset=isTip), color="blue",  geom = "label", hjust = 1, vjust=1)
inset.width <- c(1.2, 1.1, 1.1, 1)
plt <- inset(plt, node.plots, width = inset.width, height = 1, hjust = -0.75, vjust = 0.1)
plt <- inset(plt, tip.plots, width = 1, height = 1, hjust = -0.1, vjust = 0.1)
plt <- plt + xlim_tree(xlim = c(1:max(hc.tree$edge.length)+2)) 
plt <- plt + geom_tiplab(offset = 0.5)

## Export final plot
destination <- file.path(outputDirectory, "simpleINV_evolTree.pdf")
ggsave(filename = destination, plot = plt, width = 12, height = 7, useDingbats=FALSE)

## Make unrooted NJ tree
#plot(nj(dist(overlaps.m)))

## Check inversion unique in NA19240
#human.specific.regions.trueSet <- get(load("/home/porubsky/WORK/Great_apes/Human_specific_events/human.specific.regions.trueSet.RData"))
#only.human <- greatApes.inv.overlaps[greatApes.inv.overlaps$sub.group %in% unlist(vennPartitions$..values..[24])]
#only.human <- subsetByOverlaps(only.human, human.specific.regions.trueSet, invert = T)

message("DONE!!!")

