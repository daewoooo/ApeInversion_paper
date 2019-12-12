## This script analyze regions shared between human and other great apes ##
###########################################################################

## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(UpSetR) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(reshape2) )

outputDirectory <- "/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites"

message("Searching for polymorphic inversions ...")

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load Strand-seq simple INV calls ##
######################################
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']

## Create a non-redundant set from NHP simple inversions
overlaps.gr <- getDisjointOverlapsWeighted(all.SimpleInversion.calls.filt, percTh = 50)
NHP.nonred.gr <- collapseBins(overlaps.gr[order(overlaps.gr$sub.group)], id.field = 11)
NHP.nonred.gr <- sort(NHP.nonred.gr)
NHP.nonred.gr$ID <- 'NHP'

## Load ape's composite files ##
################################
chimp.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
chimp.data$ID <- 'chimpanzee'
bonobo.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
bonobo.data$ID <- 'bonobo'
gorilla.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
gorilla.data$ID <- 'gorilla'
orangutan.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
orangutan.data$ID <- 'orangutan'

## Polymorphic inversions between great apes and humans ##
##########################################################
# ## Load HGSVC collapsed inverions
# hgsvc.inv.regions.reduced <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/InversionRegions_reduced_v2_HGSVC.bed", header = TRUE, stringsAsFactors = FALSE)
# hgsvc.inv.regions.reduced.gr <- GRanges(seqnames = hgsvc.inv.regions.reduced$seqnames, ranges=IRanges(start=hgsvc.inv.regions.reduced$start, end=hgsvc.inv.regions.reduced$end))
# hgsvc.inv.regions.reduced.gr$gen <- ''
# hgsvc.inv.regions.reduced.gr$ID <- 'HGSVC'
# hgsvc.inv.regions.reduced.gr$SVclass <- 'INV'
# hgsvc.inv.regions.reduced.gr <- getSegDupOverlaps(query.gr = hgsvc.inv.regions.reduced.gr, subject.gr = seg.dup.gr)
# ## Load Ashley's GR paper inversions
# ashley.invertom <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/Ashley_GR/Supplemental_Tables_S4_polymorphic_inversions.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
# ashley.invertom.gr <- GRanges(seqnames = ashley.invertom$Chr, ranges=IRanges(start=ashley.invertom$Start, end=ashley.invertom$End))
# ashley.invertom.gr$gen <- ''
# ashley.invertom.gr$ID <- 'AshleyGR'
# ashley.invertom.gr$SVclass <- 'INV'
# ashley.invertom.gr <- getSegDupOverlaps(query.gr = ashley.invertom.gr, subject.gr = seg.dup.gr)
# ## Merge Inversion callsets from HGSVC and Ashley's GR paper
# all.human.inv <- suppressWarnings( c(hgsvc.inv.regions.reduced.gr, ashley.invertom.gr) )
# all.human.inv.reduced <- reduce(all.human.inv) # Take only nonverlapping ranges!!!
# all.human.inv.reduced$gen <- ''
# all.human.inv.reduced$ID <- 'human'
# all.human.inv.reduced$SVclass <- 'INV'
# all.human.inv.reduced$SDTrackPercOverlap <- 0
# all.human.inv.reduced$TotalUniqueBases <- 0
# all.human.inv.reduced$LongestUniqueBases <- 0

## Load HGSCV simple inversions
hgsvc.simple.inv <- read.table("/home/porubsky/WORK/Great_apes/Ashley_inversions/HGSVC_data/HGSVC_simple_inversions.csv", sep = ",", header=TRUE, stringsAsFactors = FALSE)
## Construct GRanges object
hgsvc.simple.inv.gr <- GRanges(seqnames=hgsvc.simple.inv$seqnames, 
                                       ranges=IRanges(start = hgsvc.simple.inv$innerBP_start, end = hgsvc.simple.inv$innerBP_end)
                               )
## Add genotype information
mcols(hgsvc.simple.inv.gr) <- hgsvc.simple.inv[,c(10:ncol(hgsvc.simple.inv ))]
## Filter inversion with at least 1000kb of unique sequence
hgsvc.simple.inv.gr <- getSegDupOverlaps(query.gr = hgsvc.simple.inv.gr, subject.gr = seg.dup.gr)
hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[hgsvc.simple.inv.gr$TotalUniqueBases >= 1000]

## Export HGSVC and NHP nonredundant simple inversion datasets
destination <- "/home/porubsky/WORK/Great_apes/Final_INV_calls/HGSVC.nonred.filt.RData"
save(hgsvc.simple.inv.gr, file = destination)
destination <- "/home/porubsky/WORK/Great_apes/Final_INV_calls/NHP.nonred.filt.RData"
save(NHP.nonred.gr, file = destination)

## Make Venn diagram of overlaps between HGSVC and NHP non-redundant datasets ##
################################################################################
## Remove human SD regions flanking HGSVC and NHP inversion calls
hgsvc.calls <- subtractRegions(gr = hgsvc.simple.inv.gr[,0], remove.gr = reduce(seg.dup.gr), mode = 'flanks')
hgsvc.calls$ID <- 'HGSVC'
nhp.calls <- subtractRegions(gr = NHP.nonred.gr[,0], remove.gr = reduce(seg.dup.gr), mode = 'flanks')
nhp.calls$ID <- 'NHP'
hgsvc.nhp.calls <- sort(c(hgsvc.calls, nhp.calls))
#hgsvc.nhp.overlaps <- getReciprocalOverlaps(query = hgsvc.calls, subject = nhp.calls, report = 'both')
hgsvc.nhp.overlaps <- getDisjointOverlapsWeighted(gr = hgsvc.nhp.calls, percTh = 50)
## Split overlaps by individual
overlaps.list <- split(hgsvc.nhp.overlaps$sub.group, hgsvc.nhp.overlaps$ID)
## To get correct counts of overlapping ranges within a group ID assign them unique idx
unique.NHP <- !overlaps.list[['NHP']] %in% overlaps.list[['HGSVC']]
overlaps.list[['NHP']][unique.NHP] <- seq(10000, by=1, length=length(unique.NHP[unique.NHP == TRUE]))
## Prepare venn diagram
library(Vennerable)
venn.data <- Venn(overlaps.list)
destination <- file.path(outputDirectory, "venn_HGSV2NHP.pdf")
pdf(destination, useDingbats=FALSE)
plot(venn.data)
dev.off()
## Get HGSVC inversions overlapping with NHP
library(VennDiagram)
vennPartitions <- get.venn.partitions(overlaps.list)
shared.idx <- unlist(vennPartitions$..values..[1])
shared.inv <- hgsvc.nhp.overlaps[hgsvc.nhp.overlaps$ID == 'HGSVC' & hgsvc.nhp.overlaps$sub.group %in% shared.idx]
## Evaluate reciprocal overlaps between HGSVC and NHP
#hgsvc.simple.inv.gr.percOverlap <- getReciprocalOverlaps(query = hgsvc.calls, subject = nhp.calls, report = 'query')
#percOverlap.df <- as.data.frame(hgsvc.simple.inv.gr.percOverlap[,'perc.overlap'])
#percOverlap.df$roi.num <- 1:nrow(percOverlap.df)

## Regenotype human inversions shared with non-human primates (50% reciprocal overlap)
## Genotype regions from NA19240 ##
###################################
## Set ranges to be re-genotyped
regions.to.genotype <- shared.inv[,0]
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

## Extract genotypes for all apes and merge with human's genotypes
hgsvc.simple.inv.gr <- subsetByOverlaps(hgsvc.simple.inv.gr, shared.inv) ## Select only shared inversions
hgsvc.simple.inv.gr <- hgsvc.simple.inv.gr[,c(1:9)]
hgsvc.genotypes <- mcols(hgsvc.simple.inv.gr)
nhp.genotypes <- mcols(genotyped.regions)[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]

## Get number of sites for which there is at least one HET
hets <- apply(nhp.genotypes, 1, function(x) any(x == 'HET'))
hets <- apply(hgsvc.genotypes, 1, function(x) any(x == 'HET'))

## Recode genotypes to number in order to calculate allele frequencies
hgsvc.alleles.inv <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
hgsvc.alleles.dir <- apply(hgsvc.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
nhp.alleles.inv <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 2, 'HET' = 1, 'REF' = 0, 'lowReads' = 0, 'AMB' = 0))
nhp.alleles.dir <- apply(nhp.genotypes, 2, function(x) dplyr::recode(x, 'HOM' = 0, 'HET' = 1, 'REF' = 2, 'lowReads' = 2, 'AMB' = 0))
## Calculate allele freaquencies for HGSVC and NHP genotypes
hgsvc.total.inv.alleles <- rowSums(hgsvc.alleles.inv)
hgsvc.total.dir.alleles <- rowSums(hgsvc.alleles.dir)
nhp.total.inv.alleles <- rowSums(nhp.alleles.inv)
nhp.total.dir.alleles <- rowSums(nhp.alleles.dir)
hgsvc.inv.allele.freq <- hgsvc.total.inv.alleles / (hgsvc.total.inv.alleles + hgsvc.total.dir.alleles)
nhp.inv.allele.freq <- nhp.total.inv.alleles / (nhp.total.inv.alleles + nhp.total.dir.alleles)

## Plot percentage overlap per site
plt.df <- as.data.frame(shared.inv)
plt.df$region.ID <- factor(as.character(shared.inv[,0]), levels = rev(as.character(shared.inv[,0])))

plt1 <- ggplot() + geom_col(data=plt.df, aes(y=perc.overlap, x=region.ID)) +
  coord_flip() +
  geom_hline(yintercept = c(25, 50, 75), linetype='dashed', color='white') +
  xlab('Genomic position') +
  ylab('% overlap')

## Plot shared inverted loci frequencies
plt.df <- data.frame(pos = as.character(hgsvc.simple.inv.gr), hgsvc.inv.allele.freq = hgsvc.inv.allele.freq, nhp.inv.allele.freq = nhp.inv.allele.freq * -1, stringsAsFactors = FALSE)
plt.df <- reshape2::melt(plt.df)
plt.df$pos <- factor(plt.df$pos, levels = rev(unique(plt.df$pos)))
plt2 <- ggplot(plt.df, aes(x=pos, y=value, fill=variable)) + 
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c('sienna1', 'slateblue3')) +
  ylab("Inverted loci frequency") +
  xlab("Genomic position")

## Plot number of shared INVs per chromosome
plt.df <- as.data.frame(shared.inv)
plt.df$seqnames <- factor(plt.df$seqnames, levels = rev(levels(plt.df$seqnames)))
#plt.df$seqnames <- factor(plt.df$seqnames, levels = levels(plt.df$seqnames))
plt3 <- plt.df %>% group_by(seqnames) %>% summarise(count=n()) %>%
  ggplot(aes(x=seqnames, y=count)) + geom_col() +
  geom_text(aes(label=count), color='white', vjust=0.5, hjust=1) +
  xlab("") +
  coord_flip()

## Set colors
n.colors <- length(unique(plt.df$seqnames))
colors <- rep(c('dimgray', 'gray82'), n.colors/2)
colors <- colors[1:n.colors]
plt3 <- plt.df %>% group_by(seqnames) %>% summarise(count=n()) %>%
  ggplot(aes(x=1, y=count, fill=seqnames)) + geom_bar(stat='identity') +
  geom_text(aes(label = count), size = 3, position = position_stack(vjust = 0.5)) +
  geom_text(aes(x = 1.5, label = seqnames), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors, guide='none') +
  xlab("") +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  
# plt.df <- plt.df %>% group_by(seqnames) %>% summarise(count=n()) %>% 
#   arrange(desc(seqnames)) %>%
#   mutate(lab.ypos = cumsum(count) - 0.5*count)
# 
# mycols <- brewer.pal(name = 'Dark2', n = 8)
# plt4 <- ggplot(plt.df, aes(x = "", y = count, fill = seqnames)) +
#   geom_bar(width = 1, stat = "identity", color = "white") +
#   coord_polar("y", start = 0)+
#   geom_text(aes(y = lab.ypos, label = count), color = "white") +
#   #scale_fill_manual(values = mycols) +
#   #ggtitle("Gene categories duplicated by\nInverted duplications") +
#   theme_void()

## Plot distance of inversion breakpoints to human SDs
shared.inv.SDdist <- range2rangeDistance(gr=shared.inv, userTrack=seg.dup.gr, allow.overlap = TRUE)
## To see count of polymorhpic inversion flanked by SDs within 5kbp range
#shared.inv.SDdist[shared.inv.SDdist$leftDist < 5000 & shared.inv.SDdist$rightDist < 5000]

## Prepare data for plotting
plt.df <- as.data.frame(shared.inv.SDdist)
plt.df$leftDist <- plt.df$leftDist * -1
plt.df$num <- nrow(plt.df):1 
plt.df$seqnames <- factor(plt.df$seqnames, levels=unique(plt.df$seqnames))
## Set colors for plotting
colors <- rep(c('chocolate3', 'deepskyblue4'), length(unique(plt.df$seqnames)))
colors <- colors[1:length(unique(plt.df$seqnames))]
colors <- rep(colors, table(plt.df$seqnames))
## Make a plot
plt4 <- ggplot(plt.df) + 
  geom_linerange(aes(x=num, ymin=-Inf, ymax=Inf), color='gray') +
  geom_linerange(aes(x=num, ymin=leftDist, ymax=rightDist)) +
  geom_point(aes(x=num, y=0, size=width), color=colors) +
  coord_flip() +
  scale_x_continuous(breaks = plt.df$num, labels = plt.df$seqnames, expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(-2000000, -1000000, 0, 1000000, 2000000), labels = c('2Mb','1Mb','0','1Mb','2Mb')) + 
  scale_color_manual(values=colors) +
  scale_size_continuous(name='Inv size (bp)', labels = comma, breaks = c(100000, 500000, 1000000, 2000000, 3000000)) +
  ylab("Distance to Human SDs (Mb)") +
  xlab("")

plt.df <- as.data.frame(shared.inv.SDdist)
plt.df$leftDist <- plt.df$leftDist * -1
plt.df$num <- 1:nrow(plt.df)
plt.df$seqnames <- factor(plt.df$seqnames, levels=unique(plt.df$seqnames))
plt.df$flankedSD <- FALSE
plt.df$flankedSD[plt.df$leftDist <= 5000 & plt.df$rightDist <= 5000] <- TRUE
plt5 <- ggplot(plt.df, aes(x=1, y=num, fill=flankedSD)) + 
  geom_tile() +
  coord_flip() +
  scale_fill_manual(values=c('chocolate3', 'deepskyblue4')) +
  xlab("") +
  ylab("Inversion #") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## Construct final plot
#final.plt <- plot_grid(plt2, plt1, plt4, nrow = 1, rel_widths = c(1,2,1))
final.plt <- plot_grid(plt3, plt5, ncol = 1, align = 'v', axis = 'lr')
## Export final plot
destination <- file.path(outputDirectory, "putative_polymorphic_inversions.pdf")
#ggsave(final.plt, filename = destination, width = 17, height = 8, limitsize = FALSE, useDingbats=FALSE)
ggsave(final.plt, filename = destination, width = 16, height = 3, limitsize = FALSE, useDingbats=FALSE)
## Export putative polymorphic events
write.table(as.data.frame(shared.inv[,0]), file = file.path(outputDirectory, "shared_HGSVC&NHP_sites.txt"), quote = FALSE, row.names = FALSE)
save(shared.inv, file = file.path(outputDirectory, "shared_HGSVC&NHP_sites.RData"))
## Export supplemental plot
final.plt <- plot_grid(plt1, plt2, plt4, nrow = 1, rel_widths = c(1,2,1))
destination <- file.path(outputDirectory, "putative_polymorphic_inversions_supplemental.pdf")
ggsave(final.plt, filename = destination, width = 17, height = 8, limitsize = FALSE, useDingbats=FALSE)

## Plot all genotypes
gens <- hgsvc.genotypes
hets <- apply(gens, 1, function(x) length(x[x == 'HET']))
homs <- apply(gens, 1, function(x) length(x[x == 'HOM']))
refs <- apply(gens, 1, function(x) length(x[x == 'REF']))
plt.df1 <- data.frame(HOM=homs, HET=hets, ID=1:length(homs))
plt.df1 <- melt(plt.df1, id.vars = 'ID', measure.vars = c('HOM', 'HET'))
plt.df1$facet <- 'HGSVC'
gens <- nhp.genotypes
hets <- apply(gens, 1, function(x) length(x[x == 'HET']))
homs <- apply(gens, 1, function(x) length(x[x == 'HOM']))
refs <- apply(gens, 1, function(x) length(x[x == 'REF']))
plt.df2 <- data.frame(HOM=homs, HET=hets, ID=1:length(homs))
plt.df2 <- melt(plt.df2, id.vars = 'ID', measure.vars = c('HOM', 'HET'))
plt.df2$facet <- 'NHP'
plt.df <- rbind(plt.df1, plt.df2)
plt6 <- ggplot(plt.df) + 
  geom_col(aes(x=ID, y=value, fill=variable)) +
  scale_fill_manual(values = c("coral3", "cornflowerblue"), name="") +
  scale_y_continuous(limits = c(0, 13)) +
  xlab("") +
  ylab("Count (Het|Hom)") +
  facet_grid(facet ~ .)
## Export plot
destination <- file.path(outputDirectory, "polymorphic_inversions_genotypes.pdf")
ggsave(filename = destination, plot = plt6, width = 6, height = 4, useDingbats=FALSE)

## Add genotype information
shared.inv <- shared.inv.SDdist
mcols(shared.inv) <- cbind(mcols(shared.inv), hgsvc.genotypes, nhp.genotypes)
shared.inv$HOMs <- homs
shared.inv$HETs <- hets
shared.inv$REFs <- refs
## Check what genes overlaps with these likely recurrent inversions
gencode29 <- read.table("/home/porubsky/WORK/Great_apes/Annotations/gencode.v29.annotation..processed.txt", header = TRUE)
gencode29.gr <- GRanges(seqnames=gencode29$chr, ranges=IRanges(start=gencode29$start, end=gencode29$end), gene.name=gencode29$gene_name, gene.type=gencode29$gene_type, gene.id=gencode29$gene_id)
hits <- findOverlaps(gencode29.gr, shared.inv, maxgap = 10000)
## Select only protein coding genes
mask <- which(gencode29.gr[queryHits(hits)]$gene.type == 'protein_coding')
hits <- hits[mask]
gene.overlaps <- split(gencode29.gr[queryHits(hits)]$gene.name, subjectHits(hits))
gene.overlaps.merged <- sapply(gene.overlaps, function(x) paste(x, collapse = ";"))
shared.inv$genes.10kb <- ""
shared.inv[as.numeric(names(gene.overlaps))]$genes.10kb <- gene.overlaps.merged
shared.inv.df <- as.data.frame(shared.inv)
## Export data
destination <- file.path(outputDirectory, "polymorphic_inversions_summary.tsv")
write.table(x = shared.inv.df, file = destination, quote = FALSE, row.names = FALSE)
destination <- file.path(outputDirectory, "polymorphic_inversions_geneList.txt")
gene.list <- unique(unlist(gene.overlaps))
write.table(x = gene.list , file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE)

## Compare shared inversion with caceres_biorxiv 2019 ##
########################################################
caceres.df <- read.table(file.path(outputDirectory, 'mcaceres_45polymInv_GRCh36_to_GRCh38.bed'), stringsAsFactors = FALSE) 
caceres.gr <- GRanges(seqnames=caceres.df$V1, ranges=IRanges(start=caceres.df$V2, end=caceres.df$V3), inv.id=caceres.df$V4, type=caceres.df$V5)
## Keep only events 1kb and longer
caceres.gr <- caceres.gr[width(caceres.gr) >= 1000]
## Subtract SD regions from caceres calls
caceres.gr.noSD <- subtractRegions(gr = caceres.gr, remove.gr = reduce(seg.dup.gr), mode = 'flanks')
mcols(caceres.gr.noSD) <- mcols(caceres.gr)[caceres.gr.noSD$ID,]
caceres.gr.noSD$ID <- "Caceres"
## Find overlaps with NHP inversions
#caceres.nhp.overlaps1 <- getReciprocalOverlaps(query = caceres.gr.noSD, subject = nhp.calls, report = 'query')
#caceres.nhp.overlaps2 <- getReciprocalOverlaps(query = caceres.gr.noSD, subject = hgsvc.calls, report = 'query')
#polymorph1 <- caceres.nhp.overlaps1[caceres.nhp.overlaps1$perc.overlap >= 50]
#polymorph2 <- caceres.nhp.overlaps2[caceres.nhp.overlaps2$perc.overlap >= 50]
#subsetByOverlaps(polymorph2, polymorph1, invert = T)
caceres.overlap <- getReciprocalOverlaps(query = caceres.gr.noSD, subject = shared.inv, report = 'query')
caceres.gr.noSD$shared.inv <- ''
caceres.gr.noSD$shared.inv[caceres.overlap$perc.overlap >= 50]

## Genotype regions from caceres
regions.to.genotype <- caceres.gr.noSD[,0]
## Genotype chimpanzee data
genotyped.regions <- genotypeRegions(regions = regions.to.genotype, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotyped.regions <- genotypeRegions(regions = genotyped.regions, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')

nhp.genotypes <- mcols(genotyped.regions)[,c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan')]
mcols(caceres.gr.noSD) <- cbind(mcols(caceres.gr.noSD), nhp.genotypes)

## Plot the data
plt.df <- as.data.frame(caceres.gr.noSD)
plt.df$inv.id <- factor(plt.df$inv.id, levels=plt.df$inv.id)

plt.df <- melt(data = plt.df, 
               measure.vars = c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan'), 
               id.vars = c('inv.id', 'type')
)
plt.df$level <- rep(1:4, table(plt.df$variable))

plt <- plt.df %>% ggplot() + 
  geom_tile(aes(x=inv.id, y=level, fill=value)) +
  geom_tile(aes(x=inv.id, y=5, fill=type), inherit.aes = F) +
  scale_fill_manual(values=c("cornflowerblue","coral3", "white", "cadetblue2", "cadetblue4", "gray")) +
  scale_y_continuous(breaks=1:5, labels = c('genoT_chimpanzee', 'genoT_bonobo', 'genoT_gorilla', 'genoT_orangutan', 'inv.type')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("") + ylab("")
## Export final plot
destination <- file.path(outputDirectory, "caceres_regenotyped.pdf")          
ggsave(filename = destination, plot = plt, width = 10, height = 4, useDingbats=FALSE) 

message("DONE!!!")

## TMP code ##
# ## Prepare data for plotting
# mcols(hgsvc.simple.inv.gr) <- cbind(mcols(hgsvc.simple.inv.gr), ape.genotypes)
# hgsvc.plus.apes.genotypes <- as.data.frame(hgsvc.simple.inv.gr)
# hgsvc.plus.apes.genotypes$ID <- as.character(hgsvc.simple.inv.gr[,0])
# #hgsvc.plus.apes.genotypes$chr <- sapply(hgsvc.plus.apes.genotypes$ID, function(x) strsplit(x, ":")[[1]][1])
# colnames(hgsvc.plus.apes.genotypes) <- gsub(colnames(hgsvc.plus.apes.genotypes), pattern = 'genoT_', replacement = '')
# hgsvc.plus.apes.genotypes$roi.num <- 1:nrow(hgsvc.plus.apes.genotypes)
# ## From wide to long format
# plt.df <- reshape2::melt(hgsvc.plus.apes.genotypes, id.vars = c('ID','roi.num'),
#                measure.vars = c('HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19239','NA19238','NA19240','chimpanzee','bonobo','gorilla','orangutan'))
# plt.df$variable2 <- rep(c('human','NHP'), c(9*length(hgsvc.simple.inv.gr), 4*length(hgsvc.simple.inv.gr)))
# 
# plt.df.counts <- plt.df %>% group_by(variable2, roi.num, value) %>% summarise(count=n())
# plt.df.counts$count[plt.df.counts$variable2 == 'NHP'] <- plt.df.counts$count[plt.df.counts$variable2 == 'NHP'] * -1
# plt.df.counts$value <- factor(plt.df.counts$value, levels = c("REF","AMB","lowReads", "HET","HOM"))
# plt.df.counts$color <- recode(plt.df.counts$value, "REF"="gray","AMB"="gray","lowReads"="gray", "HET"="cornflowerblue","HOM"="coral3")
# 
# ## Find likely polymorphic sites
# poly.inv <- plt.df.counts %>% group_by(roi.num, variable2) %>% mutate(polymorph = (value == 'HOM' | value == 'HET'))
# poly.inv <- poly.inv[poly.inv$polymorph == TRUE,]
# genoT.inBoth <- split(poly.inv$roi.num, poly.inv$variable2)
# shared.roi.num <- intersect(genoT.inBoth[[1]], genoT.inBoth[[2]])
# ## Select only regions that pass 50% recoprocal overlap between HGSVC regions and NHP regions
# pass.perc.overlap <- which(percOverlap.df$perc.overlap >= 50)
# shared.roi.num.pass.percOverlap <- intersect(shared.roi.num, pass.perc.overlap)
# putative.polymorph <- hgsvc.plus.apes.genotypes[shared.roi.num.pass.percOverlap,]
# putative.polymorh.gr <- hgsvc.simple.inv.gr[shared.roi.num.pass.percOverlap]
# ## Export putative polymorphic events
# write.table(putative.polymorph, file = file.path(outputDirectory, "putative_polymorphic_sites.txt"), quote = FALSE, row.names = FALSE)
# save(putative.polymorh.gr, file = file.path(outputDirectory, "putative_polymorphic_sites.RData"))
# 
# ## Prepare chromosome annotation
# chr.lines.rle <- Rle(hgsvc.plus.apes.genotypes$seqnames)
# chr.lines <- cumsum(c(0.5, runLength(chr.lines.rle)))
# chr.lines <- as.data.frame(reformat(chr.lines))
# chr.lines$midpoint <- chr.lines$V1 + (chr.lines$V2 - chr.lines$V1)/2
# chr.lines$chr <- runValue(chr.lines.rle)
# chr.lines$chr.label <- gsub(chr.lines$chr, pattern = 'chr', replacement = '')
# chr.lines$color <- rep(c('cadetblue4', 'lightblue1'), nrow(chr.lines))[1:nrow(chr.lines)]
# ## Prepare sidebar annotation
# sidebar <- data.frame(ymin=c(-4,0), ymax=c(0,9), color=c('sienna1', 'slateblue3'), ID=c('NHP', 'HGSVC'))
# sidebar$midpoint <- sidebar$ymin + (sidebar$ymax - sidebar$ymin)/2
# 
# ## Construct final plot
# plt1 <- ggplot(plt.df.counts) + 
#   geom_col(aes(x=roi.num, y=count, fill=color)) + 
#   geom_rect(data=chr.lines, aes(xmin=V1, xmax=V2, ymin=9, ymax=10, fill=color), inherit.aes = FALSE) +
#   geom_text(data=chr.lines, aes(x=midpoint, y=9.5, label=chr.label), color="white") +
#   geom_rect(data=sidebar, aes(xmin=-3, xmax=0, ymin=ymin, ymax=ymax, fill=color), inherit.aes = FALSE) +
#   geom_text(data=sidebar, aes(y=midpoint, x=-1.75, label=ID), color="white", angle=90) +
#   geom_point(data=percOverlap.df, aes(x=roi.num, y=-4.25, color=perc.overlap), shape=15, size=2, inherit.aes = FALSE) +
#   geom_point(data=putative.polymorph, aes(x=roi.num, y=-4.5), shape=24, size=2, inherit.aes = FALSE) +
#   scale_fill_identity() +
#   ylab("HOM|HET count") +
#   xlab(paste0("Regions n=", nrow(hgsvc.plus.apes.genotypes)))
#   
# ## Plot putative polymorphic sites per chromosome
# plt2 <- putative.polymorph %>% group_by(seqnames) %>% summarise(count=n()) %>%
#   ggplot(aes(x=seqnames, y=count)) + geom_col() +
#   geom_text(aes(label=count), color='white', vjust=0.5, hjust=1) +
#   xlab("") + 
#   coord_flip()
# 
# ## Plot distance of inversion breakpoints to human SDs
# putative.polymorh.inv.SDdist <- range2rangeDistance(gr=putative.polymorh.gr, userTrack=seg.dup.gr, allow.overlap = TRUE)
# 
# ## Prepare data for plotting
# plt.df <- as.data.frame(putative.polymorh.inv.SDdist)
# plt.df$leftDist <- plt.df$leftDist * -1
# plt.df$num <- nrow(plt.df):1 
# plt.df$seqnames <- factor(plt.df$seqnames, levels=unique(plt.df$seqnames))
# ## Set colors for plotting
# colors <- rep(c('chocolate3', 'deepskyblue4'), length(unique(plt.df$seqnames)))
# colors <- colors[1:length(unique(plt.df$seqnames))]
# colors <- rep(colors, table(plt.df$seqnames))
# ## Make a plot
# plt3 <- ggplot(plt.df) + 
#   geom_linerange(aes(x=num, ymin=-Inf, ymax=Inf), color='gray') +
#   geom_linerange(aes(x=num, ymin=leftDist, ymax=rightDist)) +
#   geom_point(aes(x=num, y=0, size=width), color=colors) +
#   coord_flip() +
#   scale_x_continuous(breaks = plt.df$num, labels = plt.df$seqnames) +
#   scale_y_continuous(labels = comma) + 
#   scale_color_manual(values=colors) +
#   scale_size_continuous(name='Inv size (bp)', labels = comma, breaks = c(100000, 500000, 1000000, 2000000, 3000000)) +
#   ylab("Disntance to Human SDs (bp)") +
#   xlab("")
# 
# final.plt <- plot_grid(plt1, plt2, plt3, nrow = 1, rel_widths = c(6,1,2))
