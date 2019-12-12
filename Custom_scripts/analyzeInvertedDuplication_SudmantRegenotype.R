## Load required libraries
library(tidyr)
library(primatR)
library(BSgenome.Hsapiens.UCSC.hg38)

outputdirectory <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/"

## Set to TRUE if you want to consider only events with max CN 4.
filterByCN <- TRUE

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

## Load centromere Track
cent <- read.table("/home/porubsky/WORK/Great_apes/Annotations/centromeres_GRCh38.bed.gz")
cent.gr <- GRanges(seqnames=cent$V2, ranges=IRanges(start=cent$V3, end=cent$V4))

## Create set of blacklisted regions ##
#######################################
## +/- 500kb from the centromeres
#blacklist1 <- GRanges(seqnames=cent$V2, ranges=IRanges(start=cent$V3-500000, end=cent$V4+500000))
## 5Mb from each chromosome end
#chrom.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22, 'X'))]
#blacklist2 <- GRanges(seqnames=names(chrom.len), ranges=IRanges(start = 1, end = 5000000))
#blacklist3 <- GRanges(seqnames=names(chrom.len), ranges=IRanges(start = chrom.len - 5000000, end = chrom.len))
## Merge all blacklists together
#blacklist <- reduce(c(blacklist1, blacklist2, blacklist3))
## Blacklist regions that are likely interchromosomal links
links.folder <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/Mapping_intSites/"
## Get gorilla blacklisted regions
gorilla.links <- get(load(file.path(links.folder, 'gorilla/gorilla_GRCh38align_links.RData'))) 
gorilla.links <- gorilla.links$inter.links
ranges1 <- gorilla.links[,0]
ranges2 <- gorilla.links$to.gr
blacklist.gorilla <- reduce(c(ranges1, ranges2))
## Get orangutan blacklisted regions
orangutan.links <- get(load(file.path(links.folder, 'orangutan/orangutan_GRCh38align_links.RData'))) 
orangutan.links <- orangutan.links$inter.links
ranges1 <- orangutan.links[,0]
ranges2 <- orangutan.links$to.gr
blacklist.orangutan <- reduce(c(ranges1, ranges2))
## Get chimpanzee blacklisted regions
chimpanzee.links <- get(load(file.path(links.folder, 'chimpanzee/chimpanzee_GRCh38align_links.RData'))) 
chimpanzee.links <- chimpanzee.links$inter.links
ranges1 <- chimpanzee.links[,0]
ranges2 <- chimpanzee.links$to.gr
blacklist.chimpanzee <- reduce(c(ranges1, ranges2))
## Get bonobo blacklisted regions
bonobo.links <- get(load(file.path(links.folder, 'bonobo/bonobo_GRCh38align_links.RData'))) 
bonobo.links <- bonobo.links$inter.links
ranges1 <- bonobo.links[,0]
ranges2 <- bonobo.links$to.gr
blacklist.bonobo <- reduce(c(ranges1, ranges2))

## Load Sudmant_2013 lineage-specific duplications
linspec.dups <- read.table("/home/porubsky/WORK/Great_apes/Published_data/Zev_2018/lineageSpecificFixedDuplicationsSudmant2013.csv", sep = ",", skip = 1, header = TRUE, stringsAsFactors = FALSE)
## CN genotypes
CN.gens <- linspec.dups[,5:ncol(linspec.dups)]
CN.gens <- CN.gens[,grep(colnames(CN.gens), pattern = 'Homo', ignore.case = TRUE, invert = TRUE)] ## Remove human's CN calls
## Create GRanges object
linspec.dups <- linspec.dups[,c(1,2,3)]
linspec.dups <- tidyr::separate(linspec.dups, col = region, into = c('chr', 'start', 'end'), sep = ":|-")
linspec.dups.gr <- GRanges(seqnames=linspec.dups$chr, 
                           ranges=IRanges(start=as.numeric(linspec.dups$start), 
                                          end=as.numeric(linspec.dups$end)), 
                                          lineage=linspec.dups$lineage)
## Add CN genotypes
mcols(linspec.dups.gr) <- cbind(mcols(linspec.dups.gr), CN.gens)
## Keep only regions 10kb and longer
linspec.dups.gr <- linspec.dups.gr[width(linspec.dups.gr) >= 10000]
## Keep only regions not present in human (Hsa)
linspec.dups.gr <- linspec.dups.gr[grep(linspec.dups.gr$lineage, pattern = 'Hsa|HUMAN|Hde', invert = TRUE)]
## Remove calls from the blacklisted regions
blacklist <- c(blacklist.gorilla, blacklist.orangutan, blacklist.chimpanzee, blacklist.bonobo)
linspec.dups.gr <- subsetByOverlaps(linspec.dups.gr, blacklist, invert = TRUE)
## Sort ranges
seqlevels(linspec.dups.gr) <- paste0('chr', c(1:22,'X'))
linspec.dups.gr <- sort(linspec.dups.gr)
## Remove unnecessary characters
linspec.dups.gr$lineage <- gsub(linspec.dups.gr$lineage, pattern = 'ils_clades_', replacement = '')
names <- linspec.dups.gr$lineage
## Decode orangutan
names <- gsub(names, pattern = 'Ppy|Pab|ORANGUTAN', replacement = 'orangutan')
## Decode gorilla
names <- gsub(names, pattern = 'Ggod|Ggog|Gbeg|EASTERN_GORILLA|WESTERN_GORILLA', replacement = 'gorilla')
## Decode bonobo
names <- gsub(names, pattern = 'Ppa|BONOBO', replacement = 'bonobo')
## Decode chimpanzee
names <- gsub(names, pattern = 'Ptrs|Ptrt|Ptre|Ptrv|CHIMP', replacement = 'chimpanzee')
## Expand ranges
names.list <- strsplit(names, "-")
names.list <- lapply(names.list, unique) ## Keep unique names
reps <- lengths(names.list)
linspec.dups.expand.gr <- rep(linspec.dups.gr, reps)
linspec.dups.expand.gr$ID <- unlist(names.list)
## Export table
## Liftover coordinates (GRCh36 -> GRCh38) from this table before proceeding further !!!
## Remove the range that failed to lift from GRCh36 to GRCh38 (chr2:242366422-242417043)
linspec.dups.expand.gr <- linspec.dups.expand.gr[!start(linspec.dups.expand.gr) == 242366422]	
export.df <- as.data.frame(linspec.dups.expand.gr[,'ID'])
destination <- file.path(outputdirectory, "sudmantLinSpecDups_GRCh36.bed")
write.table(export.df, file = destination, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = '\t')

## Load ape's composite files
chimp.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_chimpanzee_final.RData"))
chimp.data$ID <- 'chimpanzee'
bonobo.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_bonobo_final.RData"))
bonobo.data$ID <- 'bonobo'
gorilla.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_gorilla_final.RData"))
gorilla.data$ID <- 'gorilla'
orangutan.data <- get(load("/home/porubsky/WORK/Great_apes/Composite_files/syncReads_orangutan_final.RData"))
orangutan.data$ID <- 'orangutan'

## Load liftover data
linspec.dups.lifted <- read.table(file.path(outputdirectory, "sudmantLinSpecDups_GRCh36toGRCh38.bed"))
linspec.dups.lifted.gr <- GRanges(seqnames=linspec.dups.lifted$V1, ranges=IRanges(start=linspec.dups.lifted$V2, end=linspec.dups.lifted$V3), ID=linspec.dups.lifted$V6)

if (filterByCN) {
  ## Add mean CN info from the original(unlifted) data object
  CN.tab <- mcols(linspec.dups.expand.gr)[-c(1,ncol(mcols(linspec.dups.expand.gr)))]
  CN.tab.df <- as.data.frame(CN.tab)
  mean.CN <- rowMeans(CN.tab.df)
  linspec.dups.lifted.gr$mean.CN <- mean.CN
  ## Keep only region with mean CN < 5
  linspec.dups.lifted.gr <- linspec.dups.lifted.gr[linspec.dups.lifted.gr$mean.CN < 5]
} 

## Regenotype these regions using Strand-seq
## Genotype chimpanzee data
genotype.chimp <- linspec.dups.lifted.gr[linspec.dups.lifted.gr$ID == 'c']
genotype.chimp <- genotypeRegions(regions = genotype.chimp, directional.reads = chimp.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'chimpanzee')
## Genotype bonobo data
genotype.bonobo <- linspec.dups.lifted.gr[linspec.dups.lifted.gr$ID == 'b']
genotype.bonobo <- genotypeRegions(regions = genotype.bonobo, directional.reads = bonobo.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'bonobo')
## Genotype gorilla data
genotype.gorilla <- linspec.dups.lifted.gr[linspec.dups.lifted.gr$ID == 'g']
genotype.gorilla  <- genotypeRegions(regions = genotype.gorilla , directional.reads = gorilla.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'gorilla')
## Genotype orangutan data
genotype.orang <- linspec.dups.lifted.gr[linspec.dups.lifted.gr$ID == 'o']
genotype.orang <- genotypeRegions(regions = genotype.orang, directional.reads = orangutan.data, blacklist = seg.dup.gr, min.reads = 5, alpha = 0.05, index = 'orangutan')
## Compile all results into a table
genotype.chimp.gen <- table(genotype.chimp$genoT_chimpanzee)
chimp.df <- data.frame(ID='chimpanzee', inverted=genotype.chimp.gen['HET'] + genotype.chimp.gen['HOM'],
                       direct=genotype.chimp.gen['REF'])
genotype.bonobo.gen <- table(genotype.bonobo$genoT_bonobo)
bonobo.df <- data.frame(ID='bonobo', inverted=genotype.bonobo.gen['HET'] + genotype.bonobo.gen['HOM'],
                       direct=genotype.bonobo.gen['REF'])
genotype.gorilla.gen <- table(genotype.gorilla$genoT_gorilla)
gorilla.df <- data.frame(ID='gorilla', inverted=genotype.gorilla.gen['HET'] + genotype.gorilla.gen['HOM'],
                        direct=genotype.gorilla.gen['REF'])
genotype.orangutan.gen <- table(genotype.orang$genoT_orangutan)
orangutan.df <- data.frame(ID='orangutan', inverted=genotype.orangutan.gen['HET'] + genotype.orangutan.gen['HOM'],
                         direct=genotype.orangutan.gen['REF'])
final.tab <- rbind(chimp.df, bonobo.df, gorilla.df, orangutan.df)
## Calculate proportion of inverted versus direct duplications
sums <- rowSums(final.tab[,c('inverted','direct')])
final.tab.long <- reshape2::melt(final.tab, id.vars = 'ID', measure.vars = c('direct', 'inverted'), variable.name = 'Dir', factorsAsStrings = TRUE)
final.tab.long$perc <- (final.tab.long$value / rep(sums, 2)) * 100

## Calculate significance of a difference between count of direct vs inverted invDups
plt.df <- final.tab.long
plt.df$Dir <- factor(plt.df$Dir, levels = c('direct', 'inverted'))
signif.results.all <- list()
for (i in unique(plt.df$ID)) {
  sub.tab <- plt.df[plt.df$ID == i,]
  data.matrix <- matrix(c(sub.tab$value[sub.tab$Dir == 'inverted'],
                          sum(sub.tab$value)/2,
                          sub.tab$value[sub.tab$Dir == 'direct'],
                          sum(sub.tab$value)/2),
                        nrow = 2
  )
  ## Calculate Chi-square
  chisq.pval <- chisq.test(data.matrix)$p.value
  suppressWarnings( fish.oddsr <- fisher.test(data.matrix)$estimate )
  signif.results <- data.frame(ID=i, chisq.pval=chisq.pval, fish.oddsr=unname(fish.oddsr))
  signif.results.all[[length(signif.results.all) + 1]] <- signif.results
}
signif.results.all <- do.call(rbind, signif.results.all)
## Do multiple testing correction (Bonferroni)
signif.results.all$chisq.pval.BonfCorr <- p.adjust(p = signif.results.all$chisq.pval, method = 'bonferroni', n = length(signif.results.all$chisq.pval))
signif.results.all$pVal <- format(signif.results.all$chisq.pval.BonfCorr, scientific = TRUE, digits = 2)

plt1 <- ggplot(plt.df, aes(x=ID, y=value, group=ID)) +
  geom_col(aes(fill=Dir)) +
  geom_text(data=signif.results.all, aes(x=ID, y=Inf, label=pVal), inherit.aes = FALSE, vjust=1) +
  geom_text(aes(group=ID, label=value), position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("Counts") +
  scale_fill_manual(values = c('darkslategray4', 'indianred4'), name="Direction") +
  scale_x_discrete(labels = c('B', 'C', 'G', 'O'))

plt2 <- ggplot(plt.df) +
  geom_col(aes(x=ID, y=perc, fill=Dir)) +
  xlab("") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c('darkslategray4', 'indianred4')) +
  scale_x_discrete(labels = c('B', 'C', 'G', 'O'))
## Save final plot
final.plt <- plot_grid(plt1, plt2, nrow = 1)
#destination <- file.path(outputdirectory, "invervedVSdirect_duplications_sudmant2013_regenotyped.pdf")
destination <- file.path(outputdirectory, "invervedVSdirect_duplications_sudmant2013_regenotyped_filtCN.pdf")
ggsave(final.plt, filename = destination, width = 10, height = 5, useDingbats=FALSE)
destination <- file.path(outputdirectory, "invervedVSdirect_duplications_sudmant2013_regenotyped_filtCN.txt")
write.table(plt.df, file = destination, append = FALSE, quote = FALSE, row.names = FALSE)

message("DONE!!!")