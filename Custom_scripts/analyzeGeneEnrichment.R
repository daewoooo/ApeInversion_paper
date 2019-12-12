## This script perform gene enrichment analysis of differentially expressed genes from
## Pollen et al. against SDs that flanks NHP inversions incomparison to those that don't

## Load requried libraries
suppressPackageStartupMessages( library(regioneR) )
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg38.masked) )
suppressPackageStartupMessages( library(dplyr) )

message("Testing enrichement of Pollen et al. genes around simple inversions ...")

## Get masked parts of the genome (gaps etc.)
hg38.mask <- getMask(BSgenome.Hsapiens.UCSC.hg38.masked)

outputDirectory <- "/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/"

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/Annotations/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)
seg.dup.regions <- reduce(seg.dup.gr)

## Load all StrandS inversion calls [includes putative human specific inversions] ##
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']

## Split inversion as flanked by SDs and those that are not
all.SimpleInversion.calls.filt.dists <- range2rangeDistance(gr=all.SimpleInversion.calls.filt, userTrack=seg.dup.gr, allow.overlap = TRUE)
SDflank <- which(all.SimpleInversion.calls.filt.dists$leftDist < 5000 & all.SimpleInversion.calls.filt.dists$rightDist < 5000 & all.SimpleInversion.calls.filt.dists$leftDist >= 0 & all.SimpleInversion.calls.filt.dists$rightDist >= 0)
simpleINV.SDflank <- all.SimpleInversion.calls.filt[SDflank]
simpleINV.SDflank$SVclass <- 'SDflankINV'
simpleINV.noSDflank <- all.SimpleInversion.calls.filt[-SDflank]
simpleINV.noSDflank$SVclass <- 'noSDflankINV' 
## Get SDs that overlap with inversions flanked by SDs
inversion.breakpoints <- getRegionBoundaries(simpleINV.SDflank)
breakpointSDs <- subsetByOverlaps(seg.dup.regions, inversion.breakpoints)
## Get SDs that do not overlap with predicted inversion breakpoint
allOtherSDs <- subsetByOverlaps(seg.dup.regions, inversion.breakpoints, invert = TRUE)

## Load exported gene lists ##
##############################
## Alex Pollen
enriched.gene.list <- get(load("/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/upregul.gene.list.RData"))

## Ranges to test enrichment
test.grl <- list(breakpointSDs=breakpointSDs, allOtherSDs=allOtherSDs)

## Perform enrichment analysis ##
#################################
## The region set A is the one randomized. And the randomized A is then
## evaluated with the evaluation function, sometimes (for example, if evaluating overlaps) against teh region set B.

nperm <- 1000
permuted.l <- list()
observed.l <- list()
zscores.l <- list()
for (i in seq_along(test.grl)) {
  dataset <- names(test.grl[i])
  message("Processing dataset: ", dataset)
  
  enrich.analysis <- list()
  test.gr <- test.grl[[i]]
  for (j in seq_along(enriched.gene.list)) {
    ID <- names(enriched.gene.list[j])
    ## Skip over Concerted gene list
    if (ID == 'CONCERTED') {
      next
    }
    message("    Working on ", ID, " gene list ...")
    sub.genes <- enriched.gene.list[[j]]
    
    ## Run regioneR
    message("        Running permTEST for chimp genes ...")
    test.chimp <- permTest(
        B=test.gr, 
        A=sub.genes[['chimp']], 
        randomize.function=circularRandomizeRegions,
        #randomize.function=randomizeRegions, 
        evaluate.function=numOverlaps, 
        genome="hg38", 
        ntimes=nperm, 
        allow.overlaps=FALSE, 
        per.chromosome=TRUE, 
        mask=hg38.mask, 
        count.once=FALSE,
        mc.set.seed=FALSE,
        mc.cores=4)
      
    message("        Running permTEST for human genes ...")
    test.human <- permTest(
        B=test.gr,
        A=sub.genes[['human']],
        randomize.function=circularRandomizeRegions,
        #randomize.function=randomizeRegions, 
        evaluate.function=numOverlaps, 
        genome="hg38", 
        ntimes=nperm, 
        allow.overlaps=FALSE, 
        per.chromosome=TRUE, 
        mask=hg38.mask, 
        count.once=FALSE,
        mc.set.seed=FALSE,
        mc.cores=4)
    
    ## Get the results from regioneR
    test.chimp <- test.chimp$numOverlaps
    test.human <- test.human$numOverlaps
    ## Construct data object for plotting
    chimp.perm <- test.chimp$permuted
    human.perm <- test.human$permuted
    chimp.obs <- test.chimp$observed
    human.obs <- test.human$observed
    chimp.pval <- round(test.chimp$pval, digits = 5)
    human.pval <- round(test.human$pval, digits = 5)  
    permuted <- data.frame(perm = c(chimp.perm, human.perm),
                           sample = rep(c('chimp', 'human'), c(length(chimp.perm), length(human.perm))),
                           pval = rep(c(chimp.pval, human.pval), c(length(chimp.perm), length(human.perm))),
                           ID = ID, dataset = dataset
    )
    observed <- data.frame(obs = c(chimp.obs, human.obs), sample = c('chimp', 'human'), ID = ID, dataset = dataset)
    
    zscores <- data.frame(chimp=test.chimp$zscore, human=test.human$zscore, ID = ID, dataset = dataset)
    
    permuted.l[[length(permuted.l) + 1]] <- permuted
    observed.l[[length(observed.l) + 1]] <- observed
    zscores.l[[length(zscores.l) + 1]] <- zscores
  }
}
## Compile all results together
plt.df.permuted <- do.call(rbind, permuted.l)
plt.df.observed <- do.call(rbind, observed.l)
plt.df.zscores <- do.call(rbind, zscores.l)

## Save raw data
plt.data <- list(permuted=plt.df.permuted, observed=plt.df.observed, zscores=plt.df.zscores)
destination <- file.path(outputDirectory, "SDgeneEnrichment.RData")
save(plt.data, file = destination)

## Multiple testing correction
pvals <- plt.df.permuted %>% group_by(dataset, sample, ID) %>% summarise(pval=unique(pval))
pvals$pval.corr <-p.adjust(p = pvals$pval, method = 'bonferroni', n = nrow(pvals))

## Plot the results
plt1 <- ggplot(plt.df.permuted, aes(x = ID, y = perm, fill=factor(sample))) + 
  geom_violin(position = position_dodge()) +
  stat_summary(position = position_dodge(width = 0.9), fun.data='mean_sdl', size=0.5, geom="pointrange", color="white") +
  geom_point(data=plt.df.observed, aes(x = ID, y = obs, group=factor(sample)), position = position_dodge(width = 0.9), size=3, color='red', inherit.aes = FALSE) +
  #geom_text(data=plt.df.permuted, aes(x = ID, y = 0, group=factor(sample), label=pval), position = position_dodge(width = 0.9), inherit.aes = FALSE) +
  geom_text(data=pvals, aes(x = ID, y = 0, group=sample, label=pval.corr), position = position_dodge(width = 0.9), inherit.aes = FALSE) +
  scale_fill_manual(values = c('#31a354','#a6611a'), name="Sample") +
  xlab("") +
  ylab("Permuted overlaps (n)") +
  facet_grid(dataset ~ ., scales = 'free')

plt.df.zscores <- reshape2::melt(plt.df.zscores)
plt2 <- ggplot(plt.df.zscores) +
  geom_col(aes(x=ID, y=value, group=variable, fill=variable), position = position_dodge()) +
  scale_fill_manual(values = c('#31a354','#a6611a'), name="Sample") +
  xlab("") +
  ylab("z-score") +
  facet_grid(. ~ dataset)

final.plt <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(2,1))
## Save final plot
destination <- file.path(outputDirectory, "SDgeneEnrichment.pdf")
ggsave(filename = destination, plot = final.plt, width = 14, height = 5)

message("DONE!!!")