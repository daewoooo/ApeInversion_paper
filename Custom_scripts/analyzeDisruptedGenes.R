## Load required libraries
suppressPackageStartupMessages( library(primatR) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(RColorBrewer) )

outputDirectory <- "/home/porubsky/WORK/Great_apes/Disrupted_genes/"

message("Searching genes disrupted by inversions ...")

## Load gene lists
gencode29 <- read.table("/home/porubsky/WORK/Great_apes/Annotations/gencode.v29.annotation..processed.txt", header = TRUE)
gencode29.gr <- GRanges(seqnames=gencode29$chr, ranges=IRanges(start=gencode29$start, end=gencode29$end), gene.name=gencode29$gene_name, gene.type=gencode29$gene_type, gene.id=gencode29$gene_id)

## Load complete dataset including putative human specific inversions
all.inversion.calls.filt <- get(load("/home/porubsky/WORK/Great_apes/Final_INV_calls/all.inversion.calls.filt.RData"))
## Retain only simple inversions
all.SimpleInversion.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'INV']
## Retain only inverted duplications
all.invertedDups.calls.filt <- all.inversion.calls.filt[all.inversion.calls.filt$SVclass == 'invDup']

## Create a non-redundant set from NHP simple inversions
overlaps.gr <- getDisjointOverlapsWeighted(all.SimpleInversion.calls.filt, percTh = 50)
NHP.nonred.gr <- collapseBins(overlaps.gr[order(overlaps.gr$sub.group)], id.field = 11)
NHP.individ <- split(overlaps.gr$ID, overlaps.gr$sub.group)
NHP.nonred.gr$ID <- sapply(NHP.individ, function(x) paste(x, collapse = ","))

## Create a non-redundant set from NHP inverted duplications
overlaps.gr <- getDisjointOverlapsWeighted(all.invertedDups.calls.filt, percTh = 50)
NHP.nonred.invDup.gr <- collapseBins(overlaps.gr[order(overlaps.gr$sub.group)], id.field = 11)
NHP.individ <- split(overlaps.gr$ID, overlaps.gr$sub.group)
NHP.nonred.invDup.gr$ID <- sapply(NHP.individ, function(x) paste(x, collapse = ","))

## Get simple inversions disrupting genes  ##
#############################################
max.dist <- 1000
SimpleInversion.starts <- GRanges(seqnames=seqnames(NHP.nonred.gr), ranges=IRanges(start=start(NHP.nonred.gr), end=start(NHP.nonred.gr)+1))
mcols(SimpleInversion.starts) <- mcols(NHP.nonred.gr)[1:6]
SimpleInversion.starts$INV.id <- as.character(NHP.nonred.gr[,0])
SimpleInversion.ends <- GRanges(seqnames=seqnames(NHP.nonred.gr), ranges=IRanges(start=end(NHP.nonred.gr), end=end(NHP.nonred.gr)+1))
mcols(SimpleInversion.ends) <- mcols(NHP.nonred.gr)[1:6]
SimpleInversion.ends$INV.id <- as.character(NHP.nonred.gr[,0])
SimpleInversion.breakpoints <- c(SimpleInversion.starts, SimpleInversion.ends)

## Get overlaps between INV breakpoints and gene list
gene.hits <- findOverlaps(gencode29.gr, SimpleInversion.breakpoints, maxgap = max.dist)
## Split overlaping genes by breakpoint
broken.genes.simpINV <- gencode29.gr[queryHits(gene.hits)]
broken.genes.simpINV$ID <- SimpleInversion.breakpoints[subjectHits(gene.hits)]$ID
broken.genes.simpINV$idx <- subjectHits(gene.hits)
broken.genes.simpINV$inv.ID <- SimpleInversion.breakpoints$INV.id[subjectHits(gene.hits)]
## Keep only protein coding genes
broken.genes.simpINV.prot <- broken.genes.simpINV[broken.genes.simpINV$gene.type == 'protein_coding']
broken.genes.simpINV.prot$gene.name <- as.character(broken.genes.simpINV.prot$gene.name)
broken.genes.simpINV.prot.toPrint <- unique(broken.genes.simpINV.prot$gene.name)
#gene.IDs.toPrint <- sapply(as.character(gene.IDs.toPrint), function(x) strsplit(x, "\\.")[[1]][1])
## Export list of disrupted protein coding genes 
destination <- file.path(outputDirectory, "disrupted_proteinCoding_genes_bySimpleInv.txt")
write.table(x = broken.genes.simpINV.prot.toPrint, file = destination, append = FALSE, row.names = FALSE, quote = FALSE, col.names = FALSE)

## Get inverted duplications disrupting genes  ##
#################################################
max.dist <- 1000
invDup.breakpoints <- getRegionBoundaries(NHP.nonred.invDup.gr)
## Get overlaps between INV breakpoints and gene list
gene.hits <- findOverlaps(gencode29.gr, invDup.breakpoints, maxgap = max.dist)
## Split overlaping genes by breakpoint
broken.genes.invDup <- gencode29.gr[queryHits(gene.hits)]
broken.genes.invDup$idx <- subjectHits(gene.hits)
broken.genes.invDup$inv.ID <- invDup.breakpoints$region.ID[subjectHits(gene.hits)]
## Keep only protein coding genes
broken.genes.invDup.prot <- broken.genes.invDup[broken.genes.invDup$gene.type == 'protein_coding']
broken.genes.invDup.prot$gene.name <- as.character(broken.genes.invDup.prot$gene.name)
broken.genes.invDup.prot.toPrint <- unique(broken.genes.invDup.prot$gene.name)
## Export list of disrupted protein coding genes 
destination <- file.path(outputDirectory, "disrupted_proteinCoding_genes_byInvDups.txt")
write.table(x = broken.genes.invDup.prot.toPrint, file = destination, append = FALSE, row.names = FALSE, quote = FALSE, col.names = FALSE)
## Print all broken genes by simple inversions and inverted duplications
all.broken.genes <- union(broken.genes.simpINV.prot.toPrint, broken.genes.invDup.prot.toPrint)
double.broken.genes <- intersect(broken.genes.simpINV.prot.toPrint, broken.genes.invDup.prot.toPrint)

## Get genes duplicated within inverted duplications  ##
########################################################
gene.hits <- findOverlaps(gencode29.gr, NHP.nonred.invDup.gr, type = 'within')
## Split overlapping genes by breakpoint
dup.genes.invDup <- gencode29.gr[queryHits(gene.hits)]
dup.genes.invDup$idx <- subjectHits(gene.hits)
dup.genes.invDup$inv.ID <- as.character(NHP.nonred.invDup.gr)[subjectHits(gene.hits)]
## Keep only protein coding genes
dup.genes.invDup.prot <- dup.genes.invDup[dup.genes.invDup$gene.type == 'protein_coding']
dup.genes.invDup.prot$gene.name <- as.character(dup.genes.invDup.prot$gene.name)
dup.genes.invDup.prot.toPrint <- unique(dup.genes.invDup.prot$gene.name)
## Get gene category counts for genes embedded in invDups
dup.genes.invDup <- dup.genes.invDup[!duplicated(dup.genes.invDup$gene.name)]
gene.type.counts <- sort(table(dup.genes.invDup$gene.type), decreasing = TRUE)
other <- sum(gene.type.counts[6:length(gene.type.counts)])
names(other) <- 'other'
gene.type.counts <- c(gene.type.counts[1:5], other)
## Export list of duplicated protein coding genes within invDups
destination <- file.path(outputDirectory, "duplicated_proteinCoding_genes_byInvDups.txt")
write.table(x = dup.genes.invDup.prot.toPrint, file = destination, append = FALSE, row.names = FALSE, quote = FALSE, col.names = FALSE)

## Export all likely disrupted genes
all.broken.genes.df <- data.frame(gene=c(broken.genes.simpINV.prot.toPrint, broken.genes.invDup.prot.toPrint, dup.genes.invDup.prot.toPrint), 
                                  ID=rep(c('INV','invDup', 'insideDup'), c(length(broken.genes.simpINV.prot.toPrint), length(broken.genes.invDup.prot.toPrint), length(dup.genes.invDup.prot.toPrint))))
destination <- file.path(outputDirectory, "proteinCodingGenes_changedBySimpleINVandInvDups.txt")
write.table(x = all.broken.genes.df, file = destination, append = FALSE, row.names = FALSE, quote = FALSE, col.names = FALSE)

## Get overlaps when single inversion overlaps two different protein coding genes on either side of the inversion ##
####################################################################################################################
broken.genes.prot.uniq <- broken.genes.simpINV.prot[!duplicated(broken.genes.simpINV.prot$gene.name)]
hits.starts <- findOverlaps(SimpleInversion.starts, broken.genes.prot.uniq)
hits.ends <- findOverlaps(SimpleInversion.ends, broken.genes.prot.uniq)
## Add gene names
SimpleInversion.starts$break.gene <- ''
gene.names <- split(broken.genes.prot.uniq$gene.name[subjectHits(hits.starts)], queryHits(hits.starts))
gene.names <- sapply(gene.names, function(x) paste(x, collapse = ", "))
SimpleInversion.starts$break.gene[unique(queryHits(hits.starts))] <- gene.names

SimpleInversion.ends$break.gene <- ''
gene.names <- split(broken.genes.prot.uniq$gene.name[subjectHits(hits.ends)], queryHits(hits.ends))
gene.names <- sapply(gene.names, function(x) paste(x, collapse = ", "))
SimpleInversion.ends$break.gene[unique(queryHits(hits.ends))] <- gene.names
## Select inversion where start and end overlap with protein coding gene
hits.both <- intersect(queryHits(hits.starts), queryHits(hits.ends))
fusion.starts <- SimpleInversion.starts[hits.both]
fusion.ends <- SimpleInversion.ends[hits.both]
## Filter out inversion that are within a single gene
mask <- fusion.starts$break.gene != fusion.ends$break.gene
fusion.starts <- fusion.starts[mask]
fusion.ends <- fusion.ends[mask]

## Put together putative fusion genes
putative.fusions <- data.frame(INV.id=fusion.starts$INV.id, 
                               gen=fusion.starts$gen, 
                               ID=fusion.starts$ID,
                               left.break.gene=fusion.starts$break.gene,
                               right.break.gene=fusion.ends$break.gene)
## Export list of disrupted protein coding genes at both breakpoints (Fusions?)
destination <- file.path(outputDirectory, "disrupted_proteinCoding_genes_bothBreakpoints.txt")
write.table(x = putative.fusions, file = destination, append = FALSE, row.names = FALSE, quote = FALSE)

## Plot total number of breakpoints disrupting a protein coding gene
# plt.df <- data.frame(simpINV.breaks=length(SimpleInversion.breakpoints),
#                      invDup.breaks=length(invDup.breakpoints),
#                      simpINV.break.genes=length(unique(broken.genes.simpINV$idx)),
#                      invDup.break.genes=length(unique(broken.genes.invDup$idx)))
## Fraction of breaks disrupting a protein-coding gene
plt.df <- data.frame(simpINV.breaks = length(unique(broken.genes.simpINV$idx)) / length(SimpleInversion.breakpoints),
                     invDup.breaks = length(unique(broken.genes.invDup$idx)) / length(invDup.breakpoints)
)

suppressMessages( plt.df <- reshape2::melt(plt.df) )
plt.df$value <- round(plt.df$value, digits = 2)
plt1 <- ggplot(plt.df) + 
  geom_col(aes(x=variable, y=value, fill=variable)) +
  geom_text(aes(x=variable, y=value, label=value), vjust=-0.1) +
  scale_fill_manual(values = c('darkslateblue', 'indianred3', 'indianred1'), guide='none') +
  xlab("") +
  ylab("Fraction of inversion breakpoints likely\ndisrupting protein-coding genes") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Plot number of protein coding genes disrupted by simple inversion and inverted duplications
plt.df <- data.frame(simpINV.genes=length(broken.genes.simpINV.prot.toPrint), 
                    invDup.genes=length(broken.genes.invDup.prot.toPrint))
suppressMessages( plt.df <- reshape2::melt(plt.df) )
plt2 <- ggplot(plt.df) + 
  geom_col(aes(x=variable, y=value, fill=variable)) +
  geom_text(aes(x=variable, y=value, label=value), vjust=-0.1) +
  scale_fill_manual(values = c('darkslateblue', 'indianred3', 'indianred1'), guide='none') +
  xlab("") +
  ylab("Number of protein-coding genes likely\ndisrupted by inversion breakpoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Plot number of disrupted proteing coding genes per individual
# NHP.IDs <- unlist(sapply(broken.genes.prot$ID, function(x) strsplit(x, ",")[[1]]))
# plt.df <- as.data.frame(NHP.IDs)
# plt2 <- plt.df %>% group_by(NHP.IDs) %>% summarise(count=n()) %>% 
#   ggplot() + 
#   geom_col(aes(x=NHP.IDs, y=count, fill=NHP.IDs)) +
#   geom_text(aes(x=NHP.IDs, y=count, label=count), vjust=-0.1) +
#   scale_fill_manual(values = c('#3182bd','#31a354','#8856a7','#e6550d'), name="") +
#   xlab("") +
#   ggtitle("Disrupted protein-coding genes\nper individual") +
#   theme_bw()

## List of genes disrupted by simple inversions in all individuals
# NHP.IDs <- sapply(broken.genes.simpINV.prot$ID, function(x) strsplit(x, ",")[[1]])
# filt <- which(lengths(NHP.IDs) == 4)
# genes.disrupted.allIndivd <- broken.genes.prot[filt]
# genes.disrupted.allIndivd <- genes.disrupted.allIndivd[!duplicated(genes.disrupted.allIndivd $gene.name)]
# plt.df <- as.data.frame(genes.disrupted.allIndivd)
# plt3 <- ggplot(plt.df, aes(x=1, y=1:nrow(plt.df))) + 
#   geom_point(shape=32) +
#   geom_text(aes(label=gene.name), hjust=1) +
#   scale_x_continuous(limits = c(0.9, 1.1)) +
#   xlab("") +
#   ylab("Genes") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

## Plot gene categories completely inside invDups
plt.df <- as.data.frame(gene.type.counts)
plt.df$ID <- names(gene.type.counts)
plt.df <- plt.df %>%
  arrange(desc(ID)) %>%
  mutate(lab.ypos = cumsum(gene.type.counts) - 0.5 * gene.type.counts)

mycols <- brewer.pal(name = 'Dark2', n = 8)
plt4 <- ggplot(plt.df, aes(x = "", y = gene.type.counts, fill = ID)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = gene.type.counts), color = "white")+
  scale_fill_manual(values = mycols) +
  ggtitle("Gene categories duplicated by\nInverted duplications") +
  theme_void()
## Compile final plot and export to pdf
final.plt <- plot_grid(plt1, plt2, plt4, nrow = 1, rel_widths = c(0.5,0.5,1))
destination <- file.path(outputDirectory, "disruptedGenes_summary.pdf")
ggsave(filename = destination, plot = final.plt, width = 8, height = 4, useDingbats=FALSE)

message("DONE!!!")

#========================================================================================
# Previous code

# # List of disrupted genes
# gene.list <- as.data.frame(broken.genes[!duplicated(broken.genes$gene.name)])
# gene.list <- gene.list %>% select(gene.id, gene.name)
# gene.list$gene.id <- gsub(gene.list$gene.id, pattern = "\\.\\d+$", replacement = "") #remove suffix
# names(gene.list) <- c('ENSEMBL_GENE_ID', 'Name')
# gene.list$Species <- 'Homo Sapiens'
# write.table(gene.list, file = "/home/porubsky/WORK/Great_apes/Disrupted_genes/disrupted_genes_list.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# ## Analyze gene ontology output from DAVID ##
# #############################################
# up.tissue <- read.table("/home/porubsky/WORK/Great_apes/Disrupted_genes/DAVID_GO_uptissue.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# up.tissue <- up.tissue[order(up.tissue$Count),]
# up.tissue$Term <- factor(up.tissue$Term, levels=up.tissue$Term)
# 
# plt1 <- ggplot(up.tissue) +
#   geom_col(aes(x=Term, y=Count), fill='gray64') +
#   geom_text(aes(x=Term, y=0, label=round(PValue, digits = 5)), color="chartreuse4", hjust=0) +
#   geom_text(aes(x=Term, y=0, label=round(Bonferroni, digits = 5)), color="brown2", hjust=-2) + 
#   geom_text(aes(x=Term, y=max(Count), label=Count), hjust=0) +   
#   coord_flip() +
#   theme_bw() +
#   xlab("GO term (up.tissue)")
# 
# 
# up.tissue.detailed <- read.table("/home/porubsky/WORK/Great_apes/Disrupted_genes/DAVID_GO_UNIGENE-EST-QUARTILE_CGAP-SAGE-QUARTILE.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# up.tissue.detailed$Term <- gsub(up.tissue.detailed$Term, pattern = "_3rd", replacement = "")
# up.tissue.detailed <- up.tissue.detailed[order(up.tissue.detailed$Count),]
# up.tissue.detailed$Term <- factor(up.tissue.detailed$Term, levels=up.tissue.detailed$Term)
# 
# plt2 <- ggplot(up.tissue.detailed) +
#   geom_col(aes(x=Term, y=Count), fill='gray64') +
#   coord_flip() +
#   xlab("GO term (up.tissue)") +
#   theme_bw()
# 
# plot_grid(plt1, plt2, nrow = 1)
# 
# ## Get gene names ##
# brain.genes <- up.tissue$Genes[4]
# brain.genes <- unlist(strsplit(brain.genes, ", "))
# 
# gencode29.gr$gene.id.base <- gsub(gencode29.gr$gene.id, pattern = "\\.\\d+", replacement = "")
# brain.genes.names <- gencode29.gr[gencode29.gr$gene.id.base %in% brain.genes]
# brain.genes.names <- as.character(brain.genes.names$gene.name)
# write.table(as.data.frame(brain.genes.names), file = "/home/porubsky/WORK/Great_apes/Disrupted_genes/brain_enriched_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# ## Analyze inversions that do not overlap with any gene ##
# ##########################################################
# simpleInversion.calls.noGene <- subsetByOverlaps(all.SimpleInversion.calls.filt, vega68PLUSgencode29.gr, invert = TRUE)
# simpleInversion.calls.noGene.INVregions <- reduce(simpleInversion.calls.noGene) # inverted regions with no gene overlap
# 
# ## Count overlaps with known enhancers
# simpleInversion.calls.enhancer <- subsetByOverlaps(simpleInversion.calls.noGene.INVregions, geneHancer.gr)
# simpleInversion.calls.enhancerElite <- subsetByOverlaps(simpleInversion.calls.noGene.INVregions, geneHancer.gr[geneHancer.gr$categ == 'Elite'])
# 
# total.INVregions <- length(simpleInversion.calls.noGene.INVregions)
# enhancerElite <- length(simpleInversion.calls.enhancerElite)
# enhancer <- length(simpleInversion.calls.enhancer) - enhancerElite
# 
# enhancer.polymorhic <- sum(countOverlaps(simpleInversion.calls.enhancer, polymorphic.inversions))
# enhancer.HSinv <- sum(countOverlaps(simpleInversion.calls.enhancer, HSinv.gr))
# 
# ## Prepare data for plotting
# plt.df <- data.frame(total.regions=total.INVregions, enhancerElite=enhancerElite, enhancer=enhancer, noOverlap=total.INVregions-(enhancerElite+enhancer))
# plt.df <- melt(plt.df)
# plt.df$perc <- plt.df$value / plt.df$value[1]
# plt.df <- plt.df[-1,]
# 
# ## Prepare pie chart plot
# blank_theme <- theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     plot.title=element_text(size=14, face="bold")
#   )
# 
# plt <- ggplot(plt.df, aes(x="", y=value, fill=variable)) +
#   geom_bar(width = 1, stat = "identity") +
#   coord_polar("y", start=0) +
#   scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999")) +
#   blank_theme +
#   theme(axis.text.x=element_blank())
# 
# #plt +  geom_text(aes(y = cumsum(value) - value/2, label = value, size=5))
# 
# ## Measure distance of simple inversions to human genes ##
# ##########################################################
# ## Load functional elements annotation
# header <- read.table(gzfile("/home/porubsky/WORK/Great_apes/Annotations/RefSeq_Func_Elems_GRCh38.bed.gz"), stringsAsFactors = FALSE, sep = "\t", nrows = 1, comment.char = '&')
# header <- gsub(header, pattern = "#", replacement = "")
# funcElems.annot <- read.table(gzfile("/home/porubsky/WORK/Great_apes/Annotations/RefSeq_Func_Elems_GRCh38.bed.gz"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", skip = 1, fill = TRUE)
# colnames(funcElems.annot) <- header
# funcElems.annot.gr <- GRanges(seqnames=funcElems.annot$chrom, ranges=IRanges(start=funcElems.annot$chromStart, end=funcElems.annot$chromEnd), funcElem.type=funcElems.annot$name)
# funcElems.annot.gr <- keepStandardChromosomes(funcElems.annot.gr, pruning.mode = 'coarse')
# ## Load H3K4me1 annotation
# H3K4me1 <- data.table::fread("/home/porubsky/WORK/Great_apes/Annotations/wgEncodeBroadHistoneGm12878H3k4me1StdSig.bedgraph", header = FALSE)
# H3K4me1.gr <- GRanges(seqnames=H3K4me1$V1, ranges=IRanges(start=H3K4me1$V2, end=H3K4me1$V3), mcols=H3K4me1$V4)
# ## Load ORegAnno annotation
# ORegAnno <- read.table("/home/porubsky/WORK/Great_apes/Annotations/ORegAnno_GRCh38.bed.gz")
# ORegAnno.gr <- GRanges(seqnames=ORegAnno$V2, ranges=IRanges(start=ORegAnno$V3, end=ORegAnno$V4))
# ORegAnno.gr <- reduce(ORegAnno.gr)
# 
# # Get only inversions that do not overlaps with any gene
# SimpleInversion.noGeneOverlap <- subsetByOverlaps(all.SimpleInversion.calls.filt, vega68PLUSgencode29.gr, invert = TRUE)
# # Find distance of these inversion to the closest genes
# SimpleInversion.noGeneOverlap.distsGenes <- range2rangeDistance(gr=SimpleInversion.noGeneOverlap, userTrack=vega68PLUSgencode29.gr, allow.overlap = FALSE)
# # Select only inversion that are not further than 10kb away from a closest gene
# gene.neighbour <- which(SimpleInversion.noGeneOverlap.distsGenes$leftDist < 10000 | SimpleInversion.noGeneOverlap.distsGenes$rightDist < 10000)
# invClose.genes <- SimpleInversion.noGeneOverlap[gene.neighbour]
# 
# test.regul <- permTest(
#   A=invClose.genes, 
#   B=ORegAnno.gr, 
#   randomize.function=randomizeRegions, 
#   evaluate.function=numOverlaps, 
#   genome="hg38", 
#   ntimes=1000, 
#   allow.overlaps=FALSE, 
#   per.chromosome=TRUE, 
#   mask=hg38.mask, 
#   count.once=TRUE,
#   mc.set.seed=FALSE,
#   mc.cores=4)
# 
# 
# test.regul <- permTest(
#   A=invClose.genes, 
#   randomize.function=randomizeRegions, 
#   evaluate.function=meanInRegions,
#   x=H3K4me1.gr,
#   genome="hg38", 
#   ntimes=1000, 
#   allow.overlaps=FALSE, 
#   per.chromosome=TRUE, 
#   mask=hg38.mask, 
#   mc.set.seed=FALSE,
#   mc.cores=4)