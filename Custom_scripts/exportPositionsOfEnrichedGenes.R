## Load required libraries
suppressPackageStartupMessages( library(primatR) )

message("Searching for gene positions in gene annottaion ...")

## Load gene table (from Alex Pollen)
gene.tab <- read.table("/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/Human_Chimp_genelists.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene.tab.list <- split(gene.tab, gene.tab$ID)

## Load gene annotations
vega68 <- read.table(gzfile("/home/porubsky/WORK/Great_apes/Annotations/vega68_GRCh38_geneList_mart_export.txt.gz"), header=TRUE, sep=",", stringsAsFactors = FALSE)
gencode29 <- read.table("/home/porubsky/WORK/Great_apes/Annotations/gencode.v29.annotation..processed.txt", header = TRUE)

outputDirectory <- "/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/"

## Get gene positions from vega68 & gencode29 databases ##
##########################################################
unannot.gene.list <- list()
gene.list <- list()
for (i in seq_along(gene.tab.list)) {
  ID <- names(gene.tab.list[i])
  message("Working on ", ID, " gene list ...")
  sub.tab <- gene.tab.list[[i]]
  
  ## Search gene positions in gencode
  gencode29.select <- gencode29[gencode29$gene_name %in% sub.tab$Gene,]
  missing.genes <- !sub.tab$Gene %in% gencode29$gene_name
  vega68.select <- vega68[vega68$Gene.name %in% sub.tab$Gene[missing.genes],]
  
  unannot.genes <- sub.tab$Gene[missing.genes][!sub.tab$Gene[missing.genes] %in% vega68$Gene.name] 
  message("    ", length(unannot.genes), " genes couldn't be annotated!!! Exporting into file ...")
  unannot.genes.df <- data.frame(gene.name = unannot.genes, ID = ID)
  
  ## Transform gene list into GRanges
  gencode29.select.gr <- GRanges(seqnames=gencode29.select$chr, ranges=IRanges(start=gencode29.select$start, gencode29.select$end), gene.name=gencode29.select$gene_name, gene.type=gencode29.select$gene_type)
  vega68.select.gr <- GRanges(seqnames=paste0('chr', vega68.select$Chromosome.scaffold.name), ranges=IRanges(start=vega68.select$Gene.start..bp., vega68.select$Gene.end..bp.), gene.name=vega68.select$Gene.name, gene.type=vega68.select$Gene.type)
  suppressWarnings( all.genes.pos <- c(gencode29.select.gr, vega68.select.gr) )
  # Keep only standard chromosomes
  all.genes.pos <- keepStandardChromosomes(all.genes.pos, pruning.mode = 'coarse')
  all.genes.pos <- sort(all.genes.pos)
  # For some duplicated genes remove the duplicate !!! [CAREFUL]
  all.genes.pos <- all.genes.pos[!duplicated(all.genes.pos$gene.name)]
  
  # Add info about info if this gene is upregulated in chimp or human
  all.genes.pos$upregul.in <- sub.tab$Up.In.Species[match(all.genes.pos$gene.name, sub.tab$Gene)]
  
  ## Export GeneSets into UCSC formatted bed files
  ranges2UCSC(all.genes.pos[all.genes.pos$upregul.in == 'human'], outputDirectory = outputDirectory, index = paste0(ID,"_genes_human"), colorRGB = '255,0,0', id.field = 'gene.name')
  ranges2UCSC(all.genes.pos[all.genes.pos$upregul.in == 'chimp'], outputDirectory = outputDirectory, index = paste0(ID,"_genes_chimp"), colorRGB = '0,255,0', id.field = 'gene.name')
  
  ## Split by ID
  observed.genes <- split(all.genes.pos, all.genes.pos$upregul.in)
  
  ## Store position of observed genes
  gene.list[[ID]] <- observed.genes
  ## Store missing genes
  unannot.gene.list[[ID]] <- unannot.genes.df
}  
## Export list of genes that haven't been found
unannot.genes <- do.call(rbind,  unannot.gene.list)
destination <- file.path(outputDirectory, "unannot.gene.list.txt")
write.table(unannot.genes, file = destination, quote = FALSE, row.names = FALSE)
## Export observed gene.list in RData object
destination <- file.path(outputDirectory, 'upregul.gene.list.RData')
save(gene.list, file = destination)

## Fill in missing genes ##
###########################
## Here missing genes have been looked up manually
manualy.annot.genes <- read.table("/home/porubsky/WORK/Great_apes/Gene_enrichment_analysis/unannot.gene.list.manually.filled.csv", sep = ",", header = TRUE)
gene.list <- list()
for (i in seq_along(gene.tab.list)) {
  ID <- names(gene.tab.list[i])
  message("Working on ", ID, " gene list ...")
  sub.tab <- gene.tab.list[[i]]
  manualy.annot.genes.sub <- manualy.annot.genes[manualy.annot.genes$ID == ID,]
  manualy.annot.genes.sub <- manualy.annot.genes.sub[!is.na(manualy.annot.genes.sub$chr),]
  
  ## Search gene positions in gencode
  gencode29.select <- gencode29[gencode29$gene_name %in% sub.tab$Gene,]
  missing.genes <- !sub.tab$Gene %in% gencode29$gene_name
  vega68.select <- vega68[vega68$Gene.name %in% sub.tab$Gene[missing.genes],]
  
  ## Transform gene list into GRanges
  gencode29.select.gr <- GRanges(seqnames=gencode29.select$chr, ranges=IRanges(start=gencode29.select$start, gencode29.select$end), gene.name=gencode29.select$gene_name, gene.type=gencode29.select$gene_type)
  vega68.select.gr <- GRanges(seqnames=paste0('chr', vega68.select$Chromosome.scaffold.name), ranges=IRanges(start=vega68.select$Gene.start..bp., vega68.select$Gene.end..bp.), gene.name=vega68.select$Gene.name, gene.type=vega68.select$Gene.type)
  manualy.annot.genes.sub.gr <- GRanges(seqnames=manualy.annot.genes.sub$chr, ranges=IRanges(start=manualy.annot.genes.sub$start, end=manualy.annot.genes.sub$end), gene.name=manualy.annot.genes.sub$gene.name, gene.type='')
  suppressWarnings( all.genes.pos <- c(gencode29.select.gr, vega68.select.gr, manualy.annot.genes.sub.gr) )
  # Keep only standard chromosomes
  all.genes.pos <- keepStandardChromosomes(all.genes.pos, pruning.mode = 'coarse')
  all.genes.pos <- sort(all.genes.pos)
  # For some duplicated genes remove the duplicate !!! [CAREFUL]
  all.genes.pos <- all.genes.pos[!duplicated(all.genes.pos$gene.name)]
  
  # Add info about info if this gene is upregulated in chimp or human
  all.genes.pos$upregul.in <- sub.tab$Up.In.Species[match(all.genes.pos$gene.name, sub.tab$Gene)]
  
  ## Export GeneSets into UCSC formatted bed files
  ranges2UCSC(all.genes.pos[all.genes.pos$upregul.in == 'human'], outputDirectory = outputDirectory, index = paste0(ID,"_genes_human"), colorRGB = '255,0,0', id.field = 'gene.name')
  ranges2UCSC(all.genes.pos[all.genes.pos$upregul.in == 'chimp'], outputDirectory = outputDirectory, index = paste0(ID,"_genes_chimp"), colorRGB = '0,255,0', id.field = 'gene.name')
  
  ## Split by ID
  observed.genes <- split(all.genes.pos, all.genes.pos$upregul.in)
  
  ## Store position of observed genes
  gene.list[[ID]] <- observed.genes
} 
## Export observed gene.list in RData object
destination <- file.path(outputDirectory, 'upregul.gene.list.RData')
save(gene.list, file = destination)

message("DONE!!!")