## Script to compare haplotypes assembled using Strand-seq ##
#############################################################

## Load required libraries
suppressPackageStartupMessages( library(ape) )
suppressPackageStartupMessages( library(StrandPhaseR) )
suppressPackageStartupMessages( library(VariantAnnotation) )

## Load phased data for NHP
chimpanzee.haps <- vcf2vranges(vcfFile = "/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/StrandPhaseR_analysis/VCFfiles/chimpanzee.strandphaser.vcf", 
                               genoField = 'GT', 
                               translateBases = TRUE, 
                               genome = 'hg38')

bonobo.haps <- vcf2vranges(vcfFile = "/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/StrandPhaseR_analysis/VCFfiles/bonobo.strandphaser.vcf", 
                               genoField = 'GT', 
                               translateBases = TRUE, 
                               genome = 'hg38')

gorilla.haps <- vcf2vranges(vcfFile = "/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/StrandPhaseR_analysis/VCFfiles/gorilla.strandphaser.vcf", 
                               genoField = 'GT', 
                               translateBases = TRUE, 
                               genome = 'hg38')

orangutan.haps <- vcf2vranges(vcfFile = "/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/StrandPhaseR_analysis/VCFfiles/orangutan.strandphaser.vcf", 
                               genoField = 'GT', 
                               translateBases = TRUE, 
                               genome = 'hg38')

chromosomes <- paste0('chr', c(1:22))
hap.matrix <- list()
for (chr in chromosomes) {
  ## Extract chromosome hapltoypes
  chimpanzee.chr <- chimpanzee.haps[seqnames(chimpanzee.haps) == chr]
  bonobo.chr <- bonobo.haps[seqnames(bonobo.haps) == chr]
  gorilla.chr <- gorilla.haps[seqnames(gorilla.haps) == chr]
  orangutan.chr <- orangutan.haps[seqnames(orangutan.haps) == chr]
  ## Remove duplicated records
  chimpanzee.chr <- chimpanzee.chr[!duplicated(start(chimpanzee.chr))]
  bonobo.chr <- bonobo.chr[!duplicated(start(bonobo.chr))]
  gorilla.chr <- gorilla.chr[!duplicated(start(gorilla.chr))]
  orangutan.chr <- orangutan.chr[!duplicated(start(orangutan.chr))]
  
  snv.pos <- unique(c(start(chimpanzee.chr), start(bonobo.chr), start(gorilla.chr), start(orangutan.chr)))
  uncov <- setdiff(snv.pos, start(chimpanzee.chr))
  if (length(uncov) > 0) {
    uncov.gr <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=uncov, end=uncov))
    uncov.gr$H1.chimpanzee <- 'N'
    uncov.gr$H2.chimpanzee <- 'N'
    chimpanzee.chr <- as(chimpanzee.chr, 'GRanges')
    colnames(mcols(chimpanzee.chr)) <- c('H1.chimpanzee', 'H2.chimpanzee')
    chimpanzee.chr <- sort(c(chimpanzee.chr, uncov.gr))
  }
  
  uncov <- setdiff(snv.pos, start(bonobo.chr))
  if (length(uncov) > 0) {
    uncov.gr <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=uncov, end=uncov))
    uncov.gr$H1.bonobo <- 'N'
    uncov.gr$H2.bonobo <- 'N'
    bonobo.chr <- as(bonobo.chr, 'GRanges')
    colnames(mcols(bonobo.chr)) <- c('H1.bonobo', 'H2.bonobo')
    bonobo.chr <- sort(c(bonobo.chr, uncov.gr))
  }
  
  uncov <- setdiff(snv.pos, start(gorilla.chr))
  if (length(uncov) > 0) {
    uncov.gr <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=uncov, end=uncov))
    uncov.gr$H1.gorilla <- 'N'
    uncov.gr$H2.gorilla <- 'N'
    gorilla.chr <- as(gorilla.chr, 'GRanges')
    colnames(mcols(gorilla.chr)) <- c('H1.gorilla', 'H2.gorilla')
    gorilla.chr <- sort(c(gorilla.chr, uncov.gr))
  }
  
  uncov <- setdiff(snv.pos, start(orangutan.chr))
  if (length(uncov) > 0) {
    uncov.gr <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=uncov, end=uncov))
    uncov.gr$H1.orangutan <- 'N'
    uncov.gr$H2.orangutan <- 'N'
    orangutan.chr <- as(orangutan.chr, 'GRanges')
    colnames(mcols(orangutan.chr)) <- c('H1.orangutan', 'H2.orangutan')
    orangutan.chr <- sort(c(orangutan.chr, uncov.gr))
  }
  
  mat <- cbind(as.matrix(mcols(chimpanzee.chr)), as.matrix(mcols(bonobo.chr)), as.matrix(mcols(gorilla.chr)), as.matrix(mcols(orangutan.chr)))
  hap.matrix[[chr]] <- mat
}  
mat <- do.call(rbind, hap.matrix)
## Recode nucleotides into the numbers
mat <- apply(mat, 2, function(x) recode(as.character(x), 'A'=1, 'C'=2, 'G'=3, 'T'=4, 'N'=0, .default=0))
## Keep only site with 4 alleles
mask <- apply(mat, 1, function(x) length(x[x == 0]) <= 4)
mat <- mat[mask,]

## Prepare evolutionary tree
mat.t <- t(mat)
mat.nj <- nj(dist.gene(mat.t))
tree <- as.phylo(mat.nj)
plot(tree)

rec.mat.t <- apply(mat.t, 2, recode.fun)

## Get only position represented in at least 4 haplotypes
mask <- apply(rec.mat.t, 2, function(x) length(x[x>0]) >= 3)

rec.mat.plt <- melt(rec.mat.t[,mask]) #wide to long conversion for plotting

#colors <- c("white", brewer.pal(n=9, name="Set1"))
colors <- c("white", brewer.pal(n=9, name="Set1"))
chunks <- 100 #set number of breaks on X axis
gen.pos <- start(chimpanzee.chr1)[mask]
breaks <- floor(length(gen.pos)/chunks)
breaks <- c(0, breaks*1:chunks)
labels <- c(0, gen.pos[breaks]) #get true genomic postions as labels

#breaks <- 1:length(gen.pos[mask])
#labels <- gen.pos

#plot PSVs alleles along all sequences in MSA
plt <- ggplot(rec.mat.plt) + 
  geom_tile(aes(x = Var2, y = Var1, fill = as.factor(value))) + 
  scale_fill_manual(values = colors, name="Alleles (PSVs)") + 
  scale_x_continuous(breaks = breaks, labels = labels) +
  xlab("Genomic Position") + 
  ylab("") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))


#which(gen.pos %in% 142770740)
#which(gen.pos %in% 143981273)


recode.fun <- function(x) {
  rec.x <- x
  
  code1 <- which(x == 1)
  if (length(code1) > 0) {
    rec.x[code1] <- min(code1)
  }
  
  code2 <- which(x == 2)
  if (length(code2) > 0) {
    rec.x[code2] <- min(code2)
  }
  
  code3 <- which(x == 3)
  if (length(code3) > 0) {
    rec.x[code3] <- min(code3)
  }
  
  code4 <- which(x == 4)
  if (length(code4) > 0) {
    rec.x[code4] <- min(code4)
  }
  
  code5 <- which(x == 0)
  if (length(code5) > 0) {
    rec.x[code5] <- 0
  }
  return(rec.x)
}


## Set parameters
chromosomes <- paste0('chr', c(1:22, 'X'))
phased1 <- chimpanzee.haps
phased2 <- bonobo.haps
## Process each chromosome
#all.plots <- list()
all.comparisons <- list()
for (chr in chromosomes) {
  message("Comparing haplotypes for chromosome: ", chr)

  phased1.chr <- phased1[seqnames(phased1) == chr]
  phased2.chr <- phased2[seqnames(phased2) == chr]
  
  ## Get only shared SNVs between datasets
  shared.pos <- findOverlaps(phased1.chr, phased2.chr)
  phased1.chr <- phased1.chr[queryHits(shared.pos)]
  phased2.chr <-  phased2.chr[subjectHits(shared.pos)]
  
  ## Compare haplotypes between datasets
  phased1.chr$H1.2 <- phased2.chr$H1
  phased1.chr$H2.2 <- phased2.chr$H2
  phased1.chr$H1.comp <- 'N'
  phased1.chr$H2.comp <- 'N'
  H1.mask <- which(!is.na(phased2.chr$H1))
  H2.mask <- which(!is.na(phased2.chr$H2))
  phased1.chr$H1.comp[H1.mask] <- ifelse(phased1.chr$H1[H1.mask] == phased1.chr$H1.2[H1.mask], 'H1', 'H2')
  phased1.chr$H2.comp[H2.mask] <- ifelse(phased1.chr$H2[H2.mask] == phased1.chr$H2.2[H2.mask], 'H2', 'H1')
  
  all.comparisons[[length(all.comparisons) + 1]] <- as.data.frame(phased1.chr)
}  


## Plot haplotype comparison
plt.df <- do.call(rbind, all.comparisons)
## Sort by chromosome
plt.df$seqnames <- factor(plt.df$seqnames, levels=chromosomes)

library(scales)
plt <- ggplot() + 
    geom_linerange(data=plt.df, aes(x=start, ymin=0, ymax=1, color=H1.comp)) +
    geom_linerange(data=plt.df, aes(x=start, ymin=1, ymax=2, color=H2.comp)) +
    scale_color_manual(values = c("cadetblue4", "darkgoldenrod3", "white"), name="haps") +
    scale_x_continuous(labels = comma, expand = c(0,0)) +
    facet_grid(seqnames ~ ., switch = 'y') +
    xlab("Genomic position") +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          strip.text.y = element_text(angle = 180))  

polym.inv <- get(load("/home/porubsky/WORK/Great_apes/Polymorphic_inversion_sites/shared_HGSVC&NHP_sites.RData"))
polym.inv.df <- as.data.frame(polym.inv)

plt + geom_point(data=polym.inv.df, aes(x=start, y=2), color='red', shape=118)
