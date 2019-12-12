## This file contains an order in which individual R scripts should be run ##
#############################################################################

## Directory with a set of scripts to run
script.dir <- "/home/porubsky/WORK/Great_apes/scripts"

## Calculate some basic statistics for selected bam files and make a plot [This script needs to be RUN only once]
source(file.path(script.dir, "getBamStatSummary.R"))

## Read in and export all published calls and calls for PacBio and Illumina data [These scripts need to be RUN only once]
source(file.path(script.dir, "readPublishedInversionCalls.R"))
source(file.path(script.dir, "readBioNanoInversionCalls.R"))
source(file.path(script.dir, "readPacBioAndIlluminaInversionCalls.R"))
## Get position of enriched genes predicted by Alex Pollen
source(file.path(script.dir, "exportPositionsOfEnrichedGenes.R"))

## Preparing and plotting composite files [These scripts need to be RUN only once]
source(file.path(script.dir, "exportANDcompileCompositeFiles.R"))
## Make circular composite file plot
source("/home/porubsky/WORK/Great_apes/Plotting_scripts/plotCompositeFilesCircular.R")
## Plot examples of various SV classes from orangutan composite file
source("/home/porubsky/WORK/Great_apes/Plotting_scripts/plotInvClassExamples.R")

## [Re-run below scripts in case of changes in INV calls !!!]
## Export manually checked inversion calls
source(file.path(script.dir, "exportManuallyProcessedInversionCalls.R"))

## Remove uncertain (small) INV calls that cannot be validated based on 50% reciprocal overlap with PacBio or Illumina calls
source(file.path(script.dir, "removeUncetainInversionCalls.R"))

## Merge valid uncertain calls with certain inversion calls and keep only calls with at least 1kb of unique sequence in total
## and inverted duplications that are at least 10kb in size
source(file.path(script.dir, "mergeINVcallsAndFilter.R"))
## Summarize all INV calls in pie charts
source("/home/porubsky/WORK/Great_apes/Plotting_scripts/plotComplexPie.R")

## Annotate filtered calls by predicted copy number based on Illumina data (WSSD) and mark predicted misorients
#TODO

## Export inverted duplication calls
source(file.path(script.dir, "exportInvertedDuplicationBreakpoints.R"))

## Export final Master tables by comparing final Strand-seq calls with all the other available calls
source(file.path(script.dir, "exportMasterTables.R"))

## Prepare dotplots against de-novo assemblies as an additional layer of validation
## NOTE: Dotplots have to be made using a snakemake pipeline in a path stated below
## /net/eichler/vol27/projects/strandseq_greatape/nobackups/Ape_assemblies/validateINVbyAssembly
source(file.path(script.dir, "analyzeDotplots.R"))

## Analyze master tables and report proportions of valid vs unvalid inversion calls based on 50% reciprocal overlap
## This step requires manually curated master tables with additional info on dotplot and BioNano validation
source(file.path(script.dir, "analyzeValidatedInversionCalls.R"))

## Analyze putative human specific inversions
## Compile human inversion callset for NA19240 plus include human specific inversions
source(file.path(script.dir, "analyzeHumanSpecificEvents.R"))
## Plot example of human specific inversion
source("/home/porubsky/WORK/Great_apes/Plotting_scripts/plotHSinvExamples.R")

## Analyze polymorphic inversions
source(file.path(script.dir, "analyzePolymorphicEvents.R"))

## Analyze breakpoint hotspots
source(file.path(script.dir, "analyzeBreakpointHotspots.R"))

## Analyze evolutionary distances
source(file.path(script.dir, "analyzeSimpleInvEvolDistance.R"))

## Analyze gene enrichment
source(file.path(script.dir, "analyzeGeneEnrichment.R"))

## Analyze genes disrupted by inversions [Deprecated]
source(file.path(script.dir, "analyzeDisruptedGenes.R"))

## Analyze intra- vs inter-chromosomal links between inverted duplications
## Split-read mapping of PacBio reads has to be done beforehand!!!
source(file.path(script.dir, "analyzeInvertedDuplicationLinks.R"))

## Analyze inversions on chromosomes X
source(file.path(script.dir, "analyzeChromosomeX.R"))

## Run this code line-by-line to get some basic data summary
file.path(script.dir, "getInversionCallSummary.R")






