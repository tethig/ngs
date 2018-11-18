##### RNASeq Workflow ######

## Understanding file formats ##
# for Fastq:
https://en.wikipedia.org/wiki/FASTQ_format

# for SAM
# format specification:
https://samtools.github.io/hts-specs/
https://samtools.github.io/hts-specs/SAMv1.pdf
# explain those flags:
http://broadinstitute.github.io/picard/explain-flags.html
# optional flags in Bowtie output:
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

## Notes on software ##
# popular softwares for mapping to a genome are bowtie2 or BWA:
https://github.com/BenLangmead/bowtie2
http://bio-bwa.sourceforge.net

# TopHat uses Bowtie to determine splice boundaries:
http://ccb.jhu.edu/software/tophat/index.shtml

# Cufflinks can assemble a transcriptome and carry out differential expression analysis (script carries out stats)
https://cole-trapnell-lab.github.io/cufflinks/

# Galaxy makes life easy
https://usegalaxy.org/

# Galaxy pipeline from Tejaswini Mishra, Penn State University (2012) shows TopHat and Cufflinks in action
https://usegalaxy.org/u/tejas/p/rnaseq-workshop

# SAMtools is the always useful toolset
# install via homebrew on the mac
# also recommend fastx_toolkit, seqtk and FastQC
# fastq_detect.pl script from:
http://wiki.bits.vib.be/index.php/Identify_the_Phred_scale_of_quality_scores_used_in_fastQ

# countxpression python script for processing count data from SAM files
https://github.com/DeWitP/SFG/tree/master/scripts
## error in this script!!

## Differential expression analysis via DESeq ##
# Count data sets available here here:
http://bowtie-bio.sourceforge.net/recount/
# Use Phenotype table to deduce sample IDs

# using gilad data set in this example

### acknowledgements ###
# Dannise Ruiz, Penn State University

# Tejaswini Mishra, Penn State University Guidance (2012)

# Palumbi Lab's simple fool's guide
De Wit P, Pespeni MH, Ladner JT, Barshis DJ, Seneca F, Jaris H, Overgaard Therkildsen N, Morikawa M and Palumbi SR (2012) The simple fools guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis.  Molecular Ecology Resources 12, 1058-1067.
http://sfg.stanford.edu/expression.html

# heatmap
Anders, S. and Huber, W., 2012. Differential expression of RNA-Seq data at the gene level–the DESeq package. Heidelberg, Germany: European Molecular Biology Laboratory (EMBL).
http://www.genomatix.de/online_help/help_regionminer/DESeq_1.10.1.pdf

## DESeq workflow ##
# set working directory
setwd("/Users/bio3dickib/Desktop/working2")

# go get the file and view
ReadCounts <- read.table("gilad_count_table.txt",header=TRUE,row.names=1)
ReadCounts[1:2,]
dim(ReadCounts)

# gotta have 10 or more reads in total
big10 <- apply(ReadCounts,1,sum) >= 10
TotalReads <- ReadCounts[big10,]
dim(TotalReads)

# get DESeq from bioconductor
# may need to follow instructions here:
http://www.bioconductor.org/install/
library(DESeq) # note there is now a DESeq2

# let's label treatments (must be co-linear)
head(TotalReads)
treatments=rep(c("F","M"), each=3)

# create (and inspect) dataset
cds=newCountDataSet(TotalReads,treatments)
head(counts(cds))

# normalising by overall number of mapped reads per sample
cds=estimateSizeFactors(cds)

# showing the adjustment factors
sizeFactors(cds)

# view of normalised counts
round( head( counts( cds, normalized=TRUE ) ) )

#over-dispersion of read counts
cds=estimateDispersions(cds)

# plot
plotDispEsts<- function( cds ){
plot(
rowMeans(counts(cds,normalized=TRUE)),
fitInfo(cds)$perGeneDispEsts,
pch = '.', log="xy" )
xg <- 10^seq( -.5, 5, length.out=300 )
lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
plotDispEsts(cds)


# differential expression test!!!
result=nbinomTest(cds,"F","M")
head(result)

# does fold-change depend on read depth?
plotDE <- function( res ){
plot(
res$baseMean,
res$log2FoldChange,
log="x", pch='.', cex=.3,
col = ifelse( res$padj < .1, "red", "black" ) )
}
plotDE( result )

# distribution of p values
hist(result$pval, breaks=100, col="skyblue", border="slateblue", main="",xlab="p-values")

# pull out significant results
resultSig = result[ result$pval < 0.1, ] # maybe change to pval if vector of padj too weak
dim(resultSig)
head( resultSig[ order(resultSig$pval), ] ) # order by signficance

# top upregulated
head(resultSig[order(resultSig$foldChange,-resultSig$baseMean),])

# top downreg
head(resultSig[order(-resultSig$foldChange,- resultSig$baseMean),])

# write to file
write.csv( result, "results.csv", row.names=TRUE )

# variance stablising transform
vsdFull = varianceStabilizingTransformation( cds )

# heatmap2
library("RColorBrewer")
library("gplots")
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
