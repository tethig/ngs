### BIOCONDUCTOR
# Website
# https://www.bioconductor.org

# Source the bioconductor script
source("http://bioconductor.org/biocLite.R")
biocLite()


### BALLGOWN
# Download and load Ballgown package
biocLite("ballgown")
library(ballgown)

# Quick run through based on this:
browseVignettes("ballgown")

# Loading data
data_directory = system.file('extdata', package='ballgown') # automatically finds ballgown's installation directory
# examine data_directory:
data_directory
# create object
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')
bg

# containers: structure, expr, indexes (incl. pData), dirs, mergedDate, meas

# accessing containers
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))
exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

# plotting
plotTranscripts(gene='XLOC_000454', gown=bg, samples='sample12',
    meas='FPKM', colorby='transcript',
    main='transcripts from gene XLOC_000454: sample 12, FPKM')

# plotting contrasts
plotTranscripts('XLOC_000454', bg,
    samples=c('sample01', 'sample06', 'sample12', 'sample19'),
    meas='FPKM', colorby='transcript')

plotMeans('XLOC_000454', bg, groupvar='group', meas='FPKM', colorby='transcript')

# statistical tasks
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
View(stat_results)

# let's look at top differences
library(dplyr)
filter(stat_results, qval<0.01) -> topdiffs
filter(transcript_gene_table, t_id %in% topdiffs$id) -> toptrans

pdf("topdiffs.pdf", 14, 7)
for (i in toptrans$g_id) {
  plotMeans(i, bg, groupvar='group', meas='FPKM', colorby='transcript')
  }
dev.off()

transcriptNames(bg)


### CHROMPLOT
# website
# https://bioconductor.org/packages/release/bioc/html/chromPlot.html

# Download and load
biocLite("chromPlot")
library(chromPlot)

# Browse help
browseVignettes("chromPlot")


data(hg_gap)
View(hg_gap)
chromPlot(gaps=hg_gap)

biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GenomicFeatures)
txgr <- transcripts(txdb)
txgr
chromPlot(gaps=hg_gap, annot1=txgr)

chromPlot(gaps=hg_gap, annot1=structure(bg)$trans)

# chromHeatMap
https://bioconductor.org/packages/release/bioc/html/ChromHeatMap.html


# interesting:
https://davetang.org/muse/2013/10/03/using-gviz/
