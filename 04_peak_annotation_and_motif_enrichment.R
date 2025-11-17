
#ArchR

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(motifmatchr)
library(chromVAR)
library(GenomicRanges)
library(rhdf5)
library(HDF5Array)


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

peaks_gr <- getPeakSet(projHeme5)
peak_anno <- annotatePeak(
  peaks_gr,
  tssRegion = c(-2000, 500),  #promotor region
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
)

##FGFR3 as example
anno_df <- as.data.frame(peak_anno@anno)

fgfr3_peaks <- anno_df[
  (grepl("FGFR3", anno_df$nearestGene) | grepl("FGFR3", anno_df$SYMBOL)) & 
    anno_df$Group %in% c("PC/DC_progenitor"), 
]  #c("IBC", "IPhC", "IBC/IPhC_progenitor", "IHC", "HC_progenitor")/("PC", "DC",  "PC/DC_progenitor","OHC", "HC_progenitor")
fgfr3_peaks_gr <- GRanges(
  seqnames = fgfr3_peaks$seqnames,
  ranges = IRanges(start = fgfr3_peaks$start, end = fgfr3_peaks$end),
  strand = fgfr3_peaks$strand
)

###motif enrichment
pSet <- getPeakSet(ArchRProj = projHeme5)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
matches <- getMatches(ArchRProj = projHeme5, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

gr <- GRanges(seqnames = c("chr4"), ranges = IRanges(start = c(1758500), end = c(1758999)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))

p1 = colnames(matches)[which(assay(matches[queryHits,]))]

###GC contribution
pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["RORA_658"]]
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
ppm <- PWMatrixToProbMatrix(pwm)
library(ggseqlogo)
ggseqlogo(ppm, method = "bits")

#otherwise in Signac
library(Signac)
MotifPlot(
object = snRNA,
motifs = "MA0624.1"
)






