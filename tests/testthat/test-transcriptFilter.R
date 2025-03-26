library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plotgardener)

test_that("plotTranscripts transcritFilter works", {
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  tbg <- transcriptsBy(txdb, by="gene")
  gene <- tbg[["7168"]]
  range <- range(gene)

  par <- pgParams(
    chrom = as.character(seqnames(range)),
    chromstart = start(range)-1e4, chromend = end(range)+1e4,
    assembly = "hg38", just = c("left", "bottom")
  )
  
  pageCreate(width = 5, height = 2, showGuides = TRUE)
  
  # highlight can be a subset of filter
  hilite <- data.frame(transcript=c("ENST00000559831.6","ENST00000288398.10"), color=c("blue","red"))
  filter <- c("ENST00000559831.6","ENST00000288398.10","ENST00000651344.1","ENST00000560975.5")
 
  # should only show 4 out of 44 total txps
  plotTranscripts(
    params = par, x = 0.5, y = 1.5, width = 4, height = 1.1,
    transcriptHighlight=hilite, transcriptFilter=filter
  )
 
})
