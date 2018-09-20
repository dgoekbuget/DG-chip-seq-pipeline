library(GenomicRanges)
library("GenomicFeatures")
library(Biostrings)

#Define genome

genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

#Define functions

gr2fasta <- function(gr,outfile,length="original",genome=genome){
  if (length == "original") gr.seq <- getSeq(genome,gr)
  else {features <- gr
  wid <- width(features)
  feature.recentered <- feature.center <- features
  start(feature.center) <- start(features) + floor(wid/2)
  width(feature.center) <- 1
  start(feature.recentered) <- start(feature.center) - (length/2-1)
  end(feature.recentered) <- end(feature.center) + length/2
  gr.seq <- getSeq(genome, feature.recentered)
  }
  writeXStringSet(gr.seq,outfile)
}

