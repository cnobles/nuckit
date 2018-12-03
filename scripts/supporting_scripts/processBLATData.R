#' Converts alignment data into a GRanges object without information loss. BLAT 
#' aligns on a 0-base system, while all genomic data is read on a 1-base system.
#' 
#' @param algns data.frame of psl table, containing labeled columns. 
#' @param from character, which read is the algns object from? ("anchor" or 
#' "adrift")
#' @param refGenome BSgenome object with seqInfo of reference genome.
#' @author Christopher Nobles, Ph.D.
#' 

processBLATData <- function(algns, from, refGenome){
  stopifnot(from == "anchor" | from == "adrift")
  algns$from <- from
  algns$qtStart <- ifelse(
    algns$strand == "+",
    (algns$tStart - (algns$qStart)),
    (algns$tStart - (algns$qSize - algns$qEnd - 1)))
  algns$qtEnd <- ifelse(
    algns$strand == "+",
    (algns$tEnd + (algns$qSize - algns$qEnd - 1)),
    (algns$tEnd + (algns$qStart)))    
  
  algns.gr <- GRanges(seqnames=Rle(algns$tName),
                      ranges = IRanges(
                        start = (algns$qtStart + 1), 
                        end = (algns$qtEnd)), #Convert to 1-base
                      strand=Rle(algns$strand),
                      seqinfo=seqinfo(refGenome))
  
  mcols(algns.gr) <- algns[,c("from", "qName", "matches", "repMatches", 
                              "misMatches", "qStart", "qEnd", "qSize", 
                              "tBaseInsert")]
  rm(algns)
  algns.gr
}
