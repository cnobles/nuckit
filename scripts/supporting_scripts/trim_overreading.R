#' Trim 3' ends of sequences when the sequence has overread into
#' synthetic sequence
#' \code{trim_overreading} removes 3' ends of nucleotide sequences completely or
#' partially matching the trimSequence.
#' @param seqs ShortReadQ object of reads or unique sequences
#' @param trimSequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match the end of sequences, upon which the matching
#' portion will be trimmed from the start of the match to the end of the
#' sequence. Ambiguous nucleotides within the sequence can be used for alignment,
#' but matching sequences are not recored.
#' @param percentID numeric between 0 and 1.0 denoting the minimum percent
#' identity acceptable for an matching alignment.
#' @param maxSeqLength numeric/integer the maximum length to consider of the
#' trimSequence to use for alignments. Using the full length of sequence
#' avaliable for matching can many times be computationally intensive and
#' unnessesarily time consuming. Further, identical results can be obtained
#' using only a portion of the sequence. For example, setting the maxSeqLength
#' to 15L will only use the first 15 nucleotides of the trimSequence.
#' @author Christopher Nobles, Ph.D.

trim_overreading <- function(seqs, trimSequence,
                             percentID, maxSeqLength = NULL, 
                             minSeqLength = 3){
  suppressMessages(require(BiocGenerics))
  suppressMessages(require(Biostrings))
  stopifnot(class(seqs) %in% c("ShortReadQ", "ShortRead"))
  stopifnot(!is.null(ShortRead::id(seqs)))
  
  # Trim down trimSequence if maxSeqLength provided
  if(!is.null(maxSeqLength)){
    trimSequence <- Biostrings::DNAStringSet(
      trimSequence, start = 1L, end = min(nchar(trimSequence), maxSeqLength))
  }
  
  trimSeqs <- sapply(0:(nchar(trimSequence)-minSeqLength), function(i){
    substr(trimSequence, 1, nchar(trimSequence) - i)
  })
  
  alignments <- do.call(c, lapply(trimSeqs, function(trimSeq, seqs, percentID){
      mismatch <- round(nchar(trimSeq) - percentID*nchar(trimSeq))
      vmp <- Biostrings::vmatchPattern(
        trimSeq, ShortRead::sread(seqs), max.mismatch = mismatch, fixed = FALSE)
      idx <- which(BiocGenerics::lengths(vmp) >= 1)
      len <- BiocGenerics::lengths(vmp[idx])
      len <- len[which(len >= 1)]
      idx <- S4Vectors::Rle(values = idx, lengths = len)
      ir <- unlist(vmp)
      names(ir) <- idx
      ir
    }, seqs = seqs, percentID = percentID))
  
  seqWidths <- Biostrings::width(seqs[as.numeric(names(alignments))])
  alignments <- alignments[
    IRanges::width(alignments) == nchar(trimSequence) | 
    IRanges::end(alignments) == seqWidths]
  alignments <- split(alignments, names(alignments))
  alignments <- unlist(IRanges::reduce(alignments, min.gapwidth = 1000000))
  alignments <- alignments[order(as.numeric(names(alignments)))]
  idx <- as.numeric(names(alignments))
  
  if(any(table(idx) > 1)){
    stop("Issues with overread trimming, please adjust input parameters.")}
  
  # Trim sequences
  seqs[idx] <- IRanges::narrow(
    seqs[idx], end = IRanges::start(alignments)-1)
  return(seqs)
}
