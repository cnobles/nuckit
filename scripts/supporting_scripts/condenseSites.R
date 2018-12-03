#' Deterimine all the sites present in a GRanges Object where each row is one
#' read. Read counts ('counts') for each site will be returned in the metadata
#' column of the resulting GRanges object.
#' @param sites GRanges object
#' @param keep.cols character vector with names of metadata columns to retain 
#' after condensing the ranges. Unique values will be split upon and will 
#' delimit data found with different metadata.
#' @param list.bps logical Whether breakpoint lengths from the start site should
#' be retained in an IntegerList object in the metadata after condensing the 
#' ranges. Default is FALSE.
#' @param list.bp.counts logical Whether the read counts for each breakpoint
#' length should be included as an IntegerList along with list.bps. Default is
#' FALSE, but TRUE will force list.bps to be TRUE as well.
#' @return GRanges object with one or more metadata columns.
#' @author Christopher Nobles, Ph.D.
#' 

condenseSites <- function(sites, keep.cols = NULL, list.bps = FALSE, 
                          list.bp.counts = FALSE){
  stopifnot(class(sites) == "GRanges")
  if(list.bp.counts) list.bps <- TRUE
  if(!is.null(keep.cols)){
    keptCols <- as.data.frame(
      mcols(sites)[match(keep.cols, names(mcols(sites)))])
    splitVec <- split(1:length(sites), keptCols)
    sitesList <- lapply(splitVec, function(x) granges(sites[x]))
    colList <- lapply(splitVec, function(x) as.data.frame(keptCols[x,]))
    colList <- lapply(colList, function(x){
      names(x) <- keep.cols
      return(x)})
  }else{
    sitesList <- GRangesList(sites)
    colList <- NULL
  }
  
  unlist(GRangesList(lapply(1:length(sitesList), function(i, bps, bp.counts){
    gr <- sitesList[[i]]
    cols <- colList[[i]]
    gr.red <- reduce(
      flank(gr, -1, start = TRUE), 
      min.gapwidth = 0L, 
      with.revmap = TRUE)
  
    revmap <- gr.red$revmap
    gr.red$revmap <- NULL
    bp.df <- data.frame(
      redid = Rle(
        values = 1:length(revmap), lengths = width(revmap@partitioning)),
      oriid = unlist(revmap),
      width = width(gr[unlist(revmap)]))
    bp.df <- bp.df[order(bp.df$width),]
    bp.rle <- split(Rle(bp.df$width), bp.df$redid)
    
    if(!is.null(cols)) mcols(gr.red) <- unique(cols)      
    gr.red$counts <- width(revmap@partitioning)
    gr.red$fragLengths <- width(runValue(bp.rle)@partitioning)
    if(bps){
      gr.red$bp.widths <- runValue(bp.rle)
      if(bp.counts) gr.red$bp.counts <- runLength(bp.rle)
    }
    return(gr.red)
  },
  bps = list.bps, 
  bp.counts = list.bp.counts)))
}