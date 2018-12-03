#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Set Global options and load intiial packages ---------------------------------
options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)

code_dir <- dirname(sub("--file=", "", grep(
  "--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))

desc <- yaml::yaml.load_file(file.path(code_dir, "descriptions.yml"))


# Set up arguments and workflow of script --------------------------------------
## Argument parser =============================================================
parser <- ArgumentParser(
  description = desc$program_short_description)
parser$add_argument(
  "seqFile", nargs = "+", type = "character", help = desc$seqFile)
parser$add_argument(
  "-o", "--output", nargs = "+", type = "character", help = desc$output)
parser$add_argument(
  "-i", "--index", nargs = "*", type = "character", help = desc$index)
parser$add_argument(
  "--header", action = "store_true", help = desc$header)
parser$add_argument(
  "-n", "--negSelect", action = "store_true", help = desc$negSelect)
parser$add_argument(
  "-s", "--seq", nargs = "*", type = "character", help = desc$seq)
parser$add_argument(
  "-m", "--mismatch", nargs = "+", type = "integer", default = 0, 
  help = desc$mismatch)
parser$add_argument(
  "--any", action = "store_true", help = desc$any)
parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", default = "[\\w:-]+",
  help = desc$readNamePattern)
parser$add_argument(
  "--compress", action = "store_true", help = desc$compress)
parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", help = desc$cores)
parser$add_argument(
  "-q", "--quiet", action = "store_true", help = desc$quiet)
parser$add_argument(
  "--stat", nargs = 1, default = FALSE, type = "character", help = desc$stat)


## Parse cmd line args =========================================================
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


## Checks and balance ==========================================================
if(args$cores > 1){
  # Stop code since parallel operation has not been constructed yet
  stop("Parallel options have not yet been implemented.")
  
  if(args$cores > parallel::detectCores()){
    message(paste("Requested cores is greater than availible for system.",
      "Changing cores to max allowed."))
    args$cores <- detectCores()
  }
}else if(args$cores < 1){
  args$cores <- 1
}

if(length(args$seqFile) != length(args$output)){
  stop("The same number of input and output file names need to be provided.")
}

if(length(args$index) > 1){
  stop("Only one index file can be used at a time. Please consolidate indices.")
}

if(length(args$mismatch) != length(args$seq)){
  args$mismatch <- rep(args$mismatch[1], length(args$seq))
}

if(length(args$seq) > 0){
  args$seq <- toupper(gsub("U", "T", args$seq))
  if(any(!unlist(strsplit(paste(args$seq, collapse = ""), "")) %in% 
         names(Biostrings::IUPAC_CODE_MAP))){
    stop("Unknown nucleotides detected in input filtering sequence(s).")
  }
}

# Determine input sequence file type(s)
seqType <- unlist(strsplit(args$seqFile, "/"))
seqType <- seqType[length(seqType)]
seqType <- stringr::str_extract(seqType, ".fa[\\w]*")
if(any(!seqType %in% c(".fa", ".fq", ".fasta", ".fastq"))){
  stop(paste(
    "Unrecognized sequence file type, please convert to '*.fasta' or", 
    "'*.fastq'. Gzip compression is acceptable as well."))
}
seqType <- ifelse(seqType %in% c(".fa", ".fasta"), "fasta", "fastq")

# Determine sequence output file type(s)
if(length(args$output) > 0){
  outType <- unlist(strsplit(args$output, "/"))
  outType <- outType[length(outType)]
  outType <- stringr::str_extract(outType, ".fa[\\w]*")
  if(any(!outType %in% c(".fa", ".fq", ".fasta", ".fastq"))){
    stop(paste(
      "Unrecognized output sequence file type, please change to", 
      "'*.fasta' or '*.fastq'."))
  }
  outType <- ifelse(outType %in% c(".fa", ".fasta"), "fasta", "fastq")
}

# Identify filtering type
select_methods <- c()
if(length(args$index) == 1) select_methods <- c(select_methods, 1)
if(length(args$seqFile) > 1) select_methods <- c(select_methods, 2)
if(length(args$seq) > 0) select_methods <- c(select_methods, 3)
methods <- c(
  "input indices", "multiple file input indices", "sequence content")[
    select_methods]
filType <- paste0(
  ifelse(args$negSelect, "negative", "positive"), 
  " selection using ", 
  paste(methods, collapse = ifelse(args$any, " or ", " and ")), ".")


## Input arguments table =======================================================
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(
    c("seqFile :", "output :", "index :", "header :", "negSelect :", "seq :", 
      "mismatch :", "readNamePattern :", "compress :", "cores :"), 
    input_table$Variables),]
if(!args$quiet){
  pandoc.title("seqFiltR Inputs")
  pandoc.table(
    data.frame(input_table, row.names = NULL), 
    justify = c("left", "left"), 
    split.tables = Inf,
    caption = paste0("Filtering methods include ", filType))
}


# Load additional R-packages ---------------------------------------------------
if(args$cores > 1){
  addPacks <- c("stringr", "ShortRead", "parallel")
}else{
  addPacks <- c("stringr", "ShortRead")
}

addPacksLoaded <- suppressMessages(
  sapply(addPacks, require, character.only = TRUE))
if(!all(addPacksLoaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(addPacksLoaded), 
    "Loaded" = addPacksLoaded, 
    row.names = NULL))
  stop("Check dependancies.")
}


# Additional supporting functions ----------------------------------------------

#' Combine a list of ShortRead objects
#' 
#' @param split.seqs list of ShortRead objects
#' @author Christopher Nobles, Ph.D.
serial_append_S4 <- function(split.seqs){
  require("ShortRead")
  stopifnot(class(split.seqs) == "list")
  
  app_env <- new.env()
  app_env$seqs <- split.seqs[[1]]
  
  null <- lapply(2:(length(split.seqs)), function(i){
    app_env$seqs <- append(app_env$seqs, split.seqs[[i]])
  })
  
  return(app_env$seqs)
} 

#' Write fast(a/q) files given Biostrings and ShortRead objects.
#' 
#' @usage write_seq_file(pointer, seqs, seqType, file, compress = FALSE)
#' 
#' @param pointer ShortRead pointer generated by using ShortRead::readFastq().
#' Not required for fasta writing as sequences are already in Biostrings object.
#' @param seqs BioStrings DNAStringSet object generated by trimming from the 
#' pointer object.
#' @param seqType character either "fasta" or "fastq". The format for the 
#' output file.
#' @param file character Output file name to write sequences.
#' @param compress logical Whether to gzip compress the sequence file when 
#' written or leave it uncompressed. Default is FALSE, or uncompressed.
#' @author Christopher Nobles, Ph.D.

write_seq_files <- function(seqs, seqType, file, compress = FALSE){
  packs <- c("ShortRead")
  packsLoaded <- suppressMessages(sapply(packs, require, character.only = TRUE))
  stopifnot(all(packsLoaded))
  
  if(seqType == "fasta"){
    if(compress){
      if(grepl(".gz$", file)){
        ShortRead::writeFasta(
          seqs, file = file, width = max(width(seqs)), compress = TRUE)    
      }else{
        ShortRead::writeFasta(
          seqs, file = paste0(file, ".gz"), 
          width = max(width(seqs)), compress = TRUE)
      }
    }else{
      ShortRead::writeFasta(
        seqs, file = file, width = max(width(seqs)), compress = FALSE)
    }
  }else{
    if(compress){
      if(grepl(".gz$", file)){
        ShortRead::writeFastq(seqs, file = file, compress = TRUE)
      }else{    
        ShortRead::writeFastq(seqs, file = paste0(file, ".gz"), compress = TRUE)
      }
    }else{
      ShortRead::writeFastq(
        seqs, filepath = file, compress = FALSE)
    }
  }
}

#' Filter sequences based on input arguments
#' This function is the basis for the script.
filter_seqFile <- function(input_seqs, args){
  suppressMessages(require("stringr"))
  suppressMessages(require("ShortRead"))
  
  ## Identify sequence names matching across multiple sequence files
  if(length(input_seqs) > 1){
    multi_input_ids <- lapply(input_seqs, function(seq){
      stringr::str_extract(
        as.character(unique(id(seq))), args$readNamePattern)})
    multi_input_tbl <- table(unlist(multi_input_ids))
    if(args$negSelect){
      multi_input_names <- names(multi_input_tbl)[which(multi_input_tbl == 1)]
    }else if(args$any){
      multi_input_names <- names(multi_input_tbl)[which(multi_input_tbl > 1)]
    }else{
      multi_input_names <- names(multi_input_tbl)[
        which(multi_input_tbl == length(input_seqs))]
    }
    
    multi_filter_idx <- lapply(input_seqs, function(seqs, idx){
      ids <- stringr::str_extract(as.character(id(seqs)), args$readNamePattern)
      which(ids %in% idx)
    }, idx = multi_input_names)
  }
  
  
  ## Identify sequence names by matching to index file
  if(length(args$index) == 1){
    input_ids <- lapply(input_seqs, function(seq){
      stringr::str_extract(as.character(id(seq)), args$readNamePattern)})
    index_df <- read.delim(args$index, header = args$header)
    index <- stringr::str_extract(
      as.character(index_df[,1]), args$readNamePattern)
    
    index_filter_idx <- lapply(input_ids, function(ids, idx){
      which(ids %in% idx)
    }, idx = index)
  }
  
  
  ## Identify sequences by matching input nucleotide sequence
  if(length(args$seq) > 0){
    seq_filter_idx <- lapply(
      input_seqs, function(seqs, pattern, mismatch, neg, any){
        
        vcp <- mapply(function(pat, mis, seqs, neg){
          v <- vcountPattern(
            pat, sread(seqs), max.mismatch = mis, fixed = FALSE)
          if(neg){
            return(which(v == 0))
          }else{
            return(which(v > 0))
          }
        }, pat = pattern, mis = mismatch, 
        MoreArgs = list(seqs = seqs, neg = neg),
        SIMPLIFY = FALSE)
        
        vcp_tbl <- table(unlist(vcp))
        if(any){
          return(as.numeric(names(vcp_tbl[which(vcp_tbl >= 1)])))
        }else{
          return(as.numeric(names(vcp_tbl[which(vcp_tbl == length(pattern))])))
        }
        
      }, pattern = args$seq, mismatch = args$mismatch, 
      neg = args$negSelect, any = args$any)
    
  }
  
  
  # Consolidate indices from each method employed 
  lapply(1:length(input_seqs), function(i){
    idx <- NULL
    cnt <- 0
    if(exists("multi_filter_idx")){ 
      cnt <- cnt + 1
      idx <- c(idx, multi_filter_idx[[i]]) }
    if(exists("index_filter_idx")){ 
      cnt <- cnt + 1
      idx <- c(idx, index_filter_idx[[i]]) }
    if(exists("seq_filter_idx")){ 
      cnt <- cnt + 1
      idx <- c(idx, seq_filter_idx[[i]]) }
    if(args$any){
      return(unique(idx))
    }else{
      idx_tbl <- table(idx)
      return(as.numeric(names(idx_tbl)[idx_tbl == cnt]))
    }
  })
}


# Identify indices of input file(s) for filtering ------------------------------
input_seqs <- mapply(function(file, fileType){
  if(fileType == "fasta"){
    return(ShortRead::readFasta(file))
  }else{
    return(ShortRead::readFastq(file))
  }}, file = args$seqFile, fileType = seqType, SIMPLIFY = FALSE)

output_indices <- filter_seqFile(input_seqs, args)

output_seqs <- mapply(
  function(seqs, idx){ seqs[idx] }, 
  seqs = input_seqs, idx = output_indices, SIMPLIFY = FALSE)


# Write output files -----------------------------------------------------------
if(args$stat != FALSE){
  sampleName <- strsplit(args$output, "/", fixed = TRUE)
  sampleName <- mapply("[[", sampleName, lengths(sampleName))
  sampleName <- strsplit(sampleName, ".fa", fixed = TRUE)
  sampleName <- mapply("[[", sampleName, 1)
  write.table(
    data.frame(
      sampleName = sampleName,
      metric = "reads",
      count = lengths(output_seqs)),
    file = args$stat,
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
 
 
null <- sapply(args$output, unlink)

null <- mapply(
  write_seq_files, seqs = output_seqs, seqType = outType, file = args$output, 
  MoreArgs = list(compress = args$compress))

q()

# Notes:
## Work flow: R1xR2 --> index --> sequence
## 
## For output from each step, pass on a list of indices for the ShortRead object
## This will allow for filtering and then just returning an object that contains
## the correct indices. 
## 
## For R1xR2 type filtering: intersect / c(id()) | duplicated / table / which(bool)
## For index type filtering: which(id() %in% index) or which(!id() %in% index)
## For seq type filtering  : v[match/count]pattern(fixed = FALSE) | which(lengths() >= 1 or lengths() == 0)
