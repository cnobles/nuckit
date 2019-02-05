# Import yaml file with test paths and md5 digests
# Required libraries stringr, digest, Biostrings, data.table, yaml, pander
options(stringsAsFactors = FALSE, scipen = 99, warn = -1, window = 999)
suppressMessages(library("magrittr"))

# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Test checksums of files from a yaml input.",
  usage = "Rscript tools/rscripts/check_file_digests.R <input.file> <options>"
)

parser$add_argument(
  "file", nargs = 1, type = "character",
  help = "CSV file containing paths and checksums (md5). ie. test.digest.csv"
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = FALSE,
  help = "Output file name .csv, .tsv, or .rds format."
)

parser$add_argument(
  "-v", "--verbose", action = "store_true", 
  help = "Turns on diagnositc-based messages."
)

# Set arguments with parser ----
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

root_dir <- getwd()

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(
    seq_along(args), 
    function(i) paste(args[[i]], collapse = ", ")
  )
)

input_table <- input_table[
  match(
    c("file :", "output :", "verbose :"), 
    input_table$Variables
  ),
]

## Log inputs
if( args$verbose ){
  
  cat("List Inputs")
  print(
    data.frame(input_table),
    right = FALSE, 
    row.names = FALSE
  )
  
}


# Additional functions ----
readFile <- function(path, root){
  
  if( !file.exists(path) ){
    root_path <- file.path(root, path)
    if( !file.exists(root_path) ){
      stop("Cannot find file:", path)
    }else{
      path <- root_path
    }
  }
  
  # Read extension form path
  ext <- stringr::str_extract(path, "[\\w]+$")
  supported_ext <- c("tsv", "csv", "gz", "fasta", "fastq", "rds", "txt")
  
  stopifnot( ext %in% supported_ext )
  
  # Check additional extension if compressed
  if( ext == "gz" ){
    ext2 <- stringr::str_extract(path, "[\\w]+.gz")
    ext2 <- gsub(".gz", "", ext2)
    stopifnot( ext2 %in% supported_ext )
  }else{
    ext2 <- NA
  }
  
  exts <- c(ext, ext2)
  exts <- exts[!is.na(exts)]
  
  # Read in methods based on inputs.
  if( any(exts %in% c("tsv", "csv", "txt")) ){
    
    if( ext == "gz" ) path <- paste0("zcat ", path)
    return(data.table::fread(input = path, data.table = FALSE))
    
  }else if( any(stringr::str_detect(exts, "fast")) ){
    
    if( any(stringr::str_detect(exts, "fasta")) ){
      return(as.character(ShortRead::sread(ShortRead::readFasta(path))))
    }else{
      return(as.character(ShortRead::sread(ShortRead::readFastq(path))))
    }

  }else{
    
    return(readRDS(path))
    
  }
  
}

# Load inputs ----
digest_file <- read.csv(args$file)
paths <- digest_file$name
data_objs <- lapply(
  paths, 
  readFile, 
  root = file.path(root_dir, "etc/test_output")
)

# Check digests ----
test_digests <- sapply(data_objs, digest::digest)
check_digests <- digest_file$md5

df <- data.frame(
  "file_name" = digest_file$name,
  "md5_standard" = check_digests,
  "md5_tested" = test_digests,
  "outcome" = ifelse(test_digests == check_digests, "pass", "FAIL")
)

# Write output file if requested ----
if( args$output != FALSE ){
  if( stringr::str_detect(args$output, ".tsv$") ){
    write.table(df, file = args$output, quote = FALSE, row.names = FALSE)
  }else if( stringr::str_detect(args$output, ".csv$") ){
    write.csv(df, file = args$output, quote = FALSE, row.names = FALSE)
  }else if( stringr::str_detect(args$output, ".rds$") ){
    saveRDS(df, file = args$output)
  }else if( stringr::str_detect(args$output, ".RData$") ){
    save(df, file = args$output)
  }
}

# Log output if requested ----
if( args$verbose ){
  
  df$md5_standard <- paste0(substr(df$md5_standard, start = 1, stop = 7), "...")
  df$md5_tested <- paste0(substr(df$md5_tested, start = 1, stop = 7), "...")
  
  cat("\nList of Outcomes:\n")
  print(
    df,
    right = FALSE, 
    row.names = FALSE
  )
  
}

# Finish up and close out ----
if( all(df$outcome == "pass") ){
  q(save = "no", status = 0)
}else{
  q(save = "no", status = 1)
}

