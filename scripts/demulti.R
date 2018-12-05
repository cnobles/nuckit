#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, including data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions
# Set Global options and load intiial packages ---------------------------------
options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)

code_dir <- dirname(sub(
  "--file=", "", grep(
    "--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
))

desc <- yaml::yaml.load_file(file.path(code_dir, "demulti.desc.yml"))

# Set up and gather command line arguments -------------------------------------
## Argument parser =============================================================
parser <- ArgumentParser(description = desc$program_short_description)

parser$add_argument(
  "-m", "--manifest", type = "character", 
  help = desc$manifest
)

parser$add_argument(
  "--read1", type = "character", default = "NA", 
  help = desc$read1
)

parser$add_argument(
  "--read2", type = "character", default = "NA", 
  help = desc$read2
)

parser$add_argument(
  "--index1", type = "character", default = "NA", 
  help = desc$index1
)

parser$add_argument(
  "--index2", type = "character", default = "NA", 
  help = desc$index2
)

parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", 
  help = desc$outfolder
)

parser$add_argument(
  "-p", "--poolreps", action = "store_true", 
  help = desc$poolreps
)

parser$add_argument(
  "--singleBarcode", action = "store_true", 
  help = desc$singleBarcode
)

parser$add_argument(
  "--barcode1", nargs = 1, type = "character", default = "I1", 
  help = desc$barcode1
)

parser$add_argument(
  "--barcode2", nargs = 1, type = "character", default = "I2", 
  help = desc$barcode2
)

parser$add_argument(
  "--barcode1Length", nargs = 1, type = "integer", default = 8, 
  help = desc$barcode1Length
)

parser$add_argument(
  "--barcode2Length", nargs = 1, type = "integer", default = 8,
  help = desc$barcode2Length
)

parser$add_argument(
  "--maxMismatch", nargs = 1, type = "integer", 
  help = desc$maxMismatch
)

parser$add_argument(
  "--bc1Mismatch", nargs = 1, type = "integer", default = 0, 
  help = desc$bc1Mismatch
)

parser$add_argument(
  "--bc2Mismatch", nargs = 1, type = "integer", default = 0,
  help = desc$bc2Mismatch
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = desc$stat
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", 
  help = desc$readNamePattern
)

parser$add_argument(
  "--compress", action = "store_true", 
  help = desc$compress
)

parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", 
  help = desc$cores
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


demulti <- data.frame(
  "readType" = c("R1", "R2", "I1", "I2"),
  "path" = c(args$read1, args$read2, args$index1, args$index2)
)

demulti$barcode1 <- grepl(args$barcode1, demulti$readType)
demulti$barcode2 <- grepl(args$barcode2, demulti$readType)


if( demulti$readType[demulti$barcode1] == demulti$readType[demulti$barcode2] ){
  stop("Please select different read types for barcodes 1 and 2.")
}

if( demulti$readType[demulti$barcode1] == "NA" ){
  stop("Barcode 1 is set to a read type that is not provided.")
}

if( demulti$readType[demulti$barcode2] == "NA" ){
  stop("Barcode 2 is set to a read type that is not provided.")
}

if( args$singleBarcode ){
  demulti$barcode2 <- FALSE
}

if( !is.null(args$maxMisMatch) ){
  args$bc1Mismatch <- args$maxMismatch
  args$bc2Mismatch <- args$maxMismatch
}


input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply( 1:length(args), function(i){
    paste(args[[i]], collapse = ", ")
  } )
)

input_table <- input_table[
  match(
    c("manifest :", "index1 :", "index2 :", "read1 :", "read2 :", 
      "outfolder :", "stat :", "poolreps :", "singleBarcode :", "cores :", 
      "barcode1 :", "barcode2 :", "barcode1Length :", "barcode2Length :", 
      "bc1Mismatch :", "bc2Mismatch :", "readNamePattern :"),
    input_table$Variables
  ),
]

pandoc.title("Demultiplex Inputs")
pandoc.table(
  data.frame(input_table, row.names = NULL), 
  justify = c("left", "left"), 
  split.tables = Inf
)

# Create output directory if not currently available ---------------------------
if( !file.exists(args$outfolder) ){
  
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.")

}

# Check for required packages --------------------------------------------------
required_packs <- c("stringr", "ShortRead", "Biostrings")
present_packs <- required_packs %in% row.names(installed.packages())

if( !all(present_packs) ){
  
  pandoc.table(data.frame(
    "R-Packages" = required_packs, 
    "Installed" = present_packs, 
    row.names = NULL
  ))
  stop("Check dependancies.")

}

# Operating functions ----
parseIndexReads <- function(barcode, index.file.path, barcode.length, 
                            max.mismatch, read.name.pattern){
  
  # Load index file sequences and sequence names
  index <- ShortRead::readFastq(index.file.path)
  index <- ShortRead::narrow(index, start = 1, end = barcode.length)
  index@id <- Biostrings::BStringSet(
    stringr::str_extract(
      as.character(ShortRead::id(index)), 
      read.name.pattern)
  )
  
  # Trim barcode if necessary
  barcode <- as.character(
    Biostrings::DNAStringSet(barcode, start = 1, end = barcode.length)
  )
  
  # Identify read names with sequences above or equal to the minscore
  vmp <- Biostrings::vmatchPattern(
    barcode, ShortRead::sread(index), max.mismatch = max.mismatch, fixed = FALSE
  )
  
  return( which(lengths(vmp) == 1) )
  
}

writeDemultiplexedSequences <- function(read.file.path, type, multiplexed.data, 
                                        read.name.pattern, out.folder, compress){
  
  # Load read sequences and sequence names then write to file
  reads <- ShortRead::readFastq(read.file.path)
  reads@id <- Biostrings::BStringSet(
    stringr::str_extract(as.character(ShortRead::id(reads)), read.name.pattern)
  )
  
  reads <- reads[multiplexed.data$index]
  reads <- split(reads, multiplexed.data$sampleName)
  
  null <- lapply(
    seq_along(reads), 
    function(i, reads, type, out.folder, compress){
      
      if(compress){  
        
        filePath <- file.path(
          out.folder, paste0(names(reads[i]), ".", type, ".fastq.gz")
        )
        
      }else{
        
        filePath <- file.path(
          out.folder, paste0(names(reads[i]), ".", type, ".fastq")
        )
        
      }
      
      if( file.exists(filePath) ) unlink(filePath)
      
      ShortRead::writeFastq(reads[[i]], file = filePath, compress = compress)
      
      message(
        paste0("\nWrote ", length(reads[[i]]), " reads to:\n", filePath, "."))
      
    },
    reads = reads, type = type, out.folder = out.folder, compress = compress
  )
  
  return(list(read.file.path, type, out.folder))
  
}
# Load manifest / sample mapping file ------------------------------------------
fileExt <- unlist(strsplit(args$manifest, "\\."))
fileExt <- fileExt[length(fileExt)]

if( fileExt %in% c("yaml", "yml") ){
  
  if( !"yaml" %in% row.names(installed.packages()) ){
    stop("Package:yaml not loaded or installed.")
  }
  
  manifest <- yaml::yaml.load_file(args$manifest)
  
  if( args$singleBarcode ){
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "barcode1" = sapply( manifest$samples, function(x) x$barcode1 ),
      row.names = NULL
    )
  
  }else{
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "barcode1" = sapply( manifest$samples, function(x) x$barcode1 ),
      "barcode2" = sapply( manifest$samples, function(x) x$barcode2 ),
      row.names = NULL
    )
  
  }

}else{
  
  if( fileExt == "csv" ){
    manifest <- read.csv(args$manifest)
  }else if( fileExt == "tsv" ){
    manifest <- read.delim(args$manifest)
  }
  
  if( args$singleBarcode ){
    samples_df <- manifest[, c("sampleName", "barcode1")]
  }else{
    samples_df <- manifest[, c("sampleName", "barcode1", "barcode2")]
  }
  
}


if( !args$singleBarcode ){
  
  uniqueSamples <- nrow(samples_df[,c("barcode1", "barcode2")]) == 
    nrow(unique(samples_df[,c("barcode1", "barcode2")]))
  
  if( !uniqueSamples ) stop("Ambiguous barcoding of samples. Please correct.")

}else{

  uniqueSamples <- length(samples_df[,c("barcode1")]) == 
    length(unique(samples_df[,"barcode1"]))
  
  if( !uniqueSamples ) stop("Ambiguous barcoding of samples. Please correct.")

}

# Read in barcode sequences ----------------------------------------------------
barcode1_reads <- ShortRead::readFastq(demulti$path[demulti$barcode1])
message(paste("Reads to demultiplex : ", length(barcode1_reads)))

if( args$cores > 1 ){
  
  cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))
  
  BC1_parsed <-  parallel::parLapply(
    cluster,
    unique(samples_df$barcode1), 
    parseIndexReads,
    index.file.path = demulti$path[demulti$barcode1],
    barcode.length = args$barcode1Length,
    max.mismatch = args$bc1Mismatch,
    read.name.pattern = args$readNamePattern
  )
  
  names(BC1_parsed) <- unique(samples_df$barcode1)
  
  pandoc.title("Barcode1 breakdown:")
  pandoc.table(data.frame(
    "Barcode1" = names(BC1_parsed),
    "Read Counts" = sapply( BC1_parsed, length ),
    row.names = NULL
  ))
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- parallel::parLapply(
      cluster, 
      unique(samples_df$barcode2), 
      parseIndexReads,
      index.file.path = demulti$path[demulti$barcode2],
      barcode.length = args$barcode2Length,
      max.mismatch = args$bc2Mismatch,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
  parallel::stopCluster(cluster)
  
}else{
  
  BC1_parsed <-  lapply(
    unique(samples_df$barcode1), 
    parseIndexReads,
    index.file.path = demulti$path[demulti$barcode1],
    barcode.length = args$barcode1Length,
    max.mismatch = args$bc1Mismatch,
    read.name.pattern = args$readNamePattern
  )

  names(BC1_parsed) <- unique(samples_df$barcode1)
  
  pandoc.title("Barcode1 breakdown:")
  pandoc.table(data.frame(
    "Barcode1" = names(BC1_parsed),
    "Read Counts" = sapply( BC1_parsed, length ),
    row.names = NULL
  ))
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- lapply(
      unique(samples_df$barcode2), 
      parseIndexReads,
      index.file.path = demulti$path[demulti$barcode2],
      barcode.length = args$barcode2Length,
      max.mismatch = args$bc2Mismatch,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
}

if( !args$singleBarcode ){

  names(BC2_parsed) <- unique(samples_df$barcode2)
  pandoc.title("Barcode2 breakdown:")
  pandoc.table(data.frame(
    "Barcode2" = names(BC2_parsed),
    "Read Counts" = sapply(BC2_parsed, length),
    row.names = NULL
  ))

}

if( !args$singleBarcode ){
  
  demultiplexed_indices <- mapply(
    function(barcode1, barcode2){
      Biostrings::intersect(BC1_parsed[[barcode1]], BC2_parsed[[barcode2]])
    },
    barcode1 = samples_df$barcode1,
    barcode2 = samples_df$barcode2,
    SIMPLIFY = FALSE
  )
  
  names(demultiplexed_indices) <- paste0(
    samples_df$barcode1, samples_df$barcode2
  )
  
}else{
  
  demultiplexed_indices <- BC1_parsed
  
}

# As there is some flexibility in the barcode matching, some reads may be 
# be assigned to multiple samples. These reads are ambiguous and will be 
# removed.
ambiguous_indices <- unique(
  unlist(demultiplexed_indices)[duplicated(unlist(demultiplexed_indices))]
)

demultiplexed_indices <- lapply(demultiplexed_indices, function(x, reads){
    x[!x %in% reads]
  },
  reads = ambiguous_indices
)

# Reads by sample
samples_df$read_counts <- sapply( demultiplexed_indices, length )
pandoc.title("Read counts for each sample.")
pandoc.table(samples_df, split.tables = Inf)

# Ambiguous reads
message(paste0("Ambiguous reads: ", length(ambiguous_indices)))

# Unassigned reads
all_indices <- seq_along(barcode1_reads)

unassigned_indices <- all_indices[
  !all_indices %in% unlist(demultiplexed_indices, use.names = FALSE)
]

unassigned_indices <- unassigned_indices[
  !unassigned_indices %in% ambiguous_indices
]

message(paste0("\nUnassigned reads: ", length(unassigned_indices)))

if( args$stat != FALSE ){
  write.table(
    data.frame(
      sampleName = paste0(
        c(samples_df$sampleName, "ambiguous_reads", "unassigned_reads"), 
        ".demulti"
      ),
      metric = "reads",
      count = c(
        samples_df$read_counts, 
        length(ambiguous_indices), 
        length(unassigned_indices)
      )
    ),
    file = file.path(args$outfolder, args$stat),
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

# Create multiplex dataframe for subseting sequencing files --------------------
multiplexed_data <- data.frame(
  "sampleName" = S4Vectors::Rle(
    values = samples_df$sampleName, 
    length = sapply(demultiplexed_indices, length)
  ),
  "index" = unlist(demultiplexed_indices),
  row.names = NULL
)

ambiguous_data <- data.frame(
  "sampleName" = rep("ambiguous", length(ambiguous_indices)),
  "index" = ambiguous_indices,
  row.names = NULL
)

unassignedData <- data.frame(
  "sampleName" = rep("unassigned", length(unassigned_indices)),
  "index" = unassigned_indices,
  row.names = NULL
)

multiplexed_data <- rbind(multiplexed_data, ambiguous_data, unassignedData)
multiplexed_data$sampleName <- factor(
  multiplexed_data$sampleName,
  levels = c(samples_df$sampleName, "ambiguous", "unassigned")
)

stopifnot(nrow(multiplexed_data) == length(all_indices))

if( args$poolreps ){
  multiplexed_data$sampleName <- gsub("-\\d+$", "", multiplexed_data$sampleName)
}

message(paste0("\nReads to be written to files: ", nrow(multiplexed_data)))

# Write files to read files to outfolder directory
if( args$cores > 1 ){
  
  cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))
  
  read_list <- demulti$readType[demulti$path != "NA"]
  read_paths <- demulti$path[match(read_list, demulti$readType)]

  demultiplex <- parallel::clusterMap(
    cluster,
    writeDemultiplexedSequences,
    read.file.path = read_paths,
    type = read_list,
    MoreArgs = list(
      multiplexed.data = multiplexed_data,
      read.name.pattern = args$readNamePattern,
      out.folder = args$outfolder,
      compress = args$compress
    )
  )
  
  parallel::stopCluster(cluster)
  
}else{
  
  read_list <- demulti$readType[demulti$path != "NA"]
  read_paths <- demulti$path[match(read_list, demulti$readType)]
  
  demultiplex <- mapply(
    writeDemultiplexedSequences,
    read.file.path = read_paths,
    type = read_list,
    MoreArgs = list(
      multiplexed.data = multiplexed_data,
      read.name.pattern = args$readNamePattern,
      out.folder = args$outfolder,
      compress = args$compress
    )
  )

}

message("\nDemultiplexing complete.")
q()
