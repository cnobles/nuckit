#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, including data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions
# Set Global options and load intiial packages ----
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/demulti.desc.yml")
)

# Set up and gather command line arguments ----
## Argument parser ----
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "nuc demulti [-h/--help, -v/--version] [optional args]"
)

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
  "--idx1", type = "character", default = "NA", 
  help = desc$idx1
)

parser$add_argument(
  "--idx2", type = "character", default = "NA", 
  help = desc$idx2
)

parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", 
  help = desc$outfolder
)

parser$add_argument(
  "--bc1", nargs = 1, type = "character", default = "I1", 
  help = desc$bc1
)

parser$add_argument(
  "--bc2", nargs = 1, type = "character", default = "I2", 
  help = desc$bc2
)

parser$add_argument(
  "--bc1Man", nargs = 1, type = "character", default = "barcode1", 
  help = desc$bc1Man
)

parser$add_argument(
  "--bc2Man", nargs = 1, type = "character", default = "barcode2", 
  help = desc$bc2Man
)

parser$add_argument(
  "--bc1Len", nargs = 1, type = "integer", default = 8, 
  help = desc$bc1Len
)

parser$add_argument(
  "--bc2Len", nargs = 1, type = "integer", default = 8,
  help = desc$bc2Len
)

parser$add_argument(
  "--maxMis", nargs = 1, type = "integer", 
  help = desc$maxMis
)

parser$add_argument(
  "--bc1Mis", nargs = 1, type = "integer", default = 0, 
  help = desc$bc1Mis
)

parser$add_argument(
  "--bc2Mis", nargs = 1, type = "integer", default = 0,
  help = desc$bc2Mis
)

parser$add_argument(
  "--maxN", nargs = 1, type = "integer", default = 1,
  help = desc$maxN
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = desc$stat
)

parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", 
  help = desc$cores
)

parser$add_argument(
  "--compress", action = "store_true", 
  help = desc$compress
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
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", 
  help = desc$readNamePatter
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


demulti <- data.frame(
  "readType" = c("R1", "R2", "I1", "I2"),
  "path" = c(args$read1, args$read2, args$idx1, args$idx2)
)

demulti$bc1 <- grepl(args$bc1, demulti$readType)
demulti$bc2 <- grepl(args$bc2, demulti$readType)


if( demulti$readType[demulti$bc1] == demulti$readType[demulti$bc2] ){
  stop("Please select different read types for barcodes 1 and 2.\n")
}

if( demulti$readType[demulti$bc1] == "NA" ){
  stop("Barcode 1 is set to a read type that is not provided.\n")
}

if( demulti$readType[demulti$bc2] == "NA" ){
  stop("Barcode 2 is set to a read type that is not provided.\n")
}

if( args$singleBarcode ){
  demulti$bc2 <- FALSE
}

if( !is.null(args$maxMis) ){
  args$bc1Mis <- args$maxMis
  args$bc2Mis <- args$maxMis
}


input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply( 
    seq_along(args), 
    function(i){
      paste(args[[i]], collapse = ", ")
    }
  )
)

input_table <- input_table[
  match(
    c("manifest :", "idx1 :", "idx2 :", "read1 :", "read2 :", 
      "outfolder :", "stat :", "bc1 :", "bc2 :", "bc1Man :", "bc2Man :", 
      "bc1Len :", "bc2Len :", "bc1Mis :", "bc2Mis :", "cores :", "compress :", 
      "poolreps :", "singleBarcode :",  "readNamePattern :"
    ),
    input_table$Variables
  ),
]

cat("Demultiplex Inputs:\n")
print(
  x = data.frame(input_table, row.names = NULL), 
  right = FALSE, 
  row.names = FALSE
)

# Create output directory if not currently available ----
if( !file.exists(args$outfolder) ){
  
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.\n")

}

# Check for required packages ----
required_packs <- c("stringr", "ShortRead", "Biostrings")
present_packs <- required_packs %in% row.names(installed.packages())

if( !all(present_packs) ){
  
  cat("Missing required r-packages:\n")
  print(
    data.frame(
      "R-Packages" = required_packs, 
      "Installed" = present_packs, 
      row.names = NULL
    ), right = FALSE, row.names = FALSE)

  stop("Check dependancies.\n")

}

# Operating functions ----
parseIndexReads <- function(barcode.seqs, index.file.path, barcode.length,
                                  max.mismatch, max.N.count = 1L, 
                                  read.name.pattern){
  
  # Load index file sequences and sequence names
  index <- ShortRead::readFastq(index.file.path)
  index <- ShortRead::narrow(index, start = 1, end = barcode.length)
  unique_index_seqs <- unique(ShortRead::sread(index))
  
  # Trim barcode if necessary
  barcode_seqs <- as.character(
    Biostrings::DNAStringSet(
      unique(barcode.seqs), 
      start = 1, 
      end = barcode.length
    )
  )
  
  # Identify read names with sequences above or equal to the minscore
  bc_to_unique_idxs <- lapply(
    barcode_seqs, 
    function(x){
      
      vmp <- Biostrings::vmatchPattern(
        pattern = x,
        subject = unique_index_seqs, 
        max.mismatch = max.mismatch, 
        fixed = FALSE
      )
      
      which(lengths(vmp) == 1)
      
    }
  )
  
  # Lookup frame to match barcode sequences to index sequences
  # Sequence variability accounted for and ambiguous, degenerate, and unassigned
  # sequences identified
  degenerate_idxs <- which(
    stringr::str_count(unique_index_seqs, "N") > max.N.count
  )
  
  ambiguous_idxs <- as.numeric(names(table(unlist(bc_to_unique_idxs)))[
    table(unlist(bc_to_unique_idxs)) > 1
  ])
  
  ambiguous_idxs <- ambiguous_idxs[!ambiguous_idxs %in% degenerate_idxs]
  
  unassigned_idxs <- seq_along(unique_index_seqs)[
    !seq_along(unique_index_seqs) %in% unlist(bc_to_unique_idxs)
  ]
  
  unassigned_idxs <- unassigned_idxs[!unassigned_idxs %in% degenerate_idxs]
  
  bc_to_unique_idxs <- lapply(bc_to_unique_idxs, function(x){
    x[!x %in% c(ambiguous_idxs, unassigned_idxs, degenerate_idxs)]
  })
  
  lookup_frame <- data.frame(
    bc_seqs = factor(S4Vectors::Rle(
      values = c(unique(barcode.seqs), "ambiguous", "degenerate", "unassigned"),
      lengths = c(
        lengths(bc_to_unique_idxs), length(ambiguous_idxs), 
        length(degenerate_idxs), length(unassigned_idxs)
      )
    ), 
    levels = c(unique(barcode.seqs), "ambiguous", "degenerate", "unassigned")
    ),
    index_seqs = unique_index_seqs[
      c(unlist(bc_to_unique_idxs), ambiguous_idxs, 
        degenerate_idxs, unassigned_idxs)
      ]
  )
  
  return(split(
    seq_along(index), 
    lookup_frame$bc_seqs[
      match(as.character(ShortRead::sread(index)), lookup_frame$index_seqs)
      ]
  ))
  
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
      
      cat(
        paste0("Wrote ", length(reads[[i]]), " reads to:\n  ", filePath, ".\n"))
      
    },
    reads = reads, type = type, out.folder = out.folder, compress = compress
  )
  
  return(list(read.file.path, type, out.folder))
  
}

# Load manifest / sample mapping file ----
file_ext <- unlist(strsplit(args$manifest, "\\."))
file_ext <- file_ext[length(file_ext)]

if( file_ext %in% c("yaml", "yml") ){
  
  if( !"yaml" %in% row.names(installed.packages()) ){
    stop("Package:yaml not loaded or installed.\n")
  }
  
  manifest <- yaml::yaml.load_file(args$manifest)
  
  if( args$singleBarcode ){
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      row.names = NULL
    )
    
  }else{
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      "bc2" = sapply( manifest$samples, function(x) x[args$bc2Man] ),
      row.names = NULL
    )
    
  }

}else{
  
  if( file_ext == "csv" ){
    manifest <- read.csv(args$manifest)
  }else if( file_ext == "tsv" ){
    manifest <- read.delim(args$manifest)
  }
  
  if( args$singleBarcode ){
    samples_df <- manifest[, c("sampleName", args$bc1Man)]
    names(samples_df) <- c("sampleName", "bc1")
  }else{
    samples_df <- manifest[, c("sampleName", args$bc1Man, args$bc2Man)]
    names(samples_df) <- c("sampleName", "bc1", "bc2")
  }
  
}


if( !args$singleBarcode ){
  
  unique_samples <- nrow(samples_df[,c("bc1", "bc2")]) == 
    nrow(unique(samples_df[,c("bc1", "bc2")]))
  
  if( !unique_samples ) stop("Ambiguous barcoding of samples. Please correct.\n")

}else{

  unique_samples <- length(samples_df[,c("bc1")]) == 
    length(unique(samples_df[,"bc1"]))
  
  if( !unique_samples ) stop("Ambiguous barcoding of samples. Please correct.\n")

}

# Read in barcode sequences ----
bc1_reads <- ShortRead::readFastq(demulti$path[demulti$bc1])
all_indices <- seq_along(bc1_reads)
cat(paste("\nReads to demultiplex : ", length(bc1_reads), "\n"))

if( args$cores > 1 ){
  
  cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))
  
  BC1_parsed <-  parallel::parLapply(
    cluster,
    unique(samples_df$bc1), 
    parseIndexReads,
    index.file.path = demulti$path[demulti$bc1],
    barcode.length = args$bc1Len,
    max.mismatch = args$bc1Mis,
    max.N.count = args$maxN,
    read.name.pattern = args$readNamePattern
  )
  
  names(BC1_parsed) <- unique(samples_df$bc1)
  
  cat("\nbc1 breakdown:\n")
  print(
    data.frame(
      "bc1" = names(BC1_parsed),
      "Read Counts" = sapply( BC1_parsed, length )
    ),
    right = TRUE, 
    row.names = FALSE
  )
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- parallel::parLapply(
      cluster, 
      unique(samples_df$bc2), 
      parseIndexReads,
      index.file.path = demulti$path[demulti$bc2],
      barcode.length = args$bc2Len,
      max.mismatch = args$bc2Mis,
      max.N.count = args$maxN,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
  parallel::stopCluster(cluster)
  
}else{
  
  BC1_parsed <-  parseIndexReads(
    barcode.seqs = samples_df$bc1, 
    index.file.path = demulti$path[demulti$bc1],
    barcode.length = args$bc1Len,
    max.mismatch = args$bc1Mis,
    max.N.count = args$maxN,
    read.name.pattern = args$readNamePattern
  )

  cat("\nbc1 breakdown:\n")
  print(
    data.frame(
      "bc1" = names(BC1_parsed),
      "Read Counts" = lengths(BC1_parsed)
    ),
    right = TRUE, 
    row.names = FALSE
  )
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- parseIndexReads(
      barcode.seqs = samples_df$bc2, 
      index.file.path = demulti$path[demulti$bc2],
      barcode.length = args$bc2Len,
      max.mismatch = args$bc2Mis,
      max.N.count = args$maxN,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
}

if( !args$singleBarcode ){

  cat("\nbc2 breakdown:\n")
  print(
    data.frame(
      "bc2" = names(BC2_parsed),
      "Read Counts" = lengths(BC2_parsed)
    ),
    right = TRUE,
    row.names = FALSE
  )

}

if( !args$singleBarcode ){
  
  demultiplexed_indices <- mapply(
    function(bc1, bc2){
      Biostrings::intersect(BC1_parsed[[bc1]], BC2_parsed[[bc2]])
    },
    bc1 = samples_df$bc1,
    bc2 = samples_df$bc2,
    SIMPLIFY = FALSE
  )
  
  names(demultiplexed_indices) <- paste0(
    samples_df$bc1, samples_df$bc2
  )
  
}else{
  
  demultiplexed_indices <- BC1_parsed[samples_df$bc1]
  
}

# As there is some flexibility in the barcode matching, some reads may be 
# be assigned to multiple samples (ambiguous). Additionally, uncalled bases can
# lead to degenerate sequences (a cause of ambiguous matching), or many 
# sequences will be unassigned.
if( !args$singleBarcode ){
  
  degenerate_indices <- unique(c(BC1_parsed$degenerate, BC2_parsed$degenerate))

  ambiguous_indices <- unique(c(BC1_parsed$ambiguous, BC2_parsed$ambiguous))
  
  ambiguous_indices <- ambiguous_indices[
    !ambiguous_indices %in% degenerate_indices
  ]
  
  unassigned_indices <- unique(c(BC1_parsed$unassigned, BC2_parsed$unassigned))

  unassigned_indices <- unassigned_indices[
    !unassigned_indices %in% c(degenerate_indices, ambiguous_indices)
  ]
  
  demultiplexed_indices <- lapply(demultiplexed_indices, function(x){
    x[!x %in% c(unassigned_indices, ambiguous_indices, degenerate_indices)]
  })
  
  unassigned_indices <- c(unassigned_indices, all_indices[
    !all_indices %in% c(
      unlist(demultiplexed_indices), degenerate_indices, 
      ambiguous_indices, unassigned_indices
    )
  ])
  
}else{
  
  degenerate_indices <- BC1_parsed$degenerate
  ambiguous_indices <- BC1_parsed$ambiguous
  unassigned_indices <- BC1_parsed$unassigned

}

# Reads by sample
samples_df$read_counts <- lengths(demultiplexed_indices)
cat("\nRead counts for each sample.\n")
print(samples_df, split.tables = Inf)

# Ambiguous reads
cat(paste0("\nAmbiguous reads: ", length(ambiguous_indices), "\n"))

# Degenerate reads
cat(paste0("Degenerate reads: ", length(degenerate_indices), "\n"))

# Unassigned reads
cat(paste0("Unassigned reads: ", length(unassigned_indices), "\n"))

if( args$stat != FALSE ){
  write.table(
    data.frame(
      sampleName = paste0(
        c(
          samples_df$sampleName, "ambiguous_reads", 
          "degenerate_reads", "unassigned_reads"
        ), 
        ".demulti"
      ),
      metric = "reads",
      count = c(
        samples_df$read_counts, 
        length(ambiguous_indices), 
        length(degenerate_indices),
        length(unassigned_indices)
      )
    ),
    file = file.path(args$outfolder, args$stat),
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

# Create multiplex dataframe for subseting sequencing files ----
multiplexed_data <- data.frame(
  "sampleName" = S4Vectors::Rle(
    values = samples_df$sampleName, 
    length = lengths(demultiplexed_indices)
  ),
  "index" = unlist(demultiplexed_indices),
  row.names = NULL
)

ambiguous_data <- data.frame(
  "sampleName" = rep("ambiguous", length(ambiguous_indices)),
  "index" = ambiguous_indices,
  row.names = NULL
)

degenerate_data <- data.frame(
  "sampleName" = rep("degenerate", length(degenerate_indices)),
  "index" = degenerate_indices,
  row.names = NULL
)

unassigned_data <- data.frame(
  "sampleName" = rep("unassigned", length(unassigned_indices)),
  "index" = unassigned_indices,
  row.names = NULL
)

multiplexed_data <- rbind(
  multiplexed_data, ambiguous_data, degenerate_data, unassigned_data
)

multiplexed_data$sampleName <- factor(
  multiplexed_data$sampleName,
  levels = c(samples_df$sampleName, "ambiguous", "degenerate", "unassigned")
)

stopifnot( all(multiplexed_data$index %in% all_indices) )

if( args$poolreps ){
  multiplexed_data$sampleName <- gsub("-\\d+$", "", multiplexed_data$sampleName)
}

cat(paste0("Reads to be written to files: ", nrow(multiplexed_data), "\n"))

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

cat("Demultiplexing complete.\n")
q()
