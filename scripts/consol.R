#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
suppressMessages(library("stringr"))

code_dir <- dirname(sub("--file=", "", grep(
  "--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))

desc <- desc <- yaml::yaml.load_file(file.path(code_dir, "descriptions.yml"))

#' Set up and gather command line arguments
parser <- ArgumentParser(description = desc$program_short_description)
parser$add_argument(
  "seqFile", nargs = 1, type = "character", default = "NA",
  help = desc$seqFile)
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = "NA",
  help = desc$output)
parser$add_argument(
  "-k", "--keyFile", nargs = 1, type = "character", default = "NA",
  help = desc$keyFile)
parser$add_argument(
  "-l", "--seqName", nargs = 1, type = "character", default = "NA",
  help = desc$seqName)
parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, help = desc$stat)
parser$add_argument(
  "--compress", action = "store_true", help = desc$compress) 

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if(args$seqFile == "NA") stop("No sequence file specified. Please provide.")

# Check I/O file types
seqType <- unlist(strsplit(args$seqFile, "/"))
seqType <- seqType[length(seqType)]
seqType <- str_extract(seqType, "fa[\\w]*")
if(!seqType %in% c("fa", "fasta", "fastq")){
  stop(desc$unrecognized_file_type, " ", desc$compression_note)
}
seqType <- ifelse(seqType %in% c("fa", "fasta"), "fasta", "fastq")

if(args$output != "NA"){
  outType <- unlist(strsplit(args$output, "/"))
  outType <- outType[length(outType)]
  outType <- str_extract(outType, "fa[\\w]*")
  args$output <- unlist(strsplit(args$output, outType))[1]
  if(!outType %in% c("fa", "fasta", "fastq")){
    stop(desc$unrecognized_file_type)
  }
  outType <- ifelse(outType %in% c("fa", "fasta"), "fasta", "fastq")
  if(outType == "fastq"){
    message(desc$fastq_input)
    outType <- "fasta"
  }
  args$output <- paste0(args$output, outType)
}else{
  stop("No output file name given. Please provide.")
}

if(args$keyFile != "NA"){
  keyType <- str_extract(args$keyFile, "[\\w]+$")
  if(!keyType %in% c("csv", "tsv", "rds", "RData")){
    stop(desc$output_keyfile_type)
  }
}else{
  stop("No key file name given. Please provide.")
}
  
# Check sequence name lead
if(args$seqName == "NA"){
  parsedName <- unlist(strsplit(args$seqFile, "/"))[
    length(unlist(strsplit(args$seqFile, "/")))]
  args$seqName <- unlist(strsplit(parsedName, "fa[\\w]*"))[1]
}

# Print inputs to table
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("seqFile :", "output :", "keyFile :", "seqName :", "stat :"),
        input_table$Variables),]
pandoc.title("seqConsolidateR Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Load additional R-packages
addPacks <- c("ShortRead")
addPacksLoaded <- suppressMessages(
  sapply(addPacks, require, character.only = TRUE))
if(!all(addPacksLoaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(addPacksLoaded), 
    "Loaded" = addPacksLoaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Read sequence file
if(seqType == "fasta"){
  seqPointer <- ShortRead::readFasta(args$seqFile)
}else{
  seqPointer <- ShortRead::readFastq(args$seqFile)
}

seqs <- ShortRead::sread(seqPointer)
names(seqs) <- ShortRead::id(seqPointer)

# Generate blank files if inputs are empty
if(length(seqs) == 0){
  writeXStringSet(
    DNAStringSet(),
    filepath = args$output,
    format = "fasta",
    compress = args$compress)
  
  if(!is.null(args$keyFile)){
    key <- data.frame("readNames" = c(), "seqID" = c())
    if(keyType == "csv"){
      write.csv(key, file = args$keyFile, row.names = FALSE, quote = FALSE)
    }else if(keyType == "tsv"){
      write.table(
        key, file = args$keyFile, sep = "\t", row.names = FALSE, quote = FALSE)
    }else if(keyType == "rds"){
      saveRDS(key, file = args$keyFile)
    }else if(keyType == "RData"){
      save(key, file = args$keyFile)
    }
  }
  
  if(args$stat != FALSE){
    sampleName <- unlist(strsplit(args$output, "/"))
    sampleName <- unlist(
      strsplit(sampleName[length(sampleName)], ".fa", fixed = TRUE))[1]
    write.table(
      data.frame(
        sampleName = sampleName,
        metric = "reads",
        count = length(seqs)),
      file = args$stat,
      sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  q()
}
  

factorSeqs <- factor(as.character(seqs))

key <- data.frame(
  "readNames" = names(factorSeqs),
  "seqID" = paste0(args$seqName, as.integer(factorSeqs))
)

consolidatedSeqs <- DNAStringSet(levels(factorSeqs))
names(consolidatedSeqs) <- paste0(args$seqName, 1:length(levels(factorSeqs)))

if(!is.null(args$output)){
  if(args$compress & !grepl(".gz", args$output)){
    args$output <- paste0(args$output, ".gz")
  }
}

# Stats if requested
if(args$stat != FALSE){
  sampleName <- unlist(strsplit(args$output, "/"))
  sampleName <- unlist(
    strsplit(sampleName[length(sampleName)], ".fa", fixed = TRUE))[1]
  write.table(
    data.frame(
      sampleName = sampleName,
      metric = "reads",
      count = length(consolidatedSeqs)),
    file = args$stat,
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Write output and key files
# Output
if(is.null(args$output)){
  pandoc.table(
    data.frame(
      "seqID" = names(consolidatedSeqs),
      "sequence" = as.character(consolidatedSeqs),
      row.names = NULL),
    style = "simple")
}else{
  unlink(args$output)
  ShortRead::writeFasta(
    consolidatedSeqs, 
    file = args$output,
    width = max(width(consolidatedSeqs)),
    compress = args$compress)
}

# Key file
if(!is.null(args$keyFile)){
  if(keyType == "csv"){
    write.csv(key, file = args$keyFile, row.names = FALSE, quote = FALSE)
  }else if(keyType == "tsv"){
    write.table(
      key, file = args$keyFile, sep = "\t", row.names = FALSE, quote = FALSE)
  }else if(keyType == "rds"){
    saveRDS(key, file = args$keyFile)
  }else if(keyType == "RData"){
    save(key, file = args$keyFile)
  }
}
q()
