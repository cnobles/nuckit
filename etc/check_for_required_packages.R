# Check installed packages for seqFiltR dependencies
packs <- c(
  "ShortRead", "Biostrings", "GenomicRanges", "IRanges",
  "S4Vectors", "igraph", "Matrix", "data.table", "stringr", 
  "parallel", "yaml", "argparse"
)

present <- packs %in% row.names(installed.packages())
print(data.frame(row.names = packs, "Installed" = present))
stopifnot(all(present))
q()
