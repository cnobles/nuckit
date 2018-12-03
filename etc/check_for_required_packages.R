# Check installed packages for seqFiltR dependencies
packs <- c(
  "ShortRead", "data.table", "pander", 
  "yaml", "argparse", "stringr", "parallel")

present <- packs %in% row.names(installed.packages())
print(data.frame(row.names = packs, "Installed" = present))
q()