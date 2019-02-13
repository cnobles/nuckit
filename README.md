# NucKit

[![CircleCI](https://circleci.com/gh/cnobles/nuckit.svg?style=svg)](https://circleci.com/gh/cnobles/nuckit)
[![Documentation Status](https://readthedocs.org/projects/nuckit/badge/?version=latest)](http://nuckit.readthedocs.io/en/latest/?badge=latest)

Welcome to the Nucleotide sequence analysis took kit, a set of command line utilities for managing larger volumes of nucleotide sequences all based in the R-programming language. These utilities take advantange of BioConductoR packages, including Biostrings, ShortRead, and GenomicRanges.

For documentation on the utilities, please refer to [ReadTheDocs.io](https://nuckit.readthedocs.io/en/latest/index.html).

For quick installation, follow the commands below on the command line of a linux OS:

```
git clone https://github.com/cnobles/nuckit.git
cd nuckit
bash install.sh -b
```

An overview of the utilities included can be found using the **[-h]** flag:

```
conda activate nuckit
nuckit -h
usage: nuckit [-h/--help, -v/--version] <subcommand>

NucKit - Nucleotide sequence analysis took kit!

subcommands:
  demulti       Demultiplex Illumina sequence files given various barcoding schemes.
  trim          Trim leading and/or over-reading sequences from sequence files.
  filt          Filter sequences in one or more sequence files by various criteria.
  consol        Consolidate reads to unique sequences and generate a key file.
  couple        Couple independent BLAT alignments from paired-end sequences.

optional arguments:
  --nuckit_dir NUCKIT_DIR
                        Path to NucKit installation. (default:/home/nobles/nuckit)
  -v, --version         show program's version number and exit

For more help, see the docs at http://nuckit.readthedocs.io.
```

