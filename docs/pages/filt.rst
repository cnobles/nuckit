.. _filt:

.. contents::
   :depth: 3

Filtering Sequences
===================

Sequence files may sometimes contain only a few reads relevant to an analysis, or host a number of reads needing to be removed. Removal of unwanted sequences, or selection of desired sequences, can help with reducing downstream analysis and requirements. The **filt** utility is designed to do both positive and negative selection by indices (matching sequence names) for fasta and fastq files. Further, if two fasta or fastq files are provided (such as from paired-end sequencing, R1 and R2) then only the intersect (positive selection) or unique (negative selection) reads will be selected. Sequences can also be filtered based on sequence content using additional parameters.

Usage
-----

Below are a number of different filtering processes that **seqFiltR** can be used for:

Filter sequences based on sequence names or indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a file containing indices matching the sequence names within the provided fastq or fasta file, this utility can preform either positive or negative (**[-n]**) filtering and return the filtered sequences as a fastq or fasta file, based on output file extension. If an output file is not provided, the filtered sequences will be printed to the terminal in the input file format.

.. code-block:: shell

  # Select all sequence that match the sequence names in index.txt
  nuckit filt test.fastq -o test.filt.fastq -i index.txt
  
  # Select all sequences that don't match the sequence names in index.txt
  nuckit filt test.fastq -o test.filt.fastq -i index.txt -n


Filter two fast(a/q) files for intersecting or unique read indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If two fasta or fastq files are provided, then the utility will return two output files of the desired format containing the intersecting or unique (**[-n]** when two files are provided) sequences in each file. This process requires sequence names to match between both files, but can be controlled by the **[--readNamePattern]** argument, which will be a regular expression applied to sequence names to select specific parts of the sequence name (see **Arguments** section for default). 

.. code-block:: shell

  # Filter test.R1 and test.R2 files for sequence names matching between the two files
  nuckit filt test.R1.fastq test.R2.fastq -o test.R1.filt.fastq test.R2.filt.fastq
  
  # Filter test.R1 and test.R2 files for sequence names that do not appear in the other file
  nuckit filt test.R1.fastq test.R2.fastq -o test.R1.filt.fastq test.R2.filt.fastq -n 


Filter a fast(a/q) file for sequences containing a specific nucleotide sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Providing a nucleotide sequence (DNA or RNA, ambiguous nucleotide codes accepted) and a fastq or fasta file will return an output file, of the specified format, with sequences which contained or did not contain (**[-n]**) sequences matching the provided nucleotide sequence. The **[-m]** parameter can be used to specify the mismatches tolerated for matching.

.. code-block:: shell

  # Filter test.fastq sequences for sequences containing ACTC...CTAC, or with up to 3 mismatches
  nuckit filt test.fastq -o test.filt.fastq -s ACTCTACGGCATTAGGCTAC
  nuckit filt test.fastq -o test.filt.fastq -s ACTCTACGGCATTAGGCTAC -m 3
  
  # Filter test.fastq sequence for sequences not containing ACTC...CTAC, or up to 3 mismatches
  nuckit filt test.fastq -o test.filt.fastq -s ACTCTACGGCATTAGGCTAC -n 
  nuckit filt test.fastq -o test.filt.fastq -s ACTCTACGGCATTAGGCTAC -m 3 -n


Arguments
---------

**[seqFile(s)]** Sequence file(s) to filter. Can be fasta or fastq formats.

**[-h, --help]** Help information regarding input format and arguments available.

**[-o, --output]** File name(s) of output file(s). Can specify different format from input if desired. Must provide same number of output file names as input files.

**[-i, --index]** A single file containing sequence names to filter for. Each sequence name must match the input file sequence names after extraction of the regex provided by **[--readNamePattern]**.

**[--header]** Script assumes no header is given in index file by default. This option changes that behavior to include a header. 

**[-n, --negSelect]** Specify if the filtering process should be negative, positive by default. Negative selection will only return sequences that do not match to input criteria.

**[-s, --seq]** Filter reads by input nucleotide sequence. DNA, RNA, and ambiguous nucleotide sequences allowed.

**[-m, --mismatch]** Allowed number of mismatches for sequence matching (**[-s]**). Ignored if **[-s, --seq]** is not provided. Default is 0 mismatches.

**[--any]** If multiple methods of filtering should be used, sequences in output files will need to meet all criteria by default. Using this option will allow for sequences passing any of the criteria to be returned.

**[--readNamePattern]** Regex pattern applied to sequence names prior to any matching. Also applied to sequence names from index files. Default pattern: '[\w:-]+'

**[--compress]** Output fast(a/q) files are gzip compressed.

**[-c, --cores]** Number of maximum cores to parallel the processing during certain steps.

**[-q, --quiet]** Silences any log outputs. Will still return sequence output file contents if not given an output option.

**[--stat]** File name of output file for script stats. Output formats of .csv or .tsv are compatible. Stats will still appear in log output, if not silenced by [-q, --quiet].


Dependencies
------------

This utility relies on several R-packages that need to be installed prior to use:

* ShortRead
* stringr
* argparse
* yaml
* data.table
* parallel (if multicore processing is desired)
