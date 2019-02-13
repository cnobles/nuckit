.. _consol:

.. contents::
   :depth: 3


Consolidate Sequences
=====================

Consolidate nucleotide sequences to unique sequences with a key. Pipelines and analyses work much smoother if they don't have to repeat the same process over and over. Therefore, reducing the sequences within a sequence file to only the unique occurances can help improve performance without loosing any information.

Usage
-----

.. code-block:: shell

  nuckit consol test.fastq -o test.consol.fasta -k test.key.csv -l test.seq. 
  
  
This will generate a sequence file in fasta format partnered with a key file in a csv format. The key file will contain two columns, "readNames" and "seqID". The former matching the input file read identifiers and the latter indicating the unique identifier in the output file.


Arguments
---------

**[seqFile]** Sequence file to consolidate to unique sequences, either fasta or fastq format.

**[-h, --help]** Help information regarding input format and arguments available.

**[-o, --output]** Output fasta file name. Ex. sample.consolidated.fasta

**[-k, --keyFile]** Key file output name. Ex. sample.r1.csv.

**[-l, --seqName]** Name to append to unique sequences. Ex. sample.r1. . 

**[--stat]** File name for stat file output (CSV format).

**[--compress]** Output fasta file is gzip compressed.


Dependencies
------------

This utility relies on several R-packages that need to be installed prior to use:

* ShortRead
* Biostrings
* stringr
* argparse
* yaml
