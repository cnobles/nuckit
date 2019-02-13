.. _demulti:

.. contents::
   :depth: 3


Demultiplexing
==============

A tool for demultiplexing reads (paired-end by design, but not required). Indexing can be bases on a wide variety of designs, both single and double barcoding, on any read given (R1, R2, I1, I2). Bins ambiguous and unassigned reads in additional files.

Usage
-----

.. code-block:: shell

  nuckit demulti -m manifest.yml \
    --read1 READ1.fastq --read2 READ2.fastq --index1 INDEX1.fastq --index2 INDEX2.fastq \
    --outFolder ~/demultiplexed
  
  nuckit demulti -m manifest.yml \
    --read1 READ1.fastq --read2 READ2.fastq --index1 INDEX1.fastq --index2 INDEX2.fastq \
    --outFolder ~/demultiplexed --poolreps --maxMismatch 1 --barcode1Length 8 --barcode2Length 8 \
    --readNamePattern [\\w:-]+ --compress --cores 4
  

Sample info / manifest / mapping file format
--------------------------------------------

Manifests contain the barcode information for each sample. The minimal required information are sampleNames, barcode1 and barcode2. Examples show below:

.. code-block:: shell

  # This is an example sample manifest.
  samples :
      treated_sample-1 :
          barcode1 : TAAGGCGA
          barcode2 : TAGATCGC
      
      treated_sample-2 :
          barcode1 : TAAGGCGA
          barcode2 : CGATCCTA
          
      alt_treated_sample-1 :
          barcode1 : NATCGTCA
          barcode2 : NCTGGTAC
        
  # Alternatively, sample manifests can be input as csv files.
  sampleName,barcode1,barcode2
  treated_sample-1,TAAGGCGA,TAGATCGC
  treated_sample-2,TAAGGCGA,CGATCCTA
  alt_treated_sample-1,NATCGTCA,NCTGGTAC
  
  # or tsv files.
  sampleName	barcode1	barcode2
  treated_sample-1	TAAGGCGA	TAGATCGC
  treated_sample-2	TAAGGCGA	CGATCCTA
  alt_treated_sample-1	NATCGTCA	NCTGGTAC	


If replicates are included and are to be pooled, sample names should be in the above format (sampleName-#). Pooled replicates with be consolidated if the "-p" or "--poolreps" flag is used with the "-" as the delimiter between sampleNames and replicate designations. Likewise, sampleNames should not include the "-" symbol if pooling replicates is desired.

Barcode sequences can contain ambiguous nucleotides. It is recommended that the user chooses to use either ambiguous nucleotides or declare the maximal number of mismatches allowed. Using both at the same time may lead to undesired demultiplexing. The current substitution matrix used is a binary matrix based on NUC4.4 (ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4). This matrix is meant for all ambiguous nucleotides [ATGCSWRYKMBVHDN] to be matched to those from typical machine output [ATGCN]. The matrix used (termed binary ambiguous nucleotide scoring matrix or BAN Mat) is not designed for matching between various ambiguous nucletides to eachother, only to [ATGCN]. Both ambiguous nucleotide matching and maxMismatch allowance can be useful in their own way. For poor sequencing quality, ambigous nucleotide demultiplexing may be able to acquire more reads that contained low quality index sequences at specific positions, while maxMismatch may help if phasing/pre-phasing is a major issue in the sequencing file.


Arguments
---------

**[-h, --help]** Help information regarding input format and arguments available.

**[--read1, --read2, --index1, --index2]** File paths to the respective FASTQ files.

**[-o, --outFolder]** Output directory for demultiplexed files.

**[-p, --poolreps]** Pools replicates. Replicate designation delimited by "-". ie. sample-1, sample-2, ...

**[--singleBarcode]** Demultiplexes by only a single barcode, specify with **[--barcode1]**.

**[--barcode1, --barcode2]** Read type containing respecive barcode sequences, defaults are I1 and I2. Options are R1, R2, I1, or I2.

**[--barcode1Length, --barcode2Length]** Length of barcode sequences for demultiplexing. Default is 8 nt.

**[--maxMismatch]** Allowable mismatch in barcode sequences. Overrides **[--bc1Mismatch, --bc2Mismatch]**. 

**[--bc1Mismatch, --bc2Mismatch]** Allowable mismatch in specific barcode alignments.

**[--stat]** File name to be written in output directory of read couts for each sample. CSV file format. ie. test.stat.csv.

**[--readNamePattern]** Regex pattern to capture read names without read-type specific info.

**[--compress]** Output fastq files are gzipped.

**[-c, --cores]** Number of maximum cores to parallel the processing during certain steps.


Dependencies
------------

This demultiplexing utility is coded in R, and was developed on v3.2.2, though it should run with earlier versions given the appropriate dependencies. The script uses 6 additional packages:
  * argparse
  * yaml
  * ShortRead
  * Biostrings
  * stringr
  * parallel (if multicore processing is desired)
