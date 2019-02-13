.. _trim:

.. contents::
   :depth: 3


Sequence Trimming
=================

Trims 5' and 3' nucleotide sequences from paired-end reads based on designated sequences and/or quality scores.

Trimming diagram:

.. code-block:: shell

           leading trimming                    overreading trimming
           ----------CA                                    TCA-----
           ||||||||||||                                    ||||||||    
  [seq] 5' ----------CAAGTC----------------------------TCCATCA-----3'
  
  [result] AGTC----------------------------TCCA


Random ID collection diagram:

.. code-block:: shell

           leading trimming                    overreading trimming
           ----NNNN--CA                                    TCA-----
           ||||    ||||                                    ||||||||    
  [seq] 5' ----CGCT--CAAGTC----------------------------TCCATCA-----3'
  
  [result] AGTC----------------------------TCCA
  [random] CGCT


Usage
-----

.. code-block:: shell

  nuckit trim seqFile.fastq -o output.fastq -l leadingTrimSequence 
  
  nuckit trim seqFile.fastq -o output.fastq -l leadingTrimSequence \
    -r overReadingTrimSequence --phasing 8 --leadMisMatch 2 --overMisMatch 1 \
    --overMaxLength 20 --collectRandomIDs --compress --cores 10


Trimming sequence flexibility
-----------------------------

The arguments **[leadTrimSeq]** and **[overTrimSeq]** take character string inputs, such as below.

.. code-block:: shell

  # [leadTrimSeq] examples
  ATGCCGTTAGCTATGC	    #Fixed nucleotide sequence
  ATGCCGTTNNNNNNAGCTATGC	    #Random 6 nucleotide barcode embedded in sequence
  ATGCCGNNNNTTAGNNNNCTATGC    #Dual 4 nucleotide barcodes embedded in sequence
  # Note: Random nucleotides from within leading trim sequences containing random 
  #   embeded nucleotides can be collected using the [collectRandomIDs] argument.
  
  # [overTrimSeq] examples
  GCTAACGTAC			  #Short 10 nucleotide fixed sequence
  GCTAACGTACGTTTCAAGCTACGGACATGC    #Longer nucleotide fixed sequence
  GCTAACGTACGTTTCNNNNNNCGGACATGC    #Longer nucleotide sequence with 
  				  #  embedded random sequence
  # Note: Longer overreading sequences can take longer to process, alternatively
  #   using the argument [overMaxLength] will limit the number of nucleotides
  #   from the beginning to improve performance. Using this argument does not 
  #   change the allowed mismatch parameter, which will always be based on the
  #   full length sequence.

Using the argument **[collectRandomIDs]** will collect the random sequences from within the leading trimming sequence whenever a match is made. Random sequences are returned in the order in which they appear from 5' to 3' (left to right). Multiple output file names can be given to name each random sequence file captured after the **[collectRandomIDs]** flag. 


Quality trimming
----------------

By default, the trimming utility will quality trim from the left to right in the reads provided. This feature is controlled by the arguments **[badQualBases]**, **[qualSlidingWindow]**, and **[qualThreshold]**. This feature can be disabled by passing the flag "**[--noQualTrimming]**". **[badQualBases]** is an interger value which specifies how many bases need to have quality scores below the threshold, set by **[qualThreshold]**, within the window (**[qualSlidingWindow]**) before the sequence is trimmed. This process is conducted by ShortRead::trimTailw().


Arguments
---------

**[seqFile]** Sequence file to trim, either fasta or fastq (gzip compression tolerated).

**[-h, --help]** Show help message and exit.

**[-o,--output]** Output file name.

**[-l, --leadTrimSeq]** Sequence to trim from 5' end of reads, or the leading sequence. See above for sequence flexibility.

**[-r, --overTrimSeq]** Sequence to trim from 3' end of reads, or the overreading sequence. See above for sequence flexibility.

**[--phasing]** Number of nucleotides to remove from 5' end of sequence before trimming. Default = 0.

**[--maxMisMatch]** Maximum allowable mismatches in leading or overreading trim sequences.

**[--leadMisMatch]** Maximum allowable mismatches in leading trim sequence. Default = 0.

**[--overMisMatch]** Maximum allowable mismatches in overreading trim sequence. Default = 0.

**[--overMaxLength]** Maximum length to consider of the overTrimSeq to use for alignments. See above for in-depth explanation of this feature. Default 20 nts.

**[--overMinLength]** Minimum length to consider of the overTrimSeq to use for alignments. See README for in depth explanation of this feature. Default 3 nts.

**[--minSeqLength]** Minimum length of trimmed sequence. Any trimmed sequence with a length below this value will be filtered out. Default = 30.

**[--collectRandomIDs]** Option to collect random nucleotide sequences from trimmed portions. If used, provide an output file name.

**[--noFiltering]** Will not filter reads based on leadTrimSeq, the default behavior.

**[--noQualTrimming]** Will not quality trim reads, the default behavior.

**[--badQualBases]** Number of bases below threshold in sliding window before read is trimmed. Default = 5.

**[--qualSlidingWindow]** Slinding window size for which to assess quality scores below threshold. Default = 10.

**[--qualThreshold]** Quality threshold for trimming, minimum allowable score. Default = '?', Q30.

**[--stat]** File name to be written in output directory of read couts for each sample. CSV file format. ie. test.stat.csv.

**[--compress]** Output fastq/fasta files are gzipped.

**[-c, --cores]** Max cores to be used. If 0 (default), program will not utilize parallel processing.


Dependencies
------------

This utility is coded in R, and was developed on v3.2.2, though it should run with earlier versions given the appropriate dependencies. The script uses 9 additional packages:
  * argparse
  * yaml
  * stringr
  * ShortRead
  * Biostrings
  * parallel (if multicore processing is desired)
