.. _quickstart:

.. contents::
   :depth: 2


An Example Run
==============

Below is an example that uses the test dataset provided with NucKit. The following commands lead the user through each utility to go from unprocessed paired-end sequence data at the read level to identified genomic locations. This type of work flow has several uses in bioinformatics, including the monitoring of foreign DNA integrated into a host geneome. The sequence reads themselves represent the product of an amplification reaction to target the foreign DNA and the flanking host geneomic sequence.

.. code-block:: shell

  conda activate nuckit
  cd {path/to/nuckit}

  nuckit demulti -m etc/test.sampleinfo.csv \
      --read1 etc/test_data/Undetermined_S0_L001_R1_001.fastq.gz \
      --read2 etc/test_data/Undetermined_S0_L001_R2_001.fastq.gz \
      --idx1 etc/test_data/Undetermined_S0_L001_I1_001.fastq.gz \
      --idx2 etc/test_data/Undetermined_S0_L001_I2_001.fastq.gz \
      -o etc/test_output \
      --compress
  
  nuckit trim etc/test_output/testSeq-1.R2.fastq.gz \
      -o etc/test_output/testSeq-1.R2.trim.fastq.gz \
      -l ACATATGACAACTCAATTAAACGCGAGC --leadMismatch 3 \
      -r AGATCGGAAGAGCGTCGTGT --overMismatch 4 --overMaxLength 20 \
      --compress
  
  nuckit trim etc/test_output/testSeq-1.R1.fastq.gz \
      -o etc/test_output/testSeq-1.R1.trim.fastq.gz \
      -r GCTCGCGTTTAATTGAGTTGTCATATGT --overMismatch 4 --overMaxLength 20 \
      --compress
  
  nuckit filt etc/test_output/testSeq-1.R2.trim.fastq.gz etc/test_output/testSeq-1.R1.trim.fastq.gz \
      -o etc/test_output/testSeq-1.R2.filt.fastq etc/test_output/testSeq-1.R1.filt.fastq \
      --compress
      
  nuckit consol etc/test_output/testSeq-1.R2.filt.fastq.gz \
      -o etc/test_output/testSeq-1.R2.consol.fasta \
      -k etc/test_output/testSeq-1.R2.key.csv \
      -l testSeq1. \
      --compress
  
  nuckit consol etc/test_output/testSeq-1.R1.filt.fastq.gz \
      -o etc/test_output/testSeq-1.R1.consol.fasta \
      -k etc/test_output/testSeq-1.R1.key.csv \
      -l testSeq1. \
      --compress
  
  # Make sure you have 'blat' installed and a 2bit copy of the hg38 reference genome.
  blat hg38.2bit etc/test_output/testSeq-1.R2.consol.fasta etc/test_output/testSeq-1.R2.psl \
      -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 \
      -minScore=27 -dots=1000 -out=psl -noHead
  
  gzip etc/test_output/testSeq-1.R2.psl
  
  blat hg38.2bit etc/test_output/testSeq-1.R1.consol.fasta etc/test_output/testSeq-1.R1.psl\
      -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 \
      -minScore=27 -dots=1000 -out=psl -noHead
  
  gzip etc/test_output/testSeq-1.R1.psl
  
  nuckit couple etc/test_data/testSeq-1.R2.psl.gz etc/test_data/testSeq-1.R1.psl.gz \
      -k etc/test_output/testSeq-1.R2.key.csv etc/test_output/testSeq-1.R1.key.csv \
      -o etc/test_output/testSeq-1.uniq.csv \
      --condSites etc/test_output/testSeq-1.cond.csv \
      --chimera etc/test_output/testSeq-1.chimera.rds \
      --multihit etc/test_output/testSeq-1.multihit.rds \
      --refGenome hg38 

  conda deactivate