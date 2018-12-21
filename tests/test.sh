#!/usr/bin/env bash
set -ev

# Empty test directory of test output files
rm -rf tests/test_output

# Test for required packages
Rscript etc/check_for_required_packages.R


### Demulti Tests ###

# Test functionality of dualDemultiplexR
Rscript scripts/demulti.R -m tests/sampleInfo.tsv \
  --read1 tests/Data/Undetermined_S0_L001_R1_001.fastq.gz \
  --read2 tests/Data/Undetermined_S0_L001_R2_001.fastq.gz \
  --idx1 tests/Data/Undetermined_S0_L001_I1_001.fastq.gz \
  --idx2 tests/Data/Undetermined_S0_L001_I2_001.fastq.gz \
  -o tests/test_output --stat test.stat.csv --compress

# Check output for correct findings, 50 testSeq-1, 50 testSeq-2, and 0 testSeq-3
test_1_len=$(zcat tests/test_output/testSeq-1.R2.fastq.gz | sed -n '2~4p' | wc -l)
test_2_len=$(zcat tests/test_output/testSeq-2.I1.fastq.gz | sed -n '2~4p' | wc -l)
test_3_len=$(zcat tests/test_output/testSeq-3.R1.fastq.gz | sed -n '2~4p' | wc -l)

if [ ! $test_1_len = 50 ] | [ ! $test_2_len = 50 ] | [ ! $test_3_len = 0 ]; then
    exit 1
fi

# Show stat output.
cat tests/test_output/test.stat.csv


### Trim Tests ###

# Test for leading and overreading trimming on R2 test sequences
Rscript scripts/trim.R tests/test_output/testSeq-1.R2.fastq.gz \
  -o tests/test_output/testSeq-1.R2.trim.fastq.gz \
  -l ACATATGACAACTCAATTAAACGCGAGC --leadMismatch 3 \
  -r AGATCGGAAGAGCGTCGTGT --overMismatch 4 --overMaxLength 20 \
  --stat tests/test_output/test.R2.trim.stat.csv --compress

# Test for only overreading trimming on R1 test sequences
Rscript scripts/trim.R tests/test_output/testSeq-1.R1.fastq.gz \
  -o tests/test_output/testSeq-1.R1.trim.fastq.gz \
  -r GCTCGCGTTTAATTGAGTTGTCATATGT --overMismatch 4 --overMaxLength 20 \
  --stat tests/test_output/test.R1.trim.stat.csv --compress

# Check outputs for correct findings
test_R2_trim_len=$(zcat tests/test_output/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | wc -l)
test_R1_trim_len=$(zcat tests/test_output/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $test_R2_trim_len = 50 ] | [ ! $test_R1_trim_len = 50 ]; then
    exit 1
fi

# R2 test sequences
zcat tests/test_output/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | head -n 5

# R1 test sequences
zcat tests/test_output/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | head -n 5

# Concatenate the stat files
cat tests/test_output/*.trim.stat.csv


### Filt Tests ###

# Test for sequence filtering (positive and negative selection)
Rscript scripts/filt.R tests/Data/Undetermined_S0_L001_I1_001.fastq.gz \
  -o tests/test_output/seq_filt_I1.fastq \
  -s CGTACTAG --stat tests/test_output/seq_filt_I1.stat.csv --compress
filt_test_len=$(zcat tests/test_output/seq_filt_I1.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $filt_test_len = 50 ]; then
    exit 1
fi

Rscript scripts/filt.R tests/Data/Undetermined_S0_L001_I1_001.fastq.gz \
  -o tests/test_output/seq_negfilt_I1.fastq \
  -s CGTACTAG --stat tests/test_output/seq_negfilt_I1.stat.csv --compress -n
negfilt_test_len=$(zcat tests/test_output/seq_negfilt_I1.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $negfilt_test_len = 50 ]; then
    exit 1
fi

# Test for sequence filtering by overlapping indices
Rscript scripts/filt.R tests/test_output/testSeq-1.R2.trim.fastq.gz tests/test_output/testSeq-1.R1.trim.fastq.gz \
  -o tests/test_output/testSeq-1.R2.filt.fastq tests/test_output/testSeq-1.R1.filt.fastq \
  --stat tests/test_output/seq_filt_I1_compare.stat.csv --compress
ovlp_A_test_len=$(zcat tests/test_output/seq_A_overlap_I1.fastq.gz | sed -n '2~4p' | wc -l)
ovlp_B_test_len=$(zcat tests/test_output/seq_B_overlap_I1.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $ovlp_A_test_len = 0 ] | [ ! $ovlp_B_test_len = 0 ]; then
    exit 1
fi

# Test for index filtering
zcat tests/test_output/seq_filt_I1.fastq.gz | sed -n '1~4p' > tests/test_output/seq_index.txt
head tests/test_output/seq_index.txt -n 3

Rscript scripts/filt.R tests/Data/Undetermined_S0_L001_R1_001.fastq.gz -o tests/test_output/seq_filt_R1.fastq \
  -i tests/test_output/seq_index.txt --stat tests/test_output/seq_filt_R1.stat.csv --compress
idx_test_len=$(zcat tests/test_output/seq_filt_R1.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $idx_test_len = 50 ]; then
    exit 1
fi

# Stats of each step
cat tests/test_output/*filt*.stat.csv


### Consol Tests ###

# Test script functionality
Rscript scripts/consol.R tests/test_output/testSeq-1.R2.filt.fastq.gz \
    -o tests/test_output/testSeq-1.R2.consol.fasta -k tests/test_output/testSeq-1.R2.key.csv \
    -l testSeq1. --stat tests/test_output/testSeq-1.R2.consol.stat.csv --compress

Rscript scripts/consol.R tests/test_output/testSeq-1.R1.filt.fastq.gz \
    -o tests/test_output/testSeq-1.R1.consol.fasta -k tests/test_output/testSeq-1.R1.key.csv \
    -l testSeq1. --stat tests/test_output/testSeq-1.R1.consol.stat.csv --compress

# Check output for correct processing, 50 keys, 25 unique sequences
test_consol_key_len=$(cat tests/test_output/testSeq-1.R2.key.csv | sed '/readNames/d' | wc -l)
test_consol_seq_len=$(zcat tests/test_output/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | wc -l)

if [ ! $test_consol_key_len = 50 ] | [ ! $test_consol_seq_len = 25 ]; then
    exit 1
fi

# Beginning sequences from test output
zcat tests/test_output/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | head -n 5

# Beginning of key file for sequences
cat tests/test_output/testSeq-1.R2.key.csv | head -n 6

# Stat file
cat tests/test_output/*.consol.stat.csv


### Couple Tests ###

# Test script for functionality
Rscript scripts/couple.R tests/Data/testSeq-1.R2.psl.gz tests/Data/testSeq-1.R1.psl.gz \
  -k tests/test_output/testSeq-1.R2.key.csv tests/test_output/testSeq-1.R1.key.csv \
  -o tests/test_output/testSeq-1.uniq.csv --condSites tests/test_output/testSeq-1.cond.csv \
  --chimera tests/test_output/testSeq-1.chimera.rds \
  --multihit tests/test_output/testSeq-1.multihit.rds \
  --refGenome hg38 --stat tests/test_output/testSeq-1.couple.stat.csv

# Check output for correct findings, 50 reads / alignments and 5 sites
test_couple_uniq_len=$(cat tests/test_output/testSeq-1.uniq.csv | sed '/seqnames/d' | wc -l)
test_couple_cond_len=$(cat tests/test_output/testSeq-1.cond.csv | sed '/seqnames/d' | wc -l)

if [ ! $test_couple_uniq_len = 50 ] | [ ! $test_couple_cond_len = 5 ]; then
    exit 1
fi

head tests/test_output/testSeq-1.uniq.csv

head tests/test_output/testSeq-1.cond.csv

cat tests/test_output/testSeq-1.couple.stat.csv


### Clean up ###

rm -r tests/test_output

echo "Passed all tests."

