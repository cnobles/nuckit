#!/usr/bin/env bash
set -ev

# Test functionality of dualDemultiplexR
Rscript dualDemultiplex.R -m tests/sampleInfo.tsv \
  --read1 tests/Data/Undetermined_S0_L001_R1_001.fastq.gz \
  --read2 tests/Data/Undetermined_S0_L001_R2_001.fastq.gz \
  --index1 tests/Data/Undetermined_S0_L001_I1_001.fastq.gz \
  --index2 tests/Data/Undetermined_S0_L001_I2_001.fastq.gz \
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

# Cleanup test directory
rm -r tests/test_output

echo "Passed all tests."
#exit
