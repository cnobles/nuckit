#!/usr/bin/env bash

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to test. No environment by default."
  -p --path_nuckit  [arg] Location of NucKit source code. Default: this directory
  -n --no_ctrl            Test without NucKit control module.
  -c --conda        [arg] Location of Conda installation. Default: ${PREFIX}
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script tests NucKit, change inputs to modify testing conditions.
EOF


# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

function __err_report() {
    local error_code
    error_code=${?}
    error "Error in ${__file} in function ${1} on line ${2}"
    exit ${error_code}
}
trap '__err_report "${FUNCNAME:-.}" ${LINENO}' ERR

# help mode
if [[ "${arg_h:?}" = "1" ]]; then
    # Help exists with code 1
    help "Help using ${0}:"
fi

# verbose mode
if [[ "${arg_v:?}" = "1" ]]; then
    LOG_LEVEL="7"
fi

# debug mode
if [[ "${arg_d:?}" = "1" ]]; then
    set -o xtrace
    LOG_LEVEL="7"
fi

function debug_capture () {
    debug "$(echo -e "$(${@})")"
}

function testing_error () {
    error "${1} failed!"
    if [[ "${arg_v:?}" != 1 && "${arg_d:?}" != 1 ]]; then
        error "Try re-running with -v or -d, or file an issue on Github."
    fi
    exit 1
}

function activate_env () {
    set +o nounset
    source activate $__nuckit_env
    set -o nounset
}

function deactivate_env () {
    set +o nounset
    source deactivate
    set -o nounset
}


# ------------------------------------------------------------------------------

# Set variables
__conda_path="${arg_c:-${HOME}/miniconda3}"
__nuckit_env="${arg_e:-null}"
__nuckit_dir="${arg_p:-$(pwd)}"
__with_ctrl=true
__req_r_version="3.4.1"

__old_path=$PATH

PATH=$PATH:${__conda_path}/bin

cd $__nuckit_dir

# Test without control module
if [[ "${arg_n:?}" = "1" ]]; then
    __with_ctrl=false
fi

# Set command
function __nuckit_cmd () {
  if [[ $__with_ctrl = true ]]; then
      __cmd="nuckit ${1} ${@:2}"
  else
      __cmd="Rscript scripts/${1}.R ${@:2}"
  fi
  eval "$__cmd"
}


# Test functions
function __test_r_version () {
    local sem_version=$(R --version | grep 'R version' | cut -d ' ' -f 3)
    if (( ${sem_version//./} >= ${__req_r_version//./} )); then
        echo true
    else
        echo false
    fi
}

function __test_r_packages () {
    $(Rscript etc/check_for_required_packages.R > /dev/null && \
        echo true || echo false)
}

function __test_env_vars () {
    $(echo $NUCKIT_DIR | grep nuckit > /dev/null && echo set || echo unset)
}


# Start testing ----
info "Starting testing..."
if [[ $__nuckit_env != null ]]; then
    info "    Conda path    : ${__conda_path}"
    info "    Testing env   : ${__nuckit_env}"
fi
info "    NucKit src    : ${__nuckit_dir}"
info "    Control module: ${__with_ctrl}"

# Empty test directory of test output files
rm -rf etc/test_output

# Activate test environment if supplied
if [[ $__nuckit_env != null ]]; then
    activate_env
fi

# Test for required packages
if [[ $(__test_r_packages) == false ]]; then
    testing_error "Required R-packages not installed."
else
    info "Required R-packages installed."
fi

### Demulti Tests ###

info "Testing demulti ..."

# Test functionality of dualDemultiplexR
demulti_str="-m etc/test.sampleinfo.csv --read1 etc/test_data/Undetermined_S0_L001_R1_001.fastq.gz --read2 etc/test_data/Undetermined_S0_L001_R2_001.fastq.gz --idx1 etc/test_data/Undetermined_S0_L001_I1_001.fastq.gz --idx2 etc/test_data/Undetermined_S0_L001_I2_001.fastq.gz -o etc/test_output --stat test.stat.csv --compress"

debug_capture __nuckit_cmd demulti "$demulti_str" 2>&1

# Check output for correct findings, 50 testSeq-1, 50 testSeq-2, and 0 testSeq-3
test_1_len=$(zcat etc/test_output/testSeq-1.R2.fastq.gz | sed -n '2~4p' | wc -l)
test_2_len=$(zcat etc/test_output/testSeq-2.I1.fastq.gz | sed -n '2~4p' | wc -l)
test_3_len=$(zcat etc/test_output/testSeq-3.R1.fastq.gz | sed -n '2~4p' | wc -l)

if [[ $test_1_len != 50 || $test_2_len != 50 || $test_3_len != 0 ]]; then
    testing_error "    test 1/2: FAILED (output)"
else
    info "    test 1/2: passed (output)"
fi

demulti_stat=$(cat etc/test_output/test.stat.csv > /dev/null && echo true || echo false)

if [[ $demulti_stat = false ]]; then
    testing_error "    test 2/2: FAILED (stats)"
else
    info "    test 2/2: passed (stats)"
fi


### Trim Tests ###

info "Testing trim ..."

# Test for leading and overreading trimming on R2 test sequences
trim_str_1="etc/test_output/testSeq-1.R2.fastq.gz -o etc/test_output/testSeq-1.R2.trim.fastq.gz -l ACATATGACAACTCAATTAAACGCGAGC --leadMismatch 3 -r AGATCGGAAGAGCGTCGTGT --overMismatch 4 --overMaxLength 20 --stat etc/test_output/test.R2.trim.stat.csv --compress"

debug_capture __nuckit_cmd trim "$trim_str_1" 2>&1

# Test for only overreading trimming on R1 test sequences
trim_str_2="etc/test_output/testSeq-1.R1.fastq.gz -o etc/test_output/testSeq-1.R1.trim.fastq.gz -r GCTCGCGTTTAATTGAGTTGTCATATGT --overMismatch 4 --overMaxLength 20 --stat etc/test_output/test.R1.trim.stat.csv --compress"
  
debug_capture __nuckit_cmd trim "$trim_str_2" 2>&1

# Check outputs for correct findings
test_R2_trim_len=$(zcat etc/test_output/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | wc -l)
test_R1_trim_len=$(zcat etc/test_output/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | wc -l)

# R2 test sequences
trim_out_1=$(zcat etc/test_output/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | head -n 5 > \
    /dev/null && echo true || echo false)

# R1 test sequences
trim_out_2=$(zcat etc/test_output/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | head -n 5 > \
    /dev/null && echo true || echo false)

if [[ $test_R2_trim_len != 50 || trim_out_1 = false ]]; then
    testing_error "    test 1/3: FAILED (lead and overreading)"
else
    info "    test 1/3: passed (lead and overreading)"
fi
    
if [[ $test_R1_trim_len != 50 || trim_out_2 = false ]]; then
    testing_error "    test 2/3: FAILED (overreading)"
else
    info "    test 2/3: passed (overreading)"
fi

# Concatenate the stat files
trim_stat=$(cat etc/test_output/*.trim.stat.csv > /dev/null && echo true || echo false)

if [[ $trim_stat = false ]]; then
    testing_error "    test 3/3: FAILED (stats)"
else
    info "    test 3/3: passed (stats)"
fi


### Filt Tests ###

info "Testing filt ..."

# Test for sequence filtering (positive and negative selection)
filt_str_1="etc/test_data/Undetermined_S0_L001_I1_001.fastq.gz -o etc/test_output/seq_filt_I1.fastq -s CGTACTAG --stat etc/test_output/seq_filt_I1.stat.csv --compress"

debug_capture __nuckit_cmd filt "$filt_str_1" 2>&1

filt_test_len=$(zcat etc/test_output/seq_filt_I1.fastq.gz | sed -n '2~4p' | wc -l)

if [[ $filt_test_len != 50 ]]; then
    testing_error "    test 1/5: FAILED (positive selection)"
else
    info "    test 1/5: passed (positive selection)"
fi


filt_str_2="etc/test_data/Undetermined_S0_L001_I1_001.fastq.gz -o etc/test_output/seq_negfilt_I1.fastq -s CGTACTAG --stat etc/test_output/seq_negfilt_I1.stat.csv --compress -n"

debug_capture __nuckit_cmd filt "$filt_str_2" 2>&1

negfilt_test_len=$(zcat etc/test_output/seq_negfilt_I1.fastq.gz | sed -n '2~4p' | wc -l)

if [[ $negfilt_test_len != 50 ]]; then
    testing_error "    test 2/5: FAILED (negative selection)"
else
    info "    test 2/5: passed (negative selection)"
fi

# Test for sequence filtering by overlapping indices
filt_str_3="etc/test_output/testSeq-1.R2.trim.fastq.gz etc/test_output/testSeq-1.R1.trim.fastq.gz -o etc/test_output/testSeq-1.R2.filt.fastq etc/test_output/testSeq-1.R1.filt.fastq --stat etc/test_output/seq_filt_I1_compare.stat.csv --compress"
  
debug_capture __nuckit_cmd filt "$filt_str_3" 2>&1

ovlp_A_test_len=$(zcat etc/test_output/testSeq-1.R2.filt.fastq.gz | sed -n '2~4p' | wc -l)
ovlp_B_test_len=$(zcat etc/test_output/testSeq-1.R1.filt.fastq.gz | sed -n '2~4p' | wc -l)

if [[ $ovlp_A_test_len != 50 || $ovlp_B_test_len != 50 ]]; then
    testing_error "    test 3/5: FAILED (overlap selection)"
else
    info "    test 3/5: passed (overlap selection)"
fi

# Test for index filtering
zcat etc/test_output/seq_filt_I1.fastq.gz | sed -n '1~4p' > etc/test_output/seq_index.txt
create_test_file=$(head etc/test_output/seq_index.txt -n 3 > /dev/null && echo true || echo false)

if [[ $create_test_file = false ]]; then
    testing_error "    Could not create test file: seq_index.txt"
fi

filt_str_4="etc/test_data/Undetermined_S0_L001_R1_001.fastq.gz -o etc/test_output/seq_filt_R1.fastq -i etc/test_output/seq_index.txt --stat etc/test_output/seq_filt_R1.stat.csv --compress"

debug_capture __nuckit_cmd filt "$filt_str_4" 2>&1

idx_test_len=$(zcat etc/test_output/seq_filt_R1.fastq.gz | sed -n '2~4p' | wc -l)

if [[ $idx_test_len != 50 ]]; then
    testing_error "    test 4/5: FAILED (index selection)"
else
    info "    test 4/5: passed (index selection)"
fi

# Stats of each step

filt_stat=$(cat etc/test_output/*filt*.stat.csv > /dev/null && echo true || echo false)

if [[ $filt_stat = false ]]; then
    testing_error "    test 5/5: FAILED (stats)"
else
    info "    test 5/5: passed (stats)"
fi


### Consol Tests ###

info "Testing consol ..."

consol_str_1="etc/test_output/testSeq-1.R2.filt.fastq.gz -o etc/test_output/testSeq-1.R2.consol.fasta -k etc/test_output/testSeq-1.R2.key.csv -l testSeq1. --stat etc/test_output/testSeq-1.R2.consol.stat.csv --compress"

consol_str_2="etc/test_output/testSeq-1.R1.filt.fastq.gz -o etc/test_output/testSeq-1.R1.consol.fasta -k etc/test_output/testSeq-1.R1.key.csv -l testSeq1. --stat etc/test_output/testSeq-1.R1.consol.stat.csv --compress"

# Test script functionality
debug_capture __nuckit_cmd consol "$consol_str_1" 2>&1

debug_capture __nuckit_cmd consol "$consol_str_2" 2>&1

# Beginning sequences from test output
consol_seqs=$(zcat etc/test_output/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | head -n 5)
debug_capture echo $consol_seqs

# Beginning of key file for sequences
key_head=$(cat etc/test_output/testSeq-1.R2.key.csv | head -n 6)
debug_capture echo $key_head

# Check output for correct processing, 50 keys, 25 unique sequences
test_consol_key_len=$(cat etc/test_output/testSeq-1.R2.key.csv | sed '/readNames/d' | wc -l)
test_consol_seq_len=$(zcat etc/test_output/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | wc -l)

if [[ $test_consol_key_len != 50 || $test_consol_seq_len != 25 ]]; then
    testing_error "    test 1/2: FAILED (output)"
else
    info "    test 1/2: passed (output)"
fi

# Stat file
consol_stat=$(cat etc/test_output/*.consol.stat.csv > /dev/null && echo true || echo false)

if [[ $consol_stat = false ]]; then
    testing_error "    test 2/2: FAILED (stats)"
else
    info "    test 2/2: passed (stats)"
fi


### Couple Tests ###

info "Testing couple ..."

# Test script for functionality
couple_str="etc/test_data/testSeq-1.R2.psl.gz etc/test_data/testSeq-1.R1.psl.gz -k etc/test_output/testSeq-1.R2.key.csv etc/test_output/testSeq-1.R1.key.csv -o etc/test_output/testSeq-1.uniq.csv --condSites etc/test_output/testSeq-1.cond.csv --chimera etc/test_output/testSeq-1.chimera.rds --multihit etc/test_output/testSeq-1.multihit.rds --refGenome hg38 --stat etc/test_output/testSeq-1.couple.stat.csv"

debug_capture __nuckit_cmd couple "$couple_str" 2>&1

uniq_head=$(head etc/test_output/testSeq-1.uniq.csv 2>&1)
debug_capture echo $uniq_head

cond_head=$(head etc/test_output/testSeq-1.cond.csv 2>&1)
debug_capture echo $cond_head

# Check output for correct findings, 50 reads / alignments and 5 sites
test_couple_uniq_len=$(cat etc/test_output/testSeq-1.uniq.csv | sed '/seqnames/d' | wc -l)
test_couple_cond_len=$(cat etc/test_output/testSeq-1.cond.csv | sed '/seqnames/d' | wc -l)

if [[ $test_couple_uniq_len != 50 || $test_couple_cond_len != 5 ]]; then
    testing_error "    test 1/2: FAILED (output)"
else
    info "    test 1/2: passed (output)"
fi

couple_stat=$(cat etc/test_output/testSeq-1.couple.stat.csv > /dev/null && echo true || echo false)

if [[ $couple_stat = false ]]; then
    testing_error "    test 2/2: FAILED (stats)"
else
    info "    test 2/2: passed (stats)"
fi


### Clean up ###
info "Cleaning up..."
rm -r etc/test_output

if [[ $__nuckit_env != null ]]; then
    deactivate_env
fi

info "Passed all tests. Complete."
