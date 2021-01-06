#!/bin/bash
set -e # exit when any command fails

if [ "$#" -gt 4 ] || [ "$#" -lt 3 ]; then
    echo "Usage: $0 barcode_info_file all_pool_reads_file neg_pool_reads_file [output_file]"
    exit 1
fi

# extract barcode reads from all pool sequencing file to 'all_pool_reads.csv'
Rscript extract_barcode_reads.R "$2" all_pool_reads.csv

# extract barcode reads from negative pool sequencing file to 'neg_pool_reads.csv'
Rscript extract_barcode_reads.R "$3" neg_pool_reads.csv

# analyse read counts for barcode pairs used to obtain test results
if [ "$#" -eq 4 ]; then
    Rscript analyse_results.R "$1" all_pool_reads.csv neg_pool_reads.csv "$4"
else
    Rscript analyse_results.R "$1" all_pool_reads.csv neg_pool_reads.csv
fi

