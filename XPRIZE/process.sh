#!/bin/bash
set -e # exit when any command fails
Zscore_thresh=10 # change here if a different cutoff is needed

if [ "$#" -gt 3 ] || [ "$#" -lt 2 ]; then
    echo "Usage: $0 barcode_info_file all_pool_reads_file [neg_pool_reads_file]"
    exit 1
fi

# extract barcode reads from all pool sequencing file to 'all_pool_reads.csv'
Rscript extract_barcode_reads.R "$2" all_pool_reads.csv

# extract barcode reads from negative pool sequencing file to 'neg_pool_reads.csv'
if [ "$#" -eq 3 ]; then
    Rscript extract_barcode_reads.R "$3" neg_pool_reads.csv
fi

# analyse read counts for barcode pairs used to obtain test results
if [ "$#" -eq 3 ]; then
    Rscript analyse_results.R "$1" "$Zscore_thresh" all_pool_reads.csv neg_pool_reads.csv 
else
    Rscript analyse_results.R "$1" "$Zscore_thresh" all_pool_reads.csv 
fi

