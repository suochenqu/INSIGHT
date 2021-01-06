This folder contains an example script used for Stage 2 sequencing result analysis using XPRIZE data.

# File description:
* all_pool: folder containing the 'all pool' sequencing results in .fastq format, only forward reads are needed, file can be gzipped
* neg_pool: folder containing the 'negative pool' sequencing results in .fastq format, only forward reads are needed, file can be gzipped
* barcode_used.csv: information of used barcodes in all reactions, including barcode ID, left barcode sequence, right barcode sequence and stage 1 NASBA test results
* extract_barcode_reads.R: R script used to extract barcodes from sequencing results
* analyse_results.R: R script used to generate stage 2 results and combined results
* process.sh: bash script to streamline the analysis

# Example usage:
./process.sh barcode_used.csv all_pool/1_S1_L001_R1_001.fastq.gz neg_pool/1_S1_L001_R1_001.fastq.gz