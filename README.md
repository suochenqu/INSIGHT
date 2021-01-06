# INSIGHT
This repository contains the code and data used in the paper

Wu, Q., Suo, C., Brown, T., Wang, T., Teichmann, S. A. and Bassett, A. R. (2020+) INSIGHT: a population scale COVID-19 testing strategy combining
point-of-care diagnosis with centralised high-throughput sequencing. doi: 10.1101/2020.06.01.127019

Main files:
* `code/barcode_hopping_modelling.R` contains the code used to analyse sequencing results and model barcode hopping
* `code/figures.R` contains the code used to generate all figures in the paper
* `code/hamming.py` contains the code to generate barcodes separated by any given Hamming distance. See Supplementary Note 1 of the above referenced paper for more details.
* `code/LoD.R` contains the code used to compute LoD-95 of the INSIGHT testing technology
