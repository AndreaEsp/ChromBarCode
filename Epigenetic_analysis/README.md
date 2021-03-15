# Epigenetic analysis

### This contains scripts to run the epigenetic analysis of the binding domains.

The script compute_correlation.py requires ENCODE data as input. Accession numbers for the chromatin marks considered in the project can be found in Data/ENCODE_accession_numbers.txt. After downloading, data must be binned by summing the signal within each of the 5kb window listed in Data/genomic_coordinates (e.g. using bedtools map tool).