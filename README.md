# ChromBarCode

## Genome-wide analyses of the binding domains underlying chromatin folding

Repository for scripts and data of the project.

All the codes are written in Python 3.

The folder Data contains the results of the genome wide analysis and the file necessary to run the scripts.

Please consider the following dependencies:

* In the directory Architectural_analysis, the script most_contributing_domain.py generates a file named '1st_contribution.txt', which is used as input for the scripts most_abundant_domain.py and plot_most_contributing_matrix.py.

* In the directory Epigenetic_analysis, the script compute_correlation.py requires ENCODE files as input. The list of the employed ENCODE accession numbers can be found in Data/ENCODE_accession_numbers.txt.