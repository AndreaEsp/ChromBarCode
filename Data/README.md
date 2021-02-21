# Data

### prismr polymers

This folder contains the genome-wide inferred distribution of binding domain types, i.e. the sbs polymer model, for each chromosome of the human GM12878 cell line.

In each row of the files, the 31 columns are the number of binding sites of each type in a 5kb genomic window. The first column corresponds to the number of inert, non-interacting sites along the polymer model. The rows with all zero values correspond to regions with low mappability.

### Genomic coordinates

This folder contains the files with the genomic coordinates of the genome-wide 5kb windows, in the format: chromosome, start coordinate, end coordinate.

### Epigenetic classes

The file labels.txt contains the epigenetic classification of the binding domains of even numbered chromosomes. The file mnemonics.txt specifies the binding domain, its epigenetic class and the class name.

### Epigentic profiles

The file significant_correlations.csv contains the statistical significant correlations among binding domains of the even numbered chromosomes and epigenetic marks. The file control_correlations.csv.gz. is the random control model.

### ENCODE data

The accession numbers of the ENCODE data used in the genome wide analysis are listed in the file ENCODE_accession_numbers.txt.