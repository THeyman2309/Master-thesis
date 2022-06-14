# Long-read sequencing for the detection of tandem-repeats in patients with intellectual disability and congenital anomalies
This repository was made to store supplementary files regarding the master dissertation of Tobias Heyman. Supplementary files referred to in the thesis can be accessed under the following names. 

* Supplementary file 1: tandem_repeat_positions_hg38.bed

This file contains the exact genomic positions in human reference genome 38 (hg38)for tandem repeats that are known to be involved in diseases. This file was passed on to Straglr, TRiCoLOR and tandem-genotypes to characterize known tandem-repeats in patients with intellectual disability.
        
* Supplementary file 2: tandem_repeat_positions_T2T.bed

This file contains the exact genomic positions in the telomere-to-telomere (T2T) reference genome for tandem repeats that are known to be involved in diseases. This file was passed on to Straglr, TRiCoLOR and tandem-genotypes to characterize known tandem repeats in patients with intellectual disability.

* Supplementary file 3: tg_filter.py

This file contains the code that was written to filter tandem repeats that were detected in the genome-wide analysis by tandem-genotypes based on the length of the tandem repeat.

* Supplementary file 4: Formula1.py

This file contains the code that was written to calculate the difference in reported length of tandem repeats between Straglr and tandem-genotypes. This is thus the implementation of Formula 1 in the thesis. This file requires the detected tandem repeat lengths to be provided in six columns (first three straglr, last three tandem-genotypes, where the order of individuals is identical for straglr and tandem-genotypes per line). 

* Supplementary file 5: Formula2.py

This file contains the code that was written to calculate the difference in reported length of tandem repeats between human reference genome 38 (hg38) and the telomere-to-telomere reference genome (T2T). This is thus the implementation of Formula 2 in the thesis. This file requires the detected tandem repeat lengths to be provided in twelve columns (first six hg38, last six T2T, where the order of individuals is identical for hg38 and T2T per line). 

* Supplementary file 6: Code_R_thesis.pdf

This file contains the R code that was used to construct all figures that do not depict screenshots from tools in the master thesis. Additionally, R was used to construct supplementary Table 1. This pdf file is a knitted rmd file.
