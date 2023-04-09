## Title: 

May, J.A., Feng, Z., Adamowicz, S.J. (2023). The problem with missingness in mixed-type trait datasets and what we can do about it.

## Authors:

Jacqueline A. May. Department of Integrative Biology & Biodiversity Institute of Ontario, University of Guelph, Guelph, Ontario, Canada. 

Zeny Feng. Department of Mathematics & Statistics, University of Guelph, Guelph, Ontario, Canada.

Sarah J. Adamowicz. Department of Integrative Biology & Biodiversity Institute of Ontario, University of Guelph, Guelph, Ontario, Canada. 

## Study description:
  
These data files are associated with the manuscript entitled: "The problem with missingness in mixed-type trait datasets and what we can do about it" by Jacqueline A. May, Zeny Feng, and Sarah J. Adamowicz. This project entailed an evaluation of imputation method performances using real data at a large taxonomic scale. We aimed to ascertain whether there was a best-performing method across different taxa, trait types (biological, environmental), and data types (numerical, categorical). Missing data were simulated based on patterns in the original datasets, and imputation performances for four methods (mean/mode replacement, k-nearest neighbour, random forests, and multivariate imputation by chained equations) were evaluated with and without phylogeny in the form of gene trees. Included in this dataset are the sequence alignments used to build the phylogenetic trees, the phylogenetic trees in Newick format, and the corresponding sequence IDs.

## Data sources:

Cytochrome c oxidase subunit 1 (COI) sequence data for this project were obtained from the Barcode of Life Data System (BOLD):

> Ratnasingham, S., & Hebert, P. D. N. (2007). bold: The Barcode of Life Data System (http://www.barcodinglife.org). Molecular Ecology Notes, 7(3), 355–364. https://doi.org/10.1111/j.1471-8286.2007.01678.x. 

Nuclear sequence data for this project were obtained from:

> Burleigh JG, Kimball RT, Braun EL. 2015b. Data from: Building the avian tree of life using a large-scale, sparse supermatrix. Dryad Dataset. https://doi.org/10.5061/dryad.r6b87.

> Pyron, R. A., Burbrink, F. T., & Wiens, J. J. (2013). Data from: A phylogeny and revised 	classification of Squamata, including 4161 species of lizards and snakes. Dryad Dataset. 	https://doi.org/10.5061/dryad.82h0m

> Rabosky DL, Chang J, Title PO, Cowman PF, Sallan L, Friedman M, Kaschner K, Garilao C, Near TJ, Coll M, et al. 2019. Data from: An inverse latitudinal gradient in speciation rate for marine fishes. https://doi.org/10.5061/dryad.fc71cp4.

> Upham NS, Essenstyn JA, Jetz W. 2019b. Data from: Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. Dryad Dataset. https://doi.org/10.5061/dryad.tb03d03.


## File list:

*_COI.fasta (x20)

Chiroptera_COI_FinalImp.fasta

*.tre (x26)

*_ultra.tre (x26)

Chiroptera_COI_FinalImp.tre

Chiroptera_COI_FinalImp_ultra

Rodentia_IRBP_FinalImp.tre

Rodentia_IRBP_FinalImp_ultra.tre

COIProcessIDs.xlsx

GenBankIDs.xlsx


## File descriptions:

### Alignments:

Files with the extension "_COI.fasta" are multiple sequence alignment of COI sequence records in FASTA format. These are the centroid sequences that were selected to represent each species included in simulations of missingness and evaluations of imputation performance. The Chiroptera_COI_FinalImp.fasta file are the aligned COI centroid sequences used to build the phylogenetic tree for final imputation of the Chiroptera mixed-type trait dataset.

COI sequences were aligned using the R package “DECIPHER” v. 2.18.1:

> Wright ES. DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics. 2015 Oct 6;16(1):322.

> Wright ES. RNAconTest: comparing tools for noncoding RNA multiple sequence alignment based on structural consistency. RNA. 2020 May 1;26(5):531–40.



### Trees:

Files with the extension ".tre" are the gene trees built using either the centroid COI sequence alignments or sequences obtained from the aforementioned published nuclear multigene alignments. 

These trees are in Newick format and were built using RAxML v. 8:

> Stamatakis A. RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics. 2014;30:1312–3. 

Files with the extension "ultra.tre" are the ultrametric versions of the trees. They were made ultrametric using the "chronos" function in the R ape package:

> Paradis E, Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics. 35:526–528.

The files "Chiroptera_COI_FinalImp.tre" and "Rodentia_IRBP_FinalImp.tre" are the gene trees built for final imputation of the original datasets for Chiroptera and Rodentia, respectively. "Chiroptera_COI_FinalImp_ultra.tre" and "Rodentia_IRBP_FinalImp_ultra.tre" are the ultrametric versions of these trees.


### Sequence identifiers:

The "COIProcessIDs.xlsx" and "GenBankIDs.xlsx"" spreadsheets contain the BOLD COI process IDs and GenBank accession numbers that correspond to the sequence data records used to build the phylogenetic trees described above.
