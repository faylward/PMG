# PMG

This repo contains code used for the identification of phylogenetic marker genes (PMGs) in metagenomic datasets. This code is under active development. Eventually it will contain methods for the prediction of a variety of PMGs, but at the moment it is focused on benchmarking for the rpoB gene, or the RNA polymerase beta subunit. 

A workflow explaining the general methodology is available here (see workflow.png). 
The general methodology entails 1) prediction of Open Reading Frames (ORFs) using Prodigal, 2) Initial searches against PMG Hidden Markov Models using HMMER3, 3) Analysis of the gene neighborhoods of putative PMGs to identify fragmented ORFs, 4) The joining of fragmented ORFs, and 5) the final classification of the PMGs using a RandomForest machine learning algorithm currently implemented in R. The code currently available here uses a training set for the rpoB gene that has been manually-annotated from the annotations recovered from a set set of ~1800 Prokaryotic genomes. 

**Development Plans:**
1) Currently we employ the RandomForest machine learning algorith, but additional tests using Deep Learning may improve results. Several freely available packages use deep learning, such as MXNetR and deepr, and these can be integrated into the workflow here. 
2) Addition of multiple annotation steps. Currently we use the HMMs available from Sunagawa et al., Nature Methods, 2013 (DOI:10.1038/nmeth.2693), but HMMs for these PMGs are also present in the TIGRfam and Pfam databases. Addition of multiple annotation strategies may improve the results of the machine learning classification. 
3) Further work could be done to refine the gene neighborhood analyses. Currently we use this mainly for identifying fragmented ORFs, but this could be easily expanded to include genomic context information in the ML classification. For example, rpoB is often found adjacent to rpoA, and this contextual information can help inform annotations. Many ribosomal proteins (which are also good PMGs) are also syntenic, so there is reason to believe this approach could be broadly applicable. 
4) Application of these methods to publicly-available metagenomes would also be a useful benchmark for future work. 
