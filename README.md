# What is taxa-specific-substitution-matrix?
This is a collection of Python implemented bioinformatic algorithms for the generation of substitution matrices from intra-taxa genetic variation data. It allows you to infer and analyse the basic set of substitution matrices that may be genrated from tyour data. These algorithms were implemented and published on human genetic variation data, but may be used widely for any variational dataset.

# What's in here?
## deduce_substiution_matrices.py
This module can be used to convert aggregated genetic variant information into substitution matrices. You may create substitution matrices in all resolutios: nucleotide substitution matrix, codon substitution matrix, and amino-acid substitution matrix.

## calculate_ha.py
This module can be used to construct an amino-acid substitution matrix from a codon substitution matrix which is taken to the specified power.

## symmtery_analysis.py
This module can be used to perform a symmetry analysis on a substitution matrix.

# Basic usage

## deduce_substiution_matrices.py

## calculate_ha.py
You may convert a codon substitution matrix into an amino-acid substitution matrix

construct amino-acid substitution matrix from codon substitution matrix, use build_aa_matrix as follows:
build_aa_matrix(codon_matrix, codon_count_series)

## symmtery_analysis.py

# Please cite 
Evolutionary and Functional Lessons from Human-Specific Amino-Acid Substitution Matrices
Tair Shauli, Nadav Brandes, Michal Linial
bioRxiv 2020.05.09.086009; doi: https://doi.org/10.1101/2020.05.09.086009
