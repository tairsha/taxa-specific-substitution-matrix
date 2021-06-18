# What is taxa-specific-substitution-matrix?
This is a collection of Python implemented bioinformatic algorithms for the generation of substitution matrices from intra-taxa genetic variation data. It allows one to infer a basic set of substitution matrices that may be genrated from one's data. These algorithms were implemented and published on human genetic variation data, but may be used widely for any variational dataset.

### Usage examples can be found in usage_examples.ipynb

## What's in here?
### deduce_substitution_matrices.py
This module can be used to convert aggregated genetic variant information into substitution matrices. You may create substitution matrices in all resolutios: nucleotide substitution matrix, codon substitution matrix, and amino-acid substitution matrix.

### calculate_ha.py
This module can be used to construct an amino-acid substitution matrix from a codon substitution matrix which is taken to the specified power.

### project_utils.py
Utilities to support other modules


# Please cite 
Evolutionary and Functional Lessons from Human-Specific Amino-Acid Substitution Matrices
Tair Shauli, Nadav Brandes, Michal Linial
bioRxiv 2020.05.09.086009; doi: https://doi.org/10.1101/2020.05.09.086009
