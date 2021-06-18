# What is taxa-specific-substitution-matrix?
This is a collection of Python implemented bioinformatic algorithms for the generation of substitution matrices from intra-taxa genetic variation data. It allows one to infer a basic set of substitution matrices that may be genrated from one's data. These algorithms were implemented and published on human genetic variation data, but may be used widely for any variational dataset.

## What's in here?
### deduce_substiution_matrices.py
This module can be used to convert aggregated genetic variant information into substitution matrices. You may create substitution matrices in all resolutios: nucleotide substitution matrix, codon substitution matrix, and amino-acid substitution matrix.

### calculate_ha.py
This module can be used to construct an amino-acid substitution matrix from a codon substitution matrix which is taken to the specified power.

### project_utils.py
Utilities to support other modules

# Basic usage

### deduce_substiution_matrices.py
# Usage

`<from deduce_substiution_matrices import *>`

"""
  codon_counts contains counts of the number of occurances for each 
  codon in the specific genomic location at hand. For example: When calculating 
  a human exome substitution matrix, one would use the count of occurnces of each 
  codon in exomic regions of the human genome.
"""
codon_counts = pd.Series({'AAA': 1896,
                          'AAC': 1519,
                          'AAG': 1565,
                          'AAT': 1528,
                          'ACA': 1271,
                          'ACC': 1265,
                          'ACG': 1273,
                          'ACT': 1985,
                          'AGA': 1582,
                          'AGC': 1625,
                          'AGG': 1751,
                          'AGT': 1412,
                          'ATA': 1397,
                          'ATC': 1660,
                          'ATG': 1285,
                          'ATT': 1603,
                          'CAA': 1757,
                          'CAC': 1540,
                          'CAG': 1214,
                          'CAT': 1895,
                          'CCA': 1805,
                          'CCC': 1474,
                          'CCG': 1334,
                          'CCT': 1699,
                          'CGA': 1029,
                          'CGC': 1107,
                          'CGG': 1543,
                          'CGT': 1513,
                          'CTA': 1970,
                          'CTC': 1106,
                          'CTG': 1895,
                          'CTT': 1258,
                          'GAA': 1874,
                          'GAC': 1553,
                          'GAG': 1016,
                          'GAT': 1248,
                          'GCA': 1625,
                          'GCC': 1924,
                          'GCG': 1608,
                          'GCT': 1829,
                          'GGA': 1064,
                          'GGC': 1783,
                          'GGG': 1999,
                          'GGT': 1538,
                          'GTA': 1165,
                          'GTC': 1181,
                          'GTG': 1239,
                          'GTT': 1972,
                          'TAC': 1825,
                          'TAT': 1896,
                          'TCA': 1126,
                          'TCC': 1228,
                          'TCG': 1990,
                          'TCT': 1842,
                          'TGC': 1331,
                          'TGG': 1885,
                          'TGT': 1400,
                          'TTA': 1747,
                          'TTC': 1898,
                          'TTG': 1122,
                          'TTT': 1329}) 

"""
  Calculate a codon substitution matrix according to aggregated variant 
  information, and normliaze rows to embody the estimted contiodional 
  substitution probabilites
"""
calculate_hc1(variant_df, normalize=True, codon_counts=codon_counts)

"""
  Calculate a codon substitution matrix according to aggregated variant 
  information without normalization
"""
calculate_hc1(variant_df, normalize=False)

"""
Calculate a nucleotide substitution matrix according to aggregated variant 
information, and normliaze rows to embody the estimted contiodional substitution 
probabilites. codon_counts contains counts of the number of occurances for each 
codon in the specific genomic location at hand. For example: When calculating 
a human exome substitution matrix, one would use the count of occurnces of each 
codon in exomic regions of the human genome.
"""

calculate_hn1(variant_df, normalize=True, codon_counts=codon_counts)

"""
  Calculate a nucleotide substitution matrix according to aggregated variant 
  information without normalization
"""
calculate_hn1(variant_df, normalize=False)>`

### calculate_ha.py
You may convert a codon substitution matrix into an amino-acid substitution matrix

construct amino-acid substitution matrix from codon substitution matrix, use build_aa_matrix as follows:
build_aa_matrix(codon_matrix, codon_count_series)

# Please cite 
Evolutionary and Functional Lessons from Human-Specific Amino-Acid Substitution Matrices
Tair Shauli, Nadav Brandes, Michal Linial
bioRxiv 2020.05.09.086009; doi: https://doi.org/10.1101/2020.05.09.086009
