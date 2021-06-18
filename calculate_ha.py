import numpy as np
import pandas as pd

import project_utilities as utils


def matrix_power(hc1, power):
    return pd.DataFrame(np.linalg.matrix_power(hc1, power), index = hc1.index, columns = hc1.columns)


def aa_info(aa, codon_count_series, aa_sum_ref):
    """
    given an aa a, return a list of tuples: first in the tuple is a codon that codes for a,
    and second is the codon frequency wrt all other codons that code for a
    :param aa_sum_ref: a Series or dictionary of the sum of counts of codon appearances according to aa
    :param codon_count_series: a pandas Series containing the number of appearances of each codon in the refrenence genome
    :param aa: a string representing an amino-acid
    :return: dict of codon string and it's amino-acid relative frequency
    """
    # calculate each codon's frequency, within its multiplicity
    return {c: codon_count_series[c]/aa_sum_ref[aa] for c in utils.aa_to_codon_dct[aa]}


def build_aa_matrix(codon_matrix, codon_count_series):
    """
    Construct amino-acid substitution matrix from codon substitution matrix
    :param codon_matrix: codon substitution matrix pandas DataFrame
    :param codon_count_series: a pandas Series containing the number of appearances of each codon in the refrenence genome
    :return:
    """
    # init empty matrices
    aa_matrix = pd.DataFrame(0, index = utils.aa_lst, columns = utils.aa_lst, dtype = float)
    temp_matrix = pd.DataFrame(0, index = utils.aa_lst, columns = utils.codon_lst_non_stop)

    # deduce the amino-acid's codon usage distribution
    aa_sum_ref = {aa: 0 for aa in utils.aa_lst}
    for aa in utils.aa_lst:
        for codon in utils.aa_to_codon_dct[aa]:
            aa_sum_ref[aa] += codon_count_series[codon] # calculate codon usage distribution from the codon counts for normalization

    aa_codon_distribution_ref_dct = {aa: aa_info(aa, codon_count_series, aa_sum_ref) for aa in utils.aa_lst}

    # update rows of substitution matrix
    for aa in utils.aa_lst:
        for codon in utils.aa_to_codon_dct[aa]:
            temp_matrix.loc[aa, :] += codon_matrix.loc[codon, :] * aa_codon_distribution_ref_dct[aa][codon]

    # update columns of substitution matrix
    for aa in utils.aa_lst:
        for codon in utils.aa_to_codon_dct[aa]:
            aa_matrix.loc[:, aa] += temp_matrix.loc[:, codon]

    return aa_matrix


def calc_ha(codon_substitution_matrix, power, codon_count_series):
    """
    Construct amino-acid substitution matrix from codon substitution matrix which is taken to the specified power
    :param codon_substitution_matrix: pandas DataFrame, matrix should be of the first power
    :param power: power that HC1 is to be taken to before converting to amino-acid resolution
    :param codon_count_series: a pandas Series containing the number of appearances of each codon in the reference genome
    :return: amino-acid substitution matrix to the specified power in a pandas DataFrame
    """
    codon_substitution_matrix = matrix_power(codon_substitution_matrix, power)
    return build_aa_matrix(codon_substitution_matrix, codon_count_series)
