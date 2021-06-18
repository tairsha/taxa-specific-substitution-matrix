import pandas as pd
import project_utilities as utils

# string headers of variant DataFrame
REF_CODON = 'REF_CODON'
ALT_CODON = 'ALT_CODON'

REF_NUC = 'REF_NUCLEOTIDE'
ALT_NUC = 'ALT_NUCLEOTIDE'

AF = 'ALLELE_FREQUENCY'


def count_af(substitution_matrix, variant_df, subject_lst, subject_headers):
    """
    Calculate substitution matrix, not normalized.
    :param substitution_matrix: empty DataFrame which will contain summed substitution frequencies
    :param variant_df: DataFrame which contains aggregated variant information: for each genomic position,
    it should contain at least the source codon, target codon, source nucleotide, target nucleotide, and the variant's allele frequency
    :param subject_lst: list of codons\nucleotides
    :param subject_header: tuple of header of reference and alternative columns
    :return: not normalized substitution matrix
    """
    ref_str, alt_str = subject_headers  # get header string for DataFrame column access
    for ref in subject_lst:
        ref_df = variant_df[variant_df[ref_str] == ref]  # use only variant in which the reference is ref
        for alt in subject_lst:
            if alt == ref:  # skip synonymous as data does not contain synonymous substitutions
                continue
            ref_alt_df = ref_df[ref_df[alt_str] == alt]  # use only variant in which the alternative is alt
            substitution_matrix.at[ref, alt] = ref_alt_df[
                AF].sum()  # sum corresponding allele frequencies into substitution matrix

    return substitution_matrix


def get_nuc_counts(codon_counts):
    """
    Count the number of appearances of each nucleotide in the reference genome
    :param codon_counts: A pandas Series containing the number of appearances of each codon in the reference genome
    :return:  pandas Series containing the number of appearances of each nucleotide in the reference genome
    """
    nuc_counts = pd.Series(0, index=utils.nucleotides_lst)
    for codon in codon_counts.index:
        for nuc in codon:
            nuc_counts[nuc] += codon_counts.loc[codon]
    return nuc_counts


def normalise_matrix(substitution_matrix, counts):
    """
    Row-normalisation of the substitution matrix
    :param substitution_matrix: DataFrame which contains summed substitution frequencies
    :param counts: a pandas Series containing the number of appearances of each codon, for normalization
    :return: row-normalised DataFrame
    """
    for i in substitution_matrix.index:
        curr_counts = int(counts[i])
        if curr_counts != 0:
            substitution_matrix.loc[i, :] /= curr_counts
            substitution_matrix.at[i, i] = 1 - substitution_matrix.loc[i, :].sum()
        else:
            substitution_matrix.loc[i, :] = 0

    return substitution_matrix


def calculate_substitution_matrix(variant_df, subject_lst, subject_header, normalize=True, counts=None):
    """
    calculate a nucleotide/codon substitution matrix according to aggregated variant information
    :param variant_df: DataFrame which contains aggregated variant information: for each genomic position,
    it should contain at least the source codon, target codon, source nucleotide, target nucleotide, and the variant's allele frequency
    :param subject_lst: list of codons\nucleotides
    :param subject_header: tuple of header of reference and alternative columns
    :param counts: a pandas Series containing the number of appearances of each codon/nucleotide, for normalization. Defult value is None.
    Will be used only if normalize is set to True.
    :param normalize: boolean indicating whether the substitution matrix is to be normalized.
    :return: substitution matrix
    """

    substitution_matrix = utils.create_zero_matrix(subject_lst)
    substitution_matrix = count_af(substitution_matrix, variant_df, subject_lst, subject_header)

    if normalize:
        substitution_matrix = normalise_matrix(substitution_matrix, counts)

    return substitution_matrix


def calculate_hc1(variant_df, normalize=True, codon_counts=None):
    """
    calculate a codon substitution matrix according to aggregated variant information
    :param variant_df: DataFrame which contains aggregated variant information: for each genomic position,
    it should contain at least the source codon, target codon, source nucleotide, target nucleotide, and the variant's allele frequency
    :param codon_counts: a pandas Series containing the number of appearances of each codon, for normalization. Defult value is None.
    Will be used only if normalize is set to True.
    :param normalize: boolean indicating whether the substitution matrix is to be normalized.
    :return: codon substitution matrix
    """

    return calculate_substitution_matrix(variant_df,
                                         utils.codon_lst_non_stop,
                                         (REF_CODON, ALT_CODON),
                                         normalize=normalize,
                                         counts=codon_counts)


def calculate_hn1(variant_df, normalize=True, codon_counts=None):
    """
    calculate a nucleotide substitution matrix according to aggregated variant information
    :param variant_df: DataFrame which contains aggregated variant information: for each genomic position,
    it should contain at least the source codon, target codon, source nucleotide, target nucleotide, and the variant's allele frequency
    :param normalize: boolean indicating whether the substitution matrix is to be normalized.
    :param nuc_counts: a pandas Series containing the number of appearances of each nucleotide, for normalization. Defult value is None.
    Will be used only if normalize is set to True.
    :return: nucleotide substitution matrix
    """

    if normalize:
        nuc_counts = get_nuc_counts(codon_counts)
    else:
        nuc_counts = None
    return calculate_substitution_matrix(variant_df,
                                         utils.nucleotides_lst,
                                         (REF_NUC, ALT_NUC),
                                         normalize=normalize,
                                         counts=nuc_counts)
