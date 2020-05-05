import numpy as np


def ignore_synonymous(substitution_matrix):
    """
    Ignores synonymous probabilities by zeroing diagonal and applying row normalisation. This matrix would now represent the transition
    probabilities for only non-synonymous transitions
    :param substitution_matrix: pandas DataFrame of substitution matrix
    :return: substitution matrix representing the transition probabilities for only non-synonymous transitions
    """
    substitution_matrix = substitution_matrix.copy(deep = True)

    np.fill_diagonal(substitution_matrix.values, 0)
    for i, row in substitution_matrix.iterrows():
        substitution_matrix.loc[i] = row / sum(row)

    return substitution_matrix


def divide_diagonal(matrix):
    assert matrix.index.tolist() == matrix.columns.tolist()
    return matrix/matrix.transpose()


def get_symmetry_matrices(substitution_matrix):
    substitution_matrix = ignore_synonymous(substitution_matrix)
    substitution_matrix_symmetry = divide_diagonal(substitution_matrix)
    return substitution_matrix_symmetry
