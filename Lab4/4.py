# Calc Methods Lab 4 (Finding eigenvalues, Jacobi Iterative Method) by Diana Shatkovska FB-92, Variant 3
from numpy import array, identity, set_printoptions
from math import sqrt, copysign

set_printoptions(suppress=True, precision=5)

def spherical_norm(mtr):
    rows = len(mtr)
    diagonal_sum = sum([mtr[i][i]*mtr[i][i] for i in range(rows)])
    non_diagonal_sum = sum([mtr[i][j]*mtr[i][j] for i in range(rows) for j in range(rows) if not i == j])
    total_sum = sum([mtr[i][j]*mtr[i][j] for i in range(rows) for j in range(rows)])

    return diagonal_sum, non_diagonal_sum, total_sum


def get_max_nondiagonal(mtr):
    max_i, max_j = 1, 0
    rows = len(mtr)
    for i in range(rows):
        for j in range(rows):
            if not i == j and mtr[i][j] > mtr[max_i][max_j]:
                max_i, max_j = i, j
    return max_i, max_j


def get_translation_matrix(main_matrix):
    i, j = get_max_nondiagonal(main_matrix)
    return get_translation_matrix_windex(main_matrix, i, j)


def get_translation_matrix_windex(main_matrix, i, j):
    result = identity(len(main_matrix))
    mu = 2 * main_matrix[i][j] / (main_matrix[i][i] - main_matrix[j][j])
    c = sqrt(0.5 * (1 + 1 / sqrt(1 + mu * mu)))
    s = sqrt(0.5 * (1 - 1 / sqrt(1 + mu * mu))) * copysign(1, mu)
    result[i][i], result[j][j], result[i][j], result[j][i] = c, c, s, -s
    return result


a = array([    [6.29, 0.97, 1.00, 1.1],
                [0.97, 4.13, 1.30, 0.16],
                [1.00, 1.30, 5.47, 2.10],
                [1.1, 0.16, 2.10, 6.07]], float)


print(f"The given matrix is : \n{a}")
print('--------------------------------')

print (f"Початкова сферична норма: {spherical_norm(a)}")
for i in range(len(a) * len(a)):
    print(f'{i}------')
    print(f"Сферична норма: {spherical_norm(a)[2]}\tДіагональна: {spherical_norm(a)[0]}\tНедіагональна: {spherical_norm(a)[1]}")
    translation_matrix = get_translation_matrix(a)
    print(translation_matrix)
    a = translation_matrix.dot(a).dot(translation_matrix.transpose())

print(f"Result:\n{a}")
