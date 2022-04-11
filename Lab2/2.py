# Calc Methods Lab 2 (Solving Linear Systems) by Diana Shatkovska FB-92, Variant 3
from numpy import array, matmul, subtract

a = array([[8.30, 2.62, 4.10, 1.90],
                [3.92, 8.45, 8.78, 2.46],
                [3.77, 7.21, 8.04, 2.28],
                [2.21, 3.65, 1.69, 6.99]], float)

b = array([-10.65, 12.21, 15.45, -8.35], float)
n = len(b)


print(f"The given matrix is : \n{a}")
print(f"The b vector is : \n{b}")
print('--------------------------------')


# 1st step: elimination
def elimination(matrix, b_vector):
    for k in range(n-1):
        for i in range(k+1, n):
            if matrix[i, k] == 0:
                continue
            factor = matrix[k, k]/matrix[i, k]
            for j in range(k, n):
                matrix[i, j] = matrix[k, j] - matrix[i, j]*factor
            b_vector[i] = b_vector[k] - b_vector[i]*factor
            print(f"{k}-{i} step:\n{matrix}")
            print(f"b vector: \n{b_vector}\n")
    return matrix, b_vector

triangled_a, triangled_b  = elimination(a, b)

print(f"The triangled matrix is : \n{triangled_a}")
print(f"The b vector is : \n{triangled_b}")
print('--------------------------------')


# 2nd step: back-substitution
def back_substitution(triangled_matrix, b_vector):
    x_vector = [0] * n
    x_vector[n-1] = b_vector[n-1]/triangled_matrix[n-1, n-1]

    for i in range(n-2, -1, -1):
        sum_ax = 0
        for j in range(i+1, n):
            sum_ax += triangled_matrix[i, j] * x_vector[j]
        x_vector[i] = (b_vector[i] - sum_ax) / triangled_matrix[i, i]
    return x_vector


x = back_substitution(triangled_a, triangled_b)
print(f"The solution is : \n{x}")

# calculating r = |b-Ax|
r_vector = abs(subtract(b, matmul(a, x)))
# r_vector = [f'{el:.6f}' for el in r_vector]
print(f"Вектор нев'язки: {r_vector}")
