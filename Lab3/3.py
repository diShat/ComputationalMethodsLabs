# Calc Methods Lab 3 (Solving Linear Systems with iterative methods) by Diana Shatkovska FB-92, Variant 3
from numpy import array, zeros, subtract, matmul

# partially diagonal dominant
a = array([[6.47, -1.33, 0, 0.75],
                [0.15, 1.26, 0.74, 0.18],
                [3.77, 7.21, 8.04, 2.28],
                [2.06, 2.39, 0.95, 6.81]], float)

b = array([-16.35, -3.24, 15.45, -5.11], float)

rows = len(b)
cols = len(a[0])


print(f"The given matrix is : \n{a}")
print(f"The b vector is : \n{b}")
print('--------------------------------')


e = 0.0001

x0 = [b[i]/a[i][i] for i in range(cols)]
x = zeros(cols)

itr = 0
condition = True
while condition:
    itr += 1
    print(itr, x0)
    for i in range(cols):
        s1 = sum([a[i][j]*x[j] for j in range(i)])
        s2 = sum([a[i][j]*x0[j] for j in range(i+1, cols)])

        x[i] = round((b[i]-(s1+s2))/a[i][i],4)
    if abs(x[0]-x0[0]) < e: break
    x0 = x.copy()
    r_vector = abs(subtract(b, matmul(a, x)))
    print(f"Вектор нев'язки: {[round(num, 4) for num in r_vector]}")
print(f"The solution is:\n{x}")

# calculating r = |b-Ax|
r_vector = abs(subtract(b, matmul(a, x)))
print(f"Вектор нев'язки: {[round(num,4) for num in r_vector]}")
