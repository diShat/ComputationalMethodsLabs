# Calc Methods Lab 5 (Interpolation) by Diana Shatkovska FB-92, Variant 3

import numpy as np
import matplotlib.pyplot as plt
from math import exp

np.set_printoptions(suppress=True, precision=5)


def given_func(x):
    return (x**2)*exp(x)


X = np.array([round(x,1) for x in np.arange(-3,4.1,0.7)], float).round(1)
print(f"x: {X}")
Y = np.array([given_func(x) for x in X], float)
print(f"y: {Y}")
print("-"*20)


# define some operations over linear systems beforehand (functions from lab2)
def elimination(matrix, b_vector):
    n = len(b_vector)
    for k in range(n-1):
        for i in range(k+1, n):
            if matrix[i][k] == 0:
                continue
            factor = matrix[k][k]/matrix[i][k]
            for j in range(k, n):
                matrix[i][j] = matrix[k][j] - matrix[i][j]*factor
            b_vector[i] = b_vector[k] - b_vector[i]*factor
    return matrix, b_vector


def solve(triangled_matrix, b_vector):
    n = len(b_vector)
    x_vector = [0] * n
    x_vector[n-1] = b_vector[n-1]/triangled_matrix[n-1][n-1]

    for i in range(n-2, -1, -1):
        sum_ax = 0
        for j in range(i+1, n):
            sum_ax += triangled_matrix[i][j] * x_vector[j]
        x_vector[i] = (b_vector[i] - sum_ax) / triangled_matrix[i][i]
    return x_vector


# spline method

a = Y.copy()
N = len(X)
# finding c_і coefitients
c = np.zeros(N)
b = np.zeros(N)
h = (X[-1] - X[0]) / (N - 1)
for i in range(1, N - 1):
    b[i] = (3 / h) * (((a[i + 1] - a[i]) / h) - ((a[i] - a[i - 1]) / h))
linsys = np.zeros((N,N))
linsys[0][0] = 1
linsys[N - 1][N - 1] = 1
for i in range(1, N - 1):
    linsys[i][i - 1] = 1
    linsys[i][i] = 4
    linsys[i][i + 1] = 1
linsys, b = elimination(linsys,b)
c = solve(linsys, b)

# finding d_і coefitients
d = np.zeros(N-1)
for i in range(N - 1):
    d[i] = (c[i + 1] - c[i]) / (3 * h)

# finding b_і coefitients
b = np.zeros(N-1)
for i in range(1, N):
    b[i - 1] = (a[i] - a[i - 1]) / h - (2 * c[i - 1] + c[i]) * h / 3


def getSplineRes(a,b,c,d,x):
    spline = 0
    for i in range(N):
        if X[i] <= x <= X[i + 1]:
            spline = i
            break
    res = a[spline]+b[spline]*(x-X[spline])+c[spline]*(x-X[spline])**2+d[spline]*(x-X[spline])**3
    return res


print(f"Spline coefs:"
      f"a: {a}\nb: {b}\nc: {c}\nd: {d}")
res = getSplineRes(a,b,c,d, 2)
print(f"My spline res for x=2: {res}")


# compare to scipy's function
from scipy.interpolate import CubicSpline
sc_spline = CubicSpline(X,Y)
print(f"Scipy spline res for x=2: {sc_spline(2)}")


# get the calculation error
x_check = np.linspace(min(X), max(X), 100)
e = [abs(given_func(x) - getSplineRes(a,b,c,d, x)) for x in x_check]
print(f"Error min: {min(e).round(13)}, max: {max(e).round(10)}, average: {np.mean(e).round(10)}")

plt.figure(figsize=(15, 10))
plt.ylim(min(e), max(e))
plt.title('Spline error')
plt.plot(x_check, e, label="Orig", linewidth=5)


# plotting the figures
plt.figure(figsize=(15, 10))
plt.scatter(X, Y)

# given function
x_0 = np.linspace(min(X), max(X), 100)
y_0 = [given_func(x) for x in x_0]
plt.plot(x_0, y_0, label="Orig", linewidth=5)

# Spline implementation
x_1 = np.linspace(min(X), max(X), 100)
y_1 = [getSplineRes(a,b,c,d, x) for x in x_1]
plt.plot(x_1, y_1, label="Mine", linestyle='dashed', linewidth=5)

# scipy's function
x_2 = np.linspace(min(X), max(X), 100)
y_2 = [sc_spline(x) for x in x_2]
plt.plot(x_2, y_2, label="Scipy's", linestyle='dotted', linewidth=5)

plt.title('Spline demo')
plt.legend(loc="upper left")
plt.show()


