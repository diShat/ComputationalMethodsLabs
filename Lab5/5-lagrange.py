# Calc Methods Lab 5 (Interpolation) by Diana Shatkovska FB-92, Variant 3

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polymul
from math import exp

np.set_printoptions(suppress=True, precision=5)


def given_func(x):
    return (x**2)*exp(x)


X = np.array([round(x,1) for x in np.arange(-3,4.1,0.7)], float).round(1)
print(f"x: {X}")
Y = np.array([given_func(x) for x in X], float)
print(f"y: {Y}")
print("-"*20)


# Lagrange method
def getLagrangePolynomial(xs,ys):
    n = len(xs)
    coefs = np.zeros(len(xs))
    for k in range(n):
        p = 1.
        for j in range(n):
            if k != j:
                p = polymul(p,[1/(xs[k]-xs[j]),-xs[j]/(xs[k]-xs[j])])
        coefs += np.array(p)*ys[k]
    return coefs


def getLagrangeRes(coefs,x):
    result = 0
    for i,coef in enumerate(coefs[::-1]):
        result += coef*(x**i)
    return result


poly = getLagrangePolynomial(X,Y)
print(f"My poly coefs: {poly}")
res = getLagrangeRes(poly, 3)
print(f"My poly res for x=3: {res}")


# compare to scipy's function
from scipy.interpolate import lagrange
sc_poly = lagrange(X,Y)
print(sc_poly)
print(f"Scipy poly coefs: {[round(el, 5) for el in sc_poly.coef]}")
print(f"Scipy poly res for x=3: {sc_poly(3)}")


# get the calculation error
x_check = np.linspace(min(X), max(X), 100)
e = [abs(given_func(x) - getLagrangeRes(poly, x)) for x in x_check]
print(f"Error min: {min(e).round(13)}, max: {max(e).round(10)}, average: {np.mean(e).round(10)}")

plt.figure(figsize=(15, 10))
plt.ylim(min(e), max(e))
plt.title('Lagrange error')
plt.plot(x_check,e, label="Orig", linewidth=5)


# plotting the figures
plt.figure(figsize=(15, 10))
plt.scatter(X, Y)

# given function
x_0 = np.linspace(min(X), max(X), 100)
y_0 = [given_func(x) for x in x_0]
plt.plot(x_0, y_0, label="Orig", linewidth=5)

# Lagrange implementation
x_1 = np.linspace(min(X), max(X), 100)
y_1 = [getLagrangeRes(poly, x) for x in x_1]
plt.plot(x_1, y_1, label="Mine", linestyle='dashed', linewidth=5)

# scipy's function
poly_lagr = lagrange(X,Y)
x_2 = np.linspace(min(X), max(X), 100)
y_2 = [poly_lagr(x) for x in x_2]
plt.plot(x_2, y_2, label="Scipy's", linestyle='dotted', linewidth=5)

plt.title('Lagrange demo')
plt.legend(loc="upper left")
plt.show()


