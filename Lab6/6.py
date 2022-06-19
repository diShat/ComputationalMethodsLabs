import numpy as np
import matplotlib.pyplot as plt


def func(x, y):
    return (1 - x ** 2) * y + x*(2 - x + x**3)


H = 0.1
X = np.arange(0, 3, H)


def runge_kutta(f, x, h):
    y = [0]
    for i in range(0, len(x)-1):
        k1 = h * f(x[i], y[i])
        k2 = h * f(x[i] + 0.5*h, y[i] + 0.5*k1)
        k3 = h * f(x[i] + 0.5*h, y[i] + 0.5*k2)
        k4 = h * f(x[i] + h, y[i] + k3)
        y.append(y[i] + (k1+2*k2+2*k3+k4)/6)
    return y


def adams_bashfort(f, x, h):
    y = runge_kutta(f, x, h)
    for i in range(3, len(x)-1):
        y[i+1] = y[i] + (55 * f(x[i], y[i]) - 59 * f(x[i-1], y[i-1]) +
        37 * f(x[i-2], y[i-2]) - 9 * f(x[i-3], y[i-3])) * h / 24
    return y

exact_res = [x*x for x in X]
rk_res = runge_kutta(func, X, H)
ab_res = adams_bashfort(func, X, H)

plt.figure(figsize=(15, 10))
plt.title('Results')
plt.plot(X, exact_res, label="Exact res", linewidth=5)
plt.plot(X, rk_res, label="Runge-Kutta", linewidth=5, linestyle='dashed')
plt.plot(X, ab_res, label="Adams-Bashforth", linewidth=5, linestyle='dotted')
plt.legend(loc="upper left")
plt.show()


# get the calculation error
rk_e = [abs(y0 - y1) for y0,y1 in zip(exact_res, rk_res)]
print(f"Runge-Kutta error max: {max(rk_e).round(10)}, average: {np.mean(rk_e).round(10)}")

ab_e = [abs(y0 - y2) for y0,y2 in zip(exact_res, ab_res)]
print(f"Adams-Bashforth error max: {max(ab_e).round(10)}, average: {np.mean(ab_e).round(10)}")


plt.figure(figsize=(15, 10))
plt.title('Errors')
plt.plot(X, rk_e, label="Runge-Kutta", linewidth=5)
plt.plot(X, ab_e, label="Adams-Bashforth", linewidth=5, linestyle='dashed')
plt.legend(loc="upper left")
plt.show()